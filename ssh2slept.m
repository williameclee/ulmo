%% SSH2SLEPT - Fits monthly sea surface height anomalies to regional Slepian functions.
%
% Syntax
%   [sshSlept, dates, fit, GVN] = ssh2slept("NASASSHMonthly", domain, L)
%   [__] = ssh2slept(__, timelim, unit="m", truncation=J)
%   [__] = ssh2slept(__, TimeStep="midmonth")
%   [__] = ssh2slept(__, GVN={G, V, N}, dropPartialMonths=true)
%   [__] = ssh2slept(__, OutputFormat="traditional")
%
% Input arguments
%   product - "NASASSHMonthly" or SSHProduct.NASASSHMonthly, prepared by
%       processSSHMonthlyNasa.
%   domain - GeoDomain, region name, cell of GeoDomain constructor arguments,
%       or a two-column polygon of longitude and latitude in degrees.
%   L - Nonnegative integer maximum spherical harmonic degree.
%       The fit includes cells up to 90/L degrees outside the target for
%       positive L. At degree zero it uses only the target.
%   timelim (optional) - Inclusive start and end dates as finite datetimes or
%       serial date numbers in chronological order. Equal bounds select one epoch.
%       Size: two elements.
%       Default range: [] (full record).
%   TimeStep (name-value) - Requested epochs, [] or "midmonth". Both use
%       the native NASASSHMonthly epochs and share the same projection cache.
%       Default timestep: [].
%   unit (name-value) - Height unit for coefficients and residual RMS, "mm" or "m".
%       Default unit: "mm".
%   OutputFormat (name-value) - Coefficient orientation, "timefirst" or "traditional".
%       Default format: "timefirst".
%   truncation (name-value) - Number of leading basis functions to retain.
%       With a full basis, [] uses round(N) bounded by 1 and the column count.
%       With a truncated basis, [] retains all supplied columns.
%       Default count: [].
%   GVN (name-value) - Supplied basis {G, V, N}, as in prepareslepianbasis.
%       A supplied basis forces a fresh fit and disables cache writes.
%       Default basis: {} (compute or load the domain basis).
%   ForceNew (name-value) - Whether to refit from the prepared monthly maps.
%       Use this after replacing the source product.
%       Default flag: false.
%   SaveData (name-value) - Whether to save a new projection and eligible
%       basis and mask files.
%       Default flag: true.
%   dropPartialMonths (name-value) - Whether to omit partial months from
%       the returned coefficients, dates and diagnostics.
%       Default flag: false.
%
% Input file content
%   IFILES/SSH/NASASSHMonthly.mat contains the following variables.
%   sshs - Native monthly sea surface height anomalies. NaN marks missing cells.
%       Size: (nLon, nLat, nEpoch).
%       Unit: metres.
%   lon, lat - Regular longitude and strictly monotone latitude axes.
%       Either axis may descend. A duplicated longitude endpoint counts once.
%       Size: nLon and nLat elements, respectively.
%       Unit: degrees.
%   dates - Strictly increasing monthly midpoints as datetimes.
%       Size: (nEpoch, 1).
%   partialMonths - Logical flag for each month truncated by the source date span.
%       Size: (nEpoch, 1).
%
% Output arguments
%   sshSlept - Slepian coefficients in the requested height unit.
%       Size: (nDates, J) for "timefirst" or (J, nDates) for "traditional".
%   dates - Selected monthly midpoints as datetimes.
%       Size: (nDates, 1).
%   fit - Diagnostic table with one row per returned date. nObs counts fitted
%       cells in the target and buffer. areaFraction measures observed area
%       within the target. condNum is the fit condition number, residualRms
%       is the weighted residual RMS in the requested unit, and partialMonth
%       is the source flag.
%   GVN - Retained basis {G, V, N}, with harmonic coefficients G, concentration
%       eigenvalues V and the full Shannon number N.
%       Size: G is ((L+1)^2, J), V is (J, 1) and N is scalar.
%
% Output file content
%   IFILES/SSH/SlepianExpansions/ssh2slept-<product>-<domainId>-<L>-<J>-m.nc
%       stores the full record before date or partial-month selection.
%   coefficients - Fitted heights in metres, ordered by mode then time in MATLAB.
%   time - Monthly midpoints in seconds since 2000-01-01 00:00:00.
%   nObs, areaFraction, condNum, residualRms, partialMonth - One value per
%       source month, with residualRms in metres and partialMonth stored as uint8.
%   creationDate - File creation time in seconds since 1970-01-01 00:00:00 UTC.
%
% Notes
%   An empty month has NaN coefficients, condition number and residual RMS.
%   A wholly empty record or an underdetermined nonempty month raises an error.
%
% See also
%   ssh2lonlatt, processSSHMonthlyNasa, prepareslepianbasis, fitgriddedslepian
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function [sshSlept, dates, fit, GVN] = ssh2slept(product, domain, L, timelim, options)

    arguments (Input)
        product (1, 1) string {mustBeMember(product, "NASASSHMonthly")}
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        timelim {mustBeTimeRange} = []
        options.TimeStep {mustBeTimeStep} = []
        options.unit ...
            {mustBeTextScalar, mustBeMember(options.unit, ["m", "mm"])} = "mm"
        options.OutputFormat ...
            {mustBeTextScalar, mustBeMember(options.OutputFormat, ["timefirst", "traditional"])} = "timefirst"
        options.truncation double ...
            {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.GVN cell = {}
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.dropPartialMonths (1, 1) logical = false
    end

    % The projection cache represents only the native monthly epochs.
    timestep = options.TimeStep;
    assert(isempty(timestep) || ((ischar(timestep) || isstring(timestep)) && strcmpi(timestep, 'midmonth')), ...
        'ULMO:ssh2slept:TimeStep', 'NASASSHMonthly projections support only native ([]) or midmonth epochs.');

    if ~isempty(timelim)
        if isnumeric(timelim)
            timelim = datetime(timelim, 'ConvertFrom', 'datenum');
        end
        if any(~isfinite(timelim(:)))
            error('ULMO:notValidTimeRange', 'Projection time limits must be finite.');
        end
    end

    if iscell(domain)
        domain = GeoDomain(domain{:});
    elseif ischar(domain) || isstring(domain)
        domain = GeoDomain(domain);
    end

    assert((isa(domain, 'GeoDomain') && isscalar(domain)) || ...
        (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2), ...
        'ULMO:ssh2slept:Domain', 'Supply a GeoDomain, region name, or polygon.');

    if ~isempty(options.GVN)
        options.ForceNew = true;
        options.SaveData = false;
    end

    obsBuf = 0;

    if L > 0
        obsBuf = 90 / L;
    end

    cacheHit = false;

    if ~options.ForceNew
        [path, tmp, cacheHit] = getoutputfile(product, domain, L, options.truncation, options.SaveData, false);
    end

    if cacheHit
        % Filtering needs the flags even when diagnostics are not requested.
        outputCnt = max(nargout, 3 * options.dropPartialMonths);
        [sshSlept, dates, fit] = readCache(path, timelim, outputCnt);

        if nargout >= 4
            J = size(sshSlept, 2);
            GVN = prepareslepianbasis(domain, L, J, {}, SaveData = options.SaveData);
        end

    else
        [sshSlept, dates, fit, GVN] = fitRecord(product, domain, L, options);

        if options.SaveData
            % A forced fit skips lookup. A first default fit determines J here.
            if options.ForceNew || isempty(path)
                [path, tmp] = getoutputfile(product, domain, L, size(sshSlept, 2), true, true);
            end

            defaultJ = min((L + 1) ^ 2, max(1, round(GVN{3})));
            writeCache(path, tmp, product, L, obsBuf, sshSlept, dates, fit, defaultJ);
        end

        rows = selectDates(dates, timelim);
        sshSlept = sshSlept(rows, :);
        fit = fit(rows, :);
        dates = dates(rows);
    end

    if options.dropPartialMonths
        rows = ~fit.partialMonth;
        sshSlept = sshSlept(rows, :);
        dates = dates(rows);
        fit = fit(rows, :);
    end

    if strcmpi(options.unit, "mm")
        sshSlept = sshSlept * 1000;

        if nargout >= 3
            fit.residualRms = fit.residualRms * 1000;
        end

    end

    if strcmpi(options.OutputFormat, "traditional")
        sshSlept = permute(sshSlept, [2, 1]);
    end

end

function [path, tmp, cacheHit] = getoutputfile(product, domain, L, J, saveData, forceNew)
    root = getenv('IFILES');
    assert(isfolder(root), 'ULMO:ssh2slept:MissingRoot', 'Set IFILES to the data directory.');
    cacheHit = false;
    tmp = '';
    folder = fullfile(root, 'SSH', 'SlepianExpansions');

    if isa(domain, 'GeoDomain')
        domainId = domain.Id;
        % Earlier polar projections selected observations in rotated coordinates.
        if ismember(domain.Domain, ["arctic", "antarctica"])
            domainId = [domainId, '-rotb'];
        end

    else
        domainId = hash(domain, 'sha1');
    end

    stem = sprintf('ssh2slept-%s-%s-%d', product, domainId, L);

    if isempty(J)
        % Every truncation for this projection records the same default mode count.
        files = dir(fullfile(folder, [stem, '-*-m.nc']));
        names = string({files.name});
        first = find(matches(names, stem + "-" + digitsPattern + "-m.nc"), 1);

        if isempty(first)
            path = '';
            return
        end

        J = ncreadatt(fullfile(folder, names(first)), '/', 'defaultTruncation');
    end

    path = fullfile(folder, sprintf('%s-%d-m.nc', stem, J));
    cacheHit = ~forceNew && isfile(path);

    if saveData && ~cacheHit

        if ~isfolder(folder)
            mkdir(folder);
        end

        tmp = [tempname(folder), '.nc'];
    end

end

function [coeffs, dates, fit, basis] = fitRecord(product, domain, L, options)
    [values, ~, dates, lon, lat, partialMonth] = ssh2lonlatt(product, options.TimeStep, ...
        OutputFormat = 'ndgrid', Unit = 'm', BeQuiet = true);
    lon = lon(:);
    lat = lat(:);
    nTime = numel(dates);
    assert(isdatetime(dates) && ~isempty(dates) && all(isfinite(dates)) && ...
        issorted(dates, 'strictascend') && numel(partialMonth) == nTime, ...
        'ULMO:ssh2slept:SourceDates', 'Invalid monthly dates or partial-month flags.');
    assert(isequal(size(values, 1:3), [numel(lon), numel(lat), nTime]), ...
        'ULMO:ssh2slept:SourceGrid', 'Expected a longitude-by-latitude-by-time source grid.');

    try
        basis = prepareslepianbasis(domain, L, options.truncation, options.GVN, ...
            SaveData = options.SaveData);
        [coeffs, fits] = fitgriddedslepian(struct('ssh', values), dates, lon, lat, ...
            domain, L, basis, SaveData = options.SaveData, ForceNew = options.ForceNew);
    catch cause
        errId = 'ULMO:ssh2slept:FitFailed';

        if ismember(string(cause.identifier), ...
                ["ULMO:fitgriddedslepian:EmptyRecord", "ULMO:fitgriddedslepian:EmptyTarget"])
            errId = strrep(cause.identifier, 'fitgriddedslepian', 'ssh2slept');
        elseif strcmp(cause.identifier, 'ULMO:preparegridcoordinates:InvalidGrid')
            errId = 'ULMO:ssh2slept:SourceGrid';
        end

        throwAsCaller(addCause(MException(errId, '%s', cause.message), cause));
    end

    coeffs = coeffs.ssh;
    fit = fits.ssh;
    fit.partialMonth = partialMonth(:);
end

function writeCache(path, tmpPath, product, L, buf, coeffs, dates, fit, defaultJ)
    % Publish only after every variable is complete.
    cleanup = onCleanup(@() deletefileifpresent(tmpPath));
    % NetCDF dimension order here follows MATLAB; time is the last coefficient axis.
    write2nc(tmpPath, 'coefficients', coeffs.', Dimensions = {'mode', 'time'}, Attributes = {'units', 'm'});
    writeVector(tmpPath, 'time', seconds(dates - datetime(2000, 1, 1)), 'time');
    ncwriteatt(tmpPath, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
    ncwriteatt(tmpPath, 'time', 'calendar', 'proleptic_gregorian');

    for name = string(fit.Properties.VariableNames)
        value = fit.(name);

        if islogical(value)
            value = uint8(value);
        end

        writeVector(tmpPath, name, value, 'time');
    end

    ncwriteatt(tmpPath, 'residualRms', 'units', 'm');
    ncwriteatt(tmpPath, '/', 'product', char(product));
    ncwriteatt(tmpPath, '/', 'L', L);
    ncwriteatt(tmpPath, '/', 'observationBufferDegrees', buf);
    ncwriteatt(tmpPath, '/', 'defaultTruncation', defaultJ);
    nccreate(tmpPath, 'creationDate', 'Datatype', 'double');
    ncwrite(tmpPath, 'creationDate', posixtime(datetime('now', 'TimeZone', 'UTC')));
    ncwriteatt(tmpPath, 'creationDate', 'units', 'seconds since 1970-01-01 00:00:00 UTC');
    movefile(tmpPath, path, 'f');
end

function writeVector(path, name, vals, dims)
    nccreate(path, name, 'Dimensions', {dims, numel(vals)}, 'Datatype', class(vals));
    ncwrite(path, name, vals);
end

function [coeffs, dates, fit] = readCache(path, timeRange, oucCnt)
    dates = readshnctime(path);
    rows = find(selectDates(dates, timeRange));
    dates = dates(rows);
    info = ncinfo(path, 'coefficients');
    J = info.Size(1);
    nEpochs = numel(rows);
    coeffs = zeros(nEpochs, J);

    if nEpochs > 0
        coeffs = reshape(ncread(path, 'coefficients', [1, rows(1)], [J, nEpochs]), J, nEpochs).';
    end

    fit = table();

    if oucCnt >= 3

        for name = ["nObs", "areaFraction", "condNum", "residualRms", "partialMonth"]
            vals = zeros(nEpochs, 1);

            if nEpochs > 0
                vals = ncread(path, name, rows(1), nEpochs);
            end

            if name == "partialMonth"
                vals = logical(vals);
            end

            fit.(name) = vals;
        end

    end

end

function rows = selectDates(dates, timeRange)
    rows = true(size(dates));

    if ~isempty(timeRange)
        rows = dates >= timeRange(1) & dates <= timeRange(2);
    end

end
