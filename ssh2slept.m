%% SSH2SLEPT - Fit monthly sea surface height anomalies to regional Slepian functions.
%
% Syntax
%   [sshSlept, dates, fit, GVN] = ssh2slept("NASASSHMonthly", domain, L)
%   [__] = ssh2slept(__, timeRange, unit="m", truncation=J)
%   [__] = ssh2slept(__, TimeStep="midmonth")
%   [__] = ssh2slept(__, GVN={G,V,N})
%   [__] = ssh2slept(__, dropPartialMonths=true)
%   [__] = ssh2slept(__, OutputFormat="traditional")
%
% Input arguments
%   product - "NASASSHMonthly". Reads IFILES/SSH/NASASSHMonthly.mat.
%   domain - Concentration region. A GeoDomain, region name,
%       GeoDomain argument cell, or polygon in the native longitude frame.
%   L - Maximum spherical harmonic degree. The exterior observation buffer is
%       fixed at 90/L degrees for positive L and zero for L=0.
%   timeRange - Optional inclusive pair of datetimes or datenums. Default: all.
%   TimeStep (name-value) - Requested temporal selection, forwarded to SSH2LONLATT.
%       Default selection: [] (native epochs). For NASASSHMonthly, 'midmonth'
%       is equivalent to native and uses the same projection cache. Other
%       selections are not yet supported and raise an error before cache lookup.
%   unit - Output height unit, "mm" (default) or "m".
%   OutputFormat (name-value) - Orientation of returned coefficients.
%       "timefirst": dates by modes; "traditional": modes by dates.
%       Default format: "timefirst".
%   truncation - Optional retained mode count, with the defaults of xyzs2slep.
%   GVN - Optional supplied basis {G,V,N}. A nonempty basis forces a fresh
%       SSH fit and disables reading and writing the projection cache.
%   ForceNew - Rebuild after replacing the source product. Default false.
%   SaveData - Write newly fitted coefficients to NetCDF. Default true.
%   dropPartialMonths - Whether to omit partial months from returned outputs.
%       Applies after fitting or reading the cache; saved data retain all months.
%       Default value: false.
%
% Input file content
%   sshs - Native SSH anomalies in metres, ordered by longitude, latitude and date.
%   lon, lat - Rectangular coordinates in degrees. Longitude is regular and
%       latitude is strictly monotone, possibly nonuniform. Descending axes
%       and a duplicated cyclic longitude endpoint are supported.
%   dates - Native monthly midpoints in ascending order, as datetimes.
%   partialMonths - One logical flag per date, true for a partial month.
%
% Output arguments
%   sshSlept - Coefficients in the requested height unit and OutputFormat.
%       Size: (nDates, nModes) for "timefirst", (nModes, nDates) for "traditional".
%   dates - Native monthly midpoints, as a datetime column.
%   fit - One row per date: nObs (all fitted cells), areaFraction (observed fraction of the
%       target domain grid area), condNum, residualRms, and partialMonth.
%       Residual RMS is in the requested height unit.
%   GVN - Retained harmonic basis, eigenvalues, and full Shannon number.
%
% Notes
%   A cache miss fits the complete source record. Uniform latitude uses
%   cosine-latitude weights; nonuniform latitude uses spherical midpoint
%   bands. Cyclic endpoints count once. Empty months retain NaN coefficients.
%   A wholly empty record or an underdetermined nonempty month raises an error.
%   IFILES/SSH/SlepianExpansions holds one NetCDF-4 file per projection.
%   Names follow ssh2slept-<product>-<domainId>-<L>-<J>-m.nc. For XY
%   polygons the domain ID is hash(domain, 'sha1'). J is the fitted mode count.
%   Arctic and Antarctic IDs include -rotb for geographic observation masks.
%   The cache stores metres, independent of output unit, OutputFormat and requested dates.
%   Native and midmonth requests share this cache for NASASSHMonthly; neither
%   request temporally resamples that already monthly source.
%   Cache hits read the time coordinate and selected contiguous time slabs.
%   Requesting GVN on a cache hit loads the domain basis through glmalpha_new.
%
% See also
%   ssh2lonlatt, prepareslepianbasis, fitgriddedslepian, readshnctime
%
% Created by
%   2026/11/03
%
% Last modified
%   2026/11/07, En-Chi Lee (williameclee@gmail.com)

function [sshSlept, dates, fit, GVN] = ssh2slept(product, domain, L, timeRange, options)
    arguments (Input)
        product (1, 1) string {mustBeMember(product, "NASASSHMonthly")}
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        timeRange = []
        options.TimeStep {mustBeTimeStep} = []
        options.unit (1, 1) string {mustBeMember(options.unit, ["m", "mm"])} = "mm"
        options.OutputFormat (1, 1) string {mustBeMember(options.OutputFormat, ["timefirst", "traditional"])} = "timefirst"
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.GVN cell = {}
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.dropPartialMonths (1, 1) logical = false
    end
    % Only equivalent native/monthly requests can share the existing cache.
    % Validate before lookup so other selections cannot return monthly data.
    timestep = options.TimeStep;
    assert(isempty(timestep) || ((ischar(timestep) || isstring(timestep)) && ...
        isscalar(string(timestep)) && strcmpi(timestep, 'midmonth')), ...
        'ULMO:ssh2slept:TimeStep', ...
        'NASASSHMonthly projections currently support only native ([]) or midmonth epochs.');
    if ~isempty(timeRange)
        if isnumeric(timeRange)
            timeRange = datetime(timeRange, 'ConvertFrom', 'datenum');
        end
        assert(isdatetime(timeRange) && numel(timeRange) == 2 && ...
            all(isfinite(timeRange)) && timeRange(1) <= timeRange(2), ...
            'ULMO:ssh2slept:TimeRange', 'Supply increasing start and end dates.');
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
    observationBuffer = 0;
    if L > 0
        observationBuffer = 90/L;
    end
    cacheHit = false;
    if ~options.ForceNew
        [path, tmp, cacheHit] = getoutputfile(product, domain, L, options.truncation, options.SaveData, false);
    end

    if cacheHit
        % Filtering needs the flags even when diagnostics are not requested.
        outputCount = max(nargout, 3 * options.dropPartialMonths);
        [sshSlept, dates, fit] = readCache(path, timeRange, outputCount);
        if nargout >= 4
            J = size(sshSlept, 2);
            GVN = prepareslepianbasis(domain, L, J, {}, SaveData=options.SaveData);
        end
    else
        [sshSlept, dates, fit, GVN] = fitRecord(product, domain, L, options);
        if options.SaveData
            % A forced fit skips lookup. A first default fit determines J here.
            if options.ForceNew || isempty(path)
                [path, tmp] = getoutputfile(product, domain, L, size(sshSlept, 2), true, true);
            end
            defaultJ = min((L + 1)^2, max(1, round(GVN{3})));
            writeCache(path, tmp, product, L, observationBuffer, sshSlept, dates, fit, defaultJ);
        end
        rows = selectDates(dates, timeRange);
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
    if options.unit == "mm"
        sshSlept = sshSlept * 1000;
        if nargout >= 3
            fit.residualRms = fit.residualRms * 1000;
        end
    end
    if options.OutputFormat == "traditional"
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
        OutputFormat='ndgrid', Unit='m', BeQuiet=true);
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
            SaveData=options.SaveData);
        [coefficients, fits] = fitgriddedslepian(struct('ssh', values), dates, lon, lat, ...
            domain, L, basis, SaveData=options.SaveData, ForceNew=options.ForceNew);
    catch cause
        identifier = 'ULMO:ssh2slept:FitFailed';
        if ismember(string(cause.identifier), ["ULMO:fitgriddedslepian:EmptyRecord", ...
                "ULMO:fitgriddedslepian:EmptyTarget"])
            identifier = strrep(cause.identifier, 'fitgriddedslepian', 'ssh2slept');
        elseif strcmp(cause.identifier, 'ULMO:preparegridcoordinates:InvalidGrid')
            identifier = 'ULMO:ssh2slept:SourceGrid';
        end
        throwAsCaller(addCause(MException(identifier, '%s', cause.message), cause));
    end
    coeffs = coefficients.ssh;
    fit = fits.ssh;
    fit.partialMonth = partialMonth(:);
end

function writeCache(destination, path, product, L, observationBuffer, coeffs, dates, fit, defaultJ)
    % Publish only after every variable is complete.
    cleanup = onCleanup(@() removePartial(path));
    % NetCDF dimension order here follows MATLAB; time is the last coefficient axis.
    write2nc(path, 'coefficients', coeffs.', Dimensions={'mode','time'}, Attributes={'units','m'});
    writeVector(path, 'time', seconds(dates - datetime(2000,1,1)), 'time');
    ncwriteatt(path, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
    ncwriteatt(path, 'time', 'calendar', 'proleptic_gregorian');
    for name = string(fit.Properties.VariableNames)
        value = fit.(name);
        if islogical(value)
            value = uint8(value);
        end
        writeVector(path, name, value, 'time');
    end
    ncwriteatt(path, 'residualRms', 'units', 'm');
    ncwriteatt(path, '/', 'product', char(product));
    ncwriteatt(path, '/', 'L', L);
    ncwriteatt(path, '/', 'observationBufferDegrees', observationBuffer);
    ncwriteatt(path, '/', 'defaultTruncation', defaultJ);
    nccreate(path, 'creationDate', 'Datatype', 'double');
    ncwrite(path, 'creationDate', posixtime(datetime('now', 'TimeZone', 'UTC')));
    ncwriteatt(path, 'creationDate', 'units', 'seconds since 1970-01-01 00:00:00 UTC');
    movefile(path, destination, 'f');
end

function writeVector(path, name, values, dimension)
    nccreate(path, name, 'Dimensions', {dimension, numel(values)}, 'Datatype', class(values));
    ncwrite(path, name, values);
end

function [coeffs, dates, fit] = readCache(path, timeRange, outputCount)
    dates = readshnctime(path);
    rows = find(selectDates(dates, timeRange));
    dates = dates(rows);
    info = ncinfo(path, 'coefficients');
    J = info.Size(1);
    count = numel(rows);
    coeffs = zeros(count, J);
    if count > 0
        coeffs = reshape(ncread(path, 'coefficients', [1 rows(1)], [J count]), J, count).';
    end
    fit = table();
    if outputCount >= 3
        for name = ["nObs", "areaFraction", "condNum", "residualRms", "partialMonth"]
            values = zeros(count, 1);
            if count > 0
                values = ncread(path, name, rows(1), count);
            end
            if name == "partialMonth"
                values = logical(values);
            end
            fit.(name) = values;
        end
    end
end

function rows = selectDates(dates, timeRange)
    rows = true(size(dates));
    if ~isempty(timeRange)
        rows = dates >= timeRange(1) & dates <= timeRange(2);
    end
end

function removePartial(path)
    if isfile(path)
        delete(path);
    end
end
