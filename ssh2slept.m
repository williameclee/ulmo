%% SSH2SLEPT - Fit monthly sea surface height anomalies to regional Slepian functions.
%
% Syntax
%   [sshSlept, dates, fit, GVN] = ssh2slept("NASASSHMonthly", domain, L)
%   [__] = ssh2slept(__, timeRange, unit="m", truncation=J, GVN=basis)
%
% Input arguments
%   product - "NASASSHMonthly". Reads IFILES/SSH/NASASSHMonthly.mat.
%   domain - Concentration and observation region. A GeoDomain, region name,
%       GeoDomain argument cell, or polygon in the native longitude frame.
%   L - Maximum spherical harmonic degree.
%   timeRange - Optional inclusive pair of datetimes or datenums. Default: all.
%   unit - Output height unit, "mm" (default) or "m".
%   truncation, GVN - Optional retained mode count and basis, as in xyzs2slep.
%   ForceNew - Rebuild after replacing the source product. Default false.
%   SaveData - Write newly fitted coefficients to NetCDF. Default true.
%
% Output arguments
%   sshSlept - Date-by-mode coefficients in the requested height unit.
%   dates - Native monthly midpoints, as a datetime column.
%   fit - One row per date: nObs, areaFraction (observed fraction of the
%       target domain grid area), condNum, residualRms, and partialMonth.
%       Residual RMS is in the requested height unit.
%   GVN - Retained harmonic basis, eigenvalues, and full Shannon number.
%
% Notes
%   A cache miss fits the complete source record using cosine-latitude
%   weights on the native regular grid. Empty months retain NaN coefficients.
%   A wholly empty record or an underdetermined nonempty month raises an error.
%   IFILES/SSH/SlepianExpansions holds one NetCDF-4 file per projection.
%   The cache stores metres, independent of output unit and requested dates.
%   Cache hits read the time coordinate and selected contiguous time slabs.
%
% See also
%   XYZS2SLEP, SSH2LONLATT, DOMAINMASK, READSHNCTIME
%
% Created by
%   2026/11/03

function [sshSlept, dates, fit, GVN] = ssh2slept(product, domain, L, timeRange, options)
    arguments (Input)
        product (1, 1) string {mustBeMember(product, "NASASSHMonthly")}
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        timeRange = []
        options.unit (1, 1) string {mustBeMember(options.unit, ["m", "mm"])} = "mm"
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.GVN cell = {}
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
    end
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
    root = getenv('IFILES');
    assert(~isempty(root), 'ULMO:ssh2slept:MissingRoot', 'Set IFILES to the data directory.');
    folder = fullfile(root, 'SSH', 'SlepianExpansions');
    path = fullfile(folder, cacheName(product, domain, L, options.truncation, options.GVN));

    if isfile(path) && ~options.ForceNew
        [sshSlept, dates, fit, GVN] = readCache(path, timeRange, nargout);
    else
        source = fullfile(root, 'SSH', product + ".mat");
        [sshSlept, dates, fit, GVN] = fitRecord(source, domain, L, options);
        if options.SaveData
            if ~isfolder(folder)
                mkdir(folder);
            end
            % Publish only after every variable is complete.
            tmp = [tempname(folder), '.nc'];
            cleanup = onCleanup(@() removePartial(tmp));
            writeCache(tmp, product, L, sshSlept, dates, fit, GVN);
            movefile(tmp, path, 'f');
        end
        rows = selectDates(dates, timeRange);
        sshSlept = sshSlept(rows, :);
        fit = fit(rows, :);
        dates = dates(rows);
    end
    if options.unit == "mm"
        sshSlept = sshSlept * 1000;
        if nargout >= 3
            fit.residualRms = fit.residualRms * 1000;
        end
    end
end

function name = cacheName(product, domain, L, truncation, basis)
    % Include the actual supplied basis: signs and column order affect C.
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(uint8(jsonencode(domain)));
    for item = [{truncation}, reshape(basis, 1, [])]
        value = item{1};
        validateattributes(value, {'numeric'}, {'real'});
        digest.update(typecast(uint64(size(value)), 'int8'));
        if ~isempty(value)
            bytes = typecast(double(value(:)), 'int8');
            digest.update(bytes, 0, numel(bytes));
        end
    end
    key = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
    region = 'polygon';
    if isa(domain, 'GeoDomain')
        region = domain.Id;
    end
    name = sprintf('SSH2SLEPT-%s-%s-L%d-%s.nc', product, region, L, key);
end

function [coeffs, dates, fit, basis] = fitRecord(source, domain, L, options)
    data = load(source, 'sshs', 'dates', 'lon', 'lat', 'partialMonths');
    dates = data.dates(:);
    partialMonth = logical(data.partialMonths(:));
    lon = double(data.lon(:));
    lat = double(data.lat(:));
    nTime = numel(dates);
    assert(isdatetime(dates) && ~isempty(dates) && all(isfinite(dates)) && ...
        issorted(dates, 'strictascend') && numel(partialMonth) == nTime, ...
        'ULMO:ssh2slept:SourceDates', 'Invalid monthly dates or partial-month flags.');
    assert(isequal(size(data.sshs, 1:3), [numel(lon), numel(lat), nTime]) && ...
        all(diff(lon) == 0.5) && all(diff(lat) == 0.5), ...
        'ULMO:ssh2slept:SourceGrid', 'Expected an ascending native 0.5-degree lon-by-lat grid.');
    mask = domainmask(domain, lon, lat, Format='ndgrid', SaveData=false, BeQuiet=true);
    [gridLon, gridLat] = ndgrid(lon, lat);
    weights = cosd(gridLat);
    totalArea = sum(weights(mask));
    values = reshape(data.sshs, [], nTime);
    points = find(mask(:) & any(~isnan(values), 2));
    assert(~isempty(points), 'ULMO:ssh2slept:EmptyRecord', ...
        'No observed SSH cells in the target region over the source record.');
    values = values(points, :);
    clear data
    gridLon = gridLon(points);
    gridLat = gridLat(points);
    weights = weights(points);
    nObs = sum(~isnan(values), 1).';
    areaFraction = zeros(nTime, 1);
    condNum = nan(nTime, 1);
    residualRms = nan(nTime, 1);
    design = [];
    for iDate = find(nObs > 0).'
        field = values(:, iDate);
        try
            if isempty(design)
                [c, ~, ~, stats, design] = xyzs2slep(field, gridLon, gridLat, domain, L, ...
                    weight=weights, missingPolicy="omit", GVN=options.GVN, truncation=options.truncation);
                basis = design.GVN;
                coeffs = nan(nTime, numel(c));
            else
                [c, ~, ~, stats] = xyzs2slep(field, design, weight=weights, missingPolicy="omit");
            end
        catch cause
            err = MException('ULMO:ssh2slept:FitFailed', 'NASASSHMonthly fit failed for %s: %s', ...
                string(dates(iDate)), cause.message);
            throwAsCaller(addCause(err, cause));
        end
        coeffs(iDate, :) = c.';
        condNum(iDate) = stats.condNum;
        residualRms(iDate) = stats.residualRms;
        areaFraction(iDate) = sum(weights(~isnan(field))) / totalArea;
    end
    fit = table(nObs, areaFraction, condNum, residualRms, partialMonth);
end

function writeCache(path, product, L, coeffs, dates, fit, basis)
    % NetCDF dimension order here follows MATLAB; time is the last coefficient axis.
    write2nc(path, 'coefficients', coeffs.', Dimensions={'mode','time'}, Attributes={'units','m'});
    write2nc(path, 'G', basis{1}, Dimensions={'harmonic','mode'});
    writeVector(path, 'V', basis{2}, 'mode');
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
    ncwriteatt(path, '/', 'N', basis{3});
end

function writeVector(path, name, values, dimension)
    nccreate(path, name, 'Dimensions', {dimension, numel(values)}, 'Datatype', class(values));
    ncwrite(path, name, values);
end

function [coeffs, dates, fit, basis] = readCache(path, timeRange, outputCount)
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
    basis = {};
    if outputCount >= 4
        basis = {ncread(path, 'G'), ncread(path, 'V'), ncreadatt(path, '/', 'N')};
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
