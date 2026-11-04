%% SSH2SLEPT - Fit monthly sea surface height anomalies to regional Slepian functions.
%
% Syntax
%   [sshSlept, dates, fit, GVN] = ssh2slept("NASASSHMonthly", domain, L)
%   [__] = ssh2slept(__, timeRange, unit="m", truncation=J)
%   [__] = ssh2slept(__, GVN={G,V,N})
%
% Input arguments
%   product - "NASASSHMonthly". Reads IFILES/SSH/NASASSHMonthly.mat.
%   domain - Concentration region. A GeoDomain, region name,
%       GeoDomain argument cell, or polygon in the native longitude frame.
%   L - Maximum spherical harmonic degree. The exterior observation buffer is
%       fixed at 90/L degrees for positive L and zero for L=0.
%   timeRange - Optional inclusive pair of datetimes or datenums. Default: all.
%   unit - Output height unit, "mm" (default) or "m".
%   truncation - Optional retained mode count, with the defaults of xyzs2slep.
%   GVN - Optional supplied basis {G,V,N}. A nonempty basis forces a fresh
%       SSH fit and disables reading and writing the projection cache.
%   ForceNew - Rebuild after replacing the source product. Default false.
%   SaveData - Write newly fitted coefficients to NetCDF. Default true.
%
% Input file content
%   sshs - Native SSH anomalies in metres, ordered by longitude, latitude and date.
%   lon, lat - Ascending native grid coordinates in degrees at 0.5-degree spacing.
%   dates - Native monthly midpoints in ascending order, as datetimes.
%   partialMonths - One logical flag per date, true for a partial month.
%
% Output arguments
%   sshSlept - Date-by-mode coefficients in the requested height unit.
%   dates - Native monthly midpoints, as a datetime column.
%   fit - One row per date: nObs (all fitted cells), areaFraction (observed fraction of the
%       target domain grid area), condNum, residualRms, and partialMonth.
%       Residual RMS is in the requested height unit.
%   GVN - Retained harmonic basis, eigenvalues, and full Shannon number.
%
% Notes
%   A cache miss fits the complete source record using cosine-latitude
%   weights on the native regular grid. Empty months retain NaN coefficients.
%   A wholly empty record or an underdetermined nonempty month raises an error.
%   IFILES/SSH/SlepianExpansions holds one NetCDF-4 file per projection.
%   Names follow ssh2slept-<product>-<domainId>-<L>-<J>-m.nc. For XY
%   polygons the domain ID is hash(domain, 'sha1'). J is the fitted mode count.
%   Arctic and Antarctic IDs include -rotb for geographic observation masks.
%   The cache stores metres, independent of output unit and requested dates.
%   Cache hits read the time coordinate and selected contiguous time slabs.
%   Requesting GVN on a cache hit loads the domain basis through glmalpha_new.
%
% See also
%   XYZS2SLEP, SSH2LONLATT, DOMAINMASK, READSHNCTIME, GLMALPHA_NEW
%
% Created by
%   2026/11/03
%
% Last modified
%   2026/11/04

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
        [sshSlept, dates, fit] = readCache(path, timeRange, nargout);
        if nargout >= 4
            [G, V, ~, ~, N] = glmalpha_new(domain, L, 'BeQuiet', true);
            % Match the eigenvalue order and retained columns used by xyzs2slep.
            [V, order] = sort(V(:), 'descend');
            J = size(sshSlept, 2);
            GVN = {G(:, order(1:J)), V(1:J), N};
        end
    else
        [sshSlept, dates, fit, GVN] = fitRecord(product, domain, L, observationBuffer, options);
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
    if options.unit == "mm"
        sshSlept = sshSlept * 1000;
        if nargout >= 3
            fit.residualRms = fit.residualRms * 1000;
        end
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

function [coeffs, dates, fit, basis] = fitRecord(product, domain, L, observationBuffer, options)
    [values, ~, dates, lon, lat, partialMonth] = ssh2lonlatt(product, ...
        OutputFormat='ndgrid', Unit='m', BeQuiet=true);
    lon = double(lon(:));
    lat = double(lat(:));
    nTime = numel(dates);
    assert(isdatetime(dates) && ~isempty(dates) && all(isfinite(dates)) && ...
        issorted(dates, 'strictascend') && numel(partialMonth) == nTime, ...
        'ULMO:ssh2slept:SourceDates', 'Invalid monthly dates or partial-month flags.');
    assert(isequal(size(values, 1:3), [numel(lon), numel(lat), nTime]) && ...
        all(diff(lon) == 0.5) && all(diff(lat) == 0.5), ...
        'ULMO:ssh2slept:SourceGrid', 'Expected an ascending native 0.5-degree lon-by-lat grid.');
    % areaFraction measures observed coverage within the original target domain.
    mask = domainmask(domain, lon, lat, Format='ndgrid', SaveData=false, BeQuiet=true);
    observationMask = mask;
    % A positive buffer includes exterior cells in the fit and in nObs.
    if observationBuffer > 0
        observationMask = domainmask(domain, lon, lat, Format='ndgrid', ...
            exteriorBuffer=observationBuffer, SaveData=options.SaveData, ...
            ForceNew=options.ForceNew, BeQuiet=true);
    end
    [gridLon, gridLat] = ndgrid(lon, lat);
    weights = cosd(gridLat);
    totalArea = sum(weights(mask));
    values = reshape(values, [], nTime);
    points = find(observationMask(:) & any(~isnan(values), 2));
    assert(~isempty(points), 'ULMO:ssh2slept:EmptyRecord', ...
        'No observed SSH cells in the observation region over the source record.');
    values = values(points, :);
    gridLon = gridLon(points);
    gridLat = gridLat(points);
    weights = weights(points);
    observed = ~isnan(values);
    nObs = sum(observed, 1).';
    target = mask(points);
    assert(totalArea > 0, 'ULMO:ssh2slept:EmptyTarget', 'The target contains no grid cells.');
    areaFraction = min(1, (weights(target).' * observed(target, :)).' / totalArea);
    months = find(nObs > 0);
    try
        [c, ~, ~, stats, design] = xyzs2slep(values(:, months), gridLon, gridLat, domain, L, ...
            weight=weights, missingPolicy="omit", truncation=options.truncation, GVN=options.GVN);
    catch cause
        % Field numbers in solver errors index the nonempty months.
        field = regexpi(cause.message, '\bfield (\d+)\b', 'tokens', 'once');
        context = "NASASSHMonthly fit failed";
        if ~isempty(field)
            context = context + " for " + string(dates(months(str2double(field{1}))));
        end
        err = MException('ULMO:ssh2slept:FitFailed', '%s: %s', context, cause.message);
        throwAsCaller(addCause(err, cause));
    end
    basis = design.GVN;
    coeffs = nan(nTime, size(c, 1));
    condNum = nan(nTime, 1);
    residualRms = nan(nTime, 1);
    coeffs(months, :) = c.';
    condNum(months) = stats.condNum;
    residualRms(months) = stats.residualRms;
    fit = table(nObs, areaFraction, condNum, residualRms, partialMonth);
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
