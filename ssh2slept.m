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
    root = getenv('IFILES');
    assert(~isempty(root), 'ULMO:ssh2slept:MissingRoot', 'Set IFILES to the data directory.');
    folder = fullfile(root, 'SSH', 'SlepianExpansions');
    observationBuffer = 0;
    if L > 0
        observationBuffer = 90/L;
    end
    path = '';
    if isempty(options.GVN)
        if isa(domain, 'GeoDomain')
            domainId = domain.Id;
        else
            domainId = hash(domain, 'sha1');
        end
        stem = fullfile(folder, sprintf('ssh2slept-%s-%s-%d', product, domainId, L));
        if ~options.ForceNew
            path = findCache(stem, L, options.truncation, observationBuffer);
        end
    end

    if ~isempty(path)
        [sshSlept, dates, fit] = readCache(path, timeRange, nargout);
        if nargout >= 4
            [G, V, ~, ~, N] = glmalpha_new(domain, L, 'BeQuiet', true);
            % Match the eigenvalue order and retained columns used by xyzs2slep.
            [V, order] = sort(V(:), 'descend');
            J = size(sshSlept, 2);
            GVN = {G(:, order(1:J)), V(1:J), N};
        end
    else
        source = fullfile(root, 'SSH', product + ".mat");
        [sshSlept, dates, fit, GVN] = fitRecord(source, domain, L, observationBuffer, options);
        if isempty(options.GVN) && options.SaveData
            if ~isfolder(folder)
                mkdir(folder);
            end
            % Publish only after every variable is complete.
            tmp = [tempname(folder), '.nc'];
            cleanup = onCleanup(@() removePartial(tmp));
            defaultJ = min((L + 1)^2, max(1, round(GVN{3})));
            writeCache(tmp, product, L, observationBuffer, sshSlept, dates, fit, defaultJ);
            movefile(tmp, sprintf('%s-%d-m.nc', stem, size(sshSlept, 2)), 'f');
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

function path = findCache(stem, L, J, observationBuffer)
    path = '';
    if isempty(J)
        % Any truncation for this domain and degree identifies the default file.
        files = dir([stem, '-*-m.nc']);
        for file = files(:).'
            candidate = fullfile(file.folder, file.name);
            token = extractBetween(string(candidate), [stem, '-'], '-m.nc');
            if isscalar(token) && all(isstrprop(char(token), 'digit')) && strlength(token) > 0
                info = ncinfo(candidate);
                attr = info.Attributes(strcmp({info.Attributes.Name}, 'defaultTruncation'));
                if isempty(attr)
                    % Older projections stored the Shannon number instead.
                    N = ncreadatt(candidate, '/', 'N');
                    J = min((L + 1)^2, max(1, round(N)));
                else
                    J = attr.Value;
                end
                break
            end
        end
    end
    if ~isempty(J)
        candidate = sprintf('%s-%d-m.nc', stem, J);
        if isfile(candidate)
            % Legacy target-only files used the same name and may lack the width.
            info = ncinfo(candidate);
            attr = info.Attributes(strcmp({info.Attributes.Name}, 'observationBufferDegrees'));
            if (isempty(attr) && L == 0) || (~isempty(attr) && attr.Value == observationBuffer)
                path = candidate;
            end
        end
    end
end

function [coeffs, dates, fit, basis] = fitRecord(source, domain, L, observationBuffer, options)
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
    observationMask = mask;
    if observationBuffer > 0
        observationMask = domainmask(domain, lon, lat, Format='ndgrid', ...
            exteriorBuffer=observationBuffer, SaveData=options.SaveData, ...
            ForceNew=options.ForceNew, BeQuiet=true);
    end
    [gridLon, gridLat] = ndgrid(lon, lat);
    weights = cosd(gridLat);
    totalArea = sum(weights(mask));
    values = reshape(data.sshs, [], nTime);
    points = find(observationMask(:) & any(~isnan(values), 2));
    assert(~isempty(points), 'ULMO:ssh2slept:EmptyRecord', ...
        'No observed SSH cells in the observation region over the source record.');
    values = values(points, :);
    clear data
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

function writeCache(path, product, L, observationBuffer, coeffs, dates, fit, defaultJ)
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
