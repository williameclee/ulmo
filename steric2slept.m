%% STERIC2SLEPT - Fit temporally processed steric fields to regional Slepian functions.
% Returns raw coefficients and separate global-mean correction ingredients.
%
% Syntax
%   [c, dates, fit, GVN, h, constantSlept] = STERIC2SLEPT(product, domain, L)
%   [__] = STERIC2SLEPT(product, domain, L, timestep, timeRange, Type="total")
%   [__] = STERIC2SLEPT(__, OutputFormat="traditional", SaveData=false)
%
% Input arguments
%   product - Prepared StericProduct identifier.
%   domain - GeoDomain, region name, GeoDomain argument cell or polygon.
%   L - Maximum spherical harmonic degree, a nonnegative integer.
%   timestep (optional) - Temporal selection as in STERIC2LONLATT: a positive
%       duration or number of days, 'midmonth', 'GRACE', or []. Processing
%       precedes projection and uses native spatial coordinates.
%       Default step: [] (native timestamps).
%   timeRange (optional) - Inclusive datetime or datenum output limits.
%       Increasing limits, applied after processing the complete record.
%       Size: (1, 2) or (2, 1), or empty.
%       Default range: [] (complete record).
%   Type (name-value) - total, thermosteric, halosteric, non-thermosteric
%       (total minus thermosteric), or non-halosteric (total minus halosteric).
%       Default type: total.
%   Depth (name-value) - full, shallow or deep, as supplied by the product.
%       Default depth: full.
%   Interpolation (name-value) - Temporal method as in STERIC2LONLATT.
%       'mean' averages fixed-step windows; midmonth and GRACE always average.
%       Interpolation methods are validated by INTERPTEMPORAL when used.
%       Default method: linear.
%   OutputFormat (name-value) - timefirst or traditional coefficient ordering.
%       GMSL uses the same time dimension as the coefficients.
%       Default format: timefirst.
%   Unit (name-value) - Output heights in m or mm.
%       Default unit: m.
%   TimeFormat (name-value) - datetime or datenum output dates.
%       Default format: datetime.
%   truncation (name-value) - Retained mode count, as in XYZS2SLEP.
%       Default count: [] (bounded round of the Shannon number).
%   GVN (name-value) - Optional orthonormal basis {G,V,N} with (L+1)^2
%       harmonic rows. Supplied columns retain their order and signs.
%       Bypasses projection-cache reads/writes and uses all supplied columns
%       by default for an already truncated basis.
%       Default basis: {} (calculate the domain basis).
%   ForceNew (name-value) - Whether to rebuild all supported projections and
%       native/processed means for this configuration, regardless of outputs.
%       Does not regenerate source grids or other cache configurations.
%       Default flag: false.
%   SaveData (name-value) - Whether to save calculated means and caches.
%       false permits reads and in-memory rebuilding without writes.
%       Default flag: true.
%   BeQuiet (name-value) - Whether to suppress progress messages.
%       Default level: 0.5 (calculation/cache messages); 1 suppresses messages.
%
% Input file content
%   lon, lat, dates - Native coordinates and strictly increasing timestamps
%       in the MAT file resolved by STERICFILEPATH. Longitude spans one cycle;
%       latitude may be nonuniform. Square grids use (lat,lon,time) ordering.
%   <component>Sl - Supplied native grids, in metres, with time on axis three.
%       Fields and complete depth splits follow STERICFIELDSCHEMA.
%
% Output arguments
%   stericSlept - Raw selected component/depth coefficients.
%       Size: (nTime, J) for timefirst; (J, nTime) for traditional.
%       Unit: Unit.
%   dates - Returned timestamps in TimeFormat.
%       Size: (nTime, 1).
%   fit - Table of nObs, areaFraction, condNum and residualRms, one row per
%       returned epoch. Counts include exterior-buffer cells and zero-filled
%       deep cells; areaFraction measures fitted area inside the target.
%       Residual RMS uses Unit.
%   GVN - Retained harmonic basis, eigenvalues and full Shannon number.
%   halostericGmsl - Selected depth's independently processed halosteric GMSL.
%       Size: (nTime, 1) for timefirst; (1, nTime) for traditional.
%       Unit: Unit.
%   constantSlept - Harmonic projection of unit C00: GVN{1}(1,:).'.
%       Independent of data masks and output units/order; dimensionless.
%       Size: (J, 1).
%
% Output file content
%   <component>Slept - All supported raw coefficient fields, in metres.
%       MATLAB NetCDF dimensions: (mode, time).
%   <component>Gmsl - All applicable processed global means in metres, (time).
%   <component>Slept_<diagnostic> - Per-field fitting diagnostics, (time).
%   constantSlept - Unit-C00 projection, (mode), with units 1.
%   time, mode, creationDate - Epochs, retained mode indices and UTC creation time.
%       Attributes record field availability and the projection configuration.
%
% Notes
%   One complete NetCDF-4 file per product/temporal/domain/L/J configuration
%   lives in IFILES/STERIC/product/SlepianExpansions. Type, Depth, timeRange,
%   Unit and OutputFormat select outputs, not caches. Updates replace entire
%   files only after successful completion. Convention changes require manual
%   rebuilding. Basis caches and projection caches must remain consistent.
%   The bulk path shares STERIC2LONLATT's temporal operations and obtains all
%   native means through STERIC2GMSL. It processes grids in memory without
%   intermediate grid caches. Unsaved requests may recalculate native means
%   for each selected component/depth, following STERIC2LONLATT.
%   Full and shallow fields retain their finite observations. Supplied deep
%   grids use paired finite full/shallow support: missing deep values there
%   are zero contributions, and deep values outside that support are omitted.
%   This absent-layer convention is applied before temporal processing and
%   residual construction, in memory only. It assumes missing deep values on
%   observed upper columns mean an absent deep contribution, not lost data.
%   Native GMSL estimation is unchanged. Source layer inconsistencies and
%   nonlinear temporal interpolation can still prevent exact depth closure.
%   Raw grids and GMSL use independent finite samples during epoch averaging.
%   For timefirst, caller correction is c - h*constantSlept.'; traditional
%   uses c - constantSlept*h. This is harmonic C00 subtraction, not refitting
%   corrected grids with their observation masks. No product policy is applied.
%   Empty epochs retain NaN coefficients. Entirely empty fields and
%   underdetermined nonempty fits raise errors; previous caches are preserved.
%
% See also
%   steric2lonlatt, steric2gmsl, ssh2slept, grace2slept_new, xyzs2slep
%
% Created by
%   2026/11/05, En-Chi Lee

function [stericSlept, dates, fit, GVN, halostericGmsl, constantSlept] = ...
        steric2slept(product, domain, L, timestep, timeRange, options)
    arguments (Input)
        product (1, 1) StericProduct
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        timestep {mustBeTimeStep} = []
        timeRange {mustBeTimeRange} = []
        options.Type {mustBeTextScalar, mustBeMember(options.Type, ...
            {'total', 'thermosteric', 'halosteric', 'non-thermosteric', 'non-halosteric'})} = 'total'
        options.Depth {mustBeTextScalar, mustBeMember(options.Depth, {'full', 'shallow', 'deep'})} = 'full'
        options.Interpolation {mustBeTextScalar} = 'linear'
        options.OutputFormat {mustBeTextScalar, mustBeMember(options.OutputFormat, {'timefirst', 'traditional'})} = 'timefirst'
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, {'m', 'mm'})} = 'm'
        options.TimeFormat {mustBeTextScalar, mustBeMember(options.TimeFormat, {'datetime', 'datenum'})} = 'datetime'
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.GVN cell = {}
        options.ForceNew (1, 1) {mustBeNumericOrLogical} = false
        options.SaveData (1, 1) {mustBeNumericOrLogical} = true
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
    end
    options.Interpolation = lower(options.Interpolation);
    if ~isempty(timeRange) && isnumeric(timeRange)
        timeRange = datetime(timeRange, 'ConvertFrom', 'datenum');
    end
    if iscell(domain)
        domain = GeoDomain(domain{:});
    elseif ischar(domain) || isstring(domain)
        domain = GeoDomain(domain);
    end
    assert((isa(domain, 'GeoDomain') && isscalar(domain)) || ...
        (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2 && all(isfinite(domain), 'all')), ...
        'ULMO:steric2slept:Domain', 'Supply a GeoDomain, region name or finite polygon.');
    if isnumeric(timestep) && ~isempty(timestep), timestep = days(timestep); end
    if ischar(timestep) || isstring(timestep), timestep = lower(string(timestep)); end
    schema = stericfieldschema();
    iType = find(schema.types == options.Type);
    iDepth = find(schema.depths == options.Depth);
    field = replace(schema.meanNames(iType, iDepth), 'Gmsl', 'Slept');
    haloName = schema.meanNames(3, iDepth);
    supplied = ~isempty(options.GVN);
    [stem, configuration, domainId] = cacheStem(product, domain, L, timestep, options.Interpolation);
    path = '';
    if ~options.ForceNew && ~supplied
        path = findCache(stem, options.truncation);
    end
    if ~isempty(path) && isfile(path)
        if options.BeQuiet < 1, fprintf('[ULMO>steric2slept] Loading %s.\n', path); end
        [stericSlept, dates, fit, halostericGmsl, constantSlept] = ...
            readCache(path, field, haloName, timeRange, nargout);
        if nargout >= 4
            GVN = getBasis(domain, L, size(stericSlept, 2), {}, options.SaveData);
        end
    else
        if options.BeQuiet < 1, fprintf('[ULMO>steric2slept] Fitting all fields for %s.\n', char(product)); end
        metadata = whos('-file', stericfilepath(product));
        schema = stericfieldschema(metadata);
        assert(schema.supported(iType, iDepth), 'ULMO:steric2slept:UnsupportedField', ...
            'Product %s does not supply Type=%s, Depth=%s.', char(product), options.Type, options.Depth);
        if nargout >= 5
            assert(schema.supported(3, iDepth), 'ULMO:steric2slept:UnsupportedField', ...
                'Product %s does not supply halosteric GMSL at Depth=%s.', char(product), options.Depth);
        end
        % A supplied basis bypasses projection storage, not an explicit request
        % to save native diagnostics. SaveData=false is propagated everywhere.
        GVN = getBasis(domain, L, options.truncation, options.GVN, options.SaveData);
        constantSlept = GVN{1}(1, :).';
        [fields, means, dates, lon, lat] = prepareRecord(product, timestep, options, schema);
        [coefficients, fits, buffer, weighting] = fitRecord(fields, dates, lon, lat, domain, L, GVN, options);
        if options.SaveData && ~supplied
            defaultJ = min((L + 1)^2, max(1, round(GVN{3})));
            path = sprintf('%s-%d-m.nc', stem, size(GVN{1}, 2));
            writeCache(path, coefficients, means, dates, fits, constantSlept, ...
                product, domainId, L, defaultJ, buffer, weighting, configuration);
        end
        rows = selectDates(dates, timeRange);
        stericSlept = coefficients.(field)(rows, :);
        fit = fits.(field)(rows, :);
        dates = dates(rows);
        halostericGmsl = zeros(0, 1);
        if nargout >= 5, halostericGmsl = means.(haloName)(rows); end
    end
    if strcmp(options.Unit, 'mm')
        stericSlept = 1000 * stericSlept;
        halostericGmsl = 1000 * halostericGmsl;
        if nargout >= 3, fit.residualRms = 1000 * fit.residualRms; end
    end
    if strcmp(options.OutputFormat, 'traditional')
        stericSlept = stericSlept.';
        halostericGmsl = halostericGmsl.';
    end
    if strcmp(options.TimeFormat, 'datenum'), dates = datenum(dates); end %#ok<DATNM>
end

function [stem, configuration, domainId] = cacheStem(product, domain, L, timestep, method)
    root = getenv('IFILES');
    assert(isfolder(root), 'ULMO:steric2slept:MissingRoot', 'Set IFILES to the data directory.');
    if isa(domain, 'GeoDomain')
        domainId = domain.Id;
        if ismember(domain.Domain, ["arctic", "antarctica"]), domainId = [domainId, '-rotb']; end
    else
        domainId = hash(domain, 'sha1');
    end
    if isempty(timestep)
        configuration = 'native';
    elseif isduration(timestep)
        % 17 digits distinguish representable fixed steps, including fractions.
        configuration = sprintf('T%sseconds-%s', sprintf('%.17g', seconds(timestep)), method);
    else
        configuration = char(timestep); % Named epochs always average.
    end
    stem = fullfile(root, 'STERIC', char(product), 'SlepianExpansions', ...
        sprintf('steric2slept-%s-%s-%s-%d', char(product), configuration, domainId, L));
end

function path = findCache(stem, J)
    path = '';
    if isempty(J)
        files = dir([stem, '-*-m.nc']);
        % The extensionless stem can contain decimal points in the domain ID.
        [~, name, extension] = fileparts(stem);
        name = [name, extension];
        names = string({files.name});
        first = find(matches(names, name + "-" + digitsPattern + "-m.nc"), 1);
        if isempty(first), return; end
        J = ncreadatt(fullfile(files(first).folder, files(first).name), '/', 'defaultTruncation');
    end
    path = sprintf('%s-%d-m.nc', stem, J);
end

function basis = getBasis(domain, L, J, supplied, saveData)
    if isempty(supplied)
        [G, V, ~, ~, N] = glmalpha_new(domain, L, SaveData=saveData, BeQuiet=true);
        [V, order] = sort(V(:), 'descend');
        G = G(:, order);
    else
        assert(numel(supplied) == 3, 'ULMO:steric2slept:Basis', 'Supply {G,V,N}.');
        [G, V, N] = supplied{:};
    end
    validateattributes(G, {'double', 'single'}, {'real', 'finite', '2d', 'nonempty', 'nrows', (L+1)^2});
    validateattributes(V, {'numeric'}, {'real', 'finite', 'vector', 'numel', size(G, 2)});
    validateattributes(N, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});
    assert(size(G, 2) <= size(G, 1), 'ULMO:steric2slept:Basis', 'Too many basis columns.');
    if isempty(J)
        J = size(G, 2);
        if J == (L+1)^2, J = min(J, max(1, round(N))); end
    end
    assert(J <= size(G, 2), 'ULMO:steric2slept:Truncation', 'Truncation exceeds basis size.');
    precision = double(eps(class(G)));
    G = double(G(:, 1:J));
    assert(norm(G.'*G - eye(J), 'fro') <= 100 * precision * max(size(G)), ...
        'ULMO:steric2slept:Basis', 'The retained harmonic basis must be orthonormal.');
    V = V(:);
    basis = {G, V(1:J), N};
end

function [fields, means, dates, lon, lat] = prepareRecord(product, timestep, options, schema)
    % Obtain every diagnostic through the public native-GMSL interface.
    % A saved forced call refreshes all means; later calls reuse that set.
    % Unsaved forced calls must each bypass the existing native cache.
    means = struct();
    firstMean = true;
    for iMean = find(schema.supported(:)).'
        [iType, iDepth] = ind2sub(size(schema.supported), iMean);
        forceMean = options.ForceNew && (firstMean || ~options.SaveData);
        [values, nativeDates] = steric2gmsl(product, Type=schema.types(iType), ...
            Depth=schema.depths(iDepth), ForceNew=forceMean, ...
            SaveData=options.SaveData, BeQuiet=options.BeQuiet);
        means.(schema.meanNames(iMean)) = values;
        firstMean = false;
    end
    names = cellstr(schema.gridNames(schema.available));
    data = load(stericfilepath(product), 'lon', 'lat', names{:});
    lon = double(data.lon(:));
    lat = double(data.lat(:));
    data = rmfield(data, {'lon', 'lat'});
    data = prepareDeepSupport(data, schema);
    useMean = ~isempty(timestep) && (isstring(timestep) || strcmp(options.Interpolation, 'mean'));
    if useMean, data = addResiduals(data, schema); end
    [data, dates] = processsterictime(data, nativeDates, timestep, options.Interpolation, []);
    if ~useMean, data = addResiduals(data, schema); end
    fields = struct();
    fieldNames = replace(schema.meanNames, 'Gmsl', 'Sl');
    transposeGrid = numel(lat) == size(data.(names{1}), 1);
    for iField = find(schema.supported(:)).'
        values = data.(fieldNames(iField));
        if transposeGrid, values = permute(values, [2, 1, 3]); end
        fields.(replace(fieldNames(iField), 'Sl', 'Slept')) = values;
    end
    meanNames = fieldnames(means);
    for iMean = 1:numel(meanNames)
        means.(meanNames{iMean}) = reshape(means.(meanNames{iMean}), 1, 1, []);
    end
    [means, meanDates] = processsterictime(means, nativeDates, timestep, options.Interpolation, []);
    assert(isequal(dates, meanDates), 'ULMO:steric2slept:Dates', 'Grid and GMSL dates differ.');
    for iMean = 1:numel(meanNames), means.(meanNames{iMean}) = means.(meanNames{iMean})(:); end
    % Match the loader's full GMSL convention under nonlinear interpolation.
    for iType = 1:numel(schema.types)
        if all(schema.supported(iType, :))
            means.(schema.meanNames(iType, 1)) = ...
                means.(schema.meanNames(iType, 2)) + means.(schema.meanNames(iType, 3));
        end
    end
end

function data = prepareDeepSupport(data, schema)
    % A missing supplied deep layer contributes zero on an observed upper
    % column. Never infer an unsupported depth or reconstruct full-shallow.
    % Apply this before averaging: mean([NaN,d]) differs from mean([0,d]).
    for iType = find(schema.split).'
        full = data.(schema.gridNames(iType, 1));
        shallow = data.(schema.gridNames(iType, 2));
        deep = data.(schema.gridNames(iType, 3));
        support = isfinite(full) & isfinite(shallow);
        deep(~support) = NaN;
        deep(support & isnan(deep)) = 0;
        % Inf on fitted support remains invalid and is rejected by the fitter.
        data.(schema.gridNames(iType, 3)) = deep;
    end
end

function data = addResiduals(data, schema)
    names = replace(schema.meanNames, 'Gmsl', 'Sl');
    for iDepth = 1:numel(schema.depths)
        for iType = 4:5
            if schema.supported(iType, iDepth)
                data.(names(iType, iDepth)) = ...
                    data.(names(1, iDepth)) - data.(names(iType-2, iDepth));
            end
        end
    end
end

function [coefficients, fits, buffer, weighting] = fitRecord(fields, dates, lon, lat, domain, L, basis, options)
    keepLon = 1:numel(lon);
    if abs(abs(lon(end)-lon(1))-360) <= 3e-5, keepLon(end) = []; end
    [lon, orderLon] = sort(lon(keepLon));
    keepLon = keepLon(orderLon);
    [lat, orderLat] = sort(lat);
    steps = diff(lat);
    if all(abs(steps-steps(1)) <= 3e-5)
        latitudeWeights = cosd(lat);
        weighting = 'cosine-latitude';
    else
        edges = [lat(1)-steps(1)/2; (lat(1:end-1)+lat(2:end))/2; lat(end)+steps(end)/2];
        latitudeWeights = diff(sind(max(-90, min(90, edges))));
        weighting = 'spherical-midpoint-bands';
    end
    mask = domainmask(domain, lon, lat, Format='ndgrid', SaveData=false, BeQuiet=true);
    buffer = 0;
    if L > 0, buffer = 90/L; end
    observationMask = mask;
    if buffer > 0
        observationMask = domainmask(domain, lon, lat, Format='ndgrid', ...
            exteriorBuffer=buffer, SaveData=options.SaveData, ForceNew=options.ForceNew, BeQuiet=true);
    end
    [gridLon, gridLat] = ndgrid(lon, lat);
    weights = repmat(latitudeWeights.', numel(lon), 1);
    totalArea = sum(weights(mask));
    assert(totalArea > 0, 'ULMO:steric2slept:EmptyTarget', 'The target contains no positive-area grid cells.');
    points = find(observationMask & weights > 0);
    weights = weights(points);
    target = mask(points);
    J = size(basis{1}, 2);
    nTime = numel(dates);
    assert(nTime > 0, 'ULMO:steric2slept:EmptyRecord', 'Temporal selection produced no epochs.');
    design = [];
    coefficients = struct();
    fits = struct();
    for name = string(fieldnames(fields)).'
        values = fields.(name)(keepLon, orderLat, :);
        values = reshape(values, [], nTime);
        values = values(points, :);
        assert(~any(isinf(values), 'all'), 'ULMO:steric2slept:InvalidObservation', ...
            '%s contains infinite observations.', name);
        observed = isfinite(values);
        nObs = sum(observed, 1).';
        areaFraction = min(1, (weights(target).' * observed(target, :)).' / totalArea);
        epochs = find(nObs > 0);
        assert(~isempty(epochs), 'ULMO:steric2slept:EmptyRecord', 'No observed cells for %s.', name);
        try
            if isempty(design)
                [c, ~, ~, stats, design] = xyzs2slep(values(:, epochs), gridLon(points), gridLat(points), ...
                    domain, L, GVN=basis, truncation=J, weight=weights, missingPolicy="omit");
            else
                [c, ~, ~, stats] = xyzs2slep(values(:, epochs), design, weight=weights, missingPolicy="omit");
            end
        catch cause
            field = regexpi(cause.message, '\bfield (\d+)\b', 'tokens', 'once');
            context = "Fit failed for " + name;
            if ~isempty(field), context = context + " at " + string(dates(epochs(str2double(field{1})))); end
            throwAsCaller(addCause(MException('ULMO:steric2slept:FitFailed', '%s: %s', context, cause.message), cause));
        end
        coeffs = nan(nTime, J);
        condNum = nan(nTime, 1);
        residualRms = nan(nTime, 1);
        coeffs(epochs, :) = c.';
        condNum(epochs) = stats.condNum;
        residualRms(epochs) = stats.residualRms;
        coefficients.(name) = coeffs;
        fits.(name) = table(nObs, areaFraction, condNum, residualRms);
    end
end

function writeCache(destination, coefficients, means, dates, fits, constant, product, domainId, L, defaultJ, buffer, weighting, configuration)
    folder = fileparts(destination);
    if ~isfolder(folder), mkdir(folder); end
    path = [tempname(folder), '.nc'];
    cleanup = onCleanup(@() removePartial(path));
    % Define coordinate variables before data variables attach dimension scales.
    % Replacing implicit scales after many fields can corrupt HDF5 references.
    writeVector(path, 'time', seconds(dates-datetime(2000, 1, 1)), 'time');
    ncwriteatt(path, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
    ncwriteatt(path, 'time', 'calendar', 'proleptic_gregorian');
    writeVector(path, 'mode', (1:numel(constant)).', 'mode');
    names = string(fieldnames(coefficients));
    for name = names.'
        write2nc(path, name, coefficients.(name).', Dimensions={'mode', 'time'}, Attributes={'units', 'm'});
        for diagnostic = string(fits.(name).Properties.VariableNames)
            variable = name + "_" + diagnostic;
            writeVector(path, variable, fits.(name).(diagnostic), 'time');
            unit = '1';
            if diagnostic == "residualRms", unit = 'm'; end
            ncwriteatt(path, variable, 'units', unit);
        end
    end
    for name = string(fieldnames(means)).'
        writeVector(path, name, means.(name), 'time');
        ncwriteatt(path, name, 'units', 'm');
    end
    writeVector(path, 'constantSlept', constant, 'mode');
    ncwriteatt(path, 'constantSlept', 'units', '1');
    ncwriteatt(path, 'constantSlept', 'definition', 'G(1,:).'' for retained orthonormal 4*pi real harmonic basis; unit C00');
    ncwriteatt(path, '/', 'product', char(product));
    ncwriteatt(path, '/', 'domainId', domainId);
    ncwriteatt(path, '/', 'temporalConfiguration', configuration);
    ncwriteatt(path, '/', 'availableFields', char(join(names, ',')));
    ncwriteatt(path, '/', 'availableMeans', strjoin(fieldnames(means), ','));
    ncwriteatt(path, '/', 'L', L);
    ncwriteatt(path, '/', 'truncation', numel(constant));
    ncwriteatt(path, '/', 'defaultTruncation', defaultJ);
    ncwriteatt(path, '/', 'observationBufferDegrees', buffer);
    ncwriteatt(path, '/', 'weighting', weighting);
    ncwriteatt(path, '/', 'deepSupport', ...
        'Paired finite native full/shallow support; missing deep values are zero before temporal processing');
    nccreate(path, 'creationDate', 'Datatype', 'double');
    ncwrite(path, 'creationDate', posixtime(datetime('now', 'TimeZone', 'UTC')));
    ncwriteatt(path, 'creationDate', 'units', 'seconds since 1970-01-01 00:00:00 UTC');
    movefile(path, destination, 'f');
end

function writeVector(path, name, values, dimension)
    nccreate(path, name, 'Dimensions', {dimension, numel(values)}, 'Datatype', 'double', 'Format', 'netcdf4');
    ncwrite(path, name, double(values));
end

function [coefficients, dates, fit, halo, constant] = readCache(path, field, haloName, limits, outputCount)
    available = split(string(ncreadatt(path, '/', 'availableFields')), ',');
    assert(ismember(field, available), 'ULMO:steric2slept:UnsupportedField', 'Cache does not supply %s.', field);
    if outputCount >= 5
        availableMeans = split(string(ncreadatt(path, '/', 'availableMeans')), ',');
        assert(ismember(haloName, availableMeans), 'ULMO:steric2slept:UnsupportedField', 'Cache does not supply %s.', haloName);
    end
    dates = readshnctime(path);
    rows = find(selectDates(dates, limits));
    dates = dates(rows);
    J = ncreadatt(path, '/', 'truncation');
    count = numel(rows);
    coefficients = zeros(count, J);
    if count > 0, coefficients = reshape(ncread(path, field, [1 rows(1)], [J count]), J, count).'; end
    fit = table();
    if outputCount >= 3
        for name = ["nObs", "areaFraction", "condNum", "residualRms"]
            values = zeros(count, 1);
            if count > 0, values = reshape(ncread(path, field + "_" + name, rows(1), count), [], 1); end
            fit.(name) = values;
        end
    end
    halo = zeros(0, 1);
    if outputCount >= 5
        halo = zeros(count, 1);
        if count > 0, halo = reshape(ncread(path, haloName, rows(1), count), [], 1); end
    end
    constant = [];
    if outputCount >= 6, constant = reshape(ncread(path, 'constantSlept'), [], 1); end
end

function rows = selectDates(dates, limits)
    rows = true(size(dates));
    if ~isempty(limits), rows = dates >= limits(1) & dates <= limits(2); end
end

function removePartial(path)
    if isfile(path), delete(path); end
end
