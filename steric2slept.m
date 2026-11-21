%% STERIC2SLEPT - Fits steric fields to regional Slepian functions.
%
% Syntax
%   [c, dates, GVN, h, constSlept] = STERIC2SLEPT(product, domain, L)
%   [__] = STERIC2SLEPT(product, domain, L, timestep, timelim, Component="total")
%   [__] = STERIC2SLEPT(__, OutputFormat="traditional", SaveData=false)
%
% Input arguments
%   product - Prepared StericProduct identifier.
%   domain - GeoDomain, region name, GeoDomain argument cell or polygon.
%       Numeric polygons contain longitude/latitude vertices.
%       Size: (nVertex, 2) for a polygon.
%       Unit: degrees east and north for polygon coordinates.
%   L - Maximum spherical harmonic degree, a nonnegative integer.
%   timestep (optional) - Temporal interpolation step or averaging convention.
%       - Positive numeric scalar or duration for fixed steps
%       - 'midmonth' for calendar-month means, or 'GRACE' for CSR RL06 degree-60
%         solution windows. Both conventions average native samples and return
%         window midpoints.
%       Unit: days for numeric values.
%       Default step: [] (native timestamps).
%   timelim (optional) - Inclusive datetime or datenum output limits.
%       Increasing limits, applied after processing the complete record.
%       Size: (1, 2) or (2, 1).
%       Unit: days for numeric values.
%       Default range: [] (complete record).
%   Component (name-value) - total, thermosteric, halosteric, non-thermosteric
%       (total minus thermosteric), or non-halosteric (total minus halosteric).
%       Default type: total.
%   Depth (name-value) - native, full-depth, shallow or deep.
%       Native resolves to the prepared product coverage. Upper-ocean products
%       reject full-depth and deep requests. Default depth: native.
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
%   lon - Regular longitude coordinates spanning one complete cycle.
%       A duplicated cyclic endpoint is allowed and counted only once.
%       Size: vector with nLon elements.
%       Unit: degrees east.
%   lat - Strictly monotone latitude coordinates; spacing may be nonuniform.
%       Size: vector with nLat elements.
%       Unit: degrees north.
%   dates - Strictly increasing datetime or datenum timestamps.
%       Size: vector with nNativeTime elements.
%       Unit: days for datenum values.
%   <component>Sl - Supplied native grids in STERICFILEPATH(product).
%       Fields and complete depth splits follow STERICFIELDSCHEMA. All grids
%       share one orientation; square grids use (lat, lon, time) ordering.
%       Size: (nLat, nLon, nNativeTime) or (nLon, nLat, nNativeTime).
%       Unit: m.
%
% Output arguments
%   stericSlept - Raw selected component/depth coefficients.
%       Size: (nTime, J) for timefirst; (J, nTime) for traditional.
%       Unit: Unit.
%   dates - Returned timestamps in TimeFormat.
%       Size: (nTime, 1).
%       Unit: days for datenum output.
%   GVN - Cell array {G, V, N}: retained harmonic basis, concentration
%       eigenvalues and full Shannon number. G uses the real 4*pi-normalized
%       harmonic convention.
%       Size: (1, 3); G is ((L+1)^2, J), V is (J, 1), N is scalar.
%   haloGmsl - Selected depth's temporally processed halosteric GMSL.
%       Independent of Component. For split products, full GMSL is the sum of the
%       processed shallow/deep means. Unavailable halosteric means raise an
%       error only when this output is requested.
%       Size: (nTime, 1) for timefirst; (1, nTime) for traditional.
%       Unit: Unit.
%   constSlept - Harmonic projection of unit C00: GVN{1}(1,:).'.
%       Independent of data masks and output units/order; dimensionless.
%       Size: (J, 1).
%
% Output file content
%   <component>Slept - All supported raw coefficient fields, in metres.
%       MATLAB NetCDF dimensions: (mode, time).
%   <component>Gmsl - All applicable processed global means in metres, (time).
%   <component>Slept_<diagnostic> - Per-field fitting diagnostics, (time).
%   constantSlept - Unit-C00 projection, (mode), with units 1.
%   time - Processed epoch timestamps.
%       Size: (nTime, 1).
%       Unit: seconds since 2000-01-01 00:00:00, interpreted as UTC.
%   mode - Retained mode indices, 1:J.
%       Size: (J, 1).
%   creationDate - UTC time when the projection cache was written.
%       Size: scalar.
%       Unit: seconds since 1970-01-01 00:00:00 UTC.
%
% Notes
%   For timefirst, caller correction is c - h * constSlept.'; traditional
%   uses c - constSlept * h.
%
% See also
%   steric2lonlatt, steric2gmsl
%
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function [stericSlept, dates, GVN, haloGmsl, constSlept] = steric2slept( ...
        product, domain, L, timestep, timelim, options)

    arguments (Input)
        product (1, 1) StericProduct
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        timestep {mustBeTimeStep} = []
        timelim {mustBeTimeRange} = []
        options.Component (1, 1) string {mustBeStericComponent} = "total"
        options.Depth (1, 1) string {mustBeStericDepth} = "native"
        options.Interpolation (1, 1) string = "linear"
        options.OutputFormat (1, 1) string {mustBeMember(options.OutputFormat, ["timefirst", "traditional"])} = "timefirst"
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["m", "mm"])} = "m"
        options.TimeFormat DateFormat = "datetime"
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.GVN cell = {}
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
    end

    [timestep, options.Interpolation] = normalisetimeconfig(timestep, options.Interpolation);

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    if iscell(domain)
        domain = GeoDomain(domain{:});
    elseif ischar(domain) || isstring(domain)
        domain = GeoDomain(domain);
    end

    assert((isa(domain, 'GeoDomain') && isscalar(domain)) || ...
        (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2 && all(isfinite(domain), 'all')), ...
        'ULMO:steric2slept:Domain', 'Supply a GeoDomain, region name or finite polygon.');
    nativeDepth = product.nativeDepth();

    if string(options.Depth) == "native"
        options.Depth = nativeDepth;
    end

    assert(~ismember(string(options.Depth), ["full-depth", "deep"]) || product.supportsFullDepth(), ...
        'ULMO:steric2slept:UnsupportedField', ...
        'Product %s does not supply Depth=%s.', char(product), options.Depth);
    schema = stericfieldschema(struct('name', {}, 'size', {}), nativeDepth);
    iType = find(schema.components == options.Component);
    iDepth = find(schema.depths == options.Depth);
    field = replace(schema.meanNames(iType, iDepth), 'Gmsl', 'Slept');
    haloName = schema.meanNames(3, iDepth);
    supplied = ~isempty(options.GVN);
    [cache.stem, cache.config, cache.domainId] = ...
        cacheStem(product, domain, L, timestep, options.Interpolation);
    path = '';

    if ~options.ForceNew && ~supplied
        path = findCache(cache.stem, options.truncation);
    end

    if ~isempty(path) && isfile(path)

        if options.BeQuiet < 1
            fprintf('[ULMO>%s] Loading %s.\n', mfilename, path);
        end

        [stericSlept, dates, haloGmsl, constSlept] = loadCache( ...
            path, field, haloName, timelim, nargout);

        if nargout >= 3
            GVN = prepareslepianbasis(domain, L, size(stericSlept, 2), {}, ...
                RequireOrthonormal = true, SaveData = options.SaveData);
        end

    else
        [stericSlept, dates, GVN, haloGmsl, constSlept] = steric2slept_core( ...
            product, domain, L, timestep, timelim, options, cache, nargout >= 4);
    end

    if options.Unit == "mm"
        stericSlept = 1000 * stericSlept;
        haloGmsl = 1000 * haloGmsl;

    end

    if options.OutputFormat == "traditional"
        stericSlept = stericSlept.';
        haloGmsl = haloGmsl.';
    end

    if options.TimeFormat == "datenum"
        dates = datenum(dates); %#ok<DATNM>
    end

end

%% STERIC2SLEPT_CORE - Computes and optionally caches all supported steric projections.
% Returns the requested coefficients and GMSL in metres, with time first.
% Output dates remain datetime values until the public function formats them.
% needHaloGmsl records whether the caller requested the halosteric output.
function [stericSlept, dates, GVN, haloGmsl, constSlept] = steric2slept_core( ...
        product, domain, L, timestep, timeRange, options, cache, needHaloGmsl)

    if options.BeQuiet < 1
        fprintf('[ULMO>%s] Fitting all fields for %s.\n', mfilename, char(product));
    end

    metadata = whos('-file', stericfilepath(product));
    schema = stericfieldschema(metadata, product.nativeDepth());
    iComp = find(schema.components == options.Component);
    iDepth = find(schema.depths == options.Depth);
    field = replace(schema.meanNames(iComp, iDepth), 'Gmsl', 'Slept');
    haloName = schema.meanNames(3, iDepth);
    assert(schema.supported(iComp, iDepth), 'ULMO:steric2slept:UnsupportedField', ...
        'Product %s does not supply Component=%s, Depth=%s.', ...
        char(product), options.Component, options.Depth);

    if needHaloGmsl
        assert(schema.supported(3, iDepth), 'ULMO:steric2slept:UnsupportedField', ...
            'Product %s does not supply halosteric GMSL at Depth=%s.', ...
            char(product), options.Depth);
    end

    % Prepare the GMSLs and gridded fields on the same time axis.
    [means, meanDates, nativDates] = prepareGmsl(product, timestep, options, schema);

    % A supplied basis bypasses projection storage, not an explicit request
    % to save native diagnostics. SaveData=false is propagated everywhere.
    GVN = prepareslepianbasis(domain, L, options.truncation, options.GVN, ...
        RequireOrthonormal = true, SaveData = options.SaveData);
    constSlept = GVN{1}(1, :).';

    [fields, dates, lon, lat] = prepareGriddedFields( ...
        product, nativDates, timestep, options, schema);

    assert(isequal(dates, meanDates), 'ULMO:steric2slept:Dates', 'Grid and GMSL dates differ.');

    % Fit every supported field before selecting the requested output.
    [coeffs, fits, buffer, wghtMthd] = fitgriddedslepian( ...
        fields, dates, lon, lat, domain, L, GVN, ...
        RequireGlobalLongitude = true, ...
        SaveData = options.SaveData, ForceNew = options.ForceNew);

    % Save the complete record before selecting the requested dates.
    if options.SaveData && isempty(options.GVN)
        defaultJ = min((L + 1) ^ 2, max(1, round(GVN{3})));
        path = sprintf('%s-%d-m.nc', cache.stem, size(GVN{1}, 2));
        saveCache(path, coeffs, means, dates, fits, constSlept, ...
            product, cache.domainId, L, defaultJ, buffer, wghtMthd, cache.config);
    end

    rows = selectDates(dates, timeRange);
    stericSlept = coeffs.(field)(rows, :);
    dates = dates(rows);
    haloGmsl = zeros(0, 1);

    if needHaloGmsl
        haloGmsl = means.(haloName)(rows);
    end

end

%% cacheStem
function [stem, config, domainId] = cacheStem(product, domain, L, timestep, mthd)
    root = getenv('IFILES');
    assert(isfolder(root), 'ULMO:steric2slept:MissingRoot', 'Set IFILES to the data directory.');

    if isa(domain, 'GeoDomain')
        domainId = domain.Id;

        if ismember(domain.Domain, ["arctic", "antarctica"])
            domainId = [domainId, '-rotb'];
        end

    else
        domainId = hash(domain, 'sha1');
    end

    [~, ~, ~, config] = normalisetimeconfig(timestep, mthd);
    stem = fullfile(root, 'STERIC', char(product), 'SlepianExpansions', ...
        sprintf('steric2slept-%s-%s-%s-%d', char(product), config, domainId, L));
end

%% findCache
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

%% PREPAREGMSL - Loads and temporally processes the native GMSLs.
% Obtain every diagnostic through the public native-GMSL interface.
% A saved forced call refreshes all means; later calls reuse that set.
% Unsaved forced calls must each bypass the existing native cache.
function [means, meanDates, nativeDates] = prepareGmsl( ...
        product, timestep, options, schema)
    means = struct();
    firstMean = true;
    fieldSupported = schema.supported;
    fieldSupported(4:5, :) = false;

    % Get the GMSLs
    for iField = find(fieldSupported(:)).'
        [iComp, iDepth] = ind2sub(size(schema.supported), iField);
        % Temporary patch: We have to query STERIC2GMSL once per field, but
        % they are all recalculated and saved when the first time series is
        % requested
        forceNewGmsl = options.ForceNew && (firstMean || ~options.SaveData);
        [nativGmsl, nativeDates] = steric2gmsl(product, ...
            schema.components(iComp), schema.depths(iDepth), ...
            ForceNew = forceNewGmsl, SaveData = options.SaveData, BeQuiet = options.BeQuiet);
        means.(schema.meanNames(iField)) = nativGmsl;
        firstMean = false;
    end

    % Interpolate the GMSLs
    meanNames = fieldnames(means);

    for iField = 1:numel(meanNames)
        means.(meanNames{iField}) = reshape(means.(meanNames{iField}), 1, 1, []);
    end

    [means, meanDates] = processfieldtime(means, nativeDates, timestep, options.Interpolation, []);

    for iField = 1:numel(meanNames)
        means.(meanNames{iField}) = means.(meanNames{iField})(:);
    end

    % Make sure the shallow and and deep GMSLs sum to the full-depth GMSL
    for iComp = 1:3

        if all(schema.supported(iComp, :))
            means.(schema.meanNames(iComp, 1)) = ...
                means.(schema.meanNames(iComp, 2)) + means.(schema.meanNames(iComp, 3));
        end

    end

end

%% PREPAREGRIDDEDFIELDS - Loads and processes steric grids for Slepian fitting.
function [fields, dates, lon, lat] = prepareGriddedFields( ...
        product, nativeDates, timestep, options, schema)
    % Interpolate the fields temporally
    fieldNames = cellstr(schema.gridNames(schema.available));
    data = load(stericfilepath(product), 'lon', 'lat', fieldNames{:});
    lon = data.lon(:);
    lat = data.lat(:);
    data = rmfield(data, {'lon', 'lat'}); % Remove since they can't be interpolated
    data = fillAbsentDeepSteric(data, schema);
    [data, dates] = processfieldtime(data, nativeDates, timestep, options.Interpolation, []);

    % Compute the residual fields
    processed = struct();
    fieldNames = replace(schema.meanNames, 'Gmsl', 'Sl');

    for iField = find(schema.available(:)).'
        [iComp, iDepth] = ind2sub(size(schema.available), iField);
        processed.(fieldNames(iComp, iDepth)) = data.(schema.gridNames(iField));
    end

    processed = addResidualComponents(processed, schema);
    data = processed;

    % Load the gridded steric sea level fields
    fields = struct();
    transGrid = numel(lat) == size(data.(fieldNames(find(schema.supported, 1))), 1);

    for iField = find(schema.supported(:)).'
        stericSl = data.(fieldNames(iField));

        if transGrid
            stericSl = permute(stericSl, [2, 1, 3]);
        end

        fields.(replace(fieldNames(iField), 'Sl', 'Slept')) = stericSl;
    end

end

%% FILLABSENTDEEPSTERIC - Fills absent deep-ocean steric contributions with zero.
% A grid cell within the ocean domain can have the deep steric sea level
% marked as invalid because the cell is in shallow water. However, these
% cells should be interpreted as having 0 contribution to deep steric sea
% level instead of being marked as invalid.
function data = fillAbsentDeepSteric(data, schema)
    % Loop through all components that supports splitting into the upper-
    % and deep-ocean fields
    for iComp = find(schema.split).'
        full = data.(schema.gridNames(iComp, 1));
        upper = data.(schema.gridNames(iComp, 2));
        deep = data.(schema.gridNames(iComp, 3));
        % If a cell has valid values in both the full-depth and upper-ocean
        % field but not in the deep-ocean field, this cell is in shallow
        % water and contributes 0 to the deep-ocean steric sea levels.
        validMask = isfinite(full) & isfinite(upper);
        deep(~validMask) = NaN;
        deep(validMask & isnan(deep)) = 0;
        data.(schema.gridNames(iComp, 3)) = deep;
    end

end

%% ADDRESIDUALCOMPONENTS - Computes the residual fields from differencing
function data = addResidualComponents(data, schema)
    names = replace(schema.meanNames, 'Gmsl', 'Sl');

    for iDepth = 1:numel(schema.depths)

        for iType = 4:5 % for the non-thermosteric and non-halosteric fields

            if schema.supported(iType, iDepth)
                data.(names(iType, iDepth)) = ...
                    data.(names(1, iDepth)) - data.(names(iType - 2, iDepth));
            end

        end

    end

end

%% SAVECACHE - Saves the computation results into the cache file.
function saveCache(outPath, coeffs, means, dates, fits, const, ...
        product, domainId, L, defaultJ, buf, wghtMthd, config)
    folder = fileparts(outPath);

    if ~isfolder(folder)
        mkdir(folder);
    end

    tmpPath = [tempname(folder), '.nc'];
    cleanup = onCleanup(@() deletefileifpresent(tmpPath));

    % Save the time axis
    writeVector(tmpPath, 'time', seconds(dates - datetime(2000, 1, 1)), 'time');
    ncwriteatt(tmpPath, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
    ncwriteatt(tmpPath, 'time', 'calendar', 'proleptic_gregorian');
    % Save the Slepian mode axis
    writeVector(tmpPath, 'mode', (1:numel(const)).', 'mode');

    names = string(fieldnames(coeffs));

    % Save the coefficients and their diagnostics
    for name = names.'
        write2nc(tmpPath, name, coeffs.(name).', ...
            Dimensions = {'mode', 'time'}, Attributes = {'units', 'm'});

        for diagnostic = string(fits.(name).Properties.VariableNames)
            variable = name + "_" + diagnostic;
            writeVector(tmpPath, variable, fits.(name).(diagnostic), 'time');
            unit = '1';

            if diagnostic == "residualRms"
                unit = 'm';
            end

            ncwriteatt(tmpPath, variable, 'units', unit);
        end

    end

    % Save the the GMSL time series
    for name = string(fieldnames(means)).'
        writeVector(tmpPath, name, means.(name), 'time');
        ncwriteatt(tmpPath, name, 'units', 'm');
    end

    % Save the coefficients to construct the constant box car field
    writeVector(tmpPath, 'constantSlept', const, 'mode');
    ncwriteatt(tmpPath, 'constantSlept', 'units', '1');
    ncwriteatt(tmpPath, 'constantSlept', 'definition', 'G(1,:).'' for retained orthonormal 4*pi real harmonic basis; unit C00');
    % Save metadata
    ncwriteatt(tmpPath, '/', 'product', char(product));
    ncwriteatt(tmpPath, '/', 'domainId', domainId);
    ncwriteatt(tmpPath, '/', 'temporalConfiguration', config);
    ncwriteatt(tmpPath, '/', 'availableFields', char(join(names, ',')));
    ncwriteatt(tmpPath, '/', 'availableMeans', strjoin(fieldnames(means), ','));
    ncwriteatt(tmpPath, '/', 'L', L);
    ncwriteatt(tmpPath, '/', 'truncation', numel(const));
    ncwriteatt(tmpPath, '/', 'defaultTruncation', defaultJ);
    ncwriteatt(tmpPath, '/', 'observationBufferDegrees', buf);
    ncwriteatt(tmpPath, '/', 'weighting', wghtMthd);
    nccreate(tmpPath, 'creationDate', 'Datatype', 'double');
    ncwrite(tmpPath, 'creationDate', posixtime(datetime('now', 'TimeZone', 'UTC')));
    ncwriteatt(tmpPath, 'creationDate', 'units', 'seconds since 1970-01-01 00:00:00 UTC');

    movefile(tmpPath, outPath, 'f');
end

%% WRITEVECTOR - Writes a variable to a netCDF file as a double vector.
function writeVector(path, name, vals, dim)
    nccreate(path, name, 'Dimensions', {dim, numel(vals)}, ...
        'Datatype', 'double', 'Format', 'netcdf4');
    ncwrite(path, name, double(vals));
end

%% LOADCACHE - Loads the data from the cached file.
function [coeffs, dates, halo, const] = loadCache( ...
        path, field, haloName, limits, outCnt)
    availFields = split(string(ncreadatt(path, '/', 'availableFields')), ',');
    assert(ismember(field, availFields), 'ULMO:steric2slept:UnsupportedField', ...
        'The cache file only stored these fields: %s. But ''%s'' is requested.', ...
        strjoin(availFields, ','), field);

    if outCnt >= 4
        availGmsls = split(string(ncreadatt(path, '/', 'availableMeans')), ',');
        assert(ismember(haloName, availGmsls), 'ULMO:steric2slept:UnsupportedField', ...
            'The cache file only stored these GMSLs: %s. But ''%s'' is requested.', ...
            availGmsls, haloName);
    end

    dates = readshnctime(path);
    dateIds = find(selectDates(dates, limits));
    dates = dates(dateIds);
    J = ncreadatt(path, '/', 'truncation');
    nEpochs = numel(dateIds);
    coeffs = zeros(nEpochs, J);

    if nEpochs > 0
        coeffs = reshape(ncread(path, field, [1, dateIds(1)], [J, nEpochs]), J, nEpochs).';
    end

    halo = zeros(0, 1);

    if outCnt >= 4
        halo = zeros(nEpochs, 1);

        if nEpochs > 0
            halo = reshape(ncread(path, haloName, dateIds(1), nEpochs), [], 1);
        end

    end

    const = [];

    if outCnt >= 5
        const = reshape(ncread(path, 'constantSlept'), [], 1);
    end

end

%% SELECTDATES - Finds time steps that are inside the specified time range.
function mask = selectDates(dates, timelim)
    mask = true(size(dates));

    if isempty(timelim)
        return
    end

    mask = dates >= timelim(1) & dates <= timelim(2);

end
