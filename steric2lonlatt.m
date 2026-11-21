%% STERIC2LONLATT - Loads and interpolates gridded steric sea-level fields.
%
% Syntax
%   [steric, stericSigma, dates] = STERIC2LONLATT(product)
%   [steric, stericSigma, dates] = STERIC2LONLATT(product, timestep, meshsize)
%   [steric, stericSigma, dates] = STERIC2LONLATT(product, timestep, meshsize, timelim)
%   [steric, stericSigma, dates] = STERIC2LONLATT(__, "Name", value)
%   [steric, stericSigma, dates, lon, lat, halostericGmsl] = STERIC2LONLATT(__)
%
% Input arguments
%   product (optional) - Prepared product: SIO, EN4c14, CMEMS, NOAA, ORAS5,
%       JAMSTEC, ISAS, IPRC or IAPv4.
%       Default product: EN4c14.
%   timestep (optional) - Temporal interpolation step or averaging convention.
%       - Positive numeric scalar or duration for fixed steps
%       - 'midmonth' for calendar-month means, or 'GRACE' for CSR RL06 degree-60
%         solution windows. Both conventions average native samples and return
%         window midpoints.
%       Unit: days for numeric values.
%       Default step: [] (native timestamps).
%   meshsize (optional) - Spatial interpolation spacing.
%       A two-element [spacing, longitudeCentre] vector also sets LonOrigin
%       when that option is empty.
%       With LonOrigin alone, the spacing is 0.5 degrees.
%       Unit: degrees.
%       Default spacing: [] (no spatial interpolation unless LonOrigin is set).
%   timelim (optional) - Inclusive output time limits as datetime or datenum.
%       Applied only to returned values; caches retain the complete processed record.
%       Size: (1, 2) or (2, 1), or empty.
%       Unit: days for datenum values.
%       Default limits: [] (full record).
%   LonOrigin (name-value) - Longitude centre of the interpolated grid.
%       With meshsize alone, the grid is centred on 180 degrees east.
%       Unit: degrees east.
%       Default centre: [] (no spatial interpolation unless meshsize is set).
%   Interpolation (name-value) - Temporal and spatial interpolation method.
%       An INTERP1 method, or 'mean' for fixed-step averaging in half-open
%       windows centred on output dates. For 'midmonth' and 'GRACE', this
%       option controls spatial interpolation only. Spatial interpolation
%       uses 'linear' when 'mean' is selected.
%       Default method: linear.
%   OutputFormat (name-value) - meshgrid or ndgrid spatial ordering.
%       Default format: meshgrid.
%   TimeFormat (name-value) - datetime or datenum output dates.
%       Default format: datetime.
%   Unit (name-value) - m or mm output sea level units.
%       Default unit: m.
%   Depth (name-value) - native, full-depth, shallow, or deep.
%       Native resolves to the prepared product's depth range. Full-depth
%       means integration from the sea surface to the seafloor. Upper-ocean
%       products reject full-depth and deep requests.
%       Default depth: native.
%   Component (name-value) - Selected steric component.
%       total, thermosteric, halosteric, non-thermosteric, or non-halosteric.
%       Default component: total.
%   ForceNew (name-value) - Whether to reprocess every supported component grid and
%       recalculate every native and processed GMSL.
%       Does not regenerate native grids or other cache configurations.
%       Default flag: false.
%   SaveData (name-value) - Whether to save native means and derived grids/means.
%       false allows cache reads and in-memory rebuilds without any writes.
%       Default flag: true.
%   BeQuiet (name-value) - Whether to suppress messages.
%       true suppresses messages; false prints all messages; 0.5 prints
%       cache-loading and saving messages.
%       Default level: 0.5.
%   CallChain (name-value) - Calling function names for diagnostic messages.
%       Default chain: {}.
%
% Input file content
%   lon, lat - Longitude and latitude coordinate vectors.
%       Size: vectors with nLon and nLat elements, respectively.
%       Unit: degrees east and north, respectively.
%   dates - Native or cached timestamps.
%       Size: vector with N elements.
%   <component>Sl - Prepared grids for the requested components/depths.
%       Components are steric, thermosteric and halosteric, with shallow/deep
%       prefixes when supplied, e.g. shallowHalostericSl. Residual fields
%       also require the corresponding total steric grid.
%       Size: (nLat, nLon, N) or (nLon, nLat, N).
%       Unit: m.
%
% Output arguments
%   steric - Selected steric sea level grid.
%       Size: (nLat, nLon, N) for meshgrid; (nLon, nLat, N) for ndgrid.
%       Unit: Unit.
%   stericSigma - Zero placeholder for product uncertainty.
%       Size: same as steric.
%       Unit: Unit.
%   dates - Returned timestamps in TimeFormat.
%       Size: (N, 1).
%       Unit: days for datenum output.
%   lon - Returned longitudes.
%       Size: (1, nLon) for meshgrid; (nLon, 1) for ndgrid.
%       Unit: degrees east.
%   lat - Returned latitudes.
%       Size: (nLat, 1) for meshgrid; (1, nLat) for ndgrid.
%       Unit: degrees north.
%
%   halostericGmsl - Native global halosteric mean for the requested depth,
%       independently resampled to the returned dates.
%       Size: (N, 1).
%       Unit: Unit.
%
% Output file content
%   lon, lat, dates - Processed coordinates and timestamps.
%       Same conventions as the input file; saved before output formatting.
%   <component>Sl - Supplied component/depth grids in one configuration cache.
%       Includes stericSl, thermostericSl and halostericSl, with shallow/deep
%       prefixes where supported. Shallow-native products use shallow prefixes.
%       Size: as in the input file, with processed spatial/time dimensions.
%       Unit: m.
%
%   <component>Gmsl - Processed total, thermosteric and halosteric means when
%       requested or forced, with the names defined by STERIC2GMSL.
%       Difference-component means are never saved.
%       Size: (N, 1), aligned with saved dates.
%       Unit: m.
%
% Notes
%   Difference components are derived after processing when requested.
%   Ordinary grid-only requests do not calculate means. GMSL requests reuse
%   STERIC2GMSL caches and process all applicable series separately from grids.
%   Full-depth means of split components are sums of the processed layer means.
%
% See also
%   steric2gmsl
%
% Created by
%   2025/07/22, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function [steric, stericSigma, dates, lon, lat, halostericGmsl] = steric2lonlatt( ...
        product, timestep, meshsize, timelim, options)
    %% Initialisation
    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty} = 'EN4c14'
        timestep {mustBeTimeStep} = []
        meshsize = []
        timelim {mustBeTimeRange} = []
        options.LonOrigin {mustBeFinite, mustBeReal} = []
        options.Interpolation (1, 1) string = "linear"
        options.TimeFormat (1, 1) DateFormat = "datetime"
        options.OutputFormat (1, 1) MeshFormat = "meshgrid"
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["mm", "m"])} = "m"
        options.Component (1, 1) string {mustBeStericComponent} = "total"
        options.Depth (1, 1) string {mustBeStericDepth} = "native"
        options.ForceNew (1, 1) {mustBeNumericOrLogical} = false
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
        options.SaveData (1, 1) {mustBeNumericOrLogical} = true
        options.CallChain cell = {};
    end

    arguments (Output)
        steric (:, :, :) {mustBeReal}
        stericSigma (:, :, :) {mustBeReal}
        dates (:, 1) {mustBeVector}
        lon (:, :)
        lat (:, :)
        halostericGmsl (:, 1) double
    end

    lonOrigin = options.LonOrigin;
    nativeDepth = product.nativeDepth();

    if options.Depth == "native"
        options.Depth = nativeDepth;
    end

    if ismember(options.Depth, ["full-depth", "deep"]) && ~product.supportsFullDepth()
        error('ULMO:steric2lonlatt:UnsupportedField', ...
            'Product %s does not supply steric sea level at Depth=%s.', char(product), options.Depth);
    end

    intpMthd = lower(options.Interpolation);
    timeFmt = lower(options.TimeFormat);
    outFmt = lower(options.OutputFormat);
    unit = options.Unit;
    options.BeQuiet = uint8(double(options.BeQuiet) * 2);
    saveData = options.SaveData;
    callChain = [options.CallChain, {mfilename}];

    if numel(meshsize) == 2 && isempty(lonOrigin)
        lonOrigin = meshsize(2);
        meshsize = meshsize(1);
    elseif numel(meshsize) > 2
        error('ULMO:InvalidInput', ...
            'meshsize must be a scalar or a 2-element vector of (h, lonOrigin), got %d elements.', ...
            numel(meshsize));
    end

    if meshsize < 0
        error('ULMO:InvalidMeshSize', ...
            'Mesh size must be a positive real scalar, got %f.', meshsize);
    elseif isnan(meshsize)
        warning('ULMO:InvalidMeshSize', ...
        'Mesh size is NaN, argument ignored.');
        meshsize = [];
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    [timestep, intpMthd, useMean] = normalisetimeconfig(timestep, intpMthd);

    % 'default' and 'no halosteric drift' are deprecated Component values.
    schema = stericfieldschema(struct('name', {}, 'size', {}), nativeDepth);
    gridNames = schema.gridNames;
    meanNames = schema.meanNames;
    components = schema.components;
    depths = schema.depths;
    iComp = find(components == string(options.Component));
    iDepth = find(depths == string(options.Depth));
    cacheGridNames = schema.cacheGridNames;
    sourceComponents = iComp;

    if iComp > 3
        sourceComponents = [1, iComp - 2];
    end

    spatialIntpMthd = intpMthd;

    if spatialIntpMthd == "mean"
        spatialIntpMthd = "linear";
    end

    outPath = stericfilepath(product, timestep, meshsize, lonOrigin, intpMthd);
    nativePath = stericfilepath(product);
    isDerived = ~strcmp(outPath, nativePath);
    coordVars = {'lon', 'lat', 'dates'};
    needMeans = nargout >= 6 || options.ForceNew;
    halostericGmsl = zeros(0, 1);

    if isfile(nativePath)
        schema = nativefields(whos('-file', nativePath), nativeDepth);
    elseif isDerived && isfile(outPath)
        % Processed grids have explicit depth names, including shallow-native
        % products. Infer the available set when the source is offline.
        schema = stericfieldschema(whos('-file', outPath));

        if nativeDepth == "shallow"
            schema.supported(:, [1, 3]) = false;
        end

    else
        error('ULMO:LoadData:FileNotFound', 'Input file %s not found.', nativePath);
    end

    supported = schema.supported;

    if ~supported(iComp, iDepth)
        error('ULMO:steric2lonlatt:UnsupportedField', ...
            'Product %s does not supply Component=%s, Depth=%s.', ...
            char(product), options.Component, options.Depth);
    end

    cachedComponents = supported;
    cachedComponents(4:5, :) = false;
    % Only supplied components belong to this configuration's cache.
    % Completeness depends on required variable names, not a cache version.
    requiredGrids = cellstr(cacheGridNames(supported(1:3, :)));
    cacheHit = isDerived && ~options.ForceNew && isfile(outPath) && ...
        cachecontains(outPath, [coordVars, requiredGrids(:)']);

    if cacheHit
        inPath = outPath;
        dataVars = cellstr(cacheGridNames(sourceComponents, iDepth));
    else
        inPath = nativePath;

        if ~isfile(inPath)
            error('ULMO:LoadData:FileNotFound', 'Rebuilding the cache requires %s.', inPath);
        end

        dataVars = cellstr(gridNames(schema.available));

        if ~isDerived
            dataVars = cellstr(gridNames(sourceComponents, iDepth));
        end

    end

    if options.BeQuiet == 0 || (cacheHit && options.BeQuiet <= 1)
        fprintf('[ULMO>%s] Loading %s.\n', ...
            callchaintext(callChain), filehref(inPath, sprintf('%s steric product', product)));
    end

    coord = load(inPath, coordVars{:});

    if isnumeric(coord.dates)
        coord.dates = datetime(coord.dates, 'ConvertFrom', 'datenum');
    end

    coord.dates = coord.dates(:);
    data = load(inPath, dataVars{:});
    % Process the complete record. Output time limits are applied only after
    % saving or reading the shared configuration cache.
    if ~cacheHit && isDerived

        if useMean
            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
        end

        [data, coord.dates] = processfieldtime(data, coord.dates, timestep, intpMthd, []);

        if ~useMean
            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
        end

        data = componentfields(data, schema);
    end

    means = struct();
    saveMeans = false;

    if needMeans

        if nargout >= 6 && ~supported(3, iDepth)
            error('ULMO:steric2lonlatt:UnsupportedField', ...
                'Product %s has no halosteric GMSL at Depth=%s.', char(product), options.Depth);
        end

        names = cellstr(meanNames(cachedComponents));
        % Derived diagnostics can be reused as a complete set. Native-file
        % reuse still goes through steric2gmsl, which owns its cache checks.
        meansHit = cacheHit && cachecontains(outPath, names);

        if meansHit
            means = load(outPath, names{:});
        else

            if ~isfile(nativePath)
                error('ULMO:LoadData:FileNotFound', 'Rebuilding native means requires %s.', nativePath);
            end

            % The first saved forced call refreshes the complete native set.
            % Subsequent saved calls reuse that complete set. With saving
            % disabled, every forced call must bypass old native means; repeated
            % calculation is intentional and keeps the native API unchanged.
            firstMean = true;

            for iMean = find(cachedComponents(:))'
                [meanComponent, meanDepth] = ind2sub(size(supported), iMean);
                forceMean = options.ForceNew && (firstMean || ~saveData);
                [values, nativeDates] = steric2gmsl(product, components(meanComponent), depths(meanDepth), ...
                    ForceNew = forceMean, SaveData = saveData, ...
                    BeQuiet = double(options.BeQuiet) / 2);
                % Temporal helpers expect time on axis three. This temporary
                % shape does not make the scalar series a spatial grid.
                means.(meanNames(iMean)) = reshape(values, 1, 1, []);
                firstMean = false;
            end

            [means, meanDates] = processfieldtime(means, nativeDates, timestep, intpMthd, []);
            assert(isequal(meanDates(:), coord.dates(:)), ...
                'ULMO:steric2lonlatt:MeanDates', 'GMSL and grid dates must agree.');

            for iMean = 1:numel(names)
                means.(names{iMean}) = reshape(means.(names{iMean}), [], 1);
            end

            % Nonlinear interpolation need not commute with layer addition.
            % Define the processed full-depth mean from the processed layers so the
            % returned depths remain additive under methods such as pchip.
            for iMeanComponent = 1:3

                if all(supported(iMeanComponent, :))
                    means.(meanNames(iMeanComponent, 1)) = ...
                        means.(meanNames(iMeanComponent, 2)) + means.(meanNames(iMeanComponent, 3));
                end

            end

            saveMeans = true;
        end

        % Keep cached diagnostics in metres on the entire processed date axis.
        % Only the selected public vector is filtered and converted to mm.
        if nargout >= 6
            halostericGmsl = means.(meanNames(3, iDepth));

            if ~isempty(timelim)
                keep = coord.dates >= timelim(1) & coord.dates <= timelim(2);
                halostericGmsl = halostericGmsl(keep);
            end

            if unit == "mm"
                halostericGmsl = 1e3 * halostericGmsl;
            end

        end

    end

    % steric2gmsl alone appends native means; never overwrite the native grids
    % with processed data. A grid-cache hit lacking means needs only an append.
    if saveData && isDerived

        if ~cacheHit
            stored = data;
            stored.lon = coord.lon;
            stored.lat = coord.lat;
            stored.dates = coord.dates(:);
            names = fieldnames(means);

            for iMean = 1:numel(names)
                stored.(names{iMean}) = means.(names{iMean});
            end

            save(outPath, '-struct', 'stored', '-v7.3');
        elseif saveMeans
            save(outPath, '-struct', 'means', '-append');
        end

    end

    if isDerived
        steric = selectsteric(data, cacheGridNames, iComp, iDepth);
    else
        steric = selectsteric(data, gridNames, iComp, iDepth);
    end

    if nargout == 0
        timeFmt = 'datetime';
    end

    [steric, stericSigma, dates, lon, lat] = formatoutput( ...
        steric, coord.dates, coord.lon, coord.lat, ...
        timelim, timeFmt, outFmt, unit);

    if nargout > 0
        return
    end

    try
        plotsealeveltseries(dates, steric, lon, lat, ...
            product, unit, 'Global mean steric sea level');
    catch ME
        warning('ULMO:Plotting:PlotFailed', ...
            'Plotting steric sea level change failed with error:\n%s\nSkipping plotting.', ME.message);
    end

end

%% Subfunctions
function schema = nativefields(metadata, nativeDepth)
    schema = stericfieldschema(metadata, nativeDepth);

    if schema.incompleteSplit || schema.incompatibleSplits
        error('ULMO:steric2lonlatt:InvalidSchema', ...
        'A depth split requires explicit full-depth, shallow and deep grids.');
    end

end

function valid = cachecontains(path, names)
    valid = all(ismember(names, who('-file', path)));
end

function processed = componentfields(data, schema)
    % Assign explicit depth names to supplied components only.
    processed = struct();

    for iField = find(schema.available(:))'
        [iComp, iDepth] = ind2sub(size(schema.available), iField);
        processed.(schema.cacheGridNames(iField)) = ...
            selectsteric(data, schema.gridNames, iComp, iDepth);
    end

end

function [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, intpMthd, beQuiet, callChain)
    allStericVars = fieldnames(data);

    if isempty(meshsize) && isempty(lonOrigin)
        return
    end

    if beQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO>%s] Interpolating %d fields spatially, %s\n', ...
            callchaintext(callChain), numel(allStericVars), templine);
    end

    [data, lon, lat] = structfun(@(x) interpspatial( ...
        coord.lon, coord.lat, x, meshsize, lonOrigin, intpMthd, ...
        BeQuiet = true, CallChain = callChain), data, "UniformOutput", false);
    coord.lon = lon.(allStericVars{1});
    coord.lat = lat.(allStericVars{1});

    if beQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end

function steric = selectsteric(data, rasterNames, iComp, iDepth)

    if iComp > 3
        steric = data.(rasterNames(1, iDepth)) - data.(rasterNames(iComp - 2, iDepth));
    else
        steric = data.(rasterNames(iComp, iDepth));
    end

end

function [meshIntp, lonIntp, latIntp] = interpspatial( ...
        lon, lat, mesh, meshSize, lonOrigin, intpMthd, options)

    arguments (Input)
        lon (:, :) double
        lat (:, :) double
        mesh (:, :, :) double
        meshSize double {mustBePositive} = []
        lonOrigin double {mustBeFinite, mustBeReal} = []
        intpMthd (1, :) char = "linear"
        options.BeQuiet (1, 1) logical = false
        options.CallChain cell = {mfilename}
    end

    if isempty(meshSize)
        meshSize = 1/2;
    end

    if isempty(lonOrigin)
        lonOrigin = 180;
    end

    lon = lon(:)';
    lat = lat(:);

    % ogLonOrigin = (min(lon) + max(lon)) / 2;
    lonIntp = (-180:meshSize:180) + lonOrigin;
    latIntp = -90:meshSize:90;
    [lonnIntp, lattIntp] = meshgrid(lonIntp, latIntp);

    if isequal(lonIntp(:), lon(:)) && isequal(latIntp(:), lat(:))
        meshIntp = mesh;
        return
    end

    if options.BeQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO>%s] Interpolating spatially, %s\n', ...
            callchaintext(options.CallChain), templine);
    end

    isTrans = false;

    if size(mesh, 1) == length(lon)
        isTrans = true;
        mesh = permute(mesh, [2, 1, 3]); % (lon, lat, t) -> (lat, lon, t)
    end

    if lon(1) ~= 0

        if lon(end) - lon(1) == 360
            lon = lon(1:end - 1);
            mesh = mesh(:, 1:end - 1, :);
        end

        [lon, idx] = sort(mod(lon, 360));
        mesh = mesh(:, idx, :);
    end

    if lon(end) - lon(1) ~= 360
        lon = [lon(end) - 360, lon, lon(1) + 360];
        mesh = cat(2, mesh(:, end, :), mesh, mesh(:, 1, :));
    end

    [lonn, latt] = meshgrid(lon, lat);

    if size(lonn) ~= size(mesh, 1:2)
        error('ULMO:Interpolation:DimensionMismatch', ...
            'Longitude/latitude mesh size (%d,%d) does not match data size (%d,%d).', ...
            size(lonn), size(mesh, 1:2));
    end

    meshIntp = nan([length(latIntp), length(lonIntp), size(mesh, 3)], "like", mesh);

    parfor iDate = 1:size(meshIntp, 3)
        meshIntp(:, :, iDate) = ...
            interp2(lonn, latt, squeeze(mesh(:, :, iDate)), ...
            mod(lonnIntp, 360), lattIntp, intpMthd);
    end

    if isTrans
        meshIntp = permute(meshIntp, [2, 1, 3]);
    end

    if options.BeQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end

% Format output
function [steric, stericSigma, dates, lon, lat] = formatoutput( ...
        steric, dates, lon, lat, timelim, timeFmt, outFmt, unit)

    if ~isempty(timelim)
        isValidTime = (dates >= timelim(1)) & (dates <= timelim(2));
        steric = steric(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if strcmpi(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    switch outFmt
        case 'meshgrid'
            lon = lon(:)';
            lat = lat(:);

            if size(steric, 1) ~= length(lat)
                steric = permute(steric, [2, 1, 3]);
            end

        case 'ndgrid'
            lon = lon(:);
            lat = lat(:)';

            if size(steric, 1) ~= length(lon)
                steric = permute(steric, [2, 1, 3]);
            end

    end

    if strcmpi(unit, 'mm')
        steric = steric * 1e3;
    end

    % Currently no sigma data available
    stericSigma = zeros(size(steric), 'like', steric);

end
