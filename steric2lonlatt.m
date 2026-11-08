%% STERIC2LONLATT - Loads and interpolates prepared steric sea level grids.
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
%       Positive numeric scalar or duration for fixed steps, 'midmonth' for
%       calendar-month means, or 'GRACE' for CSR RL06 degree-60 solution
%       windows. Both conventions average native samples and return window
%       midpoints. Missing GRACE months use calendar windows adjusted to
%       meet adjacent solutions.
%       Unit: days for numeric values.
%       Default step: [] (native timestamps).
%   meshsize (optional) - Spatial interpolation spacing.
%       A two-element [spacing, longitudeCentre] vector also sets LonOrigin
%       when that option is empty.
%       With LonOrigin alone, the spacing is 0.5 degrees.
%       Unit: degrees.
%       Default spacing: [] (no spatial interpolation unless LonOrigin is set).
%   timelim (optional) - Inclusive output time limits as datetime or datenum.
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
%   Depth (name-value) - full, shallow, or deep, as prepared by the product.
%       Unsupported layers are not reconstructed or aliased to full.
%       Default depth: full.
%   Type (name-value) - Selected steric component.
%       total, thermosteric, halosteric, non-thermosteric (total minus
%       thermosteric), or non-halosteric (total minus halosteric).
%       Default type: total. The former 'default' and 'no halosteric drift'
%       values are no longer accepted. Callers apply correction explicitly.
%   ForceNew (name-value) - Whether to reprocess every supported component grid and
%       recalculate every native and processed GMSL, even without output six.
%       Does not regenerate native grids or other cache configurations.
%       Default flag: false.
%   TruncationBeforeInterpolation (name-value) - Whether to restrict temporal
%       processing to timelim and include the limits in the cache filename.
%       Default flag: false (process the full record before output selection).
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
%   steric - Selected component in an explicitly named epoch cache.
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
%   <component>Sl - All supported component grids for interpolation/epoch caches.
%       Size: as in the input file, with processed spatial/time dimensions.
%       Unit: m.
%   steric - Selected component for epoch-mean caches. Residuals retain
%       paired native samples before averaging.
%       Size: as in the input file, with processed spatial/time dimensions.
%       Unit: m.
%
%   <component>Gmsl - All applicable temporally processed native means when
%       requested or forced, with the names defined by STERIC2GMSL.
%       Size: (N, 1), aligned with saved dates.
%       Unit: m.
%
% Notes
%   Files are resolved through STERICFILEPATH under IFILES/STERIC/product.
%   Ordinary grid-only requests do not calculate means. GMSL requests reuse
%   STERIC2GMSL caches and process all applicable series separately from grids.
%   Full means of split components are sums of the processed layer means.
%   Grids and GMSL use their own finite samples during epoch averaging.
%   They share output dates, not necessarily per-cell temporal support.
%   ForceNew=true explicitly rebuilds the requested configuration and all
%   native means. SaveData=false keeps this calculation entirely in memory.
%   Rebuild other derived/figure caches manually after convention changes.
%   Legacy epoch caches containing halosteric subtraction are incompatible.
%   Calling without outputs attempts to plot the global mean of the grid.
%
% See also
%   steric2gmsl, stericfilepath, stericfieldschema, interptemporal, temporalmeanepochs
%
% Created by
%   2025/07/22, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/08, OpenAI Codex

function [steric, stericSigma, dates, lon, lat, halostericGmsl] = steric2lonlatt( ...
        product, timestep, meshsize, timelim, options)
    %% Initialisation
    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty} = 'EN4c14'
        timestep {mustBeTimeStep} = []
        meshsize = []
        timelim {mustBeTimeRange} = []
        options.LonOrigin {mustBeFinite, mustBeReal} = []
        options.Interpolation {mustBeTextScalar} = 'linear'
        options.TimeFormat (1, :) DateFormat = 'datetime'
        options.OutputFormat (1, :) MeshFormat = 'meshgrid'
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, {'mm', 'm'})} = 'm'
        options.Depth ...
            {mustBeTextScalar, mustBeMember(options.Depth, {'full', 'shallow', 'deep'})} = 'full'
        options.Type ...
            {mustBeTextScalar, mustBeMember(options.Type, {'total', 'thermosteric', 'halosteric', 'non-thermosteric', 'non-halosteric'})} = 'total'
        options.TruncationBeforeInterpolation (1, 1) {mustBeNumericOrLogical} = false
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

    % 'default' and 'no halosteric drift' are deprecated Type values.
    schema = stericfieldschema();
    gridNames = schema.gridNames;
    meanNames = schema.meanNames;
    types = schema.types;
    depths = schema.depths;
    iType = find(types == string(options.Type));
    iDepth = find(depths == string(options.Depth));
    stericVar = char(gridNames(min(iType, 3), iDepth));
    baseStericVar = char(gridNames(1, iDepth));
    % Residual types select total minus a supplied component at the same depth.
    % They are not aliases for the other named physical component.
    if iType > 3
        stericVar = char(gridNames(iType - 2, iDepth));
        vars = {stericVar, baseStericVar};
    else
        vars = {stericVar};
    end

    spatialIntpMthd = intpMthd;
    if strcmpi(spatialIntpMthd, 'mean'), spatialIntpMthd = 'linear'; end
    % Epoch caches carry the selected expression in their filename. Keeping
    % the total-component suffix excludes legacy files that already contain
    % a global-mean subtraction.
    component = '';
    if useMean
        component = sprintf('-%s-%s', options.Depth, options.Type);
    end
    outPath = stericfilepath(product, timestep, meshsize, lonOrigin, intpMthd, ...
        options.TruncationBeforeInterpolation, timelim, component);
    nativePath = stericfilepath(product);
    coordVars = {'lon', 'lat', 'dates'};
    % Ordinary grid-only reads stay lazy. An explicit rebuild refreshes means
    % even when the caller does not request the diagnostic output.
    needMeans = nargout >= 6 || options.ForceNew;
    halostericGmsl = zeros(0, 1);

    % Interpret a selected epoch field only on its explicit component path.
    cacheHit = ~options.ForceNew && isfile(outPath) && ...
        supportedfields(outPath, coordVars);
    selectedCache = cacheHit && useMean && supportedfields(outPath, {'steric'});
    cacheHit = cacheHit && (selectedCache || supportedfields(outPath, vars));
    if cacheHit
        inPath = outPath;
        if selectedCache, dataVars = {'steric'}; else, dataVars = vars; end
    else
        inPath = nativePath;
        if ~isfile(inPath)
            error('ULMO:LoadData:FileNotFound', 'Input file %s not found.', inPath);
        end
        for iVar = 1:numel(vars)
            if ~supportedfields(inPath, vars(iVar))
                error('ULMO:steric2lonlatt:UnsupportedField', ...
                    'Product %s does not supply %s.', char(product), vars{iVar});
            end
        end
        metadata = whos('-file', inPath);
        schema = nativefields(metadata);
        % A fresh derived cache stores every supported source component, not
        % only the selected output. Explicit names exclude cached GMSL vectors.
        dataVars = cellstr(gridNames(schema.available));
    end
    if options.BeQuiet == 0 || (cacheHit && options.BeQuiet <= 1)
        fprintf('[ULMO>%s] Loading %s steric product from %s.\n', ...
            callchaintext(callChain), product, inPath);
    end
    coord = load(inPath, coordVars{:});
    if isnumeric(coord.dates), coord.dates = datetime(coord.dates, 'ConvertFrom', 'datenum'); end
    coord.dates = coord.dates(:);
    data = load(inPath, dataVars{:});
    % Without early truncation, cache the full processed record and apply
    % timelim only to returned values. Both grids and means use this range.
    epochRange = [];
    if options.TruncationBeforeInterpolation, epochRange = timelim; end

    if ~cacheHit
        if useMean
            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
            % Preserve paired native support for residual field expressions:
            % mean(total - halo) can differ from mean(total) - mean(halo)
            % when their finite cells differ. This is component selection,
            % separate from the caller's later global-mean subtraction.
            data.steric = selectsteric(data, stericVar, baseStericVar, options.Type);
        end
        [data, coord.dates] = processfieldtime(data, coord.dates, timestep, intpMthd, epochRange);
        if ~useMean
            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
        end
    end

    means = struct();
    saveMeans = false;
    if needMeans
        if ~isfile(nativePath)
            error('ULMO:LoadData:FileNotFound', 'Native means require %s.', nativePath);
        end
        schema = nativefields(whos('-file', nativePath));
        supported = schema.supported;
        if nargout >= 6 && ~supported(3, iDepth)
            error('ULMO:steric2lonlatt:UnsupportedField', ...
                'Product %s has no halosteric GMSL at Depth=%s.', char(product), options.Depth);
        end
        names = cellstr(meanNames(supported));
        % Derived diagnostics can be reused as a complete set. Native-file
        % reuse still goes through steric2gmsl, which owns its cache checks.
        meansHit = cacheHit && ~strcmp(outPath, nativePath) && ...
            validmeans(outPath, names, coord.dates);
        if meansHit
            means = load(outPath, names{:});
        else
            % The first saved forced call refreshes the complete native set.
            % Subsequent saved calls reuse that complete set. With saving
            % disabled, every forced call must bypass old native means; repeated
            % calculation is intentional and keeps the native API unchanged.
            firstMean = true;
            for iMean = find(supported(:))'
                [meanType, meanDepth] = ind2sub(size(supported), iMean);
                forceMean = options.ForceNew && (firstMean || ~saveData);
                [values, nativeDates] = steric2gmsl(product, types(meanType), depths(meanDepth), ...
                    ForceNew=forceMean, SaveData=saveData, ...
                    BeQuiet=double(options.BeQuiet)/2);
                % Temporal helpers expect time on axis three. This temporary
                % shape does not make the scalar series a spatial grid.
                means.(meanNames(iMean)) = reshape(values, 1, 1, []);
                firstMean = false;
            end
            [means, meanDates] = processfieldtime(means, nativeDates, timestep, intpMthd, epochRange);
            assert(isequal(meanDates(:), coord.dates(:)), ...
                'ULMO:steric2lonlatt:MeanDates', 'GMSL and grid dates must agree.');
            for iMean = 1:numel(names)
                means.(names{iMean}) = reshape(means.(names{iMean}), [], 1);
            end
            % Nonlinear interpolation need not commute with layer addition.
            % Define the processed full mean from the processed layers so the
            % returned depths remain additive under methods such as pchip.
            for iMeanType = 1:numel(types)
                if all(supported(iMeanType, :))
                    means.(meanNames(iMeanType, 1)) = ...
                        means.(meanNames(iMeanType, 2)) + means.(meanNames(iMeanType, 3));
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
            if strcmpi(unit, 'mm'), halostericGmsl = 1e3 * halostericGmsl; end
        end
    end

    % steric2gmsl alone appends native means; never overwrite the native grids
    % with processed data. A grid-cache hit lacking means needs only an append.
    if saveData && ~strcmp(outPath, nativePath)
        if ~cacheHit
            stored = data;
            stored.lon = coord.lon;
            stored.lat = coord.lat;
            stored.dates = coord.dates(:);
            names = fieldnames(means);
            for iMean = 1:numel(names), stored.(names{iMean}) = means.(names{iMean}); end
            save(outPath, '-struct', 'stored', '-v7.3');
        elseif saveMeans
            save(outPath, '-struct', 'means', '-append');
        end
    end

    if useMean && isfield(data, 'steric')
        steric = data.steric;
    else
        steric = selectsteric(data, stericVar, baseStericVar, options.Type);
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
        warning('Slepian:Plotting:PlotFailed', ...
            'Plotting steric sea level change failed with error:\n%s\nSkipping plotting.', ME.message);
    end

end

%% Subfunctions
function schema = nativefields(metadata)
    schema = stericfieldschema(metadata);
    if schema.incompleteSplit || schema.incompatibleSplits
        error('ULMO:steric2lonlatt:InvalidSchema', ...
            'A depth split requires explicit full, shallow and deep grids.');
    end
end

function valid = supportedfields(path, names)
    metadata = whos('-file', path);
    valid = true;
    for iName = 1:numel(names)
        item = metadata(strcmp({metadata.name}, names{iName}));
        if isempty(item) || any(item.size == 0)
            valid = false;
            return
        end
    end
end

function valid = validmeans(path, names, dates)
    % Validate storage completeness, not source freshness. NaN is a valid
    % no-observation result; Inf and malformed vectors require recalculation.
    metadata = whos('-file', path);
    for iName = 1:numel(names)
        item = metadata(strcmp({metadata.name}, names{iName}));
        if isempty(item) || ~ismember(item.class, {'double', 'single'}) || ...
                item.complex || item.sparse || ~isequal(item.size, [numel(dates), 1])
            valid = false;
            return
        end
    end
    values = load(path, names{:});
    valid = all(structfun(@(x) ~any(isinf(x)), values));
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

function steric = selectsteric(data, stericVar, baseStericVar, type)

    if ismember(type, ["non-thermosteric", "non-halosteric"])
        steric = data.(baseStericVar) - data.(stericVar);
    else
        steric = data.(stericVar);
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
