%% STERIC2LONLATT - Loads and interpolates prepared steric sea level grids.
%
% Syntax
%   [steric, stericSigma, dates] = STERIC2LONLATT(product)
%   [steric, stericSigma, dates] = STERIC2LONLATT(product, timestep, meshsize)
%   [steric, stericSigma, dates] = STERIC2LONLATT(product, timestep, meshsize, timelim)
%   [steric, stericSigma, dates] = STERIC2LONLATT(__, "Name", value)
%   [steric, stericSigma, dates, lon, lat] = STERIC2LONLATT(__)
%
% Input arguments
%   product (optional) - Prepared product: SIO, EN4c14, CMEMS, NOAA, ORAS5,
%       JAMSTEC, ISAS, IPRC, or IAPv4.
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
%       Unavailable native fields raise UnsupportedField.
%       Default depth: full.
%   Type (name-value) - Selected steric component or correction.
%       total, thermosteric, halosteric, non-thermosteric (total minus
%       thermosteric), or non-halosteric (total minus halosteric).
%       'no halosteric drift' subtracts the global halosteric mean for the
%       selected depth from total steric sea level using native offsets. 'default' selects total
%       for SIO and 'no halosteric drift' for every other product.
%       Full and shallow offsets use shared native support. Deep is their
%       difference. Fixed-step interpolation acts separately on fields and
%       offsets. Epoch means average corrected native fields on valid cells.
%       Rebuild older corrected epoch caches with ForceNew.
%       Default type: default.
%   ForceNew (name-value) - Whether to bypass the derived cache and reprocess
%       the existing native grids. Does not regenerate native source data.
%       Default flag: false.
%   TruncationBeforeInterpolation (name-value) - Whether to restrict temporal
%       processing to timelim and include the limits in the cache filename.
%       Default flag: false (process the full record before output selection).
%   SaveData (name-value) - Whether to save processed grids to a derived cache.
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
%       prefixes when supplied, e.g. shallowHalostericSl. Residuals and drift
%       correction also require the corresponding total steric grid.
%       Size: (nLat, nLon, N) or (nLon, nLat, N).
%       Unit: m.
%   steric - Selected component in an existing epoch-mean cache.
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
% Output file content
%   lon, lat, dates - Processed coordinates and timestamps.
%       Same conventions as the input file; saved before output formatting.
%   <component>Sl - Processed component grids for interpolation caches.
%       Size: as in the input file, with processed spatial/time dimensions.
%       Unit: m.
%   steric - Selected component for epoch-mean caches.
%       Size: as in the input file, with processed spatial/time dimensions.
%       Unit: m.
%
% Notes
%   Files are resolved through STERICFILEPATH under IFILES/STERIC/product.
%   Calling without outputs attempts to plot the global mean of the grid.
%
% See also
%   steric2gmsl, stericfilepath, interptemporal, temporalmeanepochs
%
% Created by
%   2025/07/22, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function [steric, stericSigma, dates, lon, lat] = steric2lonlatt( ...
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
            {mustBeTextScalar, mustBeMember(options.Type, {'default', 'total', 'thermosteric', 'halosteric', 'non-thermosteric', 'non-halosteric', 'no halosteric drift'})} = 'default'
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

    if ~isempty(timestep) && isnumeric(timestep)
        timestep = days(timestep);
    end

    stericVar = "sl";

    if strcmpi(options.Type, "default")

        if strcmpi(product, "SIO")
            options.Type = "total";
        else
            options.Type = "no halosteric drift";
        end

    end

    switch lower(options.Type)
        case {'thermosteric', 'non-thermosteric'}
            stericVar = ["thermosteric", stericVar];
        case {'halosteric', 'non-halosteric', 'no halosteric drift'}
            stericVar = ["halosteric", stericVar];
        otherwise
            stericVar = ["steric", stericVar];
    end

    switch options.Depth
        case 'shallow'
            stericVar = ["shallow", stericVar];
        case 'deep'
            stericVar = ["deep", stericVar];
    end

    stericVar = char(strjoin([stericVar(1), arrayfun(@capitalise, stericVar(2:end), "UniformOutput", false)], ''));

    %% Check for existing file
    useMean = ~isempty(timestep) && (ischar(timestep) || isstring(timestep) || strcmpi(intpMthd, 'mean'));
    spatialIntpMthd = intpMthd;

    if strcmpi(spatialIntpMthd, 'mean')
        spatialIntpMthd = 'linear';
    end

    % Epoch means include the selected component and any native-time drift correction.
    comp = '';

    if useMean
        comp = sprintf('-%s-%s', options.Depth, strrep(options.Type, ' ', '-'));
    end

    outPath = stericfilepath(product, timestep, meshsize, lonOrigin, intpMthd, ...
        options.TruncationBeforeInterpolation, timelim, comp);

    baseStericVar = stericVar;

    if ismember(options.Type, ["non-thermosteric", "non-halosteric", "no halosteric drift"])
        baseStericVar = replace(stericVar, ["Thermosteric", "Halosteric"], "Steric");
        baseStericVar = char(replace(baseStericVar, ["thermosteric", "halosteric"], "steric"));
    end

    vars = unique({stericVar, baseStericVar});
    if strcmpi(options.Type, 'no halosteric drift')
        % The native helper can derive a missing halo layer from its counterpart.
        vars = {baseStericVar};
    end
    nativePath = stericfilepath(product);
    if isfile(nativePath)
        sourceVars = unique([{'lon', 'lat', 'dates', 'stericSl'}, vars]);
        for iVar = 1:numel(sourceVars)
            if ~supportedfields(nativePath, sourceVars(iVar))
                error('ULMO:steric2lonlatt:UnsupportedField', ...
                    'Product %s does not support field %s (Depth=%s, Type=%s).', ...
                    char(product), sourceVars{iVar}, options.Depth, options.Type);
            end
        end
    end
    coordVars = {'lon', 'lat', 'dates'};
    cacheVars = vars;

    if useMean
        cacheVars = {'steric'};
    end

    cacheHit = ~options.ForceNew && isfile(outPath) && ...
        supportedfields(outPath, [coordVars, cacheVars]);

    if cacheHit
        inPath = outPath;
        dataVars = cacheVars;
    else
        inPath = nativePath;

        if ~isfile(inPath)
            error('ULMO:LoadData:FileNotFound', 'Input file %s not found.', inPath);
        end

        available = who('-file', inPath);
        missing = setdiff([coordVars, vars], available);

        if ~isempty(missing)
            error('ULMO:LoadData:VariableNotFound', ...
                'Variables missing from %s: %s', inPath, strjoin(missing, ', '));
        end

        if useMean
            dataVars = vars;
        else
            dataVars = available(contains(lower(available), 'stericsl'));
        end

    end

    if options.BeQuiet == 0 || (cacheHit && options.BeQuiet <= 1)
        fprintf('[ULMO>%s] Loading %s steric product from %s.\n', ...
            callchaintext(callChain), product, inPath);
    end

    coord = load(inPath, coordVars{:});
    data = load(inPath, dataVars{:});
    emptyVars = dataVars(cellfun(@(name) isempty(data.(name)), dataVars));
    data = rmfield(data, emptyVars);
    supportedVars = fieldnames(data);

    if ~cacheHit
        epochRange = [];

        if options.TruncationBeforeInterpolation
            epochRange = timelim;
        end

        if useMean
            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
            offset = [];
            if strcmpi(options.Type, 'no halosteric drift')
                offset = nativehalocorrection(product, nativePath, coord.dates, [], timelim, options);
            end
            steric = selectsteric(data, stericVar, baseStericVar, options.Type, offset);
            epochs = temporalmeanepochs(coord.dates, timestep, epochRange);
            data = struct();
            [data.steric, coord.dates] = averageepochs(coord.dates, steric, epochs);
        else

            if ~isempty(timestep)

                if options.BeQuiet == 0
                    fprintf('[ULMO>%s] Interpolating %d fields temporally.\n', ...
                        callchaintext(callChain), numel(dataVars));
                end

                [data, dates] = structfun(@(x) interptemporal(coord.dates, x, timestep, intpMthd, ...
                    TimeRange = epochRange, BeQuiet = true, CallChain = callChain), data, "UniformOutput", false);
                coord.dates = dates.(supportedVars{1});
            end

            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
        end

        % Empty source fields remain unsupported in interpolated field caches.
        if ~useMean
            for iVar = 1:numel(emptyVars)
                data.(emptyVars{iVar}) = [];
            end
        end

        if saveData && ~strcmp(inPath, outPath)
            data.lon = coord.lon;
            data.lat = coord.lat;
            data.dates = coord.dates;
            save(outPath, '-struct', 'data', '-v7.3');

            if options.BeQuiet <= 1
                fprintf('[ULMO>%s] Saved %s.\n', callchaintext(callChain), filehref(outPath, 'steric product'));
            end

        end

    end

    if useMean
        steric = data.steric;
    else
        offset = [];
        if strcmpi(options.Type, 'no halosteric drift')
            offset = nativehalocorrection(product, nativePath, coord.dates, timestep, timelim, options);
        end
        steric = selectsteric(data, stericVar, baseStericVar, options.Type, offset);
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
% Apply the same native correction to fresh and cached uncorrected fields.
function offset = nativehalocorrection(product, inPath, dates, timestep, timelim, options)
    info = stericproductinfo(char(product));
    coord = load(inPath, 'dates');
    h = sterichalocorrection(inPath, (1:numel(coord.dates))', info.hasDeep, true);
    % Interpolate full and shallow together, then derive deep. This preserves
    % correction additivity even with nonlinear methods such as pchip.
    offsets = reshape([h.total, h.shallow]', 2, 1, []) / 1000;
    correctionDates = coord.dates;

    if ~isempty(timestep)
        timeRange = [];
        if options.TruncationBeforeInterpolation, timeRange = timelim; end
        [offsets, correctionDates] = interptemporal(coord.dates, offsets, timestep, ...
            lower(options.Interpolation), TimeRange = timeRange, BeQuiet = true);
    end

    [found, indices] = ismember(dates(:), correctionDates(:));
    assert(all(found), 'ULMO:steric2lonlatt:CorrectionDates', ...
        'Native halosteric correction does not cover output dates.');
    switch options.Depth
        case 'full'
            offset = offsets(1, 1, indices);
        case 'shallow'
            offset = offsets(2, 1, indices);
        case 'deep'
            offset = offsets(1, 1, indices) - offsets(2, 1, indices);
    end
    offset = reshape(offset, 1, 1, []);
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

function steric = selectsteric(data, stericVar, baseStericVar, type, offset)

    if ismember(type, ["non-thermosteric", "non-halosteric"])
        steric = data.(baseStericVar) - data.(stericVar);
    elseif strcmpi(type, 'no halosteric drift')
        steric = data.(baseStericVar) - offset;
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

% Check whether a specific field is saved in the .mat file
function tf = supportedfields(path, names)
    info = whos('-file', path);
    tf = true;

    for iName = 1:numel(names)
        item = info(strcmp({info.name}, names{iName}));

        if isempty(item) || any(item.size == 0)
            tf = false;
            return
        end

    end

end
