%% STERIC2LONLATT
% Loads steric sea level data and interpolate to regular latitude-longitude mesh.
%
% Syntax
%   [steric, sigma, dates] = STERIC2LONLATT(product)
%   [steric, sigma, dates] = STERIC2LONLATT(product, timestep, meshsize)
%   [steric, sigma, dates] = STERIC2LONLATT(product, timestep, meshsize, timelim)
%   [steric, sigma, dates] = STERIC2LONLATT(product, timestep, meshsize, timelim)
%   [steric, sigma, dates] = STERIC2LONLATT(__, "Name", Value)
%   [steric, sigma, dates, lon, lat] = STERIC2LONLATT(__)
%
% Input arguments
%   product - Name of the steric sea level product
%       - 'SIO' : Scripps Argo steric sea level, only include the upper
%           ocean
%       - 'EN4' : Met Office multi-mission EN4.2.2 steric sea level,
%           includes full depth
%       - 'CMEMS' : Copernicus Marine Environment Monitoring Service steric
%           sea level, includes full depth
%       - 'ORAS5' : ECMWF Ocean ReAnalysis System 5 steric sea level,
%           includes full depth
%       - 'NOAA' : NOAA steric sea level
%       The default product is 'EN4'.
%   timestep - Temporal interpolation time step
%       - Numeric or duration scalar, in units of days.
%       - 'midmonth': mean of native samples in each calendar month,
%         dated at the midpoint.
%       - 'GRACE': mean within CSR RL06 degree-60 GRACE solution windows.
%         Missing months use calendar boundaries inside the gap, with outer
%         bounds aligned to adjacent real solutions. Dates are window midpoints.
%         Monthly inputs are averaged by native timestamps, not reconstructed
%         at daily resolution.
%       When not specified, no temporal interpolation is performed.
%   meshsize - Spatial interpolation grid size in degrees
%       - Real scalar, in units of degrees.
%       When not specified, no spatial interpolation is performed.
%   timelim - Time range of interest
%       - Datetime or numeric vector, in units of datenum
%       When not specified, the full time range of the data is returned.
%   LonOrigin - Longitudinal centre of the interpolated output mesh
%       - Real scalar, in units of degrees.
%       When not specified, the midpoint of the input longitude range is
%       used. This option is only relevant when spatial interpolation is
%       performed.
%   Interpolation - Interpolation method, for both temporal and spatial
%       - An INTERP1 method, or 'mean' to average native samples.
%         For 'midmonth' and 'GRACE', this option controls spatial
%         interpolation only. Fixed-step means use
%         half-open windows centered on the output dates. Spatial
%         interpolation uses 'linear' when 'mean' is requested.
%       The default interpolation method is 'linear'. This option is only
%       relevant when temporal or spatial interpolation is performed.
%   OutputFormat - Format of the output mesh
%       - 'meshgrid' or 'ndgrid' for (lat, lon, t) or (lon, lat, t)
%           ordering.
%       The default format is 'meshgrid'.
%   TimeFormat - Format of the output time vector
%       - 'datetime' or 'datenum'.
%       The default format is 'datetime'.
%   Unit - Unit of the output steric sea level data
%       - 'm' (metres) or 'mm' (millimetres, equivalent to kg/m^2).
%       The default unit is 'm'.
%   Depth - Depth range of the steric sea level data
%       - 'full' : Full depth steric sea level
%       - 'shallow' : Upper ocean steric sea level
%       - 'deep' : Deep ocean steric sea level
%       The default is 'full'.
%	ForceNew - Logical flag to force reprocess of the data
%		The default option is false.
%   TruncationBeforeInterpolation - Whether to restrict interpolation to timelim.
%       The default is false, so the cache covers the full product time range.
%	SaveData - Logical flag to cache the spatially and temporally interpolated data
%		The default option is true.
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is 'soft true', only printing important
%       messages.
%
% Output arguments
%	ssh - SSH mesh
%		The dimensions and units depend on the options specified.
%       The default dimensions are (lat, lon, t) in units of metres.
%	sigma - Uncertainty or error mesh of the SSH data
%		The units and dimensions are the same as for SSH.
%	dates - The time stamps of the data
%		The data type depends on the TimeFormat option.
%	lon - Longitude of the mesh
%		The dimensions (i.e. row or column vector) depend on the
%       OutputFormat option.
%	lat - Latitude of the mesh
%		The dimensions (i.e. row or column vector) depend on the
%       OutputFormat option.
%
% Created by
%	2025/07/22, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%	2026/11/04, En-Chi Lee (williameclee@gmail.com)

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
    outputFmt = lower(options.OutputFormat);
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

    if strcmpi(options.Depth, 'shallow') && ismember(product, {'SIO', 'NOAA'})
        options.Depth = 'full';
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

    outputPath = stericfilepath(product, timestep, meshsize, lonOrigin, intpMthd, ...
        options.TruncationBeforeInterpolation, timelim, comp);

    baseStericVar = stericVar;

    if ismember(options.Type, ["non-thermosteric", "non-halosteric", "no halosteric drift"])
        baseStericVar = replace(stericVar, ["Thermosteric", "Halosteric"], "Steric");
        baseStericVar = char(replace(baseStericVar, ["thermosteric", "halosteric"], "steric"));
    end

    vars = unique({stericVar, baseStericVar});
    coordVars = {'lon', 'lat', 'dates'};
    cacheVars = vars;

    if useMean
        cacheVars = {'steric'};
    end

    cacheHit = ~options.ForceNew && isfile(outputPath) && ...
        all(ismember([coordVars, cacheVars], who('-file', outputPath)));

    if cacheHit
        inputPath = outputPath;
        dataVars = cacheVars;
    else
        inputPath = stericfilepath(product);

        if ~isfile(inputPath)
            error('ULMO:LoadData:FileNotFound', 'Input file %s not found.', inputPath);
        end

        available = who('-file', inputPath);
        missing = setdiff([coordVars, vars], available);

        if ~isempty(missing)
            error('ULMO:LoadData:VariableNotFound', ...
                'Variables missing from %s: %s', inputPath, strjoin(missing, ', '));
        end

        if useMean
            dataVars = vars;
        else
            dataVars = available(contains(lower(available), 'stericsl'));
        end

    end

    if options.BeQuiet == 0 || (cacheHit && options.BeQuiet <= 1)
        fprintf('[ULMO>%s] Loading %s steric product from %s.\n', ...
            callchaintext(callChain), product, inputPath);
    end

    coord = load(inputPath, coordVars{:});
    data = load(inputPath, dataVars{:});

    if ~cacheHit
        epochRange = [];

        if options.TruncationBeforeInterpolation
            epochRange = timelim;
        end

        if useMean
            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
            steric = selectsteric(data, coord, stericVar, baseStericVar, options.Type);
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
                coord.dates = dates.(dataVars{1});
            end

            [data, coord] = preparegrid(data, coord, meshsize, lonOrigin, spatialIntpMthd, options.BeQuiet, callChain);
        end

        if saveData && ~strcmp(inputPath, outputPath)
            data.lon = coord.lon;
            data.lat = coord.lat;
            data.dates = coord.dates;
            save(outputPath, '-struct', 'data', '-v7.3');

            if options.BeQuiet <= 1
                fprintf('[ULMO>%s] Saved %s.\n', callchaintext(callChain), filehref(outputPath, 'steric product'));
            end

        end

    end

    if useMean
        steric = data.steric;
    else
        steric = selectsteric(data, coord, stericVar, baseStericVar, options.Type);
    end

    if nargout == 0
        timeFmt = 'datetime';
    end

    [steric, stericSigma, dates, lon, lat] = formatoutput( ...
        steric, coord.dates, coord.lon, coord.lat, ...
        timelim, timeFmt, outputFmt, unit);

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

function steric = selectsteric(data, coord, stericVar, baseStericVar, type)

    if ismember(type, ["non-thermosteric", "non-halosteric"])
        steric = data.(baseStericVar) - data.(stericVar);
    elseif strcmpi(type, 'no halosteric drift')
        gridWeight = cosd(coord.lat(:));
        tol = length(coord.dates) * 0.8;
        lsMask = (sum(~isnan(data.(stericVar)), 3) >= tol) & (sum(data.(stericVar) ~= 0, 3) >= tol);
        gmSteric = sum(data.(stericVar) .* lsMask .* gridWeight, [1, 2], 'omitmissing') ...
            ./ sum(~isnan(data.(stericVar)) .* lsMask .* gridWeight, [1, 2], 'omitmissing');
        steric = data.(baseStericVar) - gmSteric;
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
        steric, dates, lon, lat, timelim, timeFmt, outputFmt, unit)

    if ~isempty(timelim)
        isValidTime = (dates >= timelim(1)) & (dates <= timelim(2));
        steric = steric(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if strcmpi(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    switch outputFmt
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
