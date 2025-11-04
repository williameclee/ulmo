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
%       - Character vector, 'midmonth' (at the middle of each month).
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
%       - A string of valid method for INTERP1.
%       The default intepolation method is 'linear'. This option is only
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
%	SaveData - Logical flag to save the data
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
% Authored by
%	2025/07/22, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%	2025/11/03, williameclee@arizona.edu (@williameclee)

function [steric, stericSigma, dates, lon, lat] = steric2lonlatt(product, timestep, meshsize, timelim, options)
    %% Initialisation
    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty} = 'EN4c14'
        timestep {mustBeTimeStep} = []
        meshsize = []
        timelim {mustBeTimeRange} = []
        options.LonOrigin {mustBeFinite, mustBeReal} = []
        options.Interpolation (1, :) char ...
            {mustBeMember(options.Interpolation, {'linear', 'nearest', 'next', 'previous', 'spline', 'pchip'})} = 'linear'
        options.TimeFormat (1, :) DateFormat = 'datetime'
        options.OutputFormat (1, :) MeshFormat = 'meshgrid'
        options.Unit (1, :) char {mustBeMember(options.Unit, {'mm', 'm'})} = 'm'
        options.Depth (1, :) char ...
            {mustBeMember(options.Depth, {'full', 'shallow', 'deep'})} = 'full'
        options.Type (1, :) char ...
            {mustBeMember(options.Type, {'total', 'thermosteric', 'halosteric', 'non-thermosteric', 'non-halosteric', 'no halosteric drift'})} = 'totalÏ'
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

    switch options.Type
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
    outputPath = ...
        outputpath(product, timestep, meshsize, lonOrigin, intpMthd);

    vars = {'lon', 'lat', 'dates', 'stericSl', stericVar};

    if ~options.ForceNew && exist(outputPath, 'file') && ...
            all(ismember(vars, who('-file', outputPath)))

        if options.BeQuiet <= 1
            t = tic;
            msg = 'this may take a while...';
            fprintf('[ULMO>%s] Loading <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s steric product</a>, %s\n', ...
                callchaintext(callChain), outputPath, outputPath, product, msg);
        end

        data = load(outputPath, vars{:});

        if options.BeQuiet <= 1
            fprintf(repmat('\b', 1, length(msg) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

        
        if strcmpi(options.Type, "non-thermosteric")
            baseStericVar = replace(replace(stericVar, "Thermosteric", "Steric"), "thermosteric", "steric");
            steric = data.(baseStericVar) - data.(stericVar);
        elseif strcmpi(options.Type, "non-halosteric")
            baseStericVar = replace(replace(stericVar, "Halosteric", "Steric"), "halosteric", "steric");
            steric = data.(baseStericVar) - data.(stericVar);
        elseif strcmpi(options.Type, 'no halosteric drift')
            gridWeight = cosd(data.lat(:));
            tol = length(data.dates) * 0.8;
            lsMask = (sum(~isnan(data.(stericVar)), 3) >= tol) & (sum(data.(stericVar) ~= 0, 3) >= tol);
            gmSteric = sum(data.(stericVar) .* lsMask .* gridWeight, [1, 2], 'omitmissing') ...
                ./ sum(~isnan(data.(stericVar)) .* lsMask .* gridWeight, [1, 2], 'omitmissing');
            baseStericVar = replace(replace(stericVar, "Halosteric", "Steric"), "halosteric", "steric");
            steric = data.(baseStericVar) - gmSteric;
        else
            steric = data.(stericVar);
        end

        [steric, stericSigma, dates, lon, lat] = formatoutput( ...
            steric, data.dates, data.lon, data.lat, ...
            timelim, timeFmt, outputFmt, unit);

        if nargout > 0
            return
        end

        plotsealeveltseries(dates, steric, lon, lat, ...
            product, unit, 'Global mean steric sea level');

        return
    elseif ~options.ForceNew && exist(outputPath, 'file') && ...
            ~all(ismember(vars, who('-file', outputPath)))
        notFoundVars = vars(~ismember(vars, who('-file', outputPath)));
        warning('ULMO:LoadData:VariableNotFound', ...
            'Some Variables not found in %s, reprocessing required:\n''%s''', outputPath, strjoin(notFoundVars, ''', '''));
    end

    %% Main
    inputPath = outputpath(product);

    if ~exist(inputPath, 'file')
        error('ULMO:LoadData:FileNotFound', ...
            'Input file %s not found.', inputPath);
    end

    if any(~ismember(vars, who('-file', inputPath)))
        notFoundVars = vars(~ismember(vars, who('-file', inputPath)));
        error('ULMO:LoadData:VariableNotFound', ...
            'Some Variables not found in %s:\n''%s''', inputPath, strjoin(notFoundVars, ''', '''));
    end

    allStericVars = who('-file', inputPath);
    allStericVars = allStericVars(contains(lower(allStericVars), 'stericsl'));

    coord = load(inputPath, 'lon', 'lat', 'dates');
    data = load(inputPath, allStericVars{:});

    if options.BeQuiet == 0
        fprintf('[ULMO>%s] Loaded uninterpolated <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s steric product</a>.\n', ...
            callchaintext(callChain), inputPath, inputPath, product);
    end

    if ~isempty(timestep)

        if options.BeQuiet == 0
            t = tic;
            templine = 'this may take a while...';
            fprintf('[ULMO>%s] Interpolating %d fields temporally, %s\n', ...
                callchaintext(callChain), numel(fields(data)), templine);
        end

        [data, dates] = structfun(@(x) interptemporal( ...
            coord.dates, x, timestep, intpMthd, ...
            BeQuiet = true, CallChain = callChain), data, "UniformOutput", false);
        coord.dates = dates.(allStericVars{1});

        if options.BeQuiet == 0
            fprintf(repmat('\b', 1, length(templine) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

    end

    if ~(isempty(meshsize) && isempty(lonOrigin))

        if options.BeQuiet == 0
            t = tic;
            templine = 'this may take a while...';
            fprintf('[ULMO>%s] Interpolating %d fields spatially, %s\n', ...
                callchaintext(callChain), numel(fields(data)), templine);
        end

        [data, lon, lat] = structfun(@(x) interpspatial( ...
            coord.lon, coord.lat, x, meshsize, lonOrigin, intpMthd, ...
            BeQuiet = true, CallChain = callChain), data, "UniformOutput", false);
        coord.lon = lon.(allStericVars{1});
        coord.lat = lat.(allStericVars{1});

        if options.BeQuiet == 0
            fprintf(repmat('\b', 1, length(templine) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

    end

    if saveData

        if options.BeQuiet <= 1
            t = tic;
            msg = 'this may take a while...';
            fprintf('[ULMO>%s] Saving <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s steric product</a>, %s\n', ...
                callchaintext(callChain), outputPath, outputPath, product, msg);
        end

        try
            save(outputPath, '-struct', 'coord', 'lon', 'lat', 'dates', '-append');
            save(outputPath, '-struct', 'data', allStericVars{:}, '-append');
        catch ME
            warning('ULMO:SaveData:SaveFailed', ...
                'Saving to %s failed with error:\n%s\nTrying to save afresh.', ...
                outputPath, ME.message);
            save(outputPath, '-struct', 'coord', 'lon', 'lat', 'dates', '-v7.3');
            save(outputPath, '-struct', 'data', allStericVars{:}, '-append');
        end

        if options.BeQuiet <= 1
            fprintf(repmat('\b', 1, length(msg) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

    end

    if strcmpi(options.Type, "non-thermosteric")
        baseStericVar = replace(replace(stericVar, "Thermosteric", "Steric"), "thermosteric", "steric");
        steric = data.(baseStericVar) - data.(stericVar);
    elseif strcmpi(options.Type, "non-halosteric")
        baseStericVar = replace(replace(stericVar, "Halosteric", "Steric"), "halosteric", "steric");
        steric = data.(baseStericVar) - data.(stericVar);
    elseif strcmpi(options.Type, 'no halosteric drift')
        gridWeight = cosd(coord.lat(:));
        tol = length(coord.dates) * 0.8;
        lsMask = (sum(~isnan(data.(stericVar)), 3) >= tol) & (sum(data.(stericVar) ~= 0, 3) >= tol);
        gmSteric = sum(data.(stericVar) .* lsMask .* gridWeight, [1, 2], 'omitmissing') ...
            ./ sum(~isnan(data.(stericVar)) .* lsMask .* gridWeight, [1, 2], 'omitmissing');
        baseStericVar = replace(replace(stericVar, "Halosteric", "Steric"), "halosteric", "steric");
        steric = data.(baseStericVar) - gmSteric;
    else
        steric = data.(stericVar);
    end

    [steric, stericSigma, dates, lon, lat] = formatoutput( ...
        steric, coord.dates, coord.lon, coord.lat, ...
        timelim, timeFmt, outputFmt, unit);

    if nargout > 0
        return
    end

    plotsealeveltseries(dates, steric, lon, lat, ...
        product, unit, 'Global mean steric sea level');

end

%% Subfunctions
% Interpolation
% function [meshIntp, datesIntp] = ...
%         interptemporal(dates, mesh, timeStep, intpMthd, beQuiet, callChain)

%     if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
%         datesIntp = midmonth([dates(1), dates(end)]);
%     else

%         if mean(diff(dates)) > timeStep
%             warning(sprintf('ULMO:%s:InterpolationStepTooSmall', upper(mfilename)), ...
%                 'The interpolation time step (%s) is smaller than the mean data resolution (%s).', timeStep, mean(diff(dates)));
%         end

%         datesIntp = dates(1):timeStep:dates(end);
%     end

%     if isequal(datesIntp, dates)
%         meshIntp = mesh;
%         return
%     end

%     if beQuiet == 0
%         t = tic;
%         templine = 'this may take a while...';
%         fprintf('[ULMO>%s] Interpolating temporally, %s\n', ...
%             callchaintext(callChain), templine);
%     end

%     meshFlat = reshape(mesh, [], size(mesh, 3))';
%     meshIntp = interp1(dates, meshFlat, datesIntp, intpMthd)';
%     meshIntp = reshape(meshIntp, [size(mesh, 1:2), length(datesIntp)]);

%     if beQuiet == 0
%         fprintf(repmat('\b', 1, length(templine) + 1));
%         fprintf('took %.1f seconds.\n', toc(t));
%     end

% end

function [meshIntp, lonIntp, latIntp] = ...
        interpspatial(lon, lat, mesh, meshSize, lonOrigin, intpMthd, options)

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
        disp(size(lonn));
        disp(size(mesh));
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
function [steric, stericSigma, dates, lon, lat] = ...
        formatoutput(steric, dates, lon, lat, timelim, timeFmt, outputFmt, unit)

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

% Find output path
function outputPath = ...
        outputpath(product, timeStep, meshSize, lonOrigin, intpMthd)

    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty}
        timeStep {mustBeTimeStep} = []
        meshSize {mustBePositive} = []
        lonOrigin {mustBeFinite, mustBeReal} = []
        intpMthd (1, :) char = []
    end

    outputFolder = fullfile(getenv("IFILES"), 'STERIC', char(product));

    timeStepStr = '';

    if ~isempty(timeStep)

        if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
            timeStepStr = '-Tmidmonth';
        elseif isduration(timeStep)
            timeStepStr = sprintf('-T%s', erase(sprintf('%s', timeStep), ' '));
        elseif isnumeric(timeStep)
            timeStepStr = sprintf('-T%d', timeStep);
        else
            error(sprintf('ULMO:%s:InvalidTimeStep', upper(mfilename)), ...
                'Unrecognised time step input for saving: %s', class(timeStep));
        end

    end

    spaceIntpStr = '';

    if ~isempty(meshSize) && ~isempty(lonOrigin)
        spaceIntpStr = ['-H', num2str(meshSize), 'O', num2str(lonOrigin)];
    elseif ~isempty(meshSize)
        spaceIntpStr = ['-H', num2str(meshSize)];
    elseif ~isempty(lonOrigin)
        spaceIntpStr = ['-O', num2str(lonOrigin)];
    end

    intpMethodStr = '';

    if ~(isempty(timeStepStr) && isempty(spaceIntpStr))
        intpMethodStr = ['-', lower(intpMthd)];
    end

    outputFile = sprintf('%s-StericSeaLevel%s%s%s.mat', ...
        product, timeStepStr, spaceIntpStr, intpMethodStr);

    outputPath = fullfile(outputFolder, outputFile);
end
