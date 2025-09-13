%% STERIC2LONLATT
% Loads steric sea level data and interpolate to regular latitude-longitude mesh.
%
% Syntax
%   [mesh, dates, lon, lat] = STERIC2LONLATT(product)
%   [mesh, dates, lon, lat] = STERIC2LONLATT(product, timestep, meshsize)
%   [mesh, dates, lon, lat] = STERIC2LONLATT(product, timestep, meshsize, timelim)
%   [mesh, dates, lon, lat] = STERIC2LONLATT(product, timestep, meshsize, timelim)
%   [mesh, dates, lon, lat] = STERIC2LONLATT(product, "Name", Value)
%   STERIC2LONLATT(___)
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
%       - Numeric or duration scalar, in units of days
%       - Character vector, 'midmonth'
%       When not specified, no temporal interpolation is performed.
%   meshsize - Spatial interpolation grid size in degrees
%       - Real scalar, in units of degrees
%       When not specified, no spatial interpolation is performed.
%   timelim - Time range of interest
%       - Datetime vector
%       - Numeric vector, in units of datenum
%       When not specified, the full time range of the data is returned.
%   LonOrigin - Longitude origin of the output mesh
%       - Real scalar, in units of degrees, e.g. 0, 180, -180
%       When not specified, the midpoint of the input longitude range is 
%       used.
%   Interpolation - Interpolation method
%       - Character vector, a valid method for INTERP1
%       The default is 'linear'.
%   OutputFormat - Format of the output mesh
%       - 'meshgrid' : Meshgrid format, lon and lat are row and column 
%           vectors
%       - 'ndgrid' : Ndgrid format, lon and lat are column and row vectors
%       The default is 'meshgrid'.
%   TimeFormat - Format of the output time vector
%       - 'datetime' : Datetime format
%       - 'datenum' : Datenum format
%       The default is 'datetime'.
%   Unit - Unit of the output steric sea level data
%       - 'm' : metres
%       - 'mm' : millimetres
%       The default is 'm'.
%   Depth - Depth range of the steric sea level data
%       - 'full' : Full depth steric sea level
%       - 'shallow' : Upper ocean steric sea level
%       - 'deep' : Deep ocean steric sea level
%       The default is 'full'.
%	ForceNew - Logical flag to force reprocess of the data
%		- true: Reprocess the data
%		- false: Only reprocess if previous output is not found
%		The default option is 'soft true', only reprocessing the data if
%       the output file is older than the input files. This option is not
%       explicitly exposed; to enforce this option, set it to 0.5.
%	SaveData - Logical flag to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%		The default option is true.
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is 'soft true', only printing important 
%       messages.
%
% Authored by
%	2025/07/22, williameclee@arizona.edu (@williameclee)
% Last modified by
%	2025/09/12, williameclee@arizona.edu (@williameclee)

function [steric, dates, lon, lat] = steric2lonlatt(product, timestep, meshsize, timelim, options)
    %% Initialisation
    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty} = 'EN4'
        timestep {mustBeTimeStep} = []
        meshsize {mustBePositive} = []
        timelim {mustBeTimeRange} = []
        options.LonOrigin {mustBeFinite, mustBeReal} = []
        options.Interpolation (1, :) char ...
            {mustBeMember(options.Interpolation, {'linear', 'nearest', 'next', 'previous', 'spline', 'pchip'})} = 'linear'
        options.TimeFormat (1, :) DateFormat = 'datetime'
        options.OutputFormat (1, :) MeshFormat = 'meshgrid'
        options.Unit (1, :) char {mustBeMember(options.Unit, {'mm', 'm'})} = 'm'
        options.Depth (1, :) char ...
            {mustBeMember(options.Depth, {'full', 'shallow', 'deep'})} = 'full'
        options.ForceNew (1, 1) {mustBeNumericOrLogical} = false
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
        options.SaveData (1, 1) {mustBeNumericOrLogical} = true
    end

    arguments (Output)
        steric (:, :, :) {mustBeReal}
        dates (:, 1)
        lon (:, :) {mustBeReal, mustBeVector}
        lat (:, :) {mustBeReal, mustBeVector}
    end

    lonOrigin = options.LonOrigin;
    intpMthd = lower(options.Interpolation);
    timeFmt = lower(options.TimeFormat);
    outputFmt = lower(options.OutputFormat);
    unit = options.Unit;
    depth = options.Depth;
    forceNew = options.ForceNew;
    beQuiet = uint8(double(options.BeQuiet) * 2);
    saveData = options.SaveData;

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    if ~isempty(timestep) && isnumeric(timestep)
        timestep = days(timestep);
    end

    switch depth
        case 'shallow'

            if ismember(product, {'SIO', 'NOAA'})
                stericVar = 'stericSl';
            else
                stericVar = 'shallowStericSl';
            end

        case 'deep'
            stericVar = 'deepStericSl';
        otherwise
            stericVar = 'stericSl';
    end

    %% Check for existing file
    outputPath = ...
        outputpath(product, timestep, meshsize, lonOrigin, intpMthd);

    if ~forceNew && exist(outputPath, 'file') && ...
            all(ismember({'lon', 'lat', 'dates', stericVar}, who('-file', outputPath)))
        data = load(outputPath, 'lon', 'lat', 'dates', stericVar);

        if beQuiet <= 1
            fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">steric product</a>.\n', ...
                mfilename("fullpath"), mfilename, outputPath, outputPath);
        end

        [steric, dates, lon, lat] = ...
            formatoutput(data.(stericVar), data.dates, data.lon, data.lat, timelim, timeFmt, outputFmt, unit);

        if nargout == 0
            plotsealeveltseries(dates, steric, lon, lat, product, unit);
        end

        return

    end

    %% Main
    inputPath = outputpath(product, [], [], [], []);

    if ~exist(inputPath, 'file')
        error('ULMO:LoadData:MethodNotImplemented', ...
            'Input file %s not found.', inputPath);
    end

    requiredVars = {'lon', 'lat', 'dates', stericVar};

    if any(~ismember({'lon', 'lat', 'dates', stericVar}, who('-file', inputPath)))
        error('ULMO:LoadData:VariableNotFound', ...
            'Some Variables not found in %s:\n''%s''', inputPath, strjoin(requiredVars, ''', '''));
    end

    data = load(inputPath, 'lon', 'lat', 'dates', stericVar);

    if beQuiet == 0
        fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Loaded uninterpolated <a href="matlab: fprintf(''%s\\n'');open(''%s'')">steric product</a>.\n', ...
            mfilename("fullpath"), mfilename, outputPath, outputPath);
    end

    if ~isempty(timestep)
        [data.(stericVar), data.dates] = ...
            interptemporal(data.dates, data.(stericVar), timestep, intpMthd, beQuiet);
    end

    if ~(isempty(meshsize) && isempty(lonOrigin))

        if isempty(meshsize)
            meshsize = 1/2;
        end

        if isempty(lonOrigin)
            lonOrigin = 180;
        end

        [data.(stericVar), data.lon, data.lat] = ...
            interpspatial(data.lon, data.lat, data.(stericVar), meshsize, lonOrigin, intpMthd, beQuiet);
    end

    if saveData

        try
            save(outputPath, '-struct', 'data', 'lon', 'lat', 'dates', stericVar, '-append');
        catch
            save(outputPath, '-struct', 'data', 'lon', 'lat', 'dates', stericVar, '-v7.3');
        end

        if beQuiet <= 1
            fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">steric product</a>.\n', ...
                mfilename("fullpath"), mfilename, outputPath, outputPath);
        end

    end

    [steric, dates, lon, lat] = ...
        formatoutput(data.(stericVar), data.dates, data.lon, data.lat, timelim, timeFmt, outputFmt, unit);

    if nargout == 0
        plotsealeveltseries(dates, steric, lon, lat, product, unit);
    end

end

%% Subfunctions
% Interpolation
function [meshIntp, datesIntp] = ...
        interptemporal(dates, mesh, timeStep, intpMthd, beQuiet)

    if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
        datesIntp = midmonth([dates(1), dates(end)]);
    else

        if mean(diff(dates)) > timeStep
            warning(sprintf('ULMO:%s:InterpolationStepTooSmall', upper(mfilename)), ...
                'The interpolation time step (%s) is smaller than the mean data resolution (%s).', timeStep, mean(diff(dates)));
        end

        datesIntp = dates(1):timeStep:dates(end);
    end

    if isequal(datesIntp, dates)
        meshIntp = mesh;
        return
    end

    if beQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Interpolating temporally, %s\n', ...
            mfilename("fullpath"), mfilename, templine);
    end

    meshFlat = reshape(mesh, [], size(mesh, 3))';
    meshIntp = interp1(dates, meshFlat, datesIntp, intpMthd)';
    meshIntp = reshape(meshIntp, size(mesh, 1:2), length(datesIntp));

    if beQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end

function [meshIntp, lonIntp, latIntp] = ...
        interpspatial(lon, lat, mesh, meshSize, lonOrigin, intpMthd, beQuiet)

    ogLonOrigin = (min(lon) + max(lon)) / 2;
    lonIntp = (-180:meshSize:180) + lonOrigin;
    latIntp = -90:meshSize:90;

    if isequal(lonIntp(:), lon(:)) && isequal(latIntp(:), lat(:))
        meshIntp = mesh;
        return
    end

    if beQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Interpolating spatially, %s\n', ...
            mfilename("fullpath"), mfilename, templine);
    end

    [lonnIntp, lattIntp] = meshgrid(lonIntp, latIntp);
    lonnIntp(lonnIntp < ogLonOrigin - 180) = lonnIntp(lonnIntp < ogLonOrigin - 180) + 360;
    lonnIntp(lonnIntp > ogLonOrigin + 180) = lonnIntp(lonnIntp > ogLonOrigin + 180) - 360;
    lonPad = [lon(end) - 360; lon(:); lon(1) + 360];
    meshPad = cat(2, mesh(:, end, :), mesh, mesh(:, 1, :));
    [lonn, latt] = meshgrid(lonPad, lat);

    meshIntp = nan([length(latIntp), length(lonIntp), size(mesh, 3)], "like", mesh);

    parfor iDate = 1:size(meshIntp, 3)
        meshIntp(:, :, iDate) = ...
            interp2(lonn, latt, squeeze(meshPad(:, :, iDate)), ...
            lonnIntp, lattIntp, intpMthd);
    end

    if beQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end

% Format output
function varargout = ...
        formatoutput(stericSl, dates, lon, lat, timelim, timeFmt, outputFmt, unit)

    if ~isempty(timelim)
        isValidTime = (dates >= timelim(1)) & (dates <= timelim(2));
        stericSl = stericSl(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if strcmpi(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    switch outputFmt
        case 'meshgrid'
            lon = lon(:)';
            lat = lat(:);

            if size(stericSl, 1) ~= length(lat)
                stericSl = permute(stericSl, [2, 1, 3]);
            end

        case 'ndgrid'
            lon = lon(:);
            lat = lat(:)';

            if size(stericSl, 1) ~= length(lon)
                stericSl = permute(stericSl, [2, 1, 3]);
            end

    end

    if strcmpi(unit, 'mm')
        stericSl = stericSl * 1e3;
    end

    varargout = {stericSl, dates, lon, lat};
end

% Find output path
function outputPath = ...
        outputpath(product, timeStep, meshSize, lonOrigin, intpMthd)
    outputFolder = fullfile(getenv("IFILES"), 'HOMaGE', char(product));

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

function plotsealeveltseries(dates, stericSl, lon, lat, product, unit)
    domain = GeoDomain('oceans', "Buffer", 0.5, "DefaultParams", true);
    domainMask = domainmask(domain, lon, lat, "BeQuiet", true);
    gridWeights = domainMask .* cosd(lat(:));
    gridWeights(~domainMask) = 0;
    gridWeights = gridWeights ./ sum(gridWeights, 'all', 'omitmissing');
    steric = squeeze(sum(stericSl .* gridWeights, [1, 2], 'omitmissing'));
    steric = steric - mean(steric, 'omitmissing');

    figure()
    set(gcf, "Name", 'Global mean steric sea level', 'NumberTitle', 'off');
    clf

    plot(dates, steric, "DisplayName", sprintf('Product: %s', char(product)))
    ylabel(sprintf('Steric sea level [%s]', unit))
    title('Steric sea level time series')
    legend('Location', 'best')

    xlim([min(dates), max(dates)])
end
