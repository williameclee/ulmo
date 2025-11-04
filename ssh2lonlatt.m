%% SSH2LONLATT
% Loads sea surface height (SSH) anomaly data and interpolate to regular
% latitude-longitude mesh.
%
% Syntax
%   [ssh, sigma, dates] = SSH2LONLATT(product)
%   [ssh, sigma, dates] = SSH2LONLATT(product, timestep, meshsize)
%   [ssh, sigma, dates] = SSH2LONLATT(product, timestep, meshsize, timelim)
%   [ssh, sigma, dates] = SSH2LONLATT(product, timestep, meshsize, timelim)
%   [ssh, sigma, dates] = SSH2LONLATT(__, "Name", Value)
%   [ssh, sigma, dates, lon, lat] = SSH2LONLATT(__)
%
% Input arguments
%	product - Name of the SSH product to load
%		- 'MEaSUREs': The MEaSUREs Gridded Sea Surface Height Anomalies
%           Version 2205 product.
%		The default product is 'MEaSUREs', which is the only option
%       available as of right now.
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
% Data sources
%	MEaSUREs Gridded Sea Surface Height Anomalies Version 2205
% 		https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205
%
% See also
%	INTERP1
%
% Authored by
%	2025/05/19, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%	2025/11/03, williameclee@arizona.edu (@williameclee)

function [ssh, sshSigma, dates, lon, lat] = ssh2lonlatt(product, timestep, meshsize, timelim, options)

    arguments (Input)
        product SSHProduct {mustBeScalarOrEmpty} = 'MEaSUREs'
        timestep {mustBeTimeStep} = []
        meshsize = []
        timelim {mustBeTimeRange} = []
        options.LonOrigin {mustBeFinite, mustBeReal} = []
        options.Interpolation ...
            {mustBeTextScalar, mustBeMember(options.Interpolation, {'linear', 'nearest', 'next', 'previous', 'spline', 'pchip'})} = 'linear'
        options.TimeFormat (1, :) DateFormat = 'datetime'
        options.OutputFormat (1, :) MeshFormat = 'meshgrid'
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, {'mm', 'm'})} = 'm'
        options.ForceNew (1, 1) {mustBeNumericOrLogical} = false
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
        options.SaveData (1, 1) {mustBeNumericOrLogical} = true
        options.CallChain cell = {}
    end

    arguments (Output)
        ssh (:, :, :) {mustBeReal}
        sshSigma (:, :, :) {mustBeReal}
        dates (:, 1) {mustBeVector}
        lon (:, :) {mustBeReal, mustBeVector}
        lat (:, :) {mustBeReal, mustBeVector}
    end

    lonOrigin = options.LonOrigin;
    intpMthd = options.Interpolation;
    timeFmt = char(options.TimeFormat);
    outputFmt = char(options.OutputFormat);
    unit = lower(options.Unit);
    forceNew = logical(options.ForceNew);
    beQuiet = uint8(double(options.BeQuiet) * 2);
    saveData = logical(options.SaveData);
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

    if ~isempty(timestep) && isnumeric(timestep)
        timestep = days(timestep);
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    if ~isempty(timelim) && timelim(1) > timelim(2)
        warning('The start time (%s) is later than the end time (%s), flipping the order.', timelim(1), timelim(2));
        timelim = fliplr(timelim);
    end

    %% Checking for existing data
    intpOutputPath = ...
        outputpath(product, intpMthd, timestep, meshsize, lonOrigin);

    vars = {'sshs', 'dates', 'lon', 'lat'};

    if ~forceNew && exist(intpOutputPath, 'file') && ...
            all(ismember(vars, who('-file', intpOutputPath)))

        if beQuiet <= 1
            t = tic;
            msg = 'this may take a while...';
            fprintf('[ULMO>%s] Loading <a href="matlab: fprintf(''%s\\n'');open(''%s'')">interpolated SSH product</a>, %s\n', ...
                callchaintext(callChain), intpOutputPath, intpOutputPath, msg);
        end

        if ismember('sshErrors', who('-file', intpOutputPath))
            data = load(intpOutputPath, vars{:}, 'sshErrors');
        else
            data = load(intpOutputPath, vars{:});
            data.sshErrors = [];
        end

        if beQuiet <= 1
            fprintf(repmat('\b', 1, length(msg) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

        if nargout > 0
            [ssh, sshSigma, dates, lon, lat] = formatoutput( ...
                data.sshs, data.sshErrors, data.dates, data.lon, data.lat, ...
                timelim, timeFmt, outputFmt, unit);
        else
            plotsealeveltseries(data.dates, data.sshs, data.lon, data.lat, ...
                product, unit, 'Global mean sea surface height');
        end

        return
    end

    %% Main
    inputPath = outputpath(product);

    if ~exist(inputPath, 'file')
        error('ULMO:InputDataNotFound', ...
            'Input data not found at %s.', inputPath);
    elseif ~all(ismember(vars, who('-file', inputPath)))
        notFoundVars = vars(~ismember(vars, who('-file', inputPath)));
        error('ULMO:InputDataIncomplete', ...
            'Input data at %s is missing variables: %s.', ...
            inputPath, strjoin(notFoundVars, ', '));
    end

    if ismember('sshErrors', who('-file', inputPath))
        data = load(inputPath, vars{:}, 'sshErrors');
    else
        data = load(inputPath, vars{:});
        data.sshErrors = [];
    end

    if ~isempty(timestep)
        [data.sshs, data.dates] = interptemporal( ...
            data.dates, data.sshs, timestep, intpMthd, BeQuiet = beQuiet, CallChain = callChain);

        if ~isempty(data.sshErrors)
            data.sshErrors = interptemporal( ...
                data.dates, data.sshErrors, timestep, intpMthd, beQuiet, BeQuiet = beQuiet, CallChain = callChain);
        end

    end

    if ~isempty(meshsize) || ~isempty(lonOrigin)

        if isempty(meshsize)
            meshsize = 1/6;
        end

        if isempty(lonOrigin)
            lonOrigin = 180;
        end

        [data.sshs, data.lon, data.lat] = interpspatial( ...
            data.lon, data.lat, data.sshs, meshsize, lonOrigin, intpMthd, beQuiet, callChain);

        if ~isempty(data.sshErrors)
            data.sshErrors = interpspatial( ...
                data.lon, data.lat, data.sshErrors, meshsize, lonOrigin, intpMthd, beQuiet, callChain);
        end

    end

    if saveData
        save(intpOutputPath, ...
            '-struct', 'data', 'sshs', 'sshErrors', 'dates', 'lon', 'lat', '-v7.3');

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">SSH product</a>.\n', ...
                callchaintext(callChain), intpOutputPath, intpOutputPath);
        end

    end

    if nargout > 0
        [ssh, sshSigma, dates, lon, lat] = formatoutput( ...
            data.sshs, data.sshErrors, data.dates, data.lon, data.lat, ...
            timelim, timeFmt, outputFmt, unit);
    else
        plotsealeveltseries(data.dates, data.sshs, data.lon, data.lat, ...
            product, unit, 'Global mean sea surface height');
    end

end

%% Subfunctions
% Interpolation
function [meshIntp, lonIntp, latIntp] = ...
        interpspatial(lon, lat, mesh, meshSize, lonOrigin, intpMthd, beQuiet, callChain)

    lonIntp = (-180:meshSize:180) + lonOrigin;
    latIntp = -90:meshSize:90;
    [lonnIntp, lattIntp] = meshgrid(lonIntp, latIntp);

    if isequal(lonIntp(:), lon(:)) && isequal(latIntp(:), lat(:))
        meshIntp = mesh;
        return
    end

    if beQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO>%s] Interpolating spatially, %s\n', ...
            callchaintext(callChain), templine);
    end

    % lonnIntp(lonnIntp < ogLonOrigin - 180) = lonnIntp(lonnIntp < ogLonOrigin - 180) + 360;
    % lonnIntp(lonnIntp > ogLonOrigin + 180) = lonnIntp(lonnIntp > ogLonOrigin + 180) - 360;

    lon = lon - 360 * floor(lon(1) / 360);

    if lon(1) ~= 0

        if lon(end) - lon(1) == 360
            lon = lon(1:end - 1);
            mesh = mesh(1:end - 1, :, :);
        end

        [lon, idx] = sort(mod(lon, 360));
        mesh = mesh(idx, :, :);
    end

    if lon(end) - lon(1) ~= 360
        lon = [lon(end) - 360; lon; lon(1) + 360];
        mesh = cat(1, mesh(end, :, :), mesh, mesh(1, :, :));
    end

    [lonn, latt] = meshgrid(lon, lat);
    mesh = permute(mesh, [2, 1, 3]);

    meshIntp = nan([length(latIntp), length(lonIntp), size(mesh, 3)], 'single');

    for iDate = 1:size(meshIntp, 3)
        meshIntp(:, :, iDate) = ...
            interp2(lonn, latt, squeeze(mesh(:, :, iDate)), ...
            mod(lonnIntp, 360), lattIntp, intpMthd);
    end

    meshIntp = permute(meshIntp, [2, 1, 3]);

    if beQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end

% % Load the aggregated raw data (i.e. data from all time stamps in the same array)
% function [sshs, sshErrors, dates, lon, lat] = loadAggregatedRawData( ...
%         product, forceNew, beQuiet, saveData, outputFolder, inputFolder, callChain)
%     vars = {'sshs', 'sshErrors', 'dates', 'lon', 'lat'};
%     % Output location
%     outputFile = sprintf('%s.mat', product);
%     outputPath = fullfile(outputFolder, outputFile);

%     %% Checking output file
%     if exist(outputPath, 'file') && ~forceNew && ...
%             all(ismember(vars, who('-file', outputPath)))

%         if beQuiet <= 1
%             t = tic;
%             msg = 'this may take a while...';
%             fprintf('[ULMO>%s] Loading <a href="matlab: fprintf(''%s\\n'');open(''%s'')">SSH product</a>, %s\n', ...
%                 callchaintext(callChain), outputPath, outputPath, msg);
%         end

%         data = load(outputPath, vars{:});
%         sshs = data.sshs;
%         sshErrors = data.sshErrors;
%         dates = data.dates;
%         lon = data.lon;
%         lat = data.lat;

%         if beQuiet <= 1
%             fprintf(repmat('\b', 1, length(msg) + 1));
%             fprintf('took %.1f seconds.\n', toc(t));
%         end

%         return

%     end

%     if ~exist(inputFolder, 'dir')
%         error('ULMO:LoadData:FileNotFound', ...
%             'Input folder not found at %s\n', inputFolder);
%     end

%     inputFileTemplate = 'ssh_grids_*.nc';
%     inputFiles = {dir(fullfile(inputFolder, inputFileTemplate)).name};

%     if isempty(inputFiles)
%         error('ULMO:LoadData:FileNotFound', ...
%             'No input files found at %s, data can be accessed from <a href="%s">PO.DAAC</a>\n', ...
%             inputFolder, 'https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205');
%     end

%     inputFiles = sort(inputFiles);

%     nDates = length(inputFiles);

%     iDate = 1;
%     inputFile = inputFiles{iDate};
%     inputPath = fullfile(inputFolder, inputFile);

%     %% Initialisation and preallocation
%     wbar = waitbar(iDate / nDates, ...
%         sprintf('Reading raw data (%d/%d)', iDate, nDates), ...
%         "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

%     [ssh, sshError, date, dateRef, lon, lat] = readRawData(inputPath);

%     sshsSize = [length(lon), length(lat), nDates];
%     dates = nan([1, nDates]);
%     sshs = nan(sshsSize, 'single');
%     sshErrors = nan(sshsSize, 'single');

%     dates(iDate) = date;
%     sshs(:, :, iDate) = ssh;
%     sshErrors(:, :, iDate) = sshError;

%     %% Looping over all files
%     for iDate = 2:nDates
%         waitbar(iDate / nDates, wbar, ...
%             sprintf('Reading raw data (%d/%d)', iDate, nDates));

%         if getappdata(wbar, 'canceling')
%             delete(wbar);
%             warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
%             'Processing cancelled');
%             return
%         end

%         % Read the data
%         [ssh, sshError, date] = ...
%             readRawData(fullfile(inputFolder, inputFiles{iDate}));

%         % Store the data
%         dates(iDate) = date;
%         sshs(:, :, iDate) = ssh;
%         sshErrors(:, :, iDate) = sshError;
%     end

%     dates = dateRef + days(dates);

%     if ~saveData
%         delete(wbar);
%         return
%     end

%     waitbar(iDate / nDates, wbar, 'Saving data (this may take a while)');

%     if getappdata(wbar, 'canceling')
%         delete(wbar);
%         warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
%         'Saving cancelled');
%         return
%     end

%     save(outputPath, vars{:}, '-v7.3');
%     delete(wbar);

%     if beQuiet <= 1
%         fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">SSH product</a>.\n', ...
%             callchaintext(callChain), outputPath, outputPath);
%     end

% end

% % Read the raw data of a sigle time stamp from the NetCDF file
% function [ssh, sshError, time, timeRef, lon, lat] = readRawData(inputPath)
%     % Read the data
%     ssh = ncread(inputPath, 'SLA');
%     sshError = ncread(inputPath, 'SLA_ERR');
%     time = ncread(inputPath, 'Time');

%     % Convert format
%     ssh = single(ssh);
%     sshError = single(sshError);

%     if nargout <= 3
%         return
%     end

%     % Read the time reference if asked
%     inputInfo = ncinfo(inputPath);
%     timeInfo = ...
%         inputInfo.Variables(strcmp({inputInfo.Variables.Name}, 'Time')).Attributes;
%     timeRef = timeInfo(strcmp({timeInfo.Name}, 'units')).Value;
%     timeRef = ...
%         datetime(strrep(timeRef, 'Days since ', ''), ...
%         "InputFormat", 'yyyy-MM-dd HH:mm:ss');

%     % Read lon/lat if asked
%     lon = ncread(inputPath, 'Longitude');
%     lat = ncread(inputPath, 'Latitude');
%     lon = single(lon);
%     lat = single(lat);
% end

% Format the output
function [sshs, sshErrors, dates, lon, lat] = ...
        formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt, outputFmt, unit)

    if ~isempty(timelim)
        isValidTime = ...
            (dates >= timelim(1)) & (dates <= timelim(2));
        sshs = sshs(:, :, isValidTime);

        if ~isempty(sshErrors)
            sshErrors = sshErrors(:, :, isValidTime);
        end

        dates = dates(isValidTime);
    end

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    switch outputFmt
        case 'meshgrid'
            lon = lon(:)';
            lat = lat(:);
            sshs = permute(sshs, [2, 1, 3]);

            if ~isempty(sshErrors)
                sshErrors = permute(sshErrors, [2, 1, 3]);
            end

        case 'ndgrid'
            lon = lon(:);
            lat = lat(:)';
    end

    switch unit
        case 'm'
        case 'mm'
            sshs = sshs * 1e3;

            if ~isempty(sshErrors)
                sshErrors = sshErrors * 1e3;
            end

    end

end

% Get output path for interpolated data
function outputPath = ...
        outputpath(product, intpMthd, timeStep, meshSize, lonOrigin)

    arguments (Input)
        product (1, :) char
        intpMthd (1, :) char = []
        timeStep {mustBeTimeStep} = []
        meshSize {mustBePositive} = []
        lonOrigin {mustBeFinite, mustBeReal} = []
    end

    arguments (Output)
        outputPath (1, :) char
    end

    outputFolder = fullfile(getenv('IFILES'), 'SSH');

    if isempty(intpMthd) && isempty(timeStep) && ...
            isempty(meshSize) && isempty(lonOrigin)
        outputPath = fullfile(outputFolder, sprintf('%s.mat', product));
        return
    end

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

    outputFile = sprintf('%s%s%s-%s.mat', ...
        product, timeStepStr, spaceIntpStr, intpMthd);
    outputPath = fullfile(outputFolder, outputFile);
end
