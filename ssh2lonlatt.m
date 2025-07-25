%% SSH2LONLATT
% Reads in sea surface height (anomaly) products and converts them to a 3D
% array.
%
% Syntax
%	ssh = SSH2LONLATT(product)
%	ssh = SSH2LONLATT(product, TimeRange, timeStep, meshSize)
%	ssh = SSH2LONLATT(__, "Name", Value)
%	[ssh, sshError, dates, lon, lat] = SSH2XYZ(__)
%
% Input arguments
%	product - The product to load
%		- 'MEaSUREs': The MEaSUREs Gridded Sea Surface Height Anomalies
%           Version 2205 product.
%		The default option is 'MEaSUREs', which is the only option
%       available as of right now.
%		Data types: char
%	timeStep - The duration of the interpolation
% 		When not specified, the default product time grid is used, which is
%       not necessarily regular (MEaSUREs: 5 days).
%		When specified (non-empty), the function will interpolate the data
%       to a regular time grid with the specified duration. The input can
%       be:
%		- A DURATION object (e.g. days(5)).
%		- A numeric value (e.g. 5) interpreted as days.
%		The default option is empty (no interpolation).
%		Data types: duration | [numeric]
%	TimeRange - The time range to load the data
%		A 2-element vector with the start and end dates of the time range
%       to load.
%		The input can be:
%		- A DATETIME object (e.g. datetime([2002, 2022], 1, 1)).
%		- A numeric value (e.g. [737791, 737821]) interpreted as DATENUM.
%		The default option is empty (no time range truncation). No matter
%       the time range, the full range of data is processed and saved.
%		Data types: datetime | [numeric]
%	InterpolationMethod - The interpolation method to use, if interpolation
%       is requested
%		The default option is 'linear'. This function does not check the
%       validity of the method, which is passed to the INTERP1 function.
%		Data types: char
%	TimeFormat - The time format to use for the time stamp output
%		- 'datetime': The time is returned as a DATETIME object.
%		- 'datenum': The time is returned as a numeric DATENUM object.
%		The default option is 'datetime'.
%		Data types: char
%	ForceNew - Whether to force reprocess of the data
%		- true: Reprocess the data
%		- false: Only reprocess if previous output is not found
%		The default option is 'soft-true', only reprocessing the data if
%       the output file is older than the input files. This option is not
%       explicitly exposed; to enforce this option, set it to 0.5.
%		Data types: logical | [numeric]
%	BeQuiet - Whether to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is 'soft-true', only printing the most important
%       messages. This option is not explicitly exposed; to enforce this
%       option, set it to 0.5.
%		Data types: logical | [numeric]
%	SaveData - Whether to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%		The default option is true.
%		Data types: logical | [numeric]
%
% Output arguments
%	ssh - The sea surface height (anomaly) data
%		The data is a single-precision 3D array with the dimensions
%       [lon, lat, time], and the units are in metres [m].
%	sshError - The sea surface height (anomaly) error data
%		The units and dimensions are the same as for SSH.
%	dates - The time stamps of the data
%		The data is a 1D array with the dimensions [1, time], and the units
%       are as specified by the TimeFormat input argument.
%	lon - The longitude of the data
%		The data is a 1D array with the dimensions [lon], and the units are
%       in degrees [°].
%	lat - The latitude of the data
%		The units and dimensions are the same as for LON.
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
% Last modified by
%	2025/07/22, williameclee@arizona.edu (@williameclee)

function varargout = ssh2lonlatt(varargin)
    [product, timeStep, meshSize, timelim, lonOrigin, intpMthd, timeFmt, unit, ...
         forceNew, beQuiet, saveData] = ...
        parseInputs(varargin{:});

    outputFolder = fullfile(getenv('IFILES'), product);
    inputFolder = fullfile(outputFolder, 'raw');

    %% No interpolation
    if isempty(timeStep) && isempty(meshSize) && isempty(lonOrigin)
        % No interpolation duration given, load the raw data
        [sshs, sshErrors, dates, lon, lat] = ...
            loadAggregatedRawData(product, ...
            forceNew, beQuiet, saveData, outputFolder, inputFolder);

        [sshs, sshErrors, dates, lon, lat] = ...
            formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt, unit);

        varargout = {sshs, sshErrors, dates, lon, lat};
        return

    end

    %% With interpolation
    intpOutputPath = outputpath(product, intpMthd, timeStep, meshSize, lonOrigin, outputFolder);

    if exist(intpOutputPath, 'file') && ~forceNew

        if beQuiet <= 1
            loadingMsg = ...
                sprintf('%s: loading %s, this may take a while...\n', ...
                upper(mfilename), intpOutputPath);
            fprintf(loadingMsg);
        end

        load(intpOutputPath, 'sshs', 'sshErrors', 'dates', 'lon', 'lat');

        if beQuiet <= 1
            fprintf(repmat('\b', 1, length(loadingMsg)));
            fprintf('%s: loaded %s\n', upper(mfilename), intpOutputPath);
        end

        [sshs, sshErrors, dates, lon, lat] = ...
            formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt, unit);
        varargout = {sshs, sshErrors, dates, lon, lat};

        return
    end

    [sshs, sshErrors, dates, lon, lat] = ...
        loadAggregatedRawData(product, ...
        forceNew, beQuiet, saveData, outputFolder, inputFolder);

    if ~isempty(timeStep)
        sshs = ...
            interptemporal(dates, sshs, timeStep, intpMthd, beQuiet);
        [sshErrors, dates] = ...
            interptemporal(dates, sshErrors, timeStep, intpMthd, beQuiet);
    end

    if ~isempty(meshSize) || ~isempty(lonOrigin)

        if isempty(meshSize)
            meshSize = 1/6;
        end

        if isempty(lonOrigin)
            lonOrigin = 180;
        end

        sshs = interpspatial(lon, lat, sshs, meshSize, lonOrigin, intpMthd, beQuiet);
        [sshErrors, lon, lat] = ...
            interpspatial(lon, lat, sshErrors, meshSize, lonOrigin, intpMthd, beQuiet);
    end

    if saveData
        save(intpOutputPath, ...
            'sshs', 'sshErrors', 'dates', 'lon', 'lat', '-v7.3');

        if beQuiet <= 1
            fprintf('%s: saved %s\n', upper(mfilename), intpOutputPath);
        end

    end

    [sshs, sshErrors, dates, lon, lat] = ...
        formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt, unit);
    varargout = {sshs, sshErrors, dates, lon, lat};

end

%% Subfunctions
% Interpolation
function [meshIntp, datesIntp] = ...
        interptemporal(dates, mesh, timeStep, intpMthd, beQuiet)

    if beQuiet <= 1
        fprintf('%s: Interpolating temporally, this may take a while...\n', ...
            upper(mfilename));
    end

    if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
        startDate = datetime(year(dates(1)), month(dates(1)), 1, 0, 0, 0);
        endDate = datetime(year(dates(end)), month(dates(end)), 1, 0, 0, 0) + calmonths(1);
        dmonths = ceil((endDate - startDate) / days(28));
        startDates = startDate + calmonths(0:dmonths - 1);
        endDates = startDate + calmonths(1:dmonths);
        datesIntp = startDates + (endDates - startDates) / 2;
        datesIntp = datesIntp(datesIntp >= dates(1) & datesIntp <= dates(end));
    else

        if mean(diff(dates)) > timeStep
            warning(sprintf('%s:InterpolationStepTooSmall', upper(mfilename)), ...
                'The interpolation time step (%s) is smaller than the mean data resolution (%s)', timeStep, mean(diff(dates)));
        end

        datesIntp = dates(1):timeStep:dates(end);
    end

    meshFlat = reshape(mesh, [prod(size(mesh, 1:2)), size(mesh, 3)])';
    meshIntp = interp1(dates, meshFlat, datesIntp, intpMthd)';
    meshIntp = reshape(meshIntp, [size(mesh, 1:2), length(datesIntp)]);
end

function [meshIntp, lonIntp, latIntp] = ...
        interpspatial(lon, lat, mesh, meshSize, lonOrigin, intpMthd, beQuiet)

    if beQuiet <= 1
        fprintf('%s: Interpolating spatially, this may take a while...\n', ...
            upper(mfilename));
    end

    ogLonOrigin = 180;
    lonIntp = (-180:meshSize:180) + lonOrigin;
    latIntp = -90:meshSize:90;
    [lonnIntp, lattIntp] = meshgrid(lonIntp, latIntp);
    lonnIntp(lonnIntp < ogLonOrigin - 180) = lonnIntp(lonnIntp < ogLonOrigin - 180) + 360;
    lonnIntp(lonnIntp > ogLonOrigin + 180) = lonnIntp(lonnIntp > ogLonOrigin + 180) - 360;
    lonPad = [lon(end) - 360; lon; lon(1) + 360];
    meshPad = cat(1, mesh(end, :, :), mesh, mesh(1, :, :));
    [lonn, latt] = meshgrid(lonPad, lat);

    meshPad = permute(meshPad, [2, 1, 3]);

    meshIntp = nan([length(latIntp), length(lonIntp), size(mesh, 3)], 'single');

    for iDate = 1:size(meshIntp, 3)
        meshIntp(:, :, iDate) = ...
            interp2(lonn, latt, squeeze(meshPad(:, :, iDate)), ...
            lonnIntp, lattIntp, intpMthd);
    end

    meshIntp = permute(meshIntp, [2, 1, 3]);
end

% Parse input arguments
function varargout = parseInputs(varargin)
    ip = inputParser;
    addOptional(ip, 'Product', 'MEaSUREs', ...
        @(x) ischar(validatestring(x, {'MEaSUREs'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || ((isnumeric(x) || isdatetime(x)) && length(x) == 2));
    addOptional(ip, 'TimeStep', [], ...
        @(x) isempty(x) || ((isduration(x) || isnumeric(x)) && isscalar(x)) || ischar(validatestring(x, {'midmonth'})));
    addOptional(ip, 'MeshSize', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addOptional(ip, 'LonOrigin', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'InterpolationMethod', 'linear', @ischar);
    addParameter(ip, 'TimeFormat', 'datetime', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'Unit', 'm', ...
        @(x) ischar(validatestring(x, {'mm', 'm'})));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));

    parse(ip, varargin{:});

    product = ip.Results.Product;
    timelim = ip.Results.TimeRange;

    timeStep = ip.Results.TimeStep;
    meshSize = ip.Results.MeshSize;
    lonOrigin = ip.Results.LonOrigin;
    intpMthd = ip.Results.InterpolationMethod;

    timeFmt = ip.Results.TimeFormat;
    unit = lower(ip.Results.Unit);
    forceNew = logical(ip.Results.ForceNew);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);
    saveData = logical(ip.Results.SaveData);

    if ~isempty(timeStep) && isnumeric(timeStep)
        timeStep = days(timeStep);
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    varargout = ...
        {product, timeStep, meshSize, timelim, lonOrigin, intpMthd, timeFmt, unit, forceNew, beQuiet, saveData};
end

% Load the aggregated raw data (i.e. data from all time stamps in the same array)
function varargout = loadAggregatedRawData( ...
        product, forceNew, beQuiet, saveData, outputFolder, inputFolder)
    % Output location
    outputFile = sprintf('%s.mat', product);
    outputPath = fullfile(outputFolder, outputFile);

    %% Checking output file
    if exist(outputPath, 'file') && ...
            (forceNew == 0 || (forceNew == 1 && isolder(inputFolder, outputPath, true)))

        if beQuiet <= 1
            loadingMsg = ...
                sprintf('%s: Loading %s, this may take a while...\n', ...
                upper(mfilename), outputPath);
            fprintf(loadingMsg);
        end

        load(outputPath, 'sshs', 'sshErrors', 'dates', 'lon', 'lat');
        varargout = {sshs, sshErrors, dates, lon, lat};

        if beQuiet <= 1
            fprintf(repmat('\b', 1, length(loadingMsg)));
            fprintf('%s: Loaded %s\n', upper(mfilename), outputPath);
        end

        return

    end

    if ~exist(inputFolder, 'dir')
        error(sprintf('%s:InputNotFound', upper(mfilename)), ...
            'Input folder not found at %s\n', inputFolder);
    end

    inputFileTemplate = 'ssh_grids_*.nc';
    inputFiles = {dir(fullfile(inputFolder, inputFileTemplate)).name};

    if isempty(inputFiles)
        error(sprintf('%s:InputNotFound', upper(mfilename)), ...
            'No input files found at %s, data can be accessed from %s\n', ...
            inputFolder, 'https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205');
    end

    inputFiles = sort(inputFiles);

    nDates = length(inputFiles);

    iDate = 1;
    inputFile = inputFiles{iDate};
    inputPath = fullfile(inputFolder, inputFile);

    %% Initialisation and preallocation
    wbar = waitbar(iDate / nDates, ...
        sprintf('Reading raw data (%d/%d)', iDate, nDates), ...
        "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

    [ssh, sshError, date, dateRef, lon, lat] = readRawData(inputPath);

    sshsSize = [length(lon), length(lat), nDates];
    dates = nan([1, nDates]);
    sshs = nan(sshsSize, 'single');
    sshErrors = nan(sshsSize, 'single');

    dates(iDate) = date;
    sshs(:, :, iDate) = ssh;
    sshErrors(:, :, iDate) = sshError;

    %% Looping over all files
    for iDate = 2:nDates
        waitbar(iDate / nDates, wbar, ...
            sprintf('Reading raw data (%d/%d)', iDate, nDates));

        if getappdata(wbar, 'canceling')
            delete(wbar);
            warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
            'Processing cancelled');
            return
        end

        % Read the data
        [ssh, sshError, date] = ...
            readRawData(fullfile(inputFolder, inputFiles{iDate}));

        % Store the data
        dates(iDate) = date;
        sshs(:, :, iDate) = ssh;
        sshErrors(:, :, iDate) = sshError;
    end

    dates = dateRef + days(dates);

    varargout = {sshs, sshErrors, dates, lon, lat};

    if ~saveData
        delete(wbar);
        return
    end

    waitbar(iDate / nDates, wbar, 'Saving data (this may take a while)');

    if getappdata(wbar, 'canceling')
        delete(wbar);
        warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
        'Saving cancelled');
        return
    end

    save(outputPath, 'sshs', 'sshErrors', 'dates', ...
        'lon', 'lat', '-v7.3');
    delete(wbar);

    if beQuiet <= 1
        fprintf('%s: Saved %s\n', upper(mfilename), outputPath);
    end

end

% Read the raw data of a sigle time stamp from the NetCDF file
function varargout = readRawData(inputPath)
    % Read the data
    ssh = ncread(inputPath, 'SLA');
    sshError = ncread(inputPath, 'SLA_ERR');
    time = ncread(inputPath, 'Time');

    % Convert format
    ssh = single(ssh);
    sshError = single(sshError);

    if nargout <= 3
        varargout = {ssh, sshError, time};
        return
    end

    % Read the time reference if asked
    inputInfo = ncinfo(inputPath);
    timeInfo = ...
        inputInfo.Variables(strcmp({inputInfo.Variables.Name}, 'Time')).Attributes;
    timeRef = timeInfo(strcmp({timeInfo.Name}, 'units')).Value;
    timeRef = ...
        datetime(strrep(timeRef, 'Days since ', ''), ...
        "InputFormat", 'yyyy-MM-dd HH:mm:ss');

    % Read lon/lat if asked
    lon = ncread(inputPath, 'Longitude');
    lat = ncread(inputPath, 'Latitude');
    lon = single(lon);
    lat = single(lat);

    varargout = {ssh, sshError, time, timeRef, lon, lat};
end

% Format the output
function varargout = formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt, unit)

    if ~isempty(timelim)
        isValidTime = ...
            (dates >= timelim(1)) & (dates <= timelim(2));
        sshs = sshs(:, :, isValidTime);
        sshErrors = sshErrors(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    switch unit
        case 'm'
        case 'mm'
            sshs = sshs * 1e3;
            sshErrors = sshErrors * 1e3;
        otherwise
            error(sprintf('%s:InvalidUnit', upper(mfilename)), ...
                'Unrecognised unit for sea surface height %s', unit);
    end

    varargout = {sshs, sshErrors, dates, lon, lat};

end

% Compare the date of two directories or files
function isolder = isolder(dir1, dir2, checkContent)

    dirDates = NaT([2, 1]);

    for iDir = 1:2

        if iDir == 1
            dirPath = dir1;
        else
            dirPath = dir2;
        end

        dirState = exist(dirPath, 'file');

        switch dirState
            case 2 % is a file
                dirDate = datetime(dir(dirPath).date, ...
                    "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
            case 7 % is a folder
                dirInfo = dir(dirPath);

                if ~checkContent
                    dirDate = datetime(dirInfo(strcmp({dirInfo.name}, '.')).date, ...
                        "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
                else
                    disDates = ...
                        datetime({dirInfo(~ismember({dirInfo.name}, {'.', '..', '.DS_Store'})).date}, ...
                        "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
                    dirDate = max(disDates);
                end

            case 0 % does not exist
                dirDate = datetime('now');
            otherwise
                error('Unknown file type for %s, EXIST output: %d', dirPath, dirState);
        end

        dirDates(iDir) = dirDate;

    end

    isolder = dirDates(1) < dirDates(2);

end

% Get output path for interpolated data
function outputPath = outputpath(product, intpMthd, timeStep, meshSize, lonOrigin, outputFolder)
    timeStepStr = '';

    if ~isempty(timeStep)

        if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
            timeStepStr = '-Tmidmonth';
        elseif isduration(timeStep)
            timeStepStr = sprintf('-T%s', erase(sprintf('%s', timeStep), ' '));
        elseif isnumeric(timeStep)
            timeStepStr = sprintf('-T%d', timeStep);
        else
            error(sprintf('%s:InvalidTimeStep', upper(mfilename)), ...
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
