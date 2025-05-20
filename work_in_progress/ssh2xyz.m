%% SSH2XYZ
% Reads in sea surface height (anomaly) products and converts them to a 3D
% array.
%
% Syntax
%	ssh = SSH2XYZ(product)
%	ssh = SSH2XYZ(product, InterpolationDuration, TimeRange)
%	ssh = SSH2XYZ(__, "Name", Value)
%	[ssh, sshError, dates, lon, lat] = SSH2XYZ(__)
%
% Input arguments
%	product - The product to load
%		- 'MEaSUREs': The MEaSUREs Gridded Sea Surface Height Anomalies
%           Version 2205 product.
%		The default option is 'MEaSUREs', which is the only option
%       available as of right now.
%		Data types: char
%	InterpolationDuration - The duration of the interpolation
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

function varargout = ssh2xyz(varargin)
    [product, intpDuration, timelim, intpMthd, timeFmt, ...
         forceNew, beQuiet, saveData] = ...
        parseInputs(varargin{:});

    outputFolder = fullfile(getenv('IFILES'), product);
    inputFolder = fullfile(outputFolder, 'raw');

    if isempty(intpDuration)
        % No interpolation duration given, load the raw data

        if ~beQuiet
            loadingMsg = ...
                sprintf('%s: loading %s (this may take a while)\n', ...
                upper(mfilename), outputFolder);
        end

        [sshs, sshErrors, dates, lon, lat] = ...
            loadAggregatedRawData(product, ...
            forceNew, beQuiet, saveData, outputFolder, inputFolder);

        [sshs, sshErrors, dates, lon, lat] = ...
            formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt);

        varargout = {sshs, sshErrors, dates, lon, lat};

        if ~beQuiet
            fprintf(repmat('\b', 1, length(loadingMsg)));
            fprintf('%s: loaded %s\n', upper(mfilename), outputFolder);
        end

        return

    end

    intpOutputFile = sprintf('%s-It%s.mat', ...
        product, erase(sprintf('%s', intpDuration), ' '));
    intpOutputPath = fullfile(outputFolder, intpOutputFile);

    if exist(intpOutputPath, 'file') && forceNew <= 1

        if beQuiet <= 1
            loadingMsg = ...
                sprintf('%s: loading %s (this may take a while)\n', ...
                upper(mfilename), intpOutputPath);
            fprintf(loadingMsg);
        end

        load(intpOutputPath, 'intpSshs', 'intpSshErrors', 'intpDates', 'lon', 'lat');

        [intpSshs, intpSshErrors, intpDates, lon, lat] = ...
            formatoutput(intpSshs, intpSshErrors, intpDates, lon, lat, timelim, timeFmt);
        varargout = {intpSshs, intpSshErrors, intpDates, lon, lat};

        if beQuiet <= 1
            fprintf(repmat('\b', 1, length(loadingMsg)));
            fprintf('%s: loaded %s\n', upper(mfilename), intpOutputPath);
        end

        return
    end

    [sshs, sshErrors, dates, lon, lat] = ...
        loadAggregatedRawData(product, ...
        forceNew, beQuiet + (beQuiet == 1), saveData, outputFolder, inputFolder);
    intpDates = dates(1):intpDuration:dates(end);

    flatSshs = reshape(sshs, [prod(size(sshs, 1:2)), size(sshs, 3)])';
    intpSshs = interp1(dates, flatSshs, intpDates, intpMthd)';
    intpSshs = reshape(intpSshs, [size(sshs, 1:2), length(intpDates)]);

    flatSshErrors = ...
        reshape(sshErrors, [prod(size(sshErrors, 1:2)), size(sshErrors, 3)])';
    intpSshErrors = interp1(dates, flatSshErrors, intpDates, intpMthd)';
    intpSshErrors = ...
        reshape(intpSshErrors, [size(sshErrors, 1:2), length(intpDates)]);

    [intpSshs, intpSshErrors, intpDates, lon, lat] = ...
        formatoutput(intpSshs, intpSshErrors, intpDates, lon, lat, timelim, timeFmt);
    varargout = {intpSshs, intpSshErrors, intpDates, lon, lat};

    if ~saveData
        return
    end

    save(fullfile(outputFolder, intpOutputFile), ...
        'intpSshs', 'intpSshErrors', 'intpDates', 'lon', 'lat', '-v7.3');

    if beQuiet <= 1
        fprintf('%s: saved %s\n', upper(mfilename), intpOutputFile);
    end

end

%% Subfunctions
% Parse input arguments
function varargout = parseInputs(varargin)
    ip = inputParser;
    addOptional(ip, 'Product', 'MEaSUREs', ...
        @(x) ischar(validatestring(x, {'MEaSUREs'})));
    addOptional(ip, 'InterpolationDuration', [], ...
        @(x) isempty(x) || ((isduration(x) || isnumeric(x)) && isscalar(x)));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || ((isnumeric(x) || isdatetime(x)) && length(x) == 2));
    addParameter(ip, 'InterpolationMethod', 'linear', @ischar);
    addParameter(ip, 'TimeFormat', 'datetime', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'ForceNew', 0.5, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));

    parse(ip, varargin{:});
    product = ip.Results.Product;
    intpDuration = ip.Results.InterpolationDuration;

    timelim = ip.Results.TimeRange;

    intpMthd = ip.Results.InterpolationMethod;

    timeFmt = ip.Results.TimeFormat;
    forceNew = uint8(double(ip.Results.ForceNew) * 2);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);
    saveData = logical(ip.Results.SaveData);

    if ~isempty(intpDuration) && isnumeric(intpDuration)
        intpDuration = days(intpDuration);
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    varargout = ...
        {product, intpDuration, timelim, intpMthd, timeFmt, forceNew, beQuiet, saveData};
end

% Load the aggregated raw data (i.e. data from all time stamps in the same array)
function varargout = loadAggregatedRawData( ...
        product, forceNew, beQuiet, saveData, outputFolder, inputFolder)
    % Output location
    outputFile = sprintf('%s.mat', product);
    outputPath = fullfile(outputFolder, outputFile);

    %% Checking output file
    if exist(outputPath, 'file') && forceNew <= 1
        % Output file exists and not forced to recompute
        outputDate = datetime(dir(outputPath).date, ...
            "InputFormat", 'dd-MMM-yyyy HH:mm:ss');

        if exist(inputFolder, 'dir')
            dirInfo = dir(inputFolder);
            inputDate = datetime(dirInfo(strcmp({dirInfo.name}, '.')).date, ...
                "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
        else
            inputDate = datetime('now');
        end

        disp(outputDate) % DEBUG
        disp(inputDate) % DEBUG

        if ~(forceNew == 1 && outputDate >= inputDate)
            % Passed or skipped the time check

            if beQuiet <= 1
                loadingMsg = ...
                    sprintf('%s: Loading %s (this may take a while)\n', ...
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
        elseif beQuiet <= 1
            fprintf('%s: %s out of date, recomputing...\n', ...
                upper(mfilename), outputPath);
        end

    end

    if ~exist(inputFolder, 'dir')
        error(sprintf('%s:InputNotFound', upper(mfilename)), ...
            '%s: Input folder not found at %s\n', ...
            upper(mfilename), inputFolder);
    end

    inputFileTemplate = 'ssh_grids_*.nc';
    inputFiles = {dir(fullfile(inputFolder, inputFileTemplate)).name};

    if isempty(inputFiles)
        error(sprintf('%s:InputNotFound', upper(mfilename)), ...
            '%s: No input files found at %s, data can be accessed from %s\n', ...
            upper(mfilename), inputFolder, 'https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205');
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
                '%s: Processing cancelled\n', upper(mfilename));
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
            '%s: Saving cancelled\n', upper(mfilename));
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
function varargout = formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt)

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

    varargout = {sshs, sshErrors, dates, lon, lat};

end
