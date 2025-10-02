%% GRACEDEG1
% This function reads and formats the degree-1 spherical harmonic
% correction from GRACE/GRACE-FO Technical Note 13.
%
% Syntax
%   [deg1Plmt, deg1StdPlmt, dates] = GRACEDEG1(Pcenter, Rlevel)
%   [deg1Plmt, deg1StdPlmt, dates] = GRACEDEG1({Pcenter, Rlevel})
%   [deg1Plmt, deg1StdPlmt, dates] = GRACEDEG1(__, 'Name', Value)
%
% Input arguments
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data centre is 'CSR'.
%       When the first argument is a cell array, it is interpreted as
%       {Pcenter, Rlevel}.
%       Data type: char
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06' (or numbers).
%       The default release level is 'RL06'.
%       Currently, only RL06 is guaranteed to work.
%       Data type: char | [numeric]
%   TimeFormat - Format of the output time
%       - 'datetime': Matlab DATETIME object
%       - 'datenum': Matlab double in DATENUM format
%       The default option is 'datetime'.
%       Data type: char
%   OutputFormat - Format of the output coefficients
%       - 'timefirst': Coefficients are ordered as [time, lmcosi]
%       - 'traditional': Coefficients are ordered as [lmcosi, time]
%       The default option is 'timefirst' (to be consistent with
%       GRACE2PLMT).
%       Data type: char
%	ForceNew - Logical flag to force reprocess of the data
%		- true: Reprocess the data
%		- false: Only reprocess if previous output is not found
%		The default option is 'soft-true', only reprocessing the data if
%       the output file is older than the input files. This option is not
%       explicitly exposed; to enforce this option, set it to 0.5.
%		Data types: logical | ([numeric])
%	SaveData - Logical flag to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%		The default option is true.
%		Data types: logical | ([numeric])
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is false.
%		Data types: logical | ([numeric])
%
% Output arguments
%   deg1Plmt - Potential coefficients for just degree l = 1
%       Data type: DOUBLE
%       Dimension: [nmonths x 2 x 4] | [2 x 4 x nmonths]
%           (depending on the OutputFormat input)
%   deg1StdPlmt - Standard deviation of the potential coefficients
%       The data type and dimension are the same as for deg1Plmt.
%   dates - Time stamps of the coefficients
%       The time stamps are the midpoints of the time intervals in the
%       techinical note.
%       Datatype: DATATIME | DOUBLE
%           (depending on the TimeFormat input)
%       Dimension: [nmonths x 1]
%
% See also
%   GRACE2PLMT (GRACE2PLMT_NEW), SOLVEDEG1
%
% Data sources
%   For CSR RL 6.3, the GRACE/GRACE-FO Technical Note 13 is available at
%       https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt
%   For GFZ RL 6.3, the GRACE/GRACE-FO Technical Note 13 is available at
%       https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_GFZ_RL0603.txt
%   For JPL RL 6.1, the GRACE/GRACE-FO Technical Note 13 is available at
%       https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_JPL_RL0601.txt
%       More recent versions of the JPL may be available after the time of this documentation.
%
% Last modified by
%   2025/05/20, williameclee@arizona.edu (@williameclee)
%   2021/11/18, charig@arizona.edu

function varargout = gracedeg1(varargin)
    %% Initialisation
    % Parse inputs
    [Pcenter, Rlevel, timeFmt, outputFmt, forceNew, saveData, beQuiet] = ...
        parseinputs(varargin{:});

    % Find the coefficient files
    [inputPath, outputPath] = getIOfiles(Pcenter, Rlevel);

    % Check if the output file already exists and is newer than the input
    if exist(outputPath, 'file') && ...
            (forceNew == 0 || ...
            (forceNew == 1 && isolder(inputPath, outputPath, false)))
        load(outputPath, 'deg1Plmt', 'deg1StdPlmt', 'dates')

        if exist('deg1Plmt', 'var') && ...
                exist('deg1StdPlmt', 'var') && ...
                exist('dates', 'var')

            if ~beQuiet
                fprintf('%s loaded %s\n', upper(mfilename), outputPath)
            end

            [deg1Plmt, deg1StdPlmt, dates] = ...
                formatoutput(deg1Plmt, deg1StdPlmt, dates, timeFmt, outputFmt);

            varargout = {deg1Plmt, deg1StdPlmt, dates};
            return
        end

    end

    %% Reading the data
    if ~exist(inputPath, 'file') || forceNew == 2
        [inputPath, outputPath] = ...
            downloadRemoteData(inputPath, outputPath, Pcenter, Rlevel, beQuiet);
    end

    data = fileread(inputPath);
    dataBarrier = 'end of header ===============================================================================';
    data = strsplit(data, dataBarrier);
    data = data{2};
    data = textscan(data, '%s%d%d%f%f%f%f%s%s%d');

    %% Recasting data
    % Load the date
    deg1data = [double(data{2}), double(data{3}), data{4}, data{5}, data{6}, data{7}];
    inputStartDates = datetime(data{8}, ...
        "InputFormat", 'yyyyMMdd.hhmm', "Format", 'yyyy/MM/dd hh:mm');
    inputEndDates = datetime(data{9}, ...
        "InputFormat", 'yyyyMMdd.hhmm', "Format", 'yyyy/MM/dd hh:mm');
    inputMeanDates = mean([inputStartDates, inputEndDates], 2);
    [dates, ~, inputDateId] = unique(inputMeanDates, 'sorted');

    % Recast the data to lmcosi format
    deg1Plmt = zeros([length(dates), 2, 4]);
    deg1StdPlmt = zeros([length(dates), 2, 4]);

    for iData = 1:length(deg1data)
        deg1Plmt(inputDateId(iData), 1 + deg1data(iData, 2), :) = ...
            deg1data(iData, [1, 2, 3, 4]);
        deg1StdPlmt(inputDateId(iData), 1 + deg1data(iData, 2), :) = ...
            deg1data(iData, [1, 2, 5, 6]);
    end

    % Create a save file
    if saveData
        save(outputPath, 'deg1Plmt', 'deg1StdPlmt', 'dates')

        if ~beQuiet
            fprintf('%s saved %s\n', upper(mfilename), outputPath)
        end

    end

    % Collect output
    [deg1Plmt, deg1StdPlmt, dates] = ...
        formatoutput(deg1Plmt, deg1StdPlmt, dates, timeFmt, outputFmt);
    varargout = {deg1Plmt, deg1StdPlmt, dates};
end

%% Subfunctions
% Parse inputs
function varargout = parseinputs(varargin)
    ip = inputParser;
    addOptional(ip, 'Pcenter', 'CSR', ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', 'RL06', ...
        @(x) (isscalar(x) && isnumeric(x) && x >= 4 && x < 7) || ...
        ischar(validatestring(x, {'RL04', 'RL05', 'RL06'})));
    addParameter(ip, 'TimeFormat', 'datetime', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'OutputFormat', 'timefirst', ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'ForceNew', 0.5, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));

    if ~isempty(varargin) && iscell(varargin{1})
        varargin = [varargin{1}{1:2}, varargin(2:end)];
    end

    parse(ip, varargin{:});
    Pcenter = ip.Results.Pcenter;
    Rlevel = ip.Results.Rlevel;
    timeFmt = ip.Results.TimeFormat;
    outputFmt = ip.Results.OutputFormat;
    forceNew = uint8(double(ip.Results.ForceNew) * 2);
    saveData = logical(ip.Results.SaveData);
    beQuiet = logical(ip.Results.BeQuiet);

    if isnumeric(Rlevel)
        Rlevel = sprintf('RL%02d', floor(Rlevel));
    end

    varargout = {Pcenter, Rlevel, timeFmt, outputFmt, forceNew, saveData, beQuiet};
end

% Get the input and output file names
function [inputPath, outputPath] = getIOfiles(Pcenter, Rlevel)
    % Where the original data files are kept
    if exist(getenv('GRACEDATA'), 'dir')
        dataFolder = fullfile(getenv('GRACEDATA'), 'Degree1');
    else
        dataFolder = fullfile(getenv('IFILES'), 'GRACE', 'Degree1');
    end

    if ~exist(dataFolder, 'dir')
        error(sprintf('%s:DataFolderNotFound', upper(mfilename)), ...
            'Data folder %s does not exist', dataFolder)
    end

    % The name of the original data files
    switch Rlevel
        case 'RL06'
            inputPattern = fullfile(dataFolder, ...
                sprintf('TN-13_GEOC_%s_%s.txt', Pcenter, 'RL06*'));

            try
                inputPath = fullfile(dir(inputPattern).folder, dir(inputPattern).name);
            catch
                inputPath = fullfile(dataFolder, 'placeholder.txt');
            end

        otherwise
            inputPath = fullfile(dataFolder, ...
                sprintf('deg1_%s_NH.txt', Rlevel));
    end

    outputPath = replace(inputPath, '.txt', '.mat');

end

% Format outputs
function [deg1Plmt, deg1StdPlmt, dates] = ...
        formatoutput(deg1Plmt, deg1StdPlmt, dates, timeFmt, outputFmt)

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    if strcmp(outputFmt, 'traditional')
        deg1Plmt = permute(deg1Plmt, [2, 3, 1]);
        deg1StdPlmt = permute(deg1StdPlmt, [2, 3, 1]);
    end

end

% Fetch remote data
function [inputPath, outputPath] = ...
        downloadRemoteData(inputPath, outputPath, Pcenter, Rlevel, beQuiet)

    if ~strcmp(Rlevel, 'RL06')

        if exist(inputPath, 'file')
            return
        end

        error(sprintf('%s:RemoteDataNotAvailable', upper(mfilename)), ...
            'Remote data for earlier releases (%s) not available', ...
            Rlevel)

    end

    switch Pcenter
        case 'CSR'
            remoteInputPath = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt';
        case 'GFZ'
            remoteInputPath = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_GFZ_RL0603.txt';
        case 'JPL'
            remoteInputPath = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_JPL_RL0601.txt';
    end

    inputFolder = fileparts(inputPath);

    if ~exist(inputFolder, 'dir')
        mkdir(inputFolder)
    end

    [~, inputFile, inputFileExt] = fileparts(remoteInputPath);
    inputFile = [inputFile, inputFileExt];
    inputPath = fullfile(inputFolder, inputFile);
    outputPath = replace(inputPath, '.txt', '.mat');

    try

        if ~beQuiet
            loadingMsg = sprintf('%s downloading %s\nMore recent releases might be available\n', upper(mfilename), remoteInputPath);
            fprintf(loadingMsg)
        end

        websave(inputPath, remoteInputPath);

        if ~beQuiet
            fprintf(repmat('\b', 1, length(loadingMsg)))
            fprintf('%s downloaded %s\nMore recent releases might be available\n', upper(mfilename), remoteInputPath)
        end

    catch

        if ~exist(inputPath, 'file')
            error(sprintf('%s:InputDataNotFound', upper(mfilename)), ...
                'Input file not found at %s\nAttempted to download from %s but failed', ...
                inputPath, remoteInputPath)
        end

    end

end
