%% GRACEDEG2
% Reads and formats the C20 and C30 spherical harmonic corrections from GRACE/GRACE-FO Technical Note 14.
%
% Syntax
%   [c2030data, c2030std, dates] = gracedeg2(Rlevel)
%   [c2030data, c2030std, dates] = gracedeg2(__, 'Name', Value)
%
% Input arguments
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06' (or numbers).
%       The default release level is 'RL06'.
%       Currently, only RL06 is guaranteed to work.
%       When the first argument is a cell array, it is interpreted as
%       {Pcenter, Rlevel}, and Pcenter is ignored.
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
% Outputs and saved variables
%   c2030data - C20 and C30 coefficients from SLR
%       Data type: DOUBLE
%       Dimension: [nmonths x 2] | [2 x nmonths]
%           (depending on the OutputFormat input)
%   c2030std - Standard deviation of the SLR coefficients
%       The data type and dimension are the same as for c2030data.
%   dates - Time stamps of the coefficients
%       The time stamps are the midpoints of the time intervals in the
%       techinical note.
%       Datatype: DATATIME | DOUBLE
%           (depending on the TimeFormat input)
%       Dimension: [nmonths x 1]
%
% Data sources
%   For release level 6, the GRACE/GRACE-FO Technical Note 14 is available at
%       https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt
%
% Last modified by
%   2025/05/20, williameclee@arizona.edu (@williameclee)
%   2021/11/11, charig@arizona.edu

function varargout = gracedeg2(varargin)
    %% Initialisation
    % Parse inputs
    [Rlevel, timeFmt, outputFmt, forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin{:});

    [inputPath, outputPath] = getIOfiles(Rlevel);

    vars = {'c2030data', 'c2030std', 'dates'};

    if exist(outputPath, 'file') && ...
            (forceNew == 0 || ...
            forceNew == 1 && isolder(inputPath, outputPath, false)) && ...
            all(ismember(vars, who('-file', outputPath)))
        data = load(outputPath, 'c2030data', 'c2030std', 'dates');

        if ~beQuiet
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">C20/C30 data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath);
        end

        [data.c2030data, data.c2030std, data.dates] = ...
            formatoutput(data.c2030data, data.c2030std, data.dates, timeFmt, outputFmt);
        varargout = {data.c2030data, data.c2030std, data.dates};
        return
    end

    %% Processing raw data
    if ~exist(inputPath, 'file') || forceNew == 2
        downloadRemoteData(inputPath, beQuiet, callChain);
    end

    data = fileread(inputPath);
    dataBarrier = 'Product:\n';
    data = strsplit(data, dataBarrier);
    data = data{2};
    data = textscan(data, repmat('%f', [1, 10]));

    startDates = datetime(data{1}, ...
        "ConvertFrom", 'modifiedjuliandate', "Format", 'yyyy/MM/dd hh:mm');
    endDates = datetime(data{9}, ...
        "ConvertFrom", 'modifiedjuliandate', "Format", 'yyyy/MM/dd hh:mm');
    dates = mean([startDates, endDates], 2);

    c2030data = [data{3}, data{6}];
    c2030std = [data{5}, data{8}] * 1e-10;

    % Create a save file
    if saveData
        save(outputPath, vars{:}, '-v7.3')

        if ~beQuiet
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">processed C20/C30 data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath);
        end

    end

    % Collect output
    [c2030data, c2030std, dates] = ...
        formatoutput(c2030data, c2030std, dates, timeFmt, outputFmt);
    varargout = {c2030data, c2030std, dates};
end

%% Subfunctions
% Parse inputs
function varargout = parseinputs(varargin)
    ip = inputParser;
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
    addParameter(ip, 'CallChain', {}, @iscell);

    if ~isempty(varargin) && iscell(varargin{1})
        varargin = [varargin{1}{2}, varargin(2:end)];
    end

    parse(ip, varargin{:});
    Rlevel = ip.Results.Rlevel;
    timeFmt = ip.Results.TimeFormat;
    outputFmt = ip.Results.OutputFormat;
    forceNew = uint8(double(ip.Results.ForceNew) * 2);
    saveData = logical(ip.Results.SaveData);
    beQuiet = logical(ip.Results.BeQuiet);
    callChain = [ip.Results.CallChain, {mfilename}];

    if isnumeric(Rlevel)
        Rlevel = sprintf('RL%02d', floor(Rlevel));
    end

    varargout = {Rlevel, timeFmt, outputFmt, forceNew, saveData, beQuiet, callChain};
end

% Get the input and output file names
function [inputPath, outputPath] = getIOfiles(Rlevel)

    if ~isempty(getenv('GRACEDATA'))
        dataFolder = fullfile(getenv('GRACEDATA'), 'Degree2');
    else
        dataFolder = fullfile(getenv('IFILES'), 'GRACE', 'Degree2');
    end

    % The name of the original data files
    switch Rlevel
        case 'RL06'
            inputPath = fullfile(dataFolder, ...
            'TN-14_C30_C20_GSFC_SLR.txt');
        case 'RL05'
            inputPath = fullfile(dataFolder, 'deg1_RL05_NH.txt');
        case 'RL04'
            inputPath = fullfile(dataFolder, 'deg1_RL04_NH.txt');
    end

    outputPath = replace(inputPath, '.txt', '.mat');

end

% Fetch remote data
function downloadRemoteData(inputPath, beQuiet, callChain)
    remoteInputPath = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt';

    [~, inputFolder] = fileparts(inputPath);

    if ~exist(inputFolder, 'dir')
        mkdir(inputFolder)
    end

    try

        if ~beQuiet
            t = tic;
            templine = 'this may take a while...';
            fprintf('[ULMO>%s] Downloading TN-13, %s\n', ...
                callchaintext(callChain), templine);
        end

        websave(inputPath, remoteInputPath);

        if ~beQuiet
            fprintf(repmat('\b', 1, length(templine) + 1));
            fprintf('took %.1f seconds.\nMore recent releases may be available.\n', toc(t));
        end

    catch

        if ~exist(inputPath, 'file')
            error(sprintf('%s:InputDataNotFound', upper(mfilename)), ...
                'Input file not found at %s\nAttempted to download from %s but failed', ...
                inputPath, remoteInputPath)
        end

    end

end

% Format outputs
function [c2030data, c2030sigma, dates] = ...
        formatoutput(c2030data, c2030sigma, dates, timeFmt, outputFmt)

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    if strcmp(outputFmt, 'traditional')
        c2030data = c2030data';
        c2030sigma = c2030sigma';
    end

end
