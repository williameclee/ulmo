%% AOD1B2PLMT
% Fetch the AOD1B modelled GRACE geoid (or surface mass density) coefficient time series
%
% Syntax
%   [plmt, stdPlmt, dates] = aod1b2plmt(Pcenter, Rlevel, product, L, unit)
%   [plmt, stdPlmt, dates] = aod1b2plmt(__, 'Name', Value)
%
% Inputs
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data centre is 'CSR'.
%       When the first argument is a cell array, it is interpreted as
%       {Pcenter, Rlevel, Ldata}.
%       Data type: char
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06' (or numbers).
%       The default release level is 'RL06'.
%       Currently, only RL06 is guaranteed to work.
%       Data type: char | [numeric]
%   Product - Type of product
%       - 'GAC'
%       - 'GAD'
%       The default product is 'GAC'.
%       Data type: char
%   L - The bandwidth of the output product
%       The default L is 60.
%   unit - Unit of the output
%       - 'POT': Geopotential field.
%       - 'SD': Surface mass density.
%       The default field is 'SD'.
%   TimeRange - Time range of the output
%       When specified, the output will be truncated to the specified time
%       range. If the input has two elements, it is interpreted as the
%       start and end of the time range; if the input has more than two
%       elements, it is interpreted as a list of time stamps.
%       The default time range is [] (all available data).
%       Data type: datetime | ([numeric])
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
%   ForceNew - Logical flag to force reprocess of the data
%       The default option is false.
%       Data type: logical | ([numeric])
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
% Outputs
%   plmt - SH coefficients of the gravity field
%       Units: m/s^2 (POT) | kg/m^2 (SD)
%       Data type: DOUBLE
%       Dimension: [nmonths x addmup(L) x 4] | [addmup(L) x 4 x nmonths]
%           (depending on the OutputFormat input)
%   stdPlmt - Standard deviation of the SH coefficients
%       The data type and dimension are the same as for plmt.
%   dates - Time stamps of the coefficients
%       The time stamps are the midpoints of the time intervals of the
%       input data.
%       Datatype: DATATIME | DOUBLE
%           (depending on the TimeFormat input)
%       Dimension: [nmonths x 1]
%
% See also
%   GRACEDEG1, GRACEDEG2, PLM2POT, GRACE2PLMT (GRACE2PLMT_NEW)
%
% Last modified by
%   2025/10/16, williameclee@arizona.edu (@williameclee)
%   2014/02/27, charig@princeton.edu
%   2011/05/17, fjsimons@alum.mit.edu

function varargout = aod1b2plmt(varargin)
    %% Initialisation
    % Parse inputs
    [Pcenter, Rlevel, product, Loutput, unit, timelim, ...
         outputFmt, timeFmt, ~, ~, forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin{:});

    % If this file already exists, load it.  Otherwise, or if we force it, make
    % a new one (e.g. you added extra months to the database).
    [inputFolder, outputPath] = getIOpaths(Rlevel, Pcenter, product, unit);

    vars = {'aod1bPlmt', 'aod1bStdPlmt', 'dates'};

    if ~forceNew && exist(outputPath, 'file') && all(ismember(vars, who('-file', outputPath)))
        % Load the peocessed data
        data = load(outputPath, 'aod1bPlmt', 'aod1bStdPlmt', 'dates');

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath, upper(product));
        end

        % Collect output
        [aod1bPlmt, aod1bStdPlmt, dates] = ...
            formatoutput(data.aod1bPlmt, data.aod1bStdPlmt, data.dates, Loutput, timelim, outputFmt, timeFmt);

        varargout = {aod1bPlmt, aod1bStdPlmt, dates};
        return

    end

    % Reload from the raw data
    [aod1bPlmt, aod1bStdPlmt, dates, equatorRadius, gravityParam] = ...
        aod1b2plmtCore(inputFolder, product, Rlevel, unit);

    % Save
    if saveData
        save(outputPath, ...
            'aod1bPlmt', 'aod1bStdPlmt', 'dates', 'equatorRadius', 'gravityParam');

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath, upper(product));
        end

    end

    % Collect output
    [aod1bPlmt, aod1bStdPlmt, dates] = ...
        formatoutput(aod1bPlmt, aod1bStdPlmt, dates, Loutput, timelim, outputFmt, timeFmt);

    varargout = {aod1bPlmt, aod1bStdPlmt, dates};
end

%% Subfunctions
% Heart of the programme
function [aod1bPlmt, aod1bStdPlmt, dates, equatorRadius, gravityParam] = ...
        aod1b2plmtCore(inputFolder, product, Rlevel, unit)
    inputFiles = ...
        ls2cell(fullfile(inputFolder, sprintf('%s-2_*_%s*', product, Rlevel(3:4))));
    Ldata = 180;

    % Preallocate
    nDates = length(inputFiles);
    dates = NaT([nDates, 1]);
    aod1bPlmt = nan([nDates, addmup(Ldata), 4]);
    aod1bStdPlmt = nan([nDates, addmup(Ldata), 4]);

    % Loop over the months
    wbar = waitbar(0, sprintf('Reading %s data', product), ...
        "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

    for iDate = 1:nDates
        waitbar(iDate / nDates, wbar, ...
            sprintf('Reading %s data (%d/%d)', product, iDate, nDates));

        if getappdata(wbar, 'canceling')
            delete(wbar);
            warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
            'Processing cancelled');
            fclose(logFid);
            return
        end

        % load gravity coefficients
        inputPath = fullfile(inputFolder, inputFiles{iDate});
        [gracePlm, graceStdPlm, date] = parsegracesourcefile(inputPath);
        aod1bPlmt(iDate, :, :) = gracePlm;
        aod1bStdPlmt(iDate, :, :) = graceStdPlm;
        dates(iDate) = date;
    end

    delete(wbar);

    %% Converting unit
    [~, ~, ~, ~, gravityParam, equatorRadius] = ...
        parsegracesourcefile(inputPath);

    if strcmp(unit, 'POT')
        return
    end

    % Convert gravity to surface geopotential
    aod1bPlmt(:, :, 3:4) = aod1bPlmt(:, :, 3:4) * equatorRadius;
    aod1bStdPlmt(:, :, 3:4) = aod1bStdPlmt(:, :, 3:4) * equatorRadius;
    aod1bPlmt = ...
        plm2pot(aod1bPlmt, equatorRadius, gravityParam, equatorRadius, 4);
    aod1bStdPlmt = ...
        plm2pot(aod1bStdPlmt, equatorRadius, gravityParam, equatorRadius, 4);

end

% Parse input arguments
function varargout = parseinputs(varargin)
    ip = inputParser;
    addOptional(ip, 'Pcenter', 'CSR', ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', 'RL06', ...
        @(x) (isnumeric(x) && isscalar(x)) || ...
        (ischar(validatestring(x, {'RL04', 'RL05', 'RL06'}))));
    addOptional(ip, 'Product', 'GAC', ...
        @(x) ischar(validatestring(x, {'GAC', 'GAD'})));
    addOptional(ip, 'L', [], ...
        @(x) (isnumeric(x) && x > 0) || isempty(x));
    addOptional(ip, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'POT', 'SD'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) (isdatetime(x) || isnumeric(x)) || isempty(x));
    addParameter(ip, 'OutputFormat', 'timefirst', ... % to be consistent with grace2plmt_new
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'TimeFormat', 'datenum', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'RemoveOceanMean', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'OceanDomain', GeoDomain('alloceans', "Buffer", 0.5), ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'CallChain', {}, @iscell);

    if iscell(varargin{1})
        varargin = [varargin{1}{1:2}, varargin(2:end)];
    end

    parse(ip, varargin{:});

    Pcenter = ip.Results.Pcenter;
    Rlevel = ip.Results.Rlevel;
    product = ip.Results.Product;
    Loutput = round(ip.Results.L);
    unit = ip.Results.Unit;
    timelim = ip.Results.TimeRange;
    outputFmt = ip.Results.OutputFormat;
    timeFmt = ip.Results.TimeFormat;
    rmOceanMean = logical(ip.Results.RemoveOceanMean);
    oceanDomain = ip.Results.OceanDomain;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);
    callChain = [ip.Results.CallChain, {mfilename}];

    if isnumeric(Rlevel)
        Rlevel = sprintf('RL%02d', floor(Rlevel));
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    varargout = ...
        {Pcenter, Rlevel, product, Loutput, unit, timelim, ...
         outputFmt, timeFmt, rmOceanMean, oceanDomain, forceNew, saveData, beQuiet, callChain};
end

% Get the input folder and output file names
function [inputFolder, outputPath] = ...
        getIOpaths(Rlevel, Pcenter, product, unit)

    if ~isempty(getenv('GRACEDATA'))
        outputFolder = fullfile(getenv('GRACEDATA'));
    else
        outputFolder = fullfile(getenv('IFILES'), 'GRACE');
    end

    % And the name of that save file
    if strcmp(unit, 'SD')
        outputPath = fullfile(outputFolder, ...
            sprintf('%s_%s_%s_SD.mat', Pcenter, Rlevel, product));
    else
        outputPath = fullfile(outputFolder, ...
            sprintf('%s_%s_%s.mat', Pcenter, Rlevel, product));
    end

    if ~isempty(getenv('ORIGINALGRACEDATA'))
        inputFolder = fullfile(getenv('ORIGINALGRACEDATA'), ...
            Rlevel, Pcenter);
    elseif ~isempty(getenv('GRACEDATA'))
        inputFolder = fullfile(getenv('GRACEDATA'), ...
            'raw', Rlevel, Pcenter);
    else
        inputFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'raw', Rlevel, Pcenter);
    end

end

% Format the output
function [aod1bPlmt, aod1bStdPlmt, dates] = ...
        formatoutput(aod1bPlmt, aod1bStdPlmt, dates, Loutput, timelim, ...
        outputFmt, timeFmt)

    if ~isempty(timelim)

        isValidTime = false(size(dates));

        if length(timelim) == 2
            isValidTime = dates >= timelim(1) & dates <= timelim(2);
        end

        if length(timelim) ~= 2 || ~any(isValidTime)
            [~, isValidTime] = ismember(timelim, dates);
            isValidTime = isValidTime > 0;
        end

        dates = dates(isValidTime);
        aod1bPlmt = aod1bPlmt(isValidTime, :, :);
        aod1bStdPlmt = aod1bStdPlmt(isValidTime, :, :);

    end

    if ~isempty(Loutput)
        aod1bPlmt = aod1bPlmt(:, 1:addmup(Loutput), :);
        aod1bStdPlmt = aod1bStdPlmt(:, 1:addmup(Loutput), :);
    end

    if strcmp(outputFmt, 'traditional')
        aod1bPlmt = permute(aod1bPlmt, [2, 3, 1]);
        aod1bStdPlmt = permute(aod1bStdPlmt, [2, 3, 1]);
    end

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

end
