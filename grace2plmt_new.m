%% GRACE2PLMT
% Reads in the Level-2 GRACE geoid products from either the CSR or GFZ data
% centres, does some processing, and saves them as a plmt matrix in a .mat
% file. In particular, the coefficients are reordered to our prefered
% lmcosi format, they are referenced to the WGS84 ellipsoid, the C20/C30
% coefficients are replaced with more accurate measurements from satellite
% laser ranging from Loomis et al, (2020), and the degree one coefficients
% are substituted with those from Sun et al., (2016).  You have the option
% of leaving them as geopotential or converting them to surface mass
% density using the method of Wahr et al. 1998, based on Love numbers.
%
% Syntax
%   [plmt, stdPlmt, dates] = grace2plmt(Pcenter, Rlevel, unit)
%   [plmt, stdPlmt, dates] = grace2plmt(__, 'Name', Value)
%
% Input arguments
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
%   Ldata - Bandwidth of the date product
%       In the case where there are more than one product from a data
%       centre (such as BA 60 or BB 96 standard L2 products) this allows
%       you to choose between them.
%       The default L is 60.
%       Data type: [numeric]
%   Unit - Unit of the output
%       - 'GRAV': Surface gravity (this is POT in SLEPIAN_ALPHA).
%       - 'POT': Geopotential field (surface gravity * equatorial radius).
%       - 'SD': Surface mass density.
%       The default field is 'SD'.
%       Data type: char
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
%   LoveNumSource - Source of load Love numbers
%       - 'ISSM': Love numbers from the Ice Sheet System Model (ISSM)
%       - 'Wahr': Love numbers as used in slepian_delta package
%       - 'ALMA3 <model>': Love numbers from the ALMA3 models (currently
%           has to be added manually)
%       The default source is 'ISSM'.
%   Deg1Correction, C20Correction, C30Correction - Logical flag to apply
%       these corrections
%       The default options are all true.
%       Data type: logical | ([numeric])
%   RecomputeDegree1 - Logical flag to recompute the degree 1 coefficients
%       - false: Do not recompute the degree 1 coefficients.
%       - true: Recompute the degree 1 coefficients using degree 60 input
%           data, SLE solved at degree 180, ICE-6G_D (VM5a) model is used
%           for GIA correction, and the ocean function has a buffer of
%           0.5°.
%       - Name of some GIA model (e.g. 'Caron18): Same as above, but using
%           the specified GIA model instead of ICE-6G_D (VM5a).
%       - Some cell array: These options are passed to the SOLVEDEGREE1
%           function.
%       The default option is false.
%       Data type: logical | char | cell array
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
%   plmt - SH coefficients of the gravity field
%       Units: m/s^2 (POT) | kg/m^2 (SD)
%       Data type: DOUBLE
%       Dimension: [ndates x addmup(L) x 4] | [addmup(L) x 4 x ndates]
%           (depending on the OutputFormat input)
%   stdPlmt - Standard deviation of the SH coefficients
%       The data type and dimension are the same as for plmt.
%   dates - Time stamps of the coefficients
%       The time stamps are the midpoints of the time intervals of the
%       input data.
%       Datatype: DATATIME | DOUBLE
%           (depending on the TimeFormat input)
%       Dimension: [ndates x 1]
%
% Data sources
%	GRACE data available from NASA PODAAC:
%	    https://podaac.jpl.nasa.gov/
%
% See also
%   GRACEDEG1, GRACEDEG2, PLM2POT, AOD1B2PLMT
%
% Notes
%   All the intermediate outputs originally printed on the screen are
%   printed to the log file instead.
%
% Last modified by
%   2026/01/29, williameclee@arizona.edu (@williameclee)
%   2022/05/18, charig@email.arizona.edu (@harig00)
%   2020/11/09, lashokkumar@arizona.edu
%   2019/03/18, mlubeck@email.arizona.edu
%   2011/05/17, fjsimons@alum.mit.edu (@fjsimons)

function [gracePlmt, graceStdPlmt, dates] = grace2plmt_new(varargin)
    %% Initialisation
    % Parse inputs
    [Pcenter, Rlevel, Ldata, unit, timelim, redoDeg1, ...
         deg1corr, c20corr, c30corr, outputFmt, timeFmt, Loutput, loveNumSrc, ...
         forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin{:});

    % Find the coefficient files
    [inputFolder, outputPath, logPath] = getIOpaths( ...
        Pcenter, Rlevel, Ldata, unit, loveNumSrc, deg1corr, c20corr, c30corr);

    % If this file already exists, load it.  Otherwise, or if we force it, make
    % a new one (e.g. you added extra months to the database).
    vars = {'gracePlmt', 'graceStdPlmt', 'dates'};

    if exist(outputPath, 'file') && ...
            (forceNew == 0 || ...
            (forceNew == 1 && isolder(inputFolder, outputPath, false))) && all(ismember(vars, who('-file', outputPath)))
        load(outputPath, vars{:})

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">GRACE data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath);
        end

    else
        % Reload from the raw data
        [gracePlmt, graceStdPlmt, dates, gravityParam, equatorRadius] = ...
            grace2plmtCore(Pcenter, Rlevel, Ldata, unit, loveNumSrc, inputFolder, ...
            deg1corr, c20corr, c30corr, beQuiet, logPath);

        if saveData
            save(outputPath, ...
                vars{:}, 'gravityParam', 'equatorRadius');

            if beQuiet <= 1
                fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">GRACE data</a>.\n', ...
                    callchaintext(callChain), outputPath, outputPath);
            end

        end

    end

    %% Post-processing and formatting
    % Recalculate degree 1 coefficients
    if iscell(redoDeg1)
        [myDeg1, myDeg1Std, myDeg1Dates] = ...
            solvedegree1(Pcenter, Rlevel, Ldata, redoDeg1{:}, "Unit", unit, ...
            "BeQuiet", beQuiet, "CallChain", callChain);

        if size(myDeg1, 1) > size(gracePlmt, 1)
            [~, isValidTime] = ismember(myDeg1Dates, dates);
            isValidTime = isValidTime > 0;
            myDeg1 = myDeg1(isValidTime, :);
            myDeg1Std = myDeg1Std(isValidTime, :);
        end

        gracePlmt(:, 2, 3) = myDeg1(:, 1);
        gracePlmt(:, 3, 3) = myDeg1(:, 2);
        gracePlmt(:, 3, 4) = myDeg1(:, 3);

        graceStdPlmt(:, 2, 3) = myDeg1Std(:, 1);
        graceStdPlmt(:, 3, 3) = myDeg1Std(:, 2);
        graceStdPlmt(:, 3, 4) = myDeg1Std(:, 3);
    end

    % Format output
    [gracePlmt, graceStdPlmt, dates] = ...
        formatoutput(gracePlmt, graceStdPlmt, dates, timelim, outputFmt, timeFmt, Loutput);
end

%% Subfunctions
% Heart of the programme
function [gracePlmt, graceStdPlmt, dates, gravityParam, equatorRadius] = ...
        grace2plmtCore(Pcenter, Rlevel, Ldata, unit, loveNumSrc, inputFolder, ...
        deg1corr, c20corr, c30corr, beQuiet, logPath)
    %% Computing the coefficients
    logFid = fopen(logPath, 'w');
    fprintf(logFid, '%s - Local time %s\n', ...
        upper(mfilename), datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
    fprintf(logFid, 'Data centre: %s, release level: %s\n', ...
        upper(Pcenter), Rlevel);

    % Find original data files if processed ones are not available
    [inputFiles, Ldata] = getinputfiles(Pcenter, Rlevel, Ldata, inputFolder);
    fprintf(logFid, '%d files to process\n', length(inputFiles));

    % C20 and C30 correction setup
    [tnC2030, tnC2030Std, tnC2030dates] = ...
        gracedeg2(Rlevel, "BeQuiet", (uint8(beQuiet == 1) + beQuiet));

    % Degree 1 correction setup
    [tnDeg1, tnDeg1Std, tnDeg1dates] = ...
        gracedeg1(Pcenter, Rlevel, "BeQuiet", (uint8(beQuiet == 1) + beQuiet));

    % Preallocation
    nDates = length(inputFiles);
    dates = NaT([nDates, 1]);
    gracePlmt = nan([nDates, addmup(Ldata), 4]); % l m cos sin
    graceStdPlmt = nan([nDates, addmup(Ldata), 4]); % l m cos sin
    deg1corrFailedDates = [];
    c20corrFailedDates = [];
    c30corrFailedDates = [];

    %% Loop over the months
    wbar = waitbar(0, 'Reading GRACE data', ...
        "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

    for iDate = 1:nDates
        waitbar(iDate / nDates, wbar, ...
            sprintf('Reading GRACE data (%d/%d)', iDate, nDates));

        if getappdata(wbar, 'canceling')
            delete(wbar);
            warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
            'Processing cancelled');
            fclose(logFid);
            return
        end

        % load gravity coefficients
        inputPath = fullfile(inputFolder, inputFiles{iDate});
        [gracePlm, graceStdPlm, date, dateRange] = parsegracesourcefile(inputPath);
        gracePlm(1, 3) = 0;
        dates(iDate) = date;

        fprintf(logFid, '%s - %s\n', date, inputFiles{iDate});

        if deg1corr
            [gracePlm, graceStdPlm, isSuccess] = ...
                replaceDeg1(gracePlm, graceStdPlm, date, dateRange, ...
                tnDeg1, tnDeg1Std, tnDeg1dates, logFid);

            if ~isSuccess
                deg1corrFailedDates = [deg1corrFailedDates; date]; %#ok<AGROW>
            end

        end

        if c20corr
            [gracePlm, graceStdPlm, isSuccess] = ...
                replaceC20(gracePlm, graceStdPlm, date, dateRange, ...
                tnC2030, tnC2030Std, tnC2030dates, logFid);

            if ~isSuccess
                c20corrFailedDates = [c20corrFailedDates; date]; %#ok<AGROW>
            end

        end

        % Now replace C3,0 if possible
        if c30corr
            [gracePlm, graceStdPlm, isSuccess] = ...
                replaceC30(gracePlm, graceStdPlm, date, dateRange, ...
                tnC2030, tnC2030Std, tnC2030dates, logFid);

            if ~isSuccess
                c30corrFailedDates = [c30corrFailedDates; date]; %#ok<AGROW>
            end

        end

        % Combine into one matrix
        gracePlmt(iDate, :, :) = gracePlm;
        graceStdPlmt(iDate, :, :) = graceStdPlm;
    end

    if ~beQuiet

        if ~isempty(deg1corrFailedDates)
            warning(sprintf('%s:NoTNReplacementAvailable', upper(mfilename)), ...
                'No degree 1 replacement for the following dates:\n%s', ...
                strjoin(string(deg1corrFailedDates), ', '));
        end

        if ~isempty(c20corrFailedDates)
            warning(sprintf('%s:NoTNReplacementAvailable', upper(mfilename)), ...
                'No C20 replacement for the following dates:\n%s', ...
                strjoin(string(c20corrFailedDates), ', '));
        end

        if ~isempty(c30corrFailedDates)
            warning(sprintf('%s:NoTNReplacementAvailable', upper(mfilename)), ...
                'No C30 replacement for the following dates (likely because C30 replacement is only available for GFO):\n%s', ...
                strjoin(string(c30corrFailedDates), ', '));
        end

    end

    % WGS84 reference setup
    % For now just hardcode the even zonal coefficients (J), later use
    % Frederik's GRS.m program, don't bother with the higher degrees
    wgs84C20 = 0.108262982131e-2 * -1 / sqrt(5); % will be row 4
    wgs84C40 = -0.237091120053e-5 * -1 / sqrt(5); % will be row 11
    gracePlmt(:, 4, 3) = gracePlmt(:, 4, 3) - wgs84C20;
    gracePlmt(:, 11, 3) = gracePlmt(:, 11, 3) - wgs84C40;

    fclose(logFid);
    delete(wbar);

    %% Converting unit
    % Use the actual parameters stored in the file instead of from FRALMANAC
    [~, ~, ~, ~, gravityParam, equatorRadius] = parsegracesourcefile(inputPath);

    % Convert gravity to the desired unit
    gracePlmt = convertgravity(gracePlmt, 'GRAV', unit, ...
        "EquatorRadius", equatorRadius, ...
        "InputFormat", 'lmcosi', "TimeDim", 'timefirst', "LoveNumSource", loveNumSrc);
    graceStdPlmt = convertgravity(graceStdPlmt, 'GRAV', unit, ...
        "EquatorRadius", equatorRadius, ...
        "InputFormat", 'lmcosi', "TimeDim", 'timefirst', "LoveNumSource", loveNumSrc);
end

% Parse input arguments
function varargout = parseinputs(varargin)
    % Have to define the default values and use CONDDEFVAL in this
    % complicated way because other functions (e.g. GRACE2SLEPT_NEW) may
    % pass in empty values to indicate that they want to use the default
    % values.
    dfOpt.Pcenter = 'CSR';
    dfOpt.Rlevel = 'RL06';
    dfOpt.Ldata = 60;
    dfOpt.Unit = 'SD';
    dfOpt.TimeRange = [];
    dfOpt.Deg1Correction = true;
    dfOpt.C20Correction = true;
    dfOpt.C30Correction = true;
    dfOpt.RecomputeDegree1 = false;
    dfOpt.OutputFormat = 'timefirst';
    dfOpt.TimeFormat = 'datenum';
    dfOpt.Loutput = [];
    dfOpt.LoveNumSrc = 'ISSM';
    dfOpt.ForceNew = 0.5;
    dfOpt.SaveData = true;
    dfOpt.BeQuiet = 0.5;

    ip = inputParser;
    addOptional(ip, 'Pcenter', dfOpt.Pcenter, ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', dfOpt.Rlevel, ...
        @(x) (isnumeric(x) && isscalar(x)) || ...
        (ischar(validatestring(x, {'RL04', 'RL05', 'RL06'}))));
    addOptional(ip, 'Ldata', dfOpt.Ldata, ...
        @(x) isnumeric(x) && x > 0);
    addOptional(ip, 'Unit', dfOpt.Unit, ...
        @(x) isempty(x) || ischar(validatestring(x, {'GRAV', 'POT', 'SD'})));
    addOptional(ip, 'TimeRange', dfOpt.TimeRange, ...
        @(x) ((isdatetime(x) || isnumeric(x)) && length(x) == 2) || isempty(x));
    addParameter(ip, 'Deg1Correction', dfOpt.Deg1Correction, ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'C20Correction', dfOpt.C20Correction, ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'C30Correction', dfOpt.C30Correction, ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'RecomputeDegree1', dfOpt.RecomputeDegree1, ...
        @(x) islogical(x) || iscell(x) || ischar(x) || isempty(x));
    addParameter(ip, 'OutputFormat', dfOpt.OutputFormat, ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'TimeFormat', dfOpt.TimeFormat, ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'Loutput', dfOpt.Loutput, ...
        @(x) (isscalar(x) && isnumeric(x) && x > 0) || isempty(x));
    addParameter(ip, 'LoveNumSource', dfOpt.LoveNumSrc, @(x) ischar(x) || isstring(x));
    addParameter(ip, 'ForceNew', dfOpt.ForceNew, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'SaveData', dfOpt.SaveData, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', dfOpt.BeQuiet, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'CallChain', {}, @iscell);

    if ~isempty(varargin) && iscell(varargin{1}) && numel(varargin{1}) == 3
        product = varargin{1};
        varargin = [cell(1, 3), varargin(2:end)];
        varargin{1} = product{1};
        varargin{2} = product{2};
        varargin{3} = product{3};
    end

    parse(ip, varargin{:});

    Pcenter = conddefval(ip.Results.Pcenter, dfOpt.Pcenter);
    Rlevel = conddefval(ip.Results.Rlevel, dfOpt.Rlevel);
    Ldata = conddefval(round(ip.Results.Ldata), dfOpt.Ldata);
    unit = conddefval(ip.Results.Unit, dfOpt.Unit);
    timelim = conddefval(ip.Results.TimeRange, dfOpt.TimeRange);
    deg1correction = ...
        conddefval(logical(ip.Results.Deg1Correction), dfOpt.Deg1Correction);
    c20correction = ...
        conddefval(logical(ip.Results.C20Correction), dfOpt.C20Correction);
    c30correction = ...
        conddefval(logical(ip.Results.C30Correction), dfOpt.C30Correction);
    redoDeg1 = ...
        conddefval(ip.Results.RecomputeDegree1, dfOpt.RecomputeDegree1);
    outputFmt = ...
        conddefval(ip.Results.OutputFormat, dfOpt.OutputFormat);
    timeFmt = conddefval(ip.Results.TimeFormat, dfOpt.TimeFormat);
    Loutput = conddefval(ip.Results.Loutput, dfOpt.Loutput);
    loveNumSource = ...
        char(conddefval(ip.Results.LoveNumSource, dfOpt.LoveNumSrc));
    forceNew = uint8(double(ip.Results.ForceNew) * 2);
    saveData = logical(ip.Results.SaveData);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);
    callChain = [ip.Results.CallChain, {mfilename}];

    if isnumeric(Rlevel)
        Rlevel = sprintf('RL%02d', floor(Rlevel));
    end

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    if islogical(redoDeg1) && redoDeg1
        redoDeg1 = {60, [], 180, 'ice6gd', GeoDomain('alloceans', "Buffer", 0.5)};
    elseif ischar(redoDeg1)
        redoDeg1 = {60, [], 180, redoDeg1, GeoDomain('alloceans', "Buffer", 0.5)};
    end

    varargout = ...
        {Pcenter, Rlevel, Ldata, unit, timelim, redoDeg1, ...
         deg1correction, c20correction, c30correction, outputFmt, timeFmt, Loutput, loveNumSource, ...
         forceNew, saveData, beQuiet, callChain};
end

% Format the output
function [gracePlmt, graceStdPlmt, dates] = ...
        formatoutput(gracePlmt, graceStdPlmt, dates, timelim, outputFmt, timeFmt, Loutput)

    if ~isempty(timelim)
        % Only keep the data within the specified time range
        isValidTime = dates >= timelim(1) & dates <= timelim(2);
        gracePlmt = gracePlmt(isValidTime, :, :);
        graceStdPlmt = graceStdPlmt(isValidTime, :, :);
        dates = dates(isValidTime);
    end

    if ~isempty(Loutput)

        if size(gracePlmt, 2) > addmup(Loutput)
            gracePlmt = gracePlmt(:, 1:addmup(Loutput), :);
            graceStdPlmt = graceStdPlmt(:, 1:addmup(Loutput), :);
        elseif size(gracePlmt, 2) < addmup(Loutput)
            [order, degree] = addmon(Loutput);
            gracePlmt(:, 1:addmup(Loutput), 1) = repmat(reshape(degree, 1, 1, []), [size(gracePlmt, 1), addmup(Loutput), 1]);
            gracePlmt(:, 1:addmup(Loutput), 2) = repmat(reshape(order, 1, 1, []), [size(gracePlmt, 1), addmup(Loutput), 1]);
        end

    end

    if strcmp(outputFmt, 'traditional')
        gracePlmt = permute(gracePlmt, [2, 3, 1]);
        graceStdPlmt = permute(graceStdPlmt, [2, 3, 1]);
    end

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

end

% Get the input folder and output file names
function [inputFolder, outputPath, logPath] = ...
        getIOpaths(Pcenter, Rlevel, Ldata, unit, loveNumSrc, deg1corr, c20corr, c30corr)

    if ~isempty(getenv('ORIGINALGRACEDATA'))
        inputFolder = fullfile(getenv('ORIGINALGRACEDATA'), ...
            Rlevel, Pcenter);
    elseif ~isempty(getenv('GRACEDATA'))
        inputFolder = fullfile(getenv('GRACEDATA'), 'raw', ...
            Rlevel, Pcenter);
    else
        inputFolder = fullfile(getenv('IFILES'), 'GRACE', 'raw', ...
            Rlevel, Pcenter);
    end

    % Where you would like to save the new .mat file
    if ~isempty(getenv('GRACEDATA'))
        outputFolder = fullfile(getenv('GRACEDATA'));
    else
        outputFolder = fullfile(getenv('IFILES'), 'GRACE');
    end

    switch unit % no otherwise case since input validity is already checked
        case 'POT'
            outputFile = sprintf('%s_%s_alldata_%s.mat', ...
                Pcenter, Rlevel, num2str(Ldata));
        case 'GRAV'
            outputFile = sprintf('%s_%s_alldata_%s_%s.mat', ...
                Pcenter, Rlevel, num2str(Ldata), unit);
        case 'SD'
            outputFile = sprintf('%s_%s_%s_alldata_%s_%s.mat', ...
                Pcenter, Rlevel, num2str(Ldata), unit, replace(loveNumSrc, " ", ""));
    end

    if ~c30corr
        outputFile = strrep(outputFile, 'alldata', 'alldata_nC30');
    end

    if ~c20corr
        outputFile = strrep(outputFile, 'alldata', 'alldata_nC20');
    end

    if ~deg1corr
        outputFile = strrep(outputFile, 'alldata', 'alldata_nDeg1');
    end

    outputPath = fullfile(outputFolder, outputFile);
    logPath = fullfile(outputFolder, sprintf('%s_log.txt', upper(mfilename)));

    if ~isfile(outputPath) == 2 && ~exist(inputFolder, 'dir') == 2
        error('The data you asked for are not currently stored\nPlease check the input folder %s', inputFolder)
    end

end

% Get the raw input files
function [dataFiles, Ldata] = ...
        getinputfiles(Pcenter, Rlevel, Ldata, inputFolder)
    % Only RL06 is supported for now
    if ~strcmp(Rlevel, 'RL06')
        error( ...
            sprintf('%s:LoadData:SolutionLoadingNotImplemented', upper(mfilename)), ...
            'Loading %s %s solutions is not currently implemented', ...
            upper(Pcenter), Rlevel);
    end

    % Get the data files
    switch Ldata
        case 60
            dataFiles = ls2cell(fullfile(inputFolder, ...
            'GSM-2_*_BA01_06*'));
        case 96
            dataFiles = ls2cell(fullfile(inputFolder, ...
            'GSM-2_*_BB01_06*'));
        otherwise
            error(sprintf('%s:LoadData:SolutionDNE', upper(mfilename)), ...
                '%s degree %d solutions not available', ...
                upper(Pcenter), Ldata);
    end

    % Make sure the files exist
    if isempty(dataFiles)
        error(sprintf('%s:LoadData:NoRawGRACEDataFound', upper(mfilename)), ...
            'No data files found in %s', inputFolder)
    end

end

% Do degree 1 correction
function [gracePlm, graceStdPlm, successFlag] = ...
        replaceDeg1(gracePlm, graceStdPlm, date, dateRange, tnDeg1, tnDeg1Std, tnDeg1dates, fid)
    successFlag = true;
    % Find a degree 1 data point that is within the month of our GRACE data
    iTn = tnDeg1dates(:) > dateRange(1) & tnDeg1dates(:) < dateRange(2) & ...
        ~any(isnan(tnDeg1(:, :, 3:4)), [2, 3]);

    if ~any(iTn)
        fprintf(fid, '  C10: (not replaced)   %+10.6e\n', ...
            gracePlm(2, 3));
        successFlag = false;
        return

    end

    if sum(iTn) > 1
        % We have more than one month. Use the closest value
        [~, iTn] = min(abs(date - tnDeg1dates));
    end

    fprintf(fid, '  C10: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(2, 3), tnDeg1(iTn, 1, 3));
    fprintf(fid, '  C11: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(3, 3), tnDeg1(iTn, 2, 3));
    fprintf(fid, '  S11: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(3, 4), tnDeg1(iTn, 2, 4));

    % Replacement
    gracePlm(2:3, 3:4) = squeeze(tnDeg1(iTn, :, 3:4));
    graceStdPlm(2:3, 3:4) = squeeze(tnDeg1Std(iTn, :, 3:4));
end

% Do C20 correction
function [gracePlm, graceStdPlm, successFlag] = ...
        replaceC20(gracePlm, graceStdPlm, date, dateRange, tnC2030, tnC2030Std, tnC2030dates, fid)
    successFlag = true;

    iTn = tnC2030dates > dateRange(1) & tnC2030dates < dateRange(2) & ...
        ~isnan(tnC2030(:, 1));

    if ~any(iTn)
        fprintf(fid, '  C20: (not replaced)   %+10.6e\n', ...
            gracePlm(4, 3));
        successFlag = false;
        return
    end

    if sum(iTn) > 1
        % We have more than one month. Use the closest value
        [~, iTn] = min(abs(date - tnC2030dates));
    end

    fprintf(fid, '  C20: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(4, 3), tnC2030(iTn, 1));

    % Replacement
    gracePlm(4, 3) = tnC2030(iTn, 1);
    graceStdPlm(4, 3) = tnC2030Std(iTn, 1);
end

% Do C20 correction
function [gracePlm, graceStdPlm, successFlag] = ...
        replaceC30(gracePlm, graceStdPlm, date, dateRange, tnC2030, tnC2030Std, tnC2030dates, fid)
    successFlag = true;
    iTn = tnC2030dates > dateRange(1) & tnC2030dates < dateRange(2) & ...
        ~isnan(tnC2030(:, 2));

    if ~any(iTn)
        fprintf(fid, '  C30: (not replaced)   %+10.6e\n', ...
            gracePlm(7, 3));
        successFlag = false;
        return
    end

    if sum(iTn) > 1
        % We have more than one month. Use the closest value
        [~, iTn] = min(abs(date - tnC2030dates));
    end

    fprintf(fid, '  C30: %s %+13.6e -> %+13.6e\n', ...
        date, gracePlm(7, 3), tnC2030(iTn, 2));

    % Replacement
    gracePlm(7, 3) = tnC2030(iTn, 2);
    graceStdPlm(7, 3) = tnC2030Std(iTn, 2);
end
