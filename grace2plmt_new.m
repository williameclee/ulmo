%% GRACE2PLMT_NEW - Reads and corrects GRACE spherical harmonic time series.
% Reads CSR, GFZ and JPL Level-2 GSM products or CSR mascon GSU harmonics.
% GSM fields use a WGS84 reference; GSU fields retain GGM05C. Processed
% coefficients are cached before optional low-degree recomputation.
%
% Syntax
%   [gracePlmt, graceStdPlmt, dates] = grace2plmt_new()
%   [gracePlmt, graceStdPlmt, dates] = grace2plmt_new(Pcenter, Rlevel, Ldata, Unit)
%   [gracePlmt, graceStdPlmt, dates] = grace2plmt_new(product, "Name", value)
%   [gracePlmt, graceStdPlmt, dates] = grace2plmt_new(__, ...
%       "RecomputeDegree1", {"IncludeC21S21", true})
%
% Input arguments
%   Pcenter (optional) - GRACE data centre.
%       Options: 'CSR', 'GFZ', 'JPL' or 'CSR mascon'. Alternatively, pass
%       the cell array product = {Pcenter, Rlevel, Ldata} first.
%       Default centre: 'CSR'.
%   Rlevel (optional) - Release name or number.
%       RL05 and RL06 are supported; RL04 loading is not implemented.
%       Default release level: 'RL06'.
%   Ldata (optional) - Product bandwidth.
%       RL05 uses CSR=60 and GFZ/JPL=90; JPL months may stop at 60.
%       CSR mascon uses 720. These product limits override Ldata.
%       Default bandwidth: 60.
%   Unit (optional) - Output coefficient unit.
%       'SD': surface mass density in kg/m^2, equivalent to mm of fresh
%       water. 'GRAV' and 'POT' follow the CONVERTGRAVITY conventions.
%       Default unit: 'SD'.
%   TimeRange (optional) - Inclusive start and end dates for the output.
%       Accepts datetime or datenum values.
%       Size: (1, 2) or (2, 1).
%       Default range: [] (all dates).
%   Deg1Correction (name-value) - Whether to apply the degree-1 replacement series.
%       GSU uses CSR TN-13. RL05 rejects explicit true.
%       Default option: true for RL06; false for RL05.
%   C20Correction (name-value) - Whether to apply the SLR C20 replacement series.
%       GSU uses updates relative to GGM05C. RL05 rejects explicit true.
%       Default option: true for RL06; false for RL05.
%   C30Correction (name-value) - Whether to apply the SLR C30 replacement series.
%       GSU uses updates relative to GGM05C. RL05 rejects explicit true.
%       Default option: true for RL06; false for RL05.
%   RecomputeDegree1 (name-value) - Low-degree recomputation options.
%       false: retain the field after the selected replacement corrections.
%       true: use the preset solver arguments; a GIA model name selects
%       that model for the preset. A cell array passes explicit arguments
%       to SOLVEDEGREE1, including IncludeC20, IncludeC30 and IncludeC21S21.
%       All returned coefficients and uncertainties are inserted by harmonic
%       location and matched date; missing dates fail.
%       Default option: false.
%   OutputFormat (name-value) - Coefficient array layout.
%       'timefirst': time precedes lmcosi; 'traditional': time follows it.
%       Default format: 'timefirst'.
%   TimeFormat (name-value) - Output date format.
%       Options: 'datenum' or 'datetime'.
%       Default format: 'datenum'.
%   Loutput (name-value) - Output bandwidth.
%       Truncates or pads coefficients and uncertainties with NaN while
%       preserving valid degree/order labels. An empty value keeps Ldata.
%       Default bandwidth: [].
%   LoveNumSource (name-value) - Source of the load Love numbers.
%       Options: 'ISSM', 'Wahr' or 'ALMA3 <model>'.
%       Default source: 'ISSM'.
%   ForceNew (name-value) - Whether to use the processed GRACE cache.
%       false: reuse a valid cache; 0.5: refresh if source files are newer;
%       true: always reprocess.
%       Default refresh option: 0.5.
%   SaveData (name-value) - Whether to save processed GRACE data.
%       Default option: true.
%   BeQuiet (name-value) - Whether to report progress.
%       false: full reporting; 0.5: reduced reporting; true: silent.
%       Default reporting level: 0.5.
%
% Output arguments
%   gracePlmt - GRACE spherical harmonic coefficients.
%       The four columns are degree, order, cosine and sine coefficients.
%       Missing coefficients are NaN; degree/order labels remain valid.
%       Size: (N, M, 4) for 'timefirst', or (M, 4, N) for 'traditional',
%       where M = addmup(Lout); N is the number of epochs.
%       Unit: selected Unit for coefficient columns; 'SD' uses kg/m^2.
%   graceStdPlmt - Formal standard deviations.
%       Degree/order columns retain their labels.
%       Size: same as gracePlmt.
%       Unit: same as gracePlmt for coefficient columns.
%   dates - Midpoints of the input intervals in TimeFormat.
%       Size: (N, 1).
%
% Notes
%   See docs/grace-input-formats.md for source layouts and corrections.
%
% Data sources
%   NASA PO.DAAC GRACE data. https://podaac.jpl.nasa.gov/
%
% See also
%   GRACEDEG1, GRACEDEG2, SOLVEDEGREE1, CONVERTGRAVITY
%
% Last modified
%   2026/10/26, En-Chi Lee (williameclee@arizona.edu)

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
            fprintf('[ULMO>%s] Loaded %s.\n', ...
                callchaintext(callChain), filehref(outputPath, 'GRACE data'));
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
                fprintf('[ULMO>%s] Saved %s.\n', ...
                    callchaintext(callChain), filehref(outputPath, 'GRACE data'));
            end

        end

    end

    %% Post-processing and formatting
    % Recalculate degree 1 coefficients
    if iscell(redoDeg1)
        [myDeg1, myDeg1Std, myDeg1Dates, ~, coeffLocs] = ...
            solvedegree1(Pcenter, Rlevel, Ldata, redoDeg1{:}, "Unit", unit, ...
            "BeQuiet", beQuiet, "CallChain", callChain);

        [matched, dateIndex] = ismember(dates, myDeg1Dates);
        assert(all(matched), 'ULMO:MissingLowDegreeDates', ...
        'Joint coefficient estimates do not cover all GRACE dates.');

        for iCoeff = 1:size(coeffLocs, 1)
            row = coeffLocs(iCoeff, 1);
            col = 2 + coeffLocs(iCoeff, 2);
            gracePlmt(:, row, col) = myDeg1(dateIndex, iCoeff);
            graceStdPlmt(:, row, col) = myDeg1Std(dateIndex, iCoeff);
        end

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
    logCleanup = onCleanup(@() fclose(logFid));
    fprintf(logFid, '%s - Local time %s\n', ...
        upper(mfilename), datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
    fprintf(logFid, 'Data centre: %s, release level: %s\n', ...
        upper(Pcenter), Rlevel);

    % Find original data files if processed ones are not available
    [inputFiles, Ldata] = getinputfiles(Pcenter, Rlevel, Ldata, inputFolder);
    fprintf(logFid, '%d files to process\n', length(inputFiles));

    isMascon = strcmp(Pcenter, 'CSR mascon');

    % C20 and C30 correction setup
    if c20corr || c30corr
        [tnC2030, tnC2030Std, tnC2030dates] = ...
            gracedeg2(Rlevel, "BeQuiet", (uint8(beQuiet == 1) + beQuiet));

        if isMascon
            % GSU is relative to GGM05C, whereas TN-14 supplies full-field
            % coefficients. Subtract the zero-tide GGM05C C20/C30 values
            % (also listed as the means in TN-14); sigmas are unchanged.
            % https://download.csr.utexas.edu/pub/grace/GGM05/GGM05C.ICGEM
            tnC2030 = tnC2030 - [-4.841694573200e-4, 9.571647583412e-7];
        end

    end

    % Degree 1 correction setup
    % Special handling for CSR mascon
    if strcmp(Pcenter, 'CSR mascon')
        deg1Pcenter = 'CSR';
    else
        deg1Pcenter = Pcenter;
    end

    if deg1corr
        [tnDeg1, tnDeg1Std, tnDeg1dates] = ...
            gracedeg1(deg1Pcenter, Rlevel, "BeQuiet", (uint8(beQuiet == 1) + beQuiet));
    end

    % Preallocation
    nDates = length(inputFiles);
    dates = NaT([nDates, 1]);
    gracePlmt = emptyplmt(nDates, Ldata);
    graceStdPlmt = gracePlmt;
    deg1corrFailedDates = [];
    c20corrFailedDates = [];
    c30corrFailedDates = [];

    %% Loop over the months
    if ~beQuiet
        wbar = waitbar(0, 'Reading GRACE data', ...
            "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');
        cleanup = onCleanup(@() deleteWaitbar(wbar));
    end

    for iDate = 1:nDates

        if ~beQuiet
            waitbar(iDate / nDates, wbar, ...
                sprintf('Reading GRACE data (%d/%d)', iDate, nDates));

            if getappdata(wbar, 'canceling')
                error('ULMO:ProcessCancelledByUser', 'GRACE input file reading cancelled');
            end

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

        if size(gracePlm, 1) > addmup(Ldata)
            error('ULMO:grace2plmt:UnexpectedDegree', ...
                'Input %s exceeds the expected degree %d', inputPath, Ldata);
        end

        % Keep missing higher-degree coefficients/sigmas as NaN, but retain
        % valid degree/order labels for conversion and downstream indexing.
        gracePlmt(iDate, 1:size(gracePlm, 1), 1:size(gracePlm, 2)) = gracePlm;
        graceStdPlmt(iDate, 1:size(graceStdPlm, 1), 1:size(graceStdPlm, 2)) = graceStdPlm;
    end

    if ~beQuiet

        if ~isempty(deg1corrFailedDates)
            warning('ULMO:grace2plmt:NoTNReplacementAvailable', ...
                'No degree 1 replacement for the following dates:\n%s', ...
                strjoin(string(deg1corrFailedDates), ', '));
        end

        if ~isempty(c20corrFailedDates)
            warning('ULMO:grace2plmt:NoTNReplacementAvailable', ...
                'No C20 replacement for the following dates:\n%s', ...
                strjoin(string(c20corrFailedDates), ', '));
        end

        if ~isempty(c30corrFailedDates)
            warning('ULMO:grace2plmt:NoTNReplacementAvailable', ...
                'No C30 replacement for the following dates (likely because C30 replacement is only available for GFO):\n%s', ...
                strjoin(string(c30corrFailedDates), ', '));
        end

    end

    % WGS84 reference setup
    % For now just hardcode the even zonal coefficients (J), later use
    % Frederik's GRS.m program, don't bother with the higher degrees
    if ~isMascon
        wgs84C20 = 0.108262982131e-2 * -1 / sqrt(5); % will be row 4
        wgs84C40 = -0.237091120053e-5 * -1 / sqrt(5); % will be row 11
        gracePlmt(:, 4, 3) = gracePlmt(:, 4, 3) - wgs84C20;
        gracePlmt(:, 11, 3) = gracePlmt(:, 11, 3) - wgs84C40;
    end

    % Mascon GSU coefficients already describe updates to GGM05C; applying
    % the full-field WGS84 subtraction would introduce static zonal offsets.

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
        @(x) ischar(validatestring(x, {'CSR', 'CSR mascon', 'GFZ', 'JPL'})));
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

    Pcenter = validatestring(Pcenter, {'CSR', 'CSR mascon', 'GFZ', 'JPL'});
    Rlevel = validatestring(Rlevel, {'RL04', 'RL05', 'RL06'});

    % Product bandwidth controls the rectangular output and cache identity.
    if strcmpi(Pcenter, 'CSR mascon') && (Ldata ~= 720)
        warning("ULMO:grace2plmt:InvalidLdata", ...
            'CSR mascon solution only has L=720, overriding input Ldata value of %d.', Ldata);
        Ldata = 720;
    elseif strcmp(Pcenter, 'CSR') && strcmp(Rlevel, 'RL05') && Ldata ~= 60
        warning('ULMO:grace2plmt:InvalidLdata', ...
            'CSR RL05 input is degree 60; overriding Ldata=%d. Use Loutput to resize output.', Ldata);
        Ldata = 60;
    elseif (strcmpi(Pcenter, 'JPL') && strcmpi(Rlevel, 'RL05')) && (Ldata ~= 90)
        warning("ULMO:grace2plmt:InvalidLdata", ...
            'JPL RL05 contains degree 60/90 files; using Ldata=90 instead of %d. Use Loutput to truncate.', Ldata);
        Ldata = 90;
    elseif (strcmpi(Pcenter, 'GFZ') && strcmpi(Rlevel, 'RL05')) && (Ldata ~= 90)
        warning("ULMO:grace2plmt:InvalidLdata", ...
            'GFZ RL05 solution only has L=90, overriding input Ldata value of %d.', Ldata);
        Ldata = 90;
    end

    if ~strcmp(Rlevel, 'RL06')
        % RL05 replacement series are not implemented. Omitted/empty flags
        % default to false; reject explicit true instead of ignoring intent.
        correctionNames = {'Deg1Correction', 'C20Correction', 'C30Correction'};

        for k = 1:numel(correctionNames)
            name = correctionNames{k};

            if ~ismember(name, ip.UsingDefaults) && ...
                    ~isempty(ip.Results.(name)) && logical(ip.Results.(name))
                error('ULMO:grace2plmt:UnsupportedCorrection', ...
                    '%s is not supported for %s. Omit it or set it to false.', name, Rlevel);
            end

        end

        deg1correction = false;
        c20correction = false;
        c30correction = false;
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
            padded = emptyplmt(size(gracePlmt, 1), Loutput);
            paddedStd = padded;
            padded(:, 1:size(gracePlmt, 2), :) = gracePlmt;
            paddedStd(:, 1:size(graceStdPlmt, 2), :) = graceStdPlmt;
            gracePlmt = padded;
            graceStdPlmt = paddedStd;
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

% Allocate missing coefficients without losing spherical-harmonic indices.
function plmt = emptyplmt(nDates, L)
    [order, degree] = addmon(L);
    plmt = nan(nDates, addmup(L), 4);
    plmt(:, :, 1) = repmat(degree(:)', nDates, 1);
    plmt(:, :, 2) = repmat(order(:)', nDates, 1);
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
                replace(Pcenter, " ", ""), Rlevel, num2str(Ldata));
        case 'GRAV'
            outputFile = sprintf('%s_%s_alldata_%s_%s.mat', ...
                replace(Pcenter, " ", ""), Rlevel, num2str(Ldata), unit);
        case 'SD'
            outputFile = sprintf('%s_%s_%s_alldata_%s_%s.mat', ...
                replace(Pcenter, " ", ""), Rlevel, num2str(Ldata), unit, replace(loveNumSrc, " ", ""));
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

    if strcmpi(Rlevel, 'RL05')

        switch Pcenter
            case 'CSR'
                dataFiles = ls2cell(fullfile(inputFolder, 'GSM-2_2*-2*_*_UTCSR_0060_0005*'));
            case 'GFZ'
                dataFiles = ls2cell(fullfile(inputFolder, 'GSM-2_2*-2*_*_EIGEN_G---_005a'));
            case 'JPL'
                dataFiles = ls2cell(fullfile(inputFolder, 'GSM-2_2*-2*_*_JPLEM_0000_0005*'));
            otherwise
                error('ULMO:grace2plmt:LoadData:SolutionDNE', ...
                    '%s %s solutions not available', ...
                    upper(Pcenter), Rlevel);
        end

    elseif strcmpi(Rlevel, 'RL06')

        if strcmp(Pcenter, 'CSR mascon')
            dataFiles = ls2cell(fullfile(inputFolder, 'GSU-2_*_B---_06*'));

            if isempty(dataFiles)
                error('ULMO:grace2plmt_new:LoadData:NoRawGRACEDataFound', ...
                    'No data files found in %s', inputFolder)
            end

            return
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
                error("ULMO:grace2plmt:LoadData:SolutionDNE", ...
                    '%s degree %d solutions not available', ...
                    upper(Pcenter), Ldata);
        end

    else
        error('ULMO:grace2plmt:LoadData:SolutionLoadingNotImplemented', ...
            'Loading %s %s solutions is not currently implemented', ...
            upper(Pcenter), Rlevel);
    end

    % Make sure the files exist
    if isempty(dataFiles)
        error("ULMO:grace2plmt:LoadData:NoRawGRACEDataFound", ...
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
