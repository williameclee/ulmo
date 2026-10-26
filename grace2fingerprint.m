%% GRACE2FINGERPRINT - Computes relative sea-level fingerprints from GRACE.
% Removes GIA from degrees >= 2, localises the land forcing, and solves
% the sea-level equation. Degree 1 and C21/S21 can be recomputed jointly.
%
% Syntax
%   [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint()
%   [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint(product, L, oceanDomain)
%   [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint(__, "Name", value)
%   [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint(__, "RecomputeC21S21", true)
%
% Input arguments
%   product (optional) - GRACE product specification.
%       Cell array {centre, release, degree}.
%       Size: (1, 3).
%       Default product: {'CSR', 'RL06', 60}.
%   L (optional) - SLE bandwidth.
%       Default bandwidth: 96.
%   oceanDomain (optional) - Ocean domain accepted by KERNELCP_NEW.
%       Default domain: GeoDomain('alloceans', "Buffer", 0.5).
%   GIA (name-value) - GIA model name as a character vector.
%       GIA correction is required; logical flags are not accepted.
%       Default model: 'ice6gd'.
%   RotationFeedback (name-value) - Include rotation feedback.
%       Default option: true.
%   RecomputeDegree1 (name-value) - Whether to re-estimate degree 1 with
%       SOLVEDEGREE1.
%       Default option: false.
%   RecomputeC21S21 (name-value) - Whether to jointly re-estimate degree 1,
%       C21 and S21.
%       Enables RecomputeDegree1 and passes IncludeC21S21 to SOLVEDEGREE1.
%       Default option: false.
%   LandDomain (name-value) - Domain used to localise the forcing.
%       NaN selects the complement of oceanDomain.
%       Default domain: NaN.
%   Frame (name-value) - Reference-frame label used in the cache filename.
%       Options: 'CF' or 'CM'. This option is not forwarded to SOLVESLE;
%       that call currently uses its default CF frame.
%       Default frame: 'CF'.
%   Truncation (name-value) - Output bandwidth, no greater than L.
%       An empty value uses L.
%       Default truncation: [].
%   Unit (name-value) - Output coefficient unit.
%       'SD': surface mass density in kg/m^2, equivalent to mm of fresh
%       water. 'GRAV' and 'POT' follow the CONVERTGRAVITY conventions.
%       Default unit: 'SD'.
%   TimeRange (name-value) - Inclusive start and end dates for the output.
%       Accepts datetime or datenum values.
%       Size: (1, 2) or (2, 1).
%       Default range: [] (all dates).
%   OutputFormat (name-value) - Coefficient array layout.
%       'timefirst': time precedes lmcosi; 'traditional': time follows it.
%       Default format: 'timefirst'.
%   TimeFormat (name-value) - Output date format.
%       Options: 'datenum' or 'datetime'.
%       Default format: 'datenum'.
%   OceanKernel (name-value) - Precomputed ocean-localisation matrix.
%       An empty value computes or loads the kernel.
%       Size: ((L+1)^2, (L+1)^2).
%       Default kernel: [].
%   ForceNew (name-value) - Whether to bypass the fingerprint cache.
%       Default option: false.
%   SaveData (name-value) - Whether to save the computed fingerprint.
%       Default option: true.
%   BeQuiet (name-value) - Whether to log messages.
%       true suppresses messages; 0.5 prints load/save messages.
%       Default reporting level: 0.5.
%
% Output arguments
%   rslLoadPlmt - Relative sea-level load coefficients.
%       The four columns are degree, order, cosine and sine coefficients.
%       Size: (N, M, 4) for 'timefirst', or (M, 4, N) for 'traditional',
%       Where M = addmup(output bandwidth); N is the number of dates.
%       Unit: selected Unit for coefficient columns; 'SD' uses kg/m^2.
%   rslLoadStdPlmt - Formal standard deviations.
%       Degree/order columns retain their labels.
%       Size: same as rslLoadPlmt.
%       Unit: same as rslLoadPlmt for coefficient columns.
%   dates - Epochs in the selected TimeFormat.
%       Size: (N, 1).
%
% Notes
%   Forcing remains GSM-like.
%   ForceNew applies to this fingerprint cache, not upstream solver caches.
%   Regenerate fingerprints when changing a forcing convention that is not
%   represented in the cache filename. GIA uncertainty is not propagated.
%
% See also
%   SOLVEDEGREE1, GRACE2PLMT_NEW, SOLVESLE, KERNELCP_NEW
%
% Created by
%   2024/11/20, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/26, En-Chi Lee (williameclee@arizona.edu)

function [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint(varargin)
    %% Parsing inputs
    [prod, ocnDom, lndDom, L, Lout, ...
         giaMdl, redoDeg1, redoCS21, timelim, doRot, frame, unit, ...
         ocnKernel, outputFmt, timeFmt, beQuiet, forceNew, saveData, callChain] = parseinputs(varargin);

    %% Data locating
    outPath = findoutputpath( ...
        prod, L, ocnDom, lndDom, frame, giaMdl, doRot, redoDeg1, redoCS21);

    vars = {'rslLoadPlmt', 'rslLoadStdPlmt', 'dates'};

    if ~forceNew && exist(outPath, 'file')
        data = load(outPath, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded %s.\n', ...
                callchaintext(callChain), filehref(outPath, 'fingerprint'));
        end

        [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
            formatoutput(data.rslLoadPlmt, data.rslLoadStdPlmt, data.dates, timelim, L, Lout, unit, ...
            outputFmt, timeFmt);

        return

    end

    %% Solving the sea level equation
    [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        grace2fingerprintCore(prod, ocnDom, lndDom, L, giaMdl, redoDeg1, redoCS21, doRot, ocnKernel, beQuiet);

    %% Saving and returning data
    if saveData
        save(outPath, vars{:}, '-v7.3');

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved %s.\n', ...
                callchaintext(callChain), filehref(outPath, 'fingerprint'));
        end

    end

    [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        formatoutput(rslLoadPlmt, rslLoadStdPlmt, dates, timelim, L, Lout, unit, ...
        outputFmt, timeFmt);

end

%% Subfunctions
% Heart of the programme
function [rslSdPlmt, rslSdStdPlmt, dates] = ...
        grace2fingerprintCore(prod, ocnDom, lndDom, Lsle, giaMdl, redoDeg1, redoCS21, doRot, ocnKernel, beQuiet)
    %% Loading data
    deg1Args = false;
    Ldata = prod{3};

    if redoDeg1
        deg1Args = {"Lsle", Lsle, "GIAModel", giaMdl, "OceanDomain", ocnDom, "Unit", 'SD', "IncludeC21S21", redoCS21};
    end

    [sdPlmt, sdStdPlmt, dates] = ...
        grace2plmt_new(prod{:}, "RecomputeDegree1", deg1Args, "Unit", 'SD', ...
        "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);
    sdPlmt = ensureplmdegree(sdPlmt, Ldata);
    sdStdPlmt = ensureplmdegree(sdStdPlmt, Ldata);

    % Load GIA model
    giaPlmt = gia2plmt(dates, giaMdl, "L", Ldata, ...
        "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
    % GSM degree-1 fields should be GIA-free (see Sun et al., 2016)
    sdPlmt(:, 4:end, 3:4) = sdPlmt(:, 4:end, 3:4) - giaPlmt(:, 4:end, 3:4);

    sdPlmt = permute(sdPlmt, [2, 3, 1]); % timefirst -> traditional
    sdStdPlmt = permute(sdStdPlmt, [2, 3, 1]); % timefirst -> traditional

    %% Data preprocessing
    if isempty(ocnKernel)
        ocnKernel = kernelcp_new(Lsle, ocnDom, ...
            "BeQuiet", beQuiet);
    end

    if ((isnumeric(lndDom) || ischar(lndDom) || islogical(lndDom)) && ...
            isscalar(lndDom) && isnan(lndDom))
        forcPlmt = localise(sdPlmt, "L", Lsle, "K", ocnKernel, ...
            "Inverse", true);
        forcStdPlmt = localise(sdStdPlmt, "L", Lsle, "K", ocnKernel, ...
            "Inverse", true, "IsError", true);
    else
        lndKernel = kernelcp_new(Ldata, lndDom, ...
            "rotb", true, "BeQuiet", beQuiet);
        forcPlmt = localise(sdPlmt, "L", Ldata, "K", lndKernel);
        forcStdPlmt = localise(sdStdPlmt, "L", Ldata, "K", lndKernel, "IsError", true);
    end

    forcPlmt = ensureplmdegree(forcPlmt, Lsle, 'traditional');
    forcStdPlmt = ensureplmdegree(forcStdPlmt, Lsle, 'traditional');

    [rslSdPlmt, rslSdStdPlmt] = ...
        solvesle(forcPlmt, forcStdPlmt, Lsle, ocnDom, ...
        "RotationFeedback", doRot, "OceanKernel", ocnKernel);
    rslSdPlmt = rslSdPlmt(:, 3:4, :);
    rslSdStdPlmt = rslSdStdPlmt(:, 3:4, :);
end

% Parse input arguments
function varargout = parseinputs(inputs)
    ip = inputParser;
    addOptional(ip, 'Product', {'CSR', 'RL06', 60}, ...
        @(x) iscell(x) && length(x) == 3);
    addOptional(ip, 'L', 96, @(x) isscalar(x) && isnumeric(x));
    addOptional(ip, 'OceanDomain', GeoDomain('alloceans', "Buffer", 0.5), ...
        @(x) isa(x, 'GeoDomain') || ischar(x) || (iscell(x) && length(x) == 2) ...
        || (isnumeric(x) && size(x, 2) == 2));
    addOptional(ip, 'GIA', 'ice6gd', @(x) ischar(x));
    addOptional(ip, 'RotationFeedback', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'RecomputeDegree1', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'RecomputeC21S21', false, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'Frame', 'CF', @(x) ischar(validatestring(x, {'CF', 'CM'})));
    addParameter(ip, 'LandDomain', NaN, ...
        @(x) isa(x, 'GeoDomain') || ischar(x) || (iscell(x) && length(x) == 2) ...
        || (isnumeric(x) && size(x, 2) == 2) || isnan(x));
    addParameter(ip, 'Truncation', [], ...
        @(x) (isnumeric(x) && isscalar(x)) || isempty(x));
    addOptional(ip, 'Unit', 'SD', ...
        @(x) isempty(x) || ischar(validatestring(x, {'GRAV', 'POT', 'SD'})));
    addParameter(ip, 'TimeRange', [], ...
        @(x) (isdatetime(x) || isnumeric(x)) && length(x) == 2 || isempty(x));
    addParameter(ip, 'OutputFormat', 'timefirst', ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'TimeFormat', 'datenum', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'OceanKernel', [], ...
        @(x) isnumeric(x) || isempty(x));
    addParameter(ip, 'CallChain', {}, @iscell);
    parse(ip, inputs{:});

    prod = ip.Results.Product;
    L = round(ip.Results.L);
    Lout = round(ip.Results.Truncation);
    ocnDom = ip.Results.OceanDomain;
    lndDom = ip.Results.LandDomain;
    giaMdl = lower(ip.Results.GIA);
    redoCS21 = ip.Results.RecomputeC21S21;
    redoDeg1 = logical(ip.Results.RecomputeDegree1) || redoCS21;
    frame = ip.Results.Frame;
    doRot = logical(ip.Results.RotationFeedback);
    timelim = ip.Results.TimeRange;
    unit = ip.Results.Unit;
    ocnKernel = ip.Results.OceanKernel;
    outputFmt = ip.Results.OutputFormat;
    timeFmt = ip.Results.TimeFormat;
    beQuiet = double(ip.Results.BeQuiet) * 2;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    callChain = [ip.Results.CallChain, {mfilename}];

    if isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    varargout = ...
        {prod, ocnDom, lndDom, L, Lout, ...
         giaMdl, redoDeg1, redoCS21, timelim, doRot, frame, unit, ...
         ocnKernel, outputFmt, timeFmt, beQuiet, forceNew, saveData, callChain};
end

% Get the output path
function outPath = ...
        findoutputpath(prod, L, ocnDom, lndDom, frame, giaMdl, doRot, redoDeg1, redoC21S21)
    prodFlag = sprintf("%s%s%d", prod{:});

    switch class(ocnDom)
        case 'GeoDomain'
            ocnFlag = ocnDom.Id;
        case 'char'
            ocnFlag = capitalise(ocnDom);
        case 'cell'
            ocnFlag = sprintf('%s%d', capitalise(ocnDom{1}), ocnDom{2});
        case 'numeric'
            ocnFlag = sprintf('%s', hash(ocnDom, 'sha1'));
    end

    if ~((isnumeric(lndDom) || ischar(lndDom) || islogical(lndDom)) && ...
            isscalar(lndDom) && isnan(lndDom))

        switch class(lndDom)
            case 'GeoDomain'
                lndFlag = lndDom.Id;
            case 'char'
                lndFlag = capitalise(lndDom);
            case 'cell'
                lndFlag = sprintf('%s%d', capitalise(lndDom{1}), lndDom{2});
            case 'numeric'
                lndFlag = sprintf('%s', hash(lndDom, 'sha1'));
        end

        lndFlag = sprintf('_%s', lndFlag);
    else
        lndFlag = '';

    end

    if strcmpi(frame, 'CF')
        frameFlag = '';
    else
        frameFlag = sprintf('-%s', frame);
    end

    if islogical(giaMdl) && ~giaMdl
        giaFlag = '';
    else
        giaFlag = sprintf('-%s', giaMdl);
    end

    if doRot
        rotFlag = '';
    else
        rotFlag = '-NRot';
    end

    if redoDeg1 && ~redoC21S21
        redoFlag = '-Rdeg1';
    elseif redoDeg1 && redoC21S21
        redoFlag = '-Rdeg1CS21';
    else
        redoFlag = '';
    end

    outFolder = fullfile(getenv('IFILES'), 'FINGERPRINT');

    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    outFile = sprintf('%s%s-Ls%d-%s%s%s%s%s.mat', ...
        prodFlag, giaFlag, L, ocnFlag, lndFlag, frameFlag, rotFlag, redoFlag);
    outPath = fullfile(outFolder, outFile);
end

% Format the output
function [rslPlmt, rslStdPlmt, dates] = ...
        formatoutput(rslPlmt, rslStdPlmt, dates, timelim, L, Lout, unit, ...
        outputFmt, timeFmt)

    rslPlmt = convertgravity(rslPlmt, 'SD', unit, ...
        "InputFormat", 'lmcosi', "TimeDim", 'timefirst');
    rslStdPlmt = convertgravity(rslStdPlmt, 'SD', unit, ...
        "InputFormat", 'lmcosi', "TimeDim", 'timefirst');

    if ~isempty(timelim)
        isValidTime = dates >= timelim(1) & dates <= timelim(2);
        rslPlmt = rslPlmt(:, :, isValidTime);
        rslStdPlmt = rslStdPlmt(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if ~isempty(Lout) && Lout < L
        rslPlmt = rslPlmt(1:addmup(Lout), :, :);
        rslStdPlmt = rslStdPlmt(1:addmup(Lout), :, :);
    elseif isempty(Lout)
        Lout = L;
    end

    [order, degree] = addmon(Lout);
    degree = repmat(degree, [1, 1, size(rslPlmt, 3)]);
    order = repmat(order, [1, 1, size(rslPlmt, 3)]);

    rslPlmt = cat(2, degree, order, rslPlmt);
    rslStdPlmt = cat(2, degree, order, rslStdPlmt);

    if strcmpi(outputFmt, 'timefirst')
        rslPlmt = permute(rslPlmt, [3, 1, 2]);
        rslStdPlmt = permute(rslStdPlmt, [3, 1, 2]);
    end

    if strcmpi(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

end
