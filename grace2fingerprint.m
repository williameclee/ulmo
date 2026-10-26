%% GRACE2FINGERPRINT - Computes the relative sea level fingerprint based on forcing from GRACE data.
%
% Syntax
%   [rslPlmt, time] = grace2fingerprint(L)
%   [rslPlmt, time] = grace2fingerprint(L, "GIA", model, "GIAFeedback", feedback)
%   [rslPlmt, time] = grace2fingerprint(__, "RotationFeedback", feedback)
%   [rslPlmt, time] = grace2fingerprint(__, "Name", Value)
%
% Input arguments
%   product - GRACE product
%       A cell array of the form {centre, release, degree}.
%       The default product is {'CSR', 'RL06', 60}
%   L - Bandwidth to compute spherical harmonic coefficients
%       The default value is 96
%   oceanDomain - Ocean domain
%       A geographic domain (GeoDomain object).
%       The default domain is all oceans with a buffer of 1 degree
%   "GIA" - Glacial isostatic adjustment model
%       The default model is 'ICE-6G_D'
%   "GIAFeedback" - Whether to include GIA feedback
%       The default value is true
%   "RotationFeedback" - Whether to include rotation feedback
%       The default value is true
%   "RecomputeC21S21" - Jointly reestimate degree 1, C21 and S21 (default false).
%       Enables RecomputeDegree1 and uses a separate fingerprint cache.
%   "Frame" - Reference frame
%       Centre of mass (CM) frame or centre of figure (CF) frame.
%       The default frame is the CF frame, which is the frame used in slepian_delta
%   "Truncation" - Truncate the output to a lower degree
%       The default value is the same as L
%
% Output arguments
%   rslPlmt - Relative sea level load spherical harmonic coefficients
%   time - Time stamps of the data, in datetime format
%
% Notes
%   Forcing remains GSM-like; GAC/GAD are not added or subtracted here.
%   Restore GAD separately downstream when constructing ocean-mass signals.
%   Regenerate existing fingerprints with ForceNew=true after changing
%   the degree-1 forcing convention; cache names do not track that change.
%   It doesn't make much sense to compute the GIA feedback or steric correction at this moment, so both arguments are disabled.
%
% See also
%   SOLVESLE
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2026/10/26, williameclee@arizona.edu (@williameclee)

function [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint(varargin)
    %% Parsing inputs
    [prod, ocnDom, lndDom, L, Lout, ...
         giaMdl, redoDeg1, redoC21S21, timelim, doRot, frame, unit, ...
         ocnKernel, outputFmt, timeFmt, beQuiet, forceNew, saveData, callChain] = parseinputs(varargin);

    %% Data locating
    filepath = findoutputpath( ...
        prod, L, ocnDom, lndDom, frame, giaMdl, doRot, redoDeg1, redoC21S21);

    vars = {'rslLoadPlmt', 'rslLoadStdPlmt', 'dates'};

    if ~forceNew && exist(filepath, 'file')
        data = load(filepath, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded %s.\n', ...
                callchaintext(callChain), filehref(filepath, 'fingerprint'));
        end

        [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
            formatoutput(data.rslLoadPlmt, data.rslLoadStdPlmt, data.dates, timelim, L, Lout, unit, ...
            outputFmt, timeFmt);

        return

    end

    %% Solving the sea level equation
    [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        grace2fingerprintCore(prod, ocnDom, lndDom, L, giaMdl, redoDeg1, redoC21S21, doRot, ocnKernel, beQuiet);

    %% Saving and returning data
    if saveData
        save(filepath, vars{:}, '-v7.3');

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved %s.\n', ...
                callchaintext(callChain), filehref(filepath, 'fingerprint'));
        end

    end

    [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        formatoutput(rslLoadPlmt, rslLoadStdPlmt, dates, timelim, L, Lout, unit, ...
        outputFmt, timeFmt);

end

%% Subfunctions
% Heart of the programme
function [rslSdPlmt, rslSdStdPlmt, dates] = ...
        grace2fingerprintCore(prod, ocnDom, lndDom, Lsle, giaMdl, redoDeg1, redoC21S21, doRotation, ocnKernel, beQuiet)
    %% Loading data
    deg1Args = false;
    Ldata = prod{3};

    if redoDeg1
        deg1Args = {"Lsle", Lsle, "GIAModel", giaMdl, "OceanDomain", ocnDom, "Unit", 'SD', "IncludeC21S21", redoC21S21};
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
        "RotationFeedback", doRotation, "OceanKernel", ocnKernel);
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
    redoC21S21 = ip.Results.RecomputeC21S21;
    redoDeg1 = logical(ip.Results.RecomputeDegree1) || redoC21S21;
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
         giaMdl, redoDeg1, redoC21S21, timelim, doRot, frame, unit, ...
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
