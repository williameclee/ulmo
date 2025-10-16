%% GRACE2FINGERPRINT
% Computes the relative sea level fingerprint based on forcing from GRACE data.
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
%   It doesn't make much sense to compute the GIA feedback or steric correction at this moment, so both arguments are disabled.
%
% See also
%   SOLVESLE
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/10/16, williameclee@arizona.edu (@williameclee)

function [rslLoadPlmt, rslLoadStdPlmt, dates] = grace2fingerprint(varargin)
    %% Parsing inputs
    [product, oceanDomain, landDomain, L, Loutput, ...
         giaModel, doGia, redoDeg1, rwGad, timelim, doRotation, frame, unit, ...
         outputFmt, timeFmt, beQuiet, forceNew, saveData, oceanKernel, callChain] = parseinputs(varargin);

    %% Data locating
    filepath = findoutputpath( ...
        product, L, oceanDomain, landDomain, frame, giaModel, doRotation, rwGad, redoDeg1);

    vars = {'rslLoadPlmt', 'rslLoadStdPlmt', 'dates'};

    if ~forceNew && exist(filepath, 'file')
        data = load(filepath, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">precomputed fingerprint</a>.\n', ...
                callchaintext(callChain), filepath, filepath);
        end

        [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
            formatoutput(data.rslLoadPlmt, data.rslLoadStdPlmt, data.dates, timelim, L, Loutput, unit, ...
            outputFmt, timeFmt);

        return

    end

    %% Solving the sea level equation
    [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        grace2fingerprintCore(product, oceanDomain, landDomain, L, giaModel, doGia, rwGad, redoDeg1, doRotation, oceanKernel, beQuiet);

    %% Saving and returning data
    if saveData
        save(filepath, vars{:}, '-v7.3');

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">fingerprint</a>.\n', ...
                callchaintext(callChain), filepath, filepath);
        end

    end

    [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        formatoutput(rslLoadPlmt, rslLoadStdPlmt, dates, timelim, L, Loutput, unit, ...
        outputFmt, timeFmt);

end

%% Subfunctions
% Heart of the programme
function [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        grace2fingerprintCore(product, oceanDomain, landDomain, Lsle, giaModel, doGia, rwGad, redoDeg1, doRotation, oceanKernel, beQuiet)
    %% Loading data
    deg1Args = false;
    Ldata = product{3};

    if redoDeg1
        deg1Args = {"Lsle", Lsle, "GIAModel", giaModel, "OceanDomain", oceanDomain, "Unit", 'SD'};
    end

    [gracePlmt, graceStdPlmt, dates] = ...
        grace2plmt_new(product{:}, "RecomputeDegree1", deg1Args, "Unit", 'SD', ...
        "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);
    gracePlmt = ensureplmdegree(gracePlmt, Ldata);
    graceStdPlmt = ensureplmdegree(graceStdPlmt, Ldata);

    % Add back GAC and remove GAD instead (Sun et al., 2016)
    % GAC/GAD products don't have uncertainties
    % Only replace for l >= 2 (see TN-13)
    if rwGad
        [gacPlmt, ~] = aod1b2plmt(product{1:2}, 'GAC', Ldata, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gacPlmt = ensureplmdegree(gacPlmt, Ldata);
        [gadPlmt, ~] = aod1b2plmt(product{1:2}, 'GAD', Ldata, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gadPlmt = ensureplmdegree(gadPlmt, Ldata);
        gracePlmt(:, 3:end, 3:4) = gracePlmt(:, 3:end, 3:4) ...
            + gacPlmt(:, 3:end, 3:4) - gadPlmt(:, 3:end, 3:4);
    end

    % Load GIA model
    if doGia
        giaPlmt = gia2plmt(dates, giaModel, "L", Ldata, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gracePlmt(:, 3:end, 3:4) = gracePlmt(:, 3:end, 3:4) - giaPlmt(:, 3:end, 3:4);
    end

    gracePlmt = permute(gracePlmt, [2, 3, 1]); % timefirst -> traditional
    graceStdPlmt = permute(graceStdPlmt, [2, 3, 1]); % timefirst -> traditional

    %% Data preprocessing
    if isempty(oceanKernel)
        oceanKernel = kernelcp_new(Lsle, oceanDomain, ...
            "BeQuiet", beQuiet);
    end

    if ((isnumeric(landDomain) || ischar(landDomain) || islogical(landDomain)) && ...
            isscalar(landDomain) && isnan(landDomain))
        forcingPlmt = localise(gracePlmt, "L", Lsle, "K", oceanKernel, ...
            "Inverse", true);
        forcingStdPlmt = localise(graceStdPlmt, "L", Lsle, "K", oceanKernel, ...
            "Inverse", true, "IsError", true);
    else
        landKernel = kernelcp_new(Ldata, landDomain, ...
            "rotb", true, "BeQuiet", beQuiet);
        forcingPlmt = localise(gracePlmt, "L", Ldata, "K", landKernel);
        forcingStdPlmt = localise(graceStdPlmt, "L", Ldata, "K", landKernel, "IsError", true);
    end

    forcingPlmt = ensureplmdegree(forcingPlmt, Lsle, 'traditional');
    forcingStdPlmt = ensureplmdegree(forcingStdPlmt, Lsle, 'traditional');

    [~, ~, ~, ~, ~, oceanFunPlm] = ...
        geoboxcap(Lsle, oceanDomain, "BeQuiet", beQuiet);

    [rslLoadPlmt, rslLoadStdPlmt] = ...
        solvesle(forcingPlmt, forcingStdPlmt, Lsle, oceanDomain, ...
        "RotationFeedback", doRotation, ...
        "OceanKernel", oceanKernel, "OceanFunction", oceanFunPlm);
    rslLoadPlmt = rslLoadPlmt(:, 3:4, :);
    rslLoadStdPlmt = rslLoadStdPlmt(:, 3:4, :);
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
    addOptional(ip, 'GIA', 'ice6gd', @(x) ischar(x) || islogical(x));
    addOptional(ip, 'RotationFeedback', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addOptional(ip, 'ReplaceWithGAD', true, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'RecomputeDegree1', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
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

    product = ip.Results.Product;
    L = round(ip.Results.L);
    Loutput = round(ip.Results.Truncation);
    oceanDomain = ip.Results.OceanDomain;
    landDomain = ip.Results.LandDomain;
    giaModel = ip.Results.GIA;
    replaceWGad = logical(ip.Results.ReplaceWithGAD);
    redoDeg1 = logical(ip.Results.RecomputeDegree1);
    frame = ip.Results.Frame;
    doRotation = logical(ip.Results.RotationFeedback);
    timelim = ip.Results.TimeRange;
    unit = ip.Results.Unit;
    outputFmt = ip.Results.OutputFormat;
    timeFmt = ip.Results.TimeFormat;
    beQuiet = double(ip.Results.BeQuiet) * 2;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    callChain = [ip.Results.CallChain, {mfilename}];

    oceanKernel = ip.Results.OceanKernel;

    if islogical(giaModel) && ~giaModel
        doGia = false;
    else
        doGia = true;

        if islogical(giaModel)
            giaModel = 'ice6gd';
        else
            giaModel = lower(giaModel);
        end

    end

    if isnumeric(timelim)
        timelim = datetime(timelim, "ConvertFrom", 'datenum');
    end

    varargout = ...
        {product, oceanDomain, landDomain, L, Loutput, ...
         giaModel, doGia, redoDeg1, replaceWGad, timelim, doRotation, frame, unit, ...
         outputFmt, timeFmt, beQuiet, forceNew, saveData, oceanKernel, callChain};
end

% Get the output path
function filepath = findoutputpath(product, L, oceanDomain, landDomain, frame, giaModel, doRotation, rwGad, redoDeg1)
    productName = sprintf("%s%s%d", product{:});

    switch class(oceanDomain)
        case 'GeoDomain'
            oceanName = oceanDomain.Id;
        case 'char'
            oceanName = capitalise(oceanDomain);
        case 'cell'
            oceanName = sprintf('%s%d', capitalise(oceanDomain{1}), oceanDomain{2});
        case 'numeric'
            oceanName = sprintf('%s', hash(oceanDomain, 'sha1'));
    end

    if ~((isnumeric(landDomain) || ischar(landDomain) || islogical(landDomain)) && ...
            isscalar(landDomain) && isnan(landDomain))

        switch class(landDomain)
            case 'GeoDomain'
                landName = landDomain.Id;
            case 'char'
                landName = capitalise(landDomain);
            case 'cell'
                landName = sprintf('%s%d', capitalise(landDomain{1}), landDomain{2});
            case 'numeric'
                landName = sprintf('%s', hash(landDomain, 'sha1'));
        end

        landName = sprintf('_%s', landName);
    else
        landName = '';

    end

    if rwGad
        gadName = '_GAD';
    else
        gadName = '';
    end

    if strcmpi(frame, 'CF')
        frameName = '';
    else
        frameName = sprintf('-%s', frame);
    end

    if islogical(giaModel) && ~giaModel
        giaName = '';
    else
        giaName = sprintf('-%s', giaModel);
    end

    if doRotation
        rotationName = '';
    else
        rotationName = '-NRot';
    end

    if redoDeg1
        deg1Name = '-Rdeg1';
    else
        deg1Name = '';
    end

    filefolder = fullfile(getenv('IFILES'), 'FINGERPRINT');

    if ~exist(filefolder, 'dir')
        mkdir(filefolder);
    end

    filename = sprintf('%s%s%s-Ls%d-%s%s%s%s%s.mat', ...
        productName, gadName, giaName, L, oceanName, landName, frameName, rotationName, deg1Name);
    filepath = fullfile(filefolder, filename);
end

% Format the output
function [rslLoadPlmt, rslLoadStdPlmt, dates] = ...
        formatoutput(rslLoadPlmt, rslLoadStdPlmt, dates, timelim, L, Loutput, unit, ...
        outputFmt, timeFmt)

    rslLoadPlmt = convertgravity(rslLoadPlmt, 'SD', unit, ...
        "InputFormat", 'lmcosi', "TimeDim", 'timefirst');
    rslLoadStdPlmt = convertgravity(rslLoadStdPlmt, 'SD', unit, ...
        "InputFormat", 'lmcosi', "TimeDim", 'timefirst');

    if ~isempty(timelim)
        isValidTime = dates >= timelim(1) & dates <= timelim(2);
        rslLoadPlmt = rslLoadPlmt(:, :, isValidTime);
        rslLoadStdPlmt = rslLoadStdPlmt(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if ~isempty(Loutput) && Loutput < L
        rslLoadPlmt = rslLoadPlmt(1:addmup(Loutput), :, :);
        rslLoadStdPlmt = rslLoadStdPlmt(1:addmup(Loutput), :, :);
    elseif isempty(Loutput)
        Loutput = L;
    end

    [order, degree] = addmon(Loutput);
    degree = repmat(degree, [1, 1, size(rslLoadPlmt, 3)]);
    order = repmat(order, [1, 1, size(rslLoadPlmt, 3)]);

    rslLoadPlmt = cat(2, degree, order, rslLoadPlmt);
    rslLoadStdPlmt = cat(2, degree, order, rslLoadStdPlmt);

    if strcmpi(outputFmt, 'timefirst')
        rslLoadPlmt = permute(rslLoadPlmt, [3, 1, 2]);
        rslLoadStdPlmt = permute(rslLoadStdPlmt, [3, 1, 2]);
    end

    if strcmpi(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

end
