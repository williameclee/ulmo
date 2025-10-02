%% SOLVEDEGREE1
% Recomputes the GRACE degree-1 surface mass density coefficients based on
% Sun et al. (2016) with the specified GIA model.
%
% Syntax
%   [coeffs, coeffStd, dates] = SOLVEDEGREE1(...
%       pcenter, rlevel, Ldata, Lsle, GIAModel, oceanDomain)
%   [__] = SOLVEDEGREE1(__, ...
%       "ReplaceWithGAD", rwGad, "IncludeC20", includeC20, "TimeRange", timelim)
%   [__] = SOLVEDEGREE1(__, ...
%       "ForceNew", forceNew, "SaveData", saveData, "BeQuiet", beQuiet)
%   [__, gracePlmt] = SOLVEDEGREE1(__)
%
% Input arguments
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data centre is 'CSR'.
%       When the first argument is a cell array, it is interpreted as
%       {Pcenter, Rlevel, Ldata}.
%       Data type: CHAR
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06' (or numbers).
%       The default release level is 'RL06'.
%       Currently, only RL06 is guaranteed to work.
%       Data type: CHAR | [NUMERIC]
%   Ldata - Bandwidth of the date product
%       In the case where there are more than one product from a data
%       centre (such as BA 60 or BB 96 standard L2 products) this allows
%       you to choose between them.
%       The default L is 60.
%       Data type: [NUMERIC]
%   Ltruncation - Bandwidth to truncate the date product
%       The default option is no truncation.
%       Data type: [NUMERIC]
%   Lsle - Bandwidth used to solve the sea level equation
%       The default Lsle is 96.
%       Data type: [NUMERIC]
%   GiaModel - Name of GIA model
%       The default GIA model is 'ice6gd'.
%       Data type: CHAR
%   OceanDomain - Ocean domain
%       The input format should be whatever KERNELCP_NEW accepts, e.g. a
%       GeoDomain object.
%       The default ocean domain is the global ocean (ALLOCEANS) with a
%       0.5° buffer.
%   IncludeC20 - Whether to also recompute C20
%       The default option is false, as C20 is usually replaced with SLR
%       values.
%       Data type: LOGICAL
%   ReplaceWithGAD - Whether to remove GAD instead of GAC
%       The default option is true (see Sun et al., 2016).
%       Data type: LOGICAL
%   Method - Method to solve the sea level equation
%       - 'fingerprint': Solve the sea level equation using the fingerprint
%           method (see SOLVESLE).
%       - 'uniform': Assume a uniform barystatic load over the ocean
%           domain.
%       The default method is 'fingerprint'.
%       Data type: CHAR
%   Unit - Unit of the output
%       - 'GRAV': Surface gravity (this is POT in SLEPIAN_ALPHA).
%       - 'POT': Geopotential field (surface gravity * equatorial radius).
%       - 'SD': Surface mass density.
%       The default field is 'SD'.
%       Data type: CHAR
%   ForceNew - Logical flag to force reprocess of the data
%       The default option is false.
%       Data type: LOGICAL
%	SaveData - Logical flag to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%		The default option is true.
%		Data types: LOGICAL
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is false.
%		Data types: LOGICAL
%
% Output arguments
%   coeffs - C10, C11, S11 (and maybe C20) surface mass density
%       coefficients
%   dates - dates of the coefficients
%       In DATETIME format.
%   gracePlmt - GRACE SH coefficients with the recomputed coefficients
%
% Example
%   coeffs = SOLVEDEGREE1('CSR', 'RL06', 60, 96, ...
%       'ice6gd', GeoDomain('alloceans', "Buffer", 2))
%   coeffs = SOLVEDEGREE1('CSR', 'RL06', 45, 180, ...
%       'ice6gd', GeoDomain('alloceans', "Buffer", 2), ...
%       "IncludeC20", true, "ReplaceWithGAD", true)
%
% See also
%   SOLVESLE, GRACE2PLMT
%
% Authored by
%   2025/03/18, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/10/01, williameclee@arizona.edu (@williameclee)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [pcenter, rlevel, Ldata, Ltruncation, Lsle, giaModel, oceanDomain, ...
         includeC20, includeC30, replaceWGad, method, timelim, unit, ...
         forceNew, saveData, beQuiet, onlyId, callChain] = ...
        parseinputs(varargin);

    [dataPath, ~, outputId] = ...
        getoutputfile(pcenter, rlevel, Ldata, replaceWGad, giaModel, ...
        includeC20, includeC30, Ltruncation, Lsle, oceanDomain, method);

    if onlyId
        varargout = {outputId};
        return
    end

    vars = {'coeffs', 'coeffStds', 'dates'};

    if ~forceNew && exist(dataPath, 'file') && ...
            all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">recomputed degree-1 data</a>.\n', ...
                callchaintext(callChain), dataPath, dataPath);
        end

        varargout = formatoutput(data.coeffs, data.coeffStds, data.dates, timelim, unit, nargout, [includeC20, includeC30], ...
            {pcenter, rlevel, Ldata}, beQuiet);

        return

    end

    [coeffs, coeffStds, dates] = ...
        solvedegree1Core(pcenter, rlevel, Ldata, Ltruncation, Lsle, replaceWGad, ...
        giaModel, includeC20, includeC30, oceanDomain, method, beQuiet);

    if saveData
        save(dataPath, 'coeffs', 'coeffStds', 'dates')

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">recomputed degree-1 data</a>.\n', ...
                callchaintext(callChain), dataPath, dataPath);
        end

    end

    varargout = formatoutput(coeffs, coeffStds, dates, timelim, unit, nargout, [includeC20, includeC30], ...
        {pcenter, rlevel, Ldata}, beQuiet);

end

%% Subfunctions
% Heart of the programme
function [coeffs, coeffStds, dates] = ...
        solvedegree1Core(pcenter, rlevel, Ldata, Ltruncation, Lsle, rwGad, giaModel, includeC20, includeC30, oceanDomain, method, beQuiet)
    %% Loading data
    wbar = waitbar(0, 'Loading GRACE data', ...
        "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

    [gracePlmt, graceStdPlmt, dates] = grace2plmt_new(pcenter, rlevel, Ldata, ...
        "Unit", 'SD', "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);
    gracePlmt = ensureplmdegree(gracePlmt, Ltruncation);
    gracePlmt = ensureplmdegree(gracePlmt, Lsle);
    graceStdPlmt = ensureplmdegree(graceStdPlmt, Ltruncation);
    graceStdPlmt = ensureplmdegree(graceStdPlmt, Lsle);

    % Add back GAC and remove GAD instead (Sun et al., 2016)
    % GAC/GAD products don't have uncertainties
    % Only replace for l >= 2 (see TN-13)
    if rwGad
        [gacPlmt, ~] = aod1b2plmt(pcenter, rlevel, 'GAC', Lsle, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gacPlmt = ensureplmdegree(gacPlmt, Lsle);
        [gadPlmt, ~] = aod1b2plmt(pcenter, rlevel, 'GAD', Lsle, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gadPlmt = ensureplmdegree(gadPlmt, Lsle);
        gracePlmt(:, 3:end, 3:4) = gracePlmt(:, 3:end, 3:4) ...
            + gacPlmt(:, 3:end, 3:4) - gadPlmt(:, 3:end, 3:4);
    end

    % Remove GIA signal for l >= 2 (Sun et al., 2016)
    giaPlmt = gia2plmt(dates, giaModel, "L", Lsle, ...
        "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
    gracePlmt(:, 3:end, 3:4) = gracePlmt(:, 3:end, 3:4) - giaPlmt(:, 3:end, 3:4);
    % Ignore STD of GIA for now

    %% Preparing/preallocating variables
    waitbar(0, wbar, 'Preparing kernels and variables');
    gracePlmt = permute(gracePlmt, [2, 3, 1]); % timefirst -> traditional
    graceStdPlmt = permute(graceStdPlmt, [2, 3, 1]); % timefirst -> traditional

    if getappdata(wbar, 'canceling')
        delete(wbar);
        error(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
        'Processing cancelled');
    end

    % Whether to also reestimate C20 and C30
    nCoeffs = 3 + includeC20 + includeC30; % C10, C11, S11, (C20), (C30)
    coeffs = nan([length(dates), nCoeffs]);

    % Precompute kernels
    coeffLocs = [2, 1; 3, 1; 3, 2];
    coeffIds = [2, 3, 4];

    if includeC20
        coeffLocs = [coeffLocs; 4, 1];
        coeffIds = [coeffIds, 5];
    end

    if includeC30
        coeffLocs = [coeffLocs; 7, 1];
        coeffIds = [coeffIds, 10];
    end

    % Note: Converting to SINGLE does not work because it is not compatible with the sparse SLE kernel
    oceanKernelSle = kernelcp_new(Lsle, oceanDomain, ...
        "BeQuiet", beQuiet);
    landKernelSle = eye(size(oceanKernelSle)) - oceanKernelSle;
    coeffsKernel = nan([nCoeffs, nCoeffs]);

    for iCoeff = 1:nCoeffs

        for jCoeff = 1:nCoeffs
            coeffsKernel(iCoeff, jCoeff) = oceanKernelSle(coeffIds(iCoeff), coeffIds(jCoeff));
        end

    end

    [~, ~, ~, ~, ~, oceanFunPlm] = ...
        geoboxcap(Lsle, oceanDomain, "BeQuiet", beQuiet);
    kernelOrder = kernelorder(Lsle);

    graceNoUnknownPlmt = putcoeffs(gracePlmt, 0, coeffLocs);
    oceanCoeffsNoUnkown = extractcoeffs( ...
        localise(graceNoUnknownPlmt, "L", Lsle, "K", oceanKernelSle), ...
        coeffLocs);

    % STD
    graceStdNoUnknownPlmt = putcoeffs(graceStdPlmt, 0, coeffLocs);
    oceanCoeffStdsNoUnkown = extractcoeffs( ...
        localise(graceStdNoUnknownPlmt, "L", Lsle, "K", oceanKernelSle), ...
        coeffLocs);

    %% Solving the degree-1 coefficients
    maxIter = 5;

    for iIter = 1:maxIter
        waitbar((iIter - 1) / maxIter, wbar, ...
            sprintf('Solving degree-1 coefficients iteratively (%d/%d)', iIter, maxIter));

        if getappdata(wbar, 'canceling')
            delete(wbar);
            error(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
            'Computation cancelled');
        end

        landPlmt = localise(gracePlmt, "L", Lsle, "K", landKernelSle);
        landStdPlmt = localise(graceStdPlmt, "L", Lsle, "K", landKernelSle, "IsError", true);

        switch method
            case 'fingerprint'
                [oceanPlmt, oceanStdPlmt] = ...
                    solvesle(landPlmt, landStdPlmt, Lsle, oceanDomain, ...
                    "OceanKernel", oceanKernelSle, "OceanFunction", oceanFunPlm, ...
                    "KernelOrder", kernelOrder, "InitialCondition", gracePlmt, ...
                    "RotationFeedback", true, "MaxIter", 5 + (5 * (iIter == maxIter)), "BeQuiet", true);
                oceanPlmt = localise(oceanPlmt, "L", Lsle, "K", oceanKernelSle);
                oceanStdPlmt = localise(oceanStdPlmt, "L", Lsle, "K", oceanKernelSle, "IsError", true);
                oceanCoeffs = extractcoeffs(oceanPlmt, coeffLocs);
                oceanCoeffStds = extractcoeffs(oceanStdPlmt, coeffLocs);
            case 'uniform'
                barystaticLoad = squeeze(-landPlmt(1, 3, :) / oceanFunPlm(1, 3));
                oceanCoeffs = barystaticLoad' .* oceanKernelSle((1:nCoeffs) + 1, 1);
                barystaticLoadStd = squeeze(-landStdPlmt(1, 3, :) / oceanFunPlm(1, 3));
                oceanCoeffStds = barystaticLoadStd' .* oceanKernelSle((1:nCoeffs) + 1, 1);
        end

        coeffs = coeffsKernel \ (oceanCoeffs - oceanCoeffsNoUnkown);
        gracePlmt = putcoeffs(gracePlmt, coeffs, coeffLocs);
        % STD
        coeffStds = sqrt((coeffsKernel \ eye(nCoeffs)) .^ 2 * ...
            (oceanCoeffStds .^ 2 + oceanCoeffStdsNoUnkown .^ 2));
        graceStdPlmt = putcoeffs(graceStdPlmt, coeffStds, coeffLocs);
    end

    coeffs = coeffs'; % traditional -> timefirst
    coeffStds = coeffStds'; % traditional -> timefirst

    %% Postprocessing and output
    waitbar(1, wbar, 'Postprocessing results');

    % Restore GAC/GAD
    if rwGad

        for iCoeff = 1:nCoeffs
            coeffs(:, iCoeff) = ...
                coeffs(:, iCoeff) ...
                - squeeze(gacPlmt(:, coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2)) ...
                - gadPlmt(:, coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2)));
        end

    end

    % Restore GIA signal for C20
    if includeC20 && includeC30
        coeffs(:, 4) = coeffs(:, 4) + ...
            squeeze(giaPlmt(:, 4, 3));
        coeffs(:, 5) = coeffs(:, 5) + ...
            squeeze(giaPlmt(:, 7, 3));
    elseif includeC20
        coeffs(:, 4) = coeffs(:, 4) + ...
            squeeze(giaPlmt(:, 4, 3));
    elseif includeC30
        coeffs(:, 5) = coeffs(:, 5) + ...
            squeeze(giaPlmt(:, 7, 3));
    end

    delete(wbar);

end

% Parse input arguments
function varargout = parseinputs(inputs)
    % Set default parameters
    dfOpt.Pcenter = 'CSR';
    dfOpt.Rlevel = 'RL06';
    dfOpt.Ldata = 60;
    dfOpt.Ltruncation = [];
    dfOpt.Lsle = 96;
    dfOpt.GiaModel = 'ice6gd';
    dfOpt.OceanDomain = GeoDomain('alloceans', "Buffer", 0.5); % Sun et al. (2016)
    dfOpt.IncludeC20 = false;
    dfOpt.IncludeC30 = false;
    dfOpt.RwGAD = true;
    dfOpt.Method = 'fingerprint';
    % Construct input parser
    ip = inputParser;
    addOptional(ip, 'Pcenter', dfOpt.Pcenter, ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', dfOpt.Rlevel, ...
        @(x) ischar(validatestring(x, {'RL04', 'RL05', 'RL06'})));
    addOptional(ip, 'Ldata', dfOpt.Ldata, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'Ltruncation', dfOpt.Ltruncation, ...
        @(x) (isscalar(x) && isnumeric(x) && x > 0) || isempty(x));
    addOptional(ip, 'Lsle', dfOpt.Lsle, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'GiaModel', dfOpt.GiaModel, ...
        @(x) ischar(x));
    addOptional(ip, 'OceanDomain', dfOpt.OceanDomain, ...
        @(x) isa(x, 'GeoDomain'));
    addOptional(ip, 'IncludeC20', dfOpt.IncludeC20, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'IncludeC30', dfOpt.IncludeC30, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'ReplaceWithGAD', dfOpt.RwGAD, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'Method', dfOpt.Method, ...
        @(x) ischar(validatestring(x, {'fingerprint', 'uniform'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || isdatetime(x) || isnumeric(x));
    addOptional(ip, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'GRAV', 'POT', 'SD'})));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'OnlyId', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'CallChain', {}, @iscell);

    if iscell(inputs{1})
        inputs = [inputs{1}{:}, inputs(2:end)];
    end

    parse(ip, inputs{:});

    pcenter = ip.Results.Pcenter;
    rlevel = ip.Results.Rlevel;
    Ldata = ip.Results.Ldata;
    Ltruncation = conddefval(ip.Results.Ltruncation, Ldata);
    Lsle = ip.Results.Lsle;
    giaModel = ip.Results.GiaModel;
    oceanDomain = ip.Results.OceanDomain;
    includeC20 = logical(ip.Results.IncludeC20);
    includeC30 = logical(ip.Results.IncludeC30);
    replaceWGad = logical(ip.Results.ReplaceWithGAD);
    method = ip.Results.Method;
    timelim = ip.Results.TimeRange;
    unit = ip.Results.Unit;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = double(ip.Results.BeQuiet) * 2;
    onlyId = logical(ip.Results.OnlyId);
    callChain = [ip.Results.CallChain, {mfilename}];

    if ~isempty(regexp(lower(giaModel), '^ice[-]?6g[_-]?d$', 'once'))
        giaModel = 'ICE6GD';
    elseif contains(lower(giaModel), 'caron')
        giaModel = replace(giaModel, 'caron', 'Caron');
    elseif contains(lower(giaModel), 'paulson')
        giaModel = replace(giaModel, 'paulson', 'Paulson');
    end

    varargout = ...
        {pcenter, rlevel, Ldata, Ltruncation, Lsle, giaModel, oceanDomain, ...
         includeC20, includeC30, replaceWGad, method, timelim, unit, ...
         forceNew, saveData, beQuiet, onlyId, callChain};
end

% Get the input and output file names
function [outputPath, outputFile, deg1Id] = ...
        getoutputfile(pcenter, rlevel, Ldata, replaceWGad, giaModel, ...
        includeC20, includeC30, Ltruncation, Lsle, oceanDomain, method)
    outputFolder = fullfile(getenv('GRACEDATA'), 'Degree1', 'new');

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
        fprintf('%s created folder %s\n', upper(mfilename), outputFolder);
    end

    gadFlag = '';

    if replaceWGad
        gadFlag = 'GAD-';
    end

    includeC20Flag = '';

    if includeC20 && includeC30
        includeC20Flag = '-WC2030';
    elseif includeC20
        includeC20Flag = '-WC20';
    elseif includeC30
        includeC20Flag = '-WC30';
    end

    switch method
        case 'fingerprint'
            methodFlag = '';
        case 'uniform'
            methodFlag = '-uniform';
    end

    productId = sprintf('%s%s%d', pcenter, rlevel, Ldata);
    deg1Id = sprintf('%s%s%s-Ld%d_Ls%d-%s%s', ...
        gadFlag, giaModel, includeC20Flag, ...
        Ltruncation, Lsle, oceanDomain.Id, methodFlag);

    outputFile = sprintf('%s_%s-SD.mat', productId, deg1Id);

    outputPath = fullfile(outputFolder, outputFile);
end

% Format the output
function output = formatoutput(coeffs, coeffStds, dates, timelim, unit, nOut, includeC2030, product, beQuiet)
    coeffs(:, 1:3) = convertgravity(coeffs(:, 1:3), 'SD', unit, ...
        "InputFormat", 'L', "L", 1);
    coeffStds(:, 1:3) = convertgravity(coeffStds(:, 1:3), 'SD', unit, ...
        "InputFormat", 'L', "L", 1);

    if size(coeffs, 2) == 4

        if includeC2030(1)
            L = 2;
        else
            L = 3;
        end

        coeffs(:, 4) = convertgravity(coeffs(:, 4), 'SD', unit, ...
            "InputFormat", 'L', "L", L);
        coeffStds(:, 4) = convertgravity(coeffStds(:, 4), 'SD', unit, ...
            "InputFormat", 'L', "L", L);
    elseif size(coeffs, 2) == 5
        coeffs(:, 4) = convertgravity(coeffs(:, 4), 'SD', unit, ...
            "InputFormat", 'L', "L", 2);
        coeffStds(:, 4) = convertgravity(coeffStds(:, 4), 'SD', unit, ...
            "InputFormat", 'L', "L", 2);
        coeffs(:, 5) = convertgravity(coeffs(:, 5), 'SD', unit, ...
            "InputFormat", 'L', "L", 3);
        coeffStds(:, 5) = convertgravity(coeffStds(:, 5), 'SD', unit, ...
            "InputFormat", 'L', "L", 3);
    end

    if ~isempty(timelim)
        isTimeRange = dates >= timelim(1) & dates <= timelim(2);
        dates = dates(isTimeRange);
        coeffs = coeffs(isTimeRange, :);
        coeffStds = coeffStds(isTimeRange, :);
    end

    if nOut <= 3
        output = {coeffs, coeffStds, dates};
        return
    end

    [gracePlmt, ~, dates] = grace2plmt_new(product, ...
        "Unit", unit, "TimeRange", timelim, "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);

    gracePlmt(:, 2, 3) = coeffs(:, 1);
    gracePlmt(:, 3, 3) = coeffs(:, 2);
    gracePlmt(:, 3, 4) = coeffs(:, 3);

    if size(coeffs, 2) == 4
        gracePlmt(:, 4, 3) = coeffs(:, 4);
    end

    output = {coeffs, coeffStds, dates, gracePlmt};

end

% Just get the coeffs
function coeffs = extractcoeffs(plm, coeffLocs)
    nCoeffs = size(coeffLocs, 1);
    coeffLen = size(plm, 3);
    coeffs = zeros(nCoeffs, coeffLen);

    for iCoeff = 1:nCoeffs
        coeffs(iCoeff, :) = plm(coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2), :);
    end

end

% Just put the coeffs back
function plm = putcoeffs(plm, coeffs, coeffLocs)
    nCoeffs = size(coeffLocs, 1);

    if isscalar(coeffs)
        coeffs = repmat(coeffs, [nCoeffs, 1]);
    end

    for iCoeff = 1:nCoeffs
        plm(coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2), :) = coeffs(iCoeff, :);
    end

end
