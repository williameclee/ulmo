%% SOLVEDEGREE1 - Solves for the geocentre motion from GRACE observations.
% Recomputes the GRACE degree-1 surface mass density coefficients based on
% Sun et al. (2016) with the specified GIA model.
% Estimates remain GSM-like: the solver does not add GAC or subtract GAD.
% Restore GAD separately downstream when constructing ocean-mass signals.
%
% Syntax
%   [coeffs, coeffStd, dates] = SOLVEDEGREE1(...
%       pcenter, rlevel, Ldata, Lsle, GIAModel, oceanDomain)
%   [__] = SOLVEDEGREE1(__, ...
%       "IncludeC20", includeC20, "TimeRange", timelim)
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
%   IncludeC20, IncludeC30 (name-value) - Whether to also recompute C20 or C30
%       Default option is false.
%   IncludeC21S21 - Jointly estimate C21 and S21 with degree 1.
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
%   Tolerance (name-value) - Absolute tolerance.
%       Estimated coefficient and standard-deviation vectors are checked
%       separately; every coefficient in every epoch must pass.
%       Set 0 to use exactly MaxIter steps.
%       Default absolute tolerance is 1e-6 kg/m^2 = 1e-6 mm.
%   MaxIter (name-value) - Maximum number of the outer iteration.
%       Default number of iterations is 100.
%   SleTolerance (name-value) - Inner SLE absolute tolerance.
%       Set 0 to use exactly SleMaxIter steps.
%       Default absolute tolerance is 1e-6 mm.
%   SleMaxIter (name-value) - Maximum number of inner SLE iterations.
%       Default number of iterations is 100.
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
%   coeffs - Columns C10, C11, S11, followed by optional C20, C30,
%       then C21 and S21. Values use the requested Unit (default SD).
%   coeffStds - Formal uncertainties in the same order and units.
%   coeffLocs - Fifth output: [lmcosi row, cosine/sine column (1/2)].
%   dates - dates of the coefficients
%       In DATETIME format.
%   gracePlmt - GRACE SH coefficients with the recomputed coefficients
%
% Example
%   coeffs = SOLVEDEGREE1('CSR', 'RL06', 60, 96, ...
%       'ice6gd', GeoDomain('alloceans', "Buffer", 2))
%   coeffs = SOLVEDEGREE1('CSR', 'RL06', 45, 180, ...
%       'ice6gd', GeoDomain('alloceans', "Buffer", 2), ...
%       "IncludeC20", true)
%
% See also
%   SOLVESLE, GRACE2PLMT
%
% References
%   Sun, Y., Ditmar, P., and Riva, R. E. M. (2015). Observed changes
%       in the Earth’s dynamic oblateness from GRACE data and
%       geophysical models. Journal of Geodesy, 90(1):81-89. doi:
%       10.1007/s00190-015-0852-y.
%   Sun, Y., Riva, R. E. M., and Ditmar, P. (2016). Optimizing
%       estimates of annual variations and trends in geocenter
%       motion and J2 from a combination of GRACE data and
%       geophysical models. Journal of Geophysical Research: Solid
%       Earth, 121(11):8352-8370. doi: 10.1002/2016JB013073.
%
% Author
%   2025/03/18, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified by
%   2026/10/26, En-Chi Lee (williameclee@arizona.edu)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [prod, Lused, Lsle, giaMdl, ocnDom, ...
         inclC20, inclC30, inclC21S21, mthd, timelim, unit, iterOpts, ...
         onlyId, forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin);

    [dataPath, ~, outputId] = ...
        getoutputfile(prod, giaMdl, inclC20, inclC30, inclC21S21, Lused, Lsle, ocnDom, mthd);

    if onlyId
        varargout = {outputId};
        return
    end

    vars = {'coeffs', 'coeffStds', 'dates'};

    if ~forceNew && exist(dataPath, 'file') && ...
            all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded %s.\n', ...
                callchaintext(callChain), filehref(dataPath, 'recomputed degree-1 data'));
        end

        varargout = formatoutput(data.coeffs, data.coeffStds, data.dates, timelim, unit, nargout, [inclC20, inclC30, inclC21S21], ...
            prod, beQuiet);

        return

    end

    [coeffs, coeffStds, dates] = ...
        solvedegree1Core(prod, Lused, Lsle, ...
        giaMdl, inclC20, inclC30, inclC21S21, ocnDom, mthd, beQuiet, iterOpts);

    if saveData
        save(dataPath, 'coeffs', 'coeffStds', 'dates')

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved %s.\n', ...
                callchaintext(callChain), filehref(dataPath, 'recomputed degree-1 data'));
        end

    end

    varargout = formatoutput(coeffs, coeffStds, dates, timelim, unit, nargout, [inclC20, inclC30, inclC21S21], ...
        prod, beQuiet);

end

%% Subfunctions
% Heart of the programme
function [coeffs, coeffStds, dates] = ...
        solvedegree1Core(prod, Lused, Lsle, giaMdl, inclC20, inclC30, inclC21S21, ocnDom, mthd, beQuiet, iterOpts)
    %% Loading data
    if ~beQuiet
        wbar = waitbar(0, 'Loading GRACE data', ...
            "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');
        cleanup = onCleanup(@() deleteWaitbar(wbar));
    end

    [sdPlmt, sdStdPlmt, dates] = grace2plmt_new(prod{:}, ...
        "Unit", 'SD', "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);
    sdPlmt = ensureplmdegree(sdPlmt, Lused);
    sdPlmt = ensureplmdegree(sdPlmt, Lsle);
    sdStdPlmt = ensureplmdegree(sdStdPlmt, Lused);
    sdStdPlmt = ensureplmdegree(sdStdPlmt, Lsle);

    % Remove GIA signal for l >= 2 (Sun et al., 2016)
    giaPlmt = gia2plmt(dates, giaMdl, "L", Lsle, ...
        "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
    sdPlmt(:, 4:end, 3:4) = sdPlmt(:, 4:end, 3:4) - giaPlmt(:, 4:end, 3:4);
    % Ignore STD of GIA for now

    %% Preparing/preallocating variables
    if ~beQuiet
        waitbar(0, wbar, 'Preparing kernels and variables');

        if getappdata(wbar, 'canceling')
            error('ULMO:ProcessCancelledByUser', 'Geocentre computation cancelled');
        end

    end

    sdPlmt = permute(sdPlmt, [2, 3, 1]); % timefirst -> traditional
    sdStdPlmt = permute(sdStdPlmt, [2, 3, 1]); % timefirst -> traditional

    % Number of simultaneously estimated low-degree coefficients
    nCoeffs = 3 + inclC20 + inclC30 + 2 * inclC21S21; % Degree 1 and optional degree-2/3 terms

    % Precompute kernels
    [coeffLocs, coeffIds] = lowdegreecoefficients(inclC20, inclC30, inclC21S21);

    % Note: Converting to SINGLE does not work because it is not compatible with the sparse SLE kernel
    ocnKernelSle = kernelcp_new(Lsle, ocnDom, ...
        "BeQuiet", beQuiet);
    lndKernelSle = eye(size(ocnKernelSle)) - ocnKernelSle;
    coeffsKernel = nan([nCoeffs, nCoeffs]);

    for iCoeff = 1:nCoeffs

        for jCoeff = 1:nCoeffs
            coeffsKernel(iCoeff, jCoeff) = ocnKernelSle(coeffIds(iCoeff), coeffIds(jCoeff));
        end

    end

    kernelOrder = kernelorder(Lsle);

    % Set unknown coefficients to 0
    sdPlmt = putcoeffs(sdPlmt, 0, coeffLocs);
    sdStdPlmt = putcoeffs(sdStdPlmt, 0, coeffLocs);
    knownOcnSdCoeffs = extractcoeffs( ...
        localise(sdPlmt, "L", Lsle, "K", ocnKernelSle), ...
        coeffLocs);
    knownOcnSdCoeffStds = extractcoeffs( ...
        localise(sdStdPlmt, "L", Lsle, "K", ocnKernelSle, "IsError", true), ...
        coeffLocs);

    %% Solving the degree-1 coefficients
    maxIter = iterOpts.MaxIter;
    tol = iterOpts.Tolerance;
    converged = false;
    coeffs = extractcoeffs(sdPlmt, coeffLocs);
    coeffStds = extractcoeffs(sdStdPlmt, coeffLocs);

    for iIter = 1:maxIter
        prevCoeffs = coeffs;
        prevStds = coeffStds;

        if ~beQuiet
            waitbar((iIter - 1) / maxIter, wbar, ...
                sprintf('Solving degree-1 coefficients iteratively (%d/%d)', iIter, maxIter));

            if getappdata(wbar, 'canceling')
                error('ULMO:ProcessCancelledByUser', 'Geocentre computation cancelled');
            end

        end

        lndPlmt = localise(sdPlmt, "L", Lsle, "K", lndKernelSle);
        lndStdPlmt = localise(sdStdPlmt, "L", Lsle, "K", lndKernelSle, "IsError", true);

        switch mthd
            case 'fingerprint'
                [ocnPlmt, ocnStdPlmt] = ...
                    solvesle(lndPlmt, lndStdPlmt, Lsle, ocnDom, ...
                    "OceanKernel", ocnKernelSle, "KernelOrder", kernelOrder, ...
                    "InitialRSL", sdPlmt / 1000, ... % kg/m^2 -> m
                    "RotationFeedback", true, "MaxIter", iterOpts.SleMaxIter, ...
                    "Tolerance", iterOpts.SleTolerance, "BeQuiet", true);
                ocnPlmt = localise(ocnPlmt, "L", Lsle, "K", ocnKernelSle);
                ocnStdPlmt = localise(ocnStdPlmt, "L", Lsle, "K", ocnKernelSle, "IsError", true);
                ocnCoeffs = extractcoeffs(ocnPlmt, coeffLocs);
                ocnCoeffStds = extractcoeffs(ocnStdPlmt, coeffLocs);
            case 'uniform'
                ocnArea = ocnKernelSle(1, 1);
                eustatic = squeeze(-lndPlmt(1, 3, :) / ocnArea);
                ocnCoeffs = eustatic' .* ocnKernelSle(coeffIds, 1);
                eustaticStd = squeeze(-lndStdPlmt(1, 3, :) / ocnArea);
                ocnCoeffStds = eustaticStd' .* ocnKernelSle(coeffIds, 1);
        end

        coeffs = coeffsKernel \ (ocnCoeffs - knownOcnSdCoeffs);
        sdPlmt = putcoeffs(sdPlmt, coeffs, coeffLocs);
        % STD
        coeffStds = sqrt((coeffsKernel \ eye(nCoeffs)) .^ 2 * ...
            (ocnCoeffStds .^ 2 + knownOcnSdCoeffStds .^ 2));
        sdStdPlmt = putcoeffs(sdStdPlmt, coeffStds, coeffLocs);

        if any(~isfinite(coeffs), 'all') || any(~isfinite(coeffStds), 'all')
            error('ULMO:solvedegree1:NonFiniteIterate', ...
            'Estimated coefficients or uncertainties became nonfinite.');
        end

        maxChange = max(absoluteiterationchange(coeffs, prevCoeffs), ...
            absoluteiterationchange(coeffStds, prevStds));

        if tol > 0 && maxChange <= tol
            converged = true;
            break
        end

    end

    if tol > 0 && ~converged
        error('ULMO:solvedegree1:NonConvergence', ...
            'Solution did not converge in %d iterations (maximum absolute change %.6g; tolerance %.6g).', ...
            maxIter, maxChange, tol);
    end

    coeffs = coeffs'; % traditional -> timefirst
    coeffStds = coeffStds'; % traditional -> timefirst

    %% Postprocessing and output
    if ~beQuiet
        waitbar(1, wbar, 'Postprocessing results');
    end

    % Restore the GIA signal for every estimated coefficient of degree >= 2.
    for iCoeff = 4:nCoeffs
        coeffs(:, iCoeff) = coeffs(:, iCoeff) + ...
            squeeze(giaPlmt(:, coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2)));
    end

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
    addParameter(ip, 'IncludeC21S21', false, @(x) islogical(x) && isscalar(x));
    addOptional(ip, 'Method', dfOpt.Method, ...
        @(x) ischar(validatestring(x, {'fingerprint', 'uniform'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || isdatetime(x) || isnumeric(x));
    addOptional(ip, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'GRAV', 'POT', 'SD'})));
    valdLimit = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == fix(x);
    valdTol = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0;
    addParameter(ip, 'MaxIter', 100, valdLimit);
    addParameter(ip, 'Tolerance', 1e-6, valdTol);
    addParameter(ip, 'SleMaxIter', 100, valdLimit);
    addParameter(ip, 'SleTolerance', 1e-9, valdTol);
    addParameter(ip, 'ForceNew', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'OnlyId', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'CallChain', {}, @iscell);

    if ~isempty(inputs) && iscell(inputs{1})
        inputs = [inputs{1}, inputs(2:end)];
    end

    parse(ip, inputs{:});

    pcenter = ip.Results.Pcenter;
    rlevel = ip.Results.Rlevel;
    Ldata = ip.Results.Ldata;
    Lused = conddefval(ip.Results.Ltruncation, Ldata);
    Lsle = ip.Results.Lsle;
    giaMdl = ip.Results.GiaModel;
    ocnDom = ip.Results.OceanDomain;
    inclC20 = logical(ip.Results.IncludeC20);
    inclC30 = logical(ip.Results.IncludeC30);
    inclC21S21 = ip.Results.IncludeC21S21;
    mthd = ip.Results.Method;
    timelim = ip.Results.TimeRange;
    unit = ip.Results.Unit;
    onlyId = logical(ip.Results.OnlyId);
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = double(ip.Results.BeQuiet) * 2;
    callChain = [ip.Results.CallChain, {mfilename}];

    iterOpts = struct('MaxIter', ip.Results.MaxIter, ...
        'Tolerance', ip.Results.Tolerance, 'SleMaxIter', ip.Results.SleMaxIter, ...
        'SleTolerance', ip.Results.SleTolerance);

    if ~isempty(regexp(lower(giaMdl), '^ice[-]?6g[_-]?d$', 'once'))
        giaMdl = 'ICE6GD';
    elseif contains(lower(giaMdl), 'caron')
        giaMdl = replace(giaMdl, 'caron', 'Caron');
    elseif contains(lower(giaMdl), 'paulson')
        giaMdl = replace(giaMdl, 'paulson', 'Paulson');
    end

    varargout = ...
        {{pcenter, rlevel, Ldata}, Lused, Lsle, giaMdl, ocnDom, ...
          inclC20, inclC30, inclC21S21, mthd, timelim, unit, iterOpts, ...
          onlyId, forceNew, saveData, beQuiet, callChain};
end

% Get the input and output file names
function [outPath, outFile, deg1Flag] = ...
        getoutputfile(prod, giaMdl, ...
        inclC20, inclC30, inclC21S21, Lused, Lsle, ocnDom, mthd)
    outFolder = fullfile(getenv('GRACEDATA'), 'Degree1', 'new');

    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    inclC20Flag = '';

    if inclC20 && inclC30
        inclC20Flag = '-WC2030';
    elseif inclC20
        inclC20Flag = '-WC20';
    elseif inclC30
        inclC20Flag = '-WC30';
    end

    if inclC21S21
        inclC20Flag = [inclC20Flag, '-WC21S21'];
    end

    switch mthd
        case 'fingerprint'
            mthdFlag = '';
        case 'uniform'
            mthdFlag = '-uniform';
    end

    prodFlag = sprintf('%s%s%d', prod{:});
    deg1Flag = sprintf('%s%s-Ld%d_Ls%d-%s%s', ...
        giaMdl, inclC20Flag, ...
        Lused, Lsle, ocnDom.Id, mthdFlag);
    outFile = sprintf('%s_%s-SD.mat', prodFlag, deg1Flag);

    outPath = fullfile(outFolder, outFile);
end

% Format the output
function out = formatoutput(coeffs, coeffStds, dates, timelim, unit, nOut, includeC2030, product, beQuiet)
    [coeffLocs, ~, degrees] = lowdegreecoefficients(includeC2030(1), includeC2030(2), includeC2030(3));

    for iCoeff = 1:numel(degrees)
        coeffs(:, iCoeff) = convertgravity(coeffs(:, iCoeff), 'SD', unit, ...
            "InputFormat", 'L', "L", degrees(iCoeff));
        coeffStds(:, iCoeff) = convertgravity(coeffStds(:, iCoeff), 'SD', unit, ...
            "InputFormat", 'L', "L", degrees(iCoeff));
    end

    if ~isempty(timelim)
        isTimeRange = dates >= timelim(1) & dates <= timelim(2);
        dates = dates(isTimeRange);
        coeffs = coeffs(isTimeRange, :);
        coeffStds = coeffStds(isTimeRange, :);
    end

    if nOut <= 3
        out = {coeffs, coeffStds, dates};
        return
    end

    [gracePlmt, ~, dates] = grace2plmt_new(product, ...
        "Unit", unit, "TimeRange", timelim, "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);

    for iCoeff = 1:size(coeffLocs, 1)
        gracePlmt(:, coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2)) = coeffs(:, iCoeff);
    end

    out = {coeffs, coeffStds, dates, gracePlmt, coeffLocs};

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
