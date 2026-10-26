%% SOLVEDEGREE1 - Jointly estimates low-degree GRACE coefficients.
% Estimates degree 1 and optional C20, C30, C21 and S21 using an ocean
% constraint and the specified GIA model, following Sun et al. (2016).
%
% Syntax
%   [coeffs, coeffStds, dates] = solvedegree1()
%   [coeffs, coeffStds, dates] = solvedegree1(product, "Name", value)
%   [coeffs, coeffStds, dates] = solvedegree1(Pcenter, Rlevel, Ldata, ...
%       Ltruncation, Lsle, GiaModel, OceanDomain)
%   [coeffs, coeffStds, dates, gracePlmt, coeffLocs] = solvedegree1(__)
%   [coeffs, coeffStds, dates] = solvedegree1(__, "IncludeC21S21", true)
%
% Input arguments
%   Pcenter (optional) - GRACE data centre.
%       Options: 'CSR', 'GFZ' or 'JPL'. Alternatively, pass the cell array
%       product = {Pcenter, Rlevel, Ldata} as the first argument.
%       Default centre: 'CSR'.
%   Rlevel (optional) - Release name.
%       Options: 'RL04', 'RL05' or 'RL06'. Available data and corrections
%       depend on the release.
%       Default release level: 'RL06'.
%   Ldata (optional) - Product bandwidth.
%       Default bandwidth: 60.
%   Ltruncation (optional) - Input truncation degree.
%       Default truncation: Ldata (no truncation).
%   Lsle (optional) - SLE bandwidth.
%       Default bandwidth: 96.
%   GiaModel (optional) - GIA model name.
%       Default model: 'ice6gd'.
%   OceanDomain (optional) - GeoDomain used for the ocean constraint.
%       Default domain: GeoDomain('alloceans', "Buffer", 0.5).
%   IncludeC20 (name-value) - Also estimate C20.
%       Default option: false.
%   IncludeC30 (name-value) - Also estimate C30.
%       Default option: false.
%   IncludeC21S21 (name-value) - Also estimate C21 and S21 jointly.
%       Scalar logical.
%       Default option: false.
%   Method (name-value) - Method used for the ocean constraint.
%       'fingerprint': use SOLVESLE; 'uniform': use a uniform ocean load.
%       Default method: 'fingerprint'.
%   TimeRange (name-value) - Inclusive start and end dates for the output.
%       Size: (1, 2) or (2, 1).
%       Default range: [] (all dates).
%   Unit (name-value) - Output coefficient unit.
%       'SD': surface mass density in kg/m^2, equivalent to mm of fresh
%       water. 'GRAV' and 'POT' follow the CONVERTGRAVITY conventions.
%       Default unit: 'SD'.
%   Tolerance (name-value) - Outer absolute tolerance.
%       Checks coefficients and standard deviations at every epoch.
%       Set 0 for exactly MaxIter steps.
%       Unit: kg/m^2.
%       Default tolerance: 1e-6.
%   MaxIter (name-value) - Maximum outer iterations.
%       Default iteration limit: 100.
%   SleTolerance (name-value) - Inner SLE absolute tolerance.
%       Set 0 for exactly SleMaxIter steps.
%       Unit: m.
%       Default tolerance: 1e-9.
%   SleMaxIter (name-value) - Maximum inner SLE iterations.
%       Default iteration limit: 100.
%   ForceNew (name-value) - Whether to bypass the coefficient cache.
%       Default option: false.
%   SaveData (name-value) - Whether to save computed coefficients.
%       Default option: true.
%   BeQuiet (name-value) - Whether to report progress.
%       false: messages and waitbar; 0.5: messages only; true: silent.
%       Default reporting level: 0.5.
%
% Output arguments
%   coeffs - Estimated coefficients.
%       Column order: C10, C11, S11, optional C20, optional C30, then
%       optional C21, S21. N is the number of dates and
%       K = 3 + IncludeC20 + IncludeC30 + 2 * IncludeC21S21.
%       Size: (N, K).
%       Unit: selected Unit; 'SD' uses kg/m^2.
%   coeffStds - Formal standard deviations.
%       Uses the same coefficient order as coeffs. GIA uncertainty is
%       not included.
%       Size: (N, K).
%       Unit: same as coeffs.
%   dates - Coefficient epochs as a datetime column vector.
%       Size: (N, 1).
%   gracePlmt - GRACE field with the estimated coefficients inserted.
%       The four columns are degree, order, cosine and sine coefficients.
%       Uses the GRACE reader's epoch order; coeffs, coeffStds and dates
%       retain the solved epoch order. Every GRACE epoch must have an estimate.
%       Size: (N, addmup(Ldata), 4).
%       Unit: selected Unit for coefficient columns.
%   coeffLocs - Locations corresponding to columns of coeffs.
%       Each row is [lmcosi row, cosine/sine selector (1/2)].
%       Size: (K, 2).
%
% Notes
%   GIA is removed for degrees >= 2 during estimation and restored for
%   estimated coefficients of those degrees. Degree 1 remains GIA-free.
%   Cached solutions are assumed converged at sufficiently small tolerances;
%   iteration settings are not part of cache filenames. Use ForceNew=true
%   for other settings and SaveData=false to preserve an existing cache.
%
% References
%   Sun, Y, Ditmar, P, and Riva, REM. (2015). Observed changes
%       in the Earth's dynamic oblateness from GRACE data and geophysical
%       models. Journal of Geodesy, 90:81-89. doi: 10.1007/s00190-015-0852-y.
%   Sun, Y, Riva, REM, and Ditmar, P. (2016). Optimizing estimates
%       of annual variations and trends in geocenter motion and J2 from a
%       combination of GRACE data and geophysical models. Journal of
%       Geophysical Research: Solid Earth, 121(11):8352-8370.
%       doi: 10.1002/2016JB013073.
%
% See also
%   GRACE2PLMT_NEW, GRACE2FINGERPRINT, SOLVESLE
%
% Created by
%   2025/03/18, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/26, En-Chi Lee (williameclee@arizona.edu)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [prod, Lused, Lsle, giaMdl, ocnDom, ...
         inclC20, inclC30, inclCS21, mthd, timelim, unit, iterOpts, ...
         forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin);

    dataPath = getoutputfile(prod, giaMdl, inclC20, inclC30, inclCS21, Lused, Lsle, ocnDom, mthd);

    vars = {'coeffs', 'coeffStds', 'dates'};

    if ~forceNew && exist(dataPath, 'file') && ...
            all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded %s.\n', ...
                callchaintext(callChain), filehref(dataPath, 'recomputed degree-1 data'));
        end

        varargout = formatoutput(data.coeffs, data.coeffStds, data.dates, ...
            timelim, unit, nargout, [inclC20, inclC30, inclCS21], ...
            prod, beQuiet);

        return

    end

    [coeffs, coeffStds, dates] = ...
        solvedegree1Core(prod, Lused, Lsle, ...
        giaMdl, inclC20, inclC30, inclCS21, ocnDom, mthd, beQuiet, iterOpts);

    if saveData
        save(dataPath, 'coeffs', 'coeffStds', 'dates')

        if beQuiet <= 1
            fprintf('[ULMO>%s] Saved %s.\n', ...
                callchaintext(callChain), filehref(dataPath, 'recomputed degree-1 data'));
        end

    end

    varargout = formatoutput(coeffs, coeffStds, dates, timelim, unit, nargout, ...
        [inclC20, inclC30, inclCS21], prod, beQuiet);

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
    % Prepare the ocean constraint rows once for signal and variance propagation.
    coeffOcnKernel = ocnKernelSle(coeffIds, :);
    coeffOcnKernelSq = coeffOcnKernel .^ 2;
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
    % Only the jointly estimated coefficients enter the ocean constraint.
    knownOcnCoeffs = coeffOcnKernel * packcoeffs(sdPlmt, kernelOrder);
    knownOcnCoeffVars = coeffOcnKernelSq * packcoeffs(sdStdPlmt, kernelOrder) .^ 2;

    % Only the estimated coefficients change between outer iterations.
    % Keep the fixed land signal and variance, plus their update columns.
    knownLndSdPlms = lndKernelSle * packcoeffs(sdPlmt, kernelOrder);
    knownLndVarPlms = (lndKernelSle .^ 2) * packcoeffs(sdStdPlmt, kernelOrder) .^ 2;
    coeffLndKernel = lndKernelSle(:, coeffIds);
    coeffLndKernelSq = coeffLndKernel .^ 2;
    clear lndKernelSle

    %% Solving the degree-1 coefficients
    maxIter = iterOpts.MaxIter;
    tol = iterOpts.Tolerance;
    converged = false;
    coeffs = extractcoeffs(sdPlmt, coeffLocs);
    coeffStds = extractcoeffs(sdStdPlmt, coeffLocs);
    clear sdStdPlmt

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

        lndPlmt = unpackcoeffs(knownLndSdPlms + coeffLndKernel * coeffs, ...
            kernelOrder, sdPlmt);
        lndStdPlmt = unpackcoeffs( ...
            sqrt(knownLndVarPlms + coeffLndKernelSq * coeffStds .^ 2), ...
            kernelOrder, sdPlmt);

        switch mthd
            case 'fingerprint'
                [ocnPlmt, ocnStdPlmt] = ...
                    solvesle(lndPlmt, lndStdPlmt, Lsle, ocnDom, ...
                    "OceanKernel", ocnKernelSle, "KernelOrder", kernelOrder, ...
                    "InitialRSL", sdPlmt / 1000, ... % kg/m^2 -> m
                    "RotationFeedback", true, "MaxIter", iterOpts.SleMaxIter, ...
                    "Tolerance", iterOpts.SleTolerance, "BeQuiet", true);
                ocnCoeffs = coeffOcnKernel * packcoeffs(ocnPlmt, kernelOrder);
                ocnCoeffVars = coeffOcnKernelSq * packcoeffs(ocnStdPlmt, kernelOrder) .^ 2;
            case 'uniform'
                ocnArea = ocnKernelSle(1, 1);
                eustatic = squeeze(-lndPlmt(1, 3, :) / ocnArea);
                ocnCoeffs = eustatic' .* ocnKernelSle(coeffIds, 1);
                eustaticStd = squeeze(-lndStdPlmt(1, 3, :) / ocnArea);
                ocnCoeffVars = (eustaticStd' .* ocnKernelSle(coeffIds, 1)) .^ 2;
        end

        coeffs = coeffsKernel \ (ocnCoeffs - knownOcnCoeffs);
        sdPlmt = putcoeffs(sdPlmt, coeffs, coeffLocs);
        % STD
        coeffStds = sqrt((coeffsKernel \ eye(nCoeffs)) .^ 2 * ...
            (ocnCoeffVars + knownOcnCoeffVars));

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
          forceNew, saveData, beQuiet, callChain};
end

% Get the input and output file names
function outPath = ...
        getoutputfile(prod, giaMdl, ...
        inclC20, inclC30, inclC21S21, Lused, Lsle, ocnDom, mthd)
    outFolder = fullfile(getenv('GRACEDATA'), 'Degree1', 'new');

    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    moreCoeffFlag = '';

    if inclC20
        moreCoeffFlag = [moreCoeffFlag, 'C20'];
    end

    if inclC30
        moreCoeffFlag = [moreCoeffFlag, 'C30'];
    end

    if inclC21S21
        moreCoeffFlag = [moreCoeffFlag, 'C21S21'];
    end

    if ~isempty(moreCoeffFlag)
        moreCoeffFlag = ['-W', moreCoeffFlag];
    end

    switch mthd
        case 'fingerprint'
            mthdFlag = '';
        case 'uniform'
            mthdFlag = '-uniform';
    end

    prodFlag = sprintf('%s%s%d', prod{:});
    deg1Flag = sprintf('%s%s-Ld%d_Ls%d-%s%s', ...
        giaMdl, moreCoeffFlag, ...
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

    [gracePlmt, ~, graceDates] = grace2plmt_new(product, ...
        "Unit", unit, "TimeRange", timelim, "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);

    [matched, dateIdx] = ismember(graceDates, dates);
    assert(all(matched), 'ULMO:MissingLowDegreeDates', ...
    'Joint coefficient estimates do not cover all GRACE dates.');

    for iCoeff = 1:size(coeffLocs, 1)
        gracePlmt(:, coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2)) = coeffs(dateIdx, iCoeff);
    end

    out = {coeffs, coeffStds, dates, gracePlmt, coeffLocs};

end

% Pack matching-bandwidth lmcosi fields into kernel order, with epochs in columns.
function coeffs = packcoeffs(plm, kernelOrder)
    coeffs = plm(:, 3:4, :);
    coeffs = reshape(coeffs, [], size(plm, 3));
    coeffs = coeffs(kernelOrder, :);
end

% Restore matching-bandwidth lmcosi coefficients, preserving degree/order columns.
function plm = unpackcoeffs(coeffs, kernelOrder, plm)
    nRows = size(plm, 1);
    nEpochs = size(plm, 3);
    vals = zeros(2 * nRows, nEpochs);
    vals(kernelOrder, :) = coeffs;
    plm(:, 3:4, :) = reshape(vals, nRows, 2, nEpochs);
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
