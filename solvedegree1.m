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
%   [coeffs, coeffStds, dates] = solvedegree1(__, RotationFeedback=false)
%   [coeffs, coeffStds, dates] = solvedegree1(__, OceanKernel=O, ...
%       EarthResponseKernel=G, ForcingKernel=F, RotationFeedback=false)
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
%       'fingerprint': use the elastic SLE; 'uniform': use a uniform ocean load.
%       Default method: 'fingerprint'.
%   RotationFeedback (name-value) - Whether to include rotational feedback.
%       Applies to Method='fingerprint'; ignored for Method='uniform'.
%       Results without rotation use the filename suffix -NRot.
%       Default option: true.
%   OceanKernel (name-value) - Precomputed ocean-localisation matrix.
%       Must match OceanDomain and Lsle, in kernelorder(Lsle) order.
%       Size: ((Lsle+1)^2, (Lsle+1)^2). Unit: dimensionless.
%       Default kernel: [] (compute or load with KERNELCP_NEW).
%   EarthResponseKernel (name-value) - Precomputed elastic Earth response.
%       Supply together with ForcingKernel, as returned by PREPARESLEKERNELS
%       for the same ocean kernel, Lsle and RotationFeedback, using its
%       default CF frame and ISSM Love numbers. Ignored for Method='uniform'.
%       Size and ordering: same as OceanKernel. Unit: dimensionless.
%       Default kernel: [] (construct both SLE kernels).
%   ForcingKernel (name-value) - Precomputed mass-constrained forcing matrix.
%       Supply together with EarthResponseKernel; see its matching requirements.
%       Ignored for Method='uniform'.
%       Size and ordering: same as OceanKernel. Unit: dimensionless.
%       Default kernel: [] (construct both SLE kernels).
%   TimeRange (name-value) - Inclusive start and end dates for the output.
%       Size: (1, 2) or (2, 1).
%       Default range: [] (all dates).
%   Unit (name-value) - Output coefficient unit.
%       'SD': surface mass density in kg/m^2, equivalent to mm of fresh
%       water. 'GRAV' and 'POT' follow the CONVERTGRAVITY conventions.
%       Default unit: 'SD'.
%   RelativeTolerance (name-value) - Relative residual for response-map solves.
%       Default tolerance: 1e-12.
%   Tolerance (name-value) - Positive absolute low-degree residual tolerance.
%       Unit: kg/m^2.
%       Default tolerance: 1e-6.
%   SleMaxIter (name-value) - Maximum iterations per response-map solve.
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
%   Explicitly supplied kernels bypass coefficient-cache reads. SaveData still
%   controls writes to the filename for the selected settings. Matrix contents
%   cannot establish their ocean or rotation provenance: the caller must supply
%   matching kernels, or use SaveData=false for custom physical settings.
%
%   Match the selected ocean coefficients of the observed load and SLE:
%       D = EU^T * O * A^(-1) * KF * (I - O),
%       A = I - KS,
%       DU = D * EU - EU^T * O * EU,
%       DnotU = EU^T * O * EnotU - D * EnotU,
%       DU * LU = DnotU * LnotU,
%       T = DU \ DnotU,
%       LU = T * LnotU.
%   where:
%     - O: Ocean localisation kernel (ocnKer)
%     - KS: RSL load's response kernel to itself, applied by SOLVESLESYSTEM
%     - KF: RSL load's response kernel to forcing (forcKer)
%     - EU: Unit vectors selecting the unknown coefficients (coeffIds)
%     - EnotU: Unit vectors selecting the known coefficients (knownIds)
%     - D: Selected ocean response, constructed by PREPARELOWDEGREEMAP
%     - DU: Unknown-coefficient constraint matrix (unknownConstraint)
%     - DnotU: Known-coefficient constraint matrix (knownConstraint)
%     - T: Map from known to unknown loads (known2unknownMap)
%     - LnotU: Known load coefficients (known(knownIds,:))
%     - LU: Estimated load coefficients (coeffs before transposing)
%   For Method='uniform', it instead uses
%       D = -(EU^T * O * e00 / O00) * e00^T * (I - O),
%   where O00 = ocnKer(1,1) and e00 selects the first (C00) coefficient.
%
%   Independent known-coefficient variances propagate as
%       V[LU] = (T.^2) * V[LnotU].
%   where:
%     - V[LnotU]: Input variances (knownStd(knownIds,:).^2)
%   The variance formula assumes C[LnotU] is diagonal.
%
%   Check the ocean constraint before restoring GIA:
%       residual = DU * LU - DnotU * LnotU.
%   Each entry must be within Tolerance, in kg/m^2.
%
%   GIA is removed for degrees >= 2 during estimation and restored for
%   estimated coefficients of those degrees. Degree 1 remains GIA-free.
%   GIA uncertainty is excluded.
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
%   GRACE2PLMT_NEW, GRACE2FINGERPRINT, SOLVESLE, PREPARESLEKERNELS, PREPARELOWDEGREEMAP
%
% Created by
%   2025/03/18, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/29, En-Chi Lee (williameclee@arizona.edu)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [prod, Lused, Lsle, giaMdl, ocnDom, ...
         inclC20, inclC30, inclCS21, mthd, timelim, unit, solveOpts, ...
         forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin);

    dataPath = getoutputfile(prod, giaMdl, inclC20, inclC30, inclCS21, ...
        Lused, Lsle, ocnDom, mthd, solveOpts.RotationFeedback);

    vars = {'coeffs', 'coeffStds', 'dates'};

    if ~forceNew && ~solveOpts.HasKernels && exist(dataPath, 'file') && ...
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

    [coeffs, coeffStds, dates] = solvedegree1Core(prod, Lused, Lsle, ...
        giaMdl, inclC20, inclC30, inclCS21, ocnDom, mthd, beQuiet, solveOpts);

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

% Heart of the programme
function [coeffs, coeffStds, dates] = solvedegree1Core( ...
        prod, Lused, Lsle, giaMdl, inclC20, inclC30, inclCS21, ocnDom, mthd, beQuiet, solveOpts)
    %% Initialisation
    [coeffLocs, coeffIds] = lowdegreecoefficients(inclC20, inclC30, inclCS21);
    nCoeffs = numel(coeffIds);
    kerOrder = kernelorder(Lsle);

    %% Data and kernel preparation
    [sdPlmt, sdStdPlmt, giaPlmt, dates] = ...
        preparedata(prod, Lused, Lsle, giaMdl, coeffLocs, beQuiet);

    ocnKer = solveOpts.OceanKernel;

    if isempty(ocnKer)
        ocnKer = kernelcp_new(Lsle, ocnDom, "BeQuiet", beQuiet);
    end

    if strcmp(mthd, 'fingerprint')
        earthRespKer = solveOpts.EarthResponseKernel;
        forcKer = solveOpts.ForcingKernel;

        if isempty(earthRespKer)
            [earthRespKer, forcKer] = prepareslekernels(Lsle, ocnKer, ...
                RotationFeedback = solveOpts.RotationFeedback);
        end
    else
        earthRespKer = [];
        forcKer = [];
    end

    %% Ocean constraint
    % T = DU \ DnotU; reuse the map for every epoch
    [kn2unknMap, knIds, ~, unknCstr, knCstr] = preparelowdegreemap( ...
        ocnKer, earthRespKer, forcKer, coeffIds, Method = mthd, ...
        RelativeTolerance = solveOpts.RelativeTolerance, MaxIter = solveOpts.SleMaxIter);
    % Pack loads and standard deviations with epochs in columns
    known = packcoeffs(sdPlmt, kerOrder);
    knownStd = packcoeffs(sdStdPlmt, kerOrder);
    % LU = T * LnotU
    coeffs = kn2unknMap * known(knIds, :);
    % sqrt(V[LU]) = sqrt((T.^2) * V[LnotU])
    coeffStds = sqrt((kn2unknMap .^ 2) * knownStd(knIds, :) .^ 2);

    if any(~isfinite(coeffs), 'all') || any(~isfinite(coeffStds), 'all')
        error('ULMO:solvedegree1:NonFiniteIterate', ...
        'Estimated coefficients or uncertainties became nonfinite.');
    end

    % Check DU * LU - DnotU * LnotU in surface mass density units.
    resid = unknCstr * coeffs - knCstr * known(knIds, :);

    if max(abs(resid), [], 'all') > solveOpts.Tolerance
        error('ULMO:solvedegree1:NonConvergence', ...
            ['The low-degree constraint residual (max: %.5g) exceeds tolerance (%.5g); ', ...
         'the relative tolerance (%.5g) in the response map solver is likely too high.'], ...
            max(abs(resid), [], 'all'), solveOpts.Tolerance, solveOpts.RelativeTolerance);
    end

    coeffs = coeffs'; % traditional -> timefirst
    coeffStds = coeffStds'; % traditional -> timefirst

    %% Postprocessing and output
    % Restore the GIA signal for every estimated coefficient of degree >= 2.
    for iCoeff = 4:nCoeffs
        coeffs(:, iCoeff) = coeffs(:, iCoeff) + ...
            squeeze(giaPlmt(:, coeffLocs(iCoeff, 1), 2 + coeffLocs(iCoeff, 2)));
    end

end

%% PREPAREDATA - Loads known GRACE fields and subtract GIA
% Signal and STD use traditional lmcosi layout at Lsle, with unknown
% coefficients zeroed. GIA retains timefirst layout for restoration.
function [sdPlmt, sdStdPlmt, giaPlmt, dates] = preparedata( ...
        prod, Lused, Lsle, giaMdl, coeffLocs, beQuiet)
    [sdPlmt, sdStdPlmt, dates] = grace2plmt_new(prod{:}, ...
        "Unit", 'SD', "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);
    % First truncate, then pad
    sdPlmt = ensureplmdegree(sdPlmt, Lused);
    sdPlmt = ensureplmdegree(sdPlmt, Lsle);
    sdStdPlmt = ensureplmdegree(sdStdPlmt, Lused);
    sdStdPlmt = ensureplmdegree(sdStdPlmt, Lsle);

    % Remove GIA signal for l >= 2, l <= 1 is not used anyway
    giaPlmt = gia2plmt(dates, giaMdl, "L", Lsle, ...
        "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
    sdPlmt(:, 4:end, 3:4) = sdPlmt(:, 4:end, 3:4) - giaPlmt(:, 4:end, 3:4);
    % Ignore STD of GIA for now

    sdPlmt = permute(sdPlmt, [2, 3, 1]); % timefirst -> traditional
    sdStdPlmt = permute(sdStdPlmt, [2, 3, 1]); % timefirst -> traditional
    % The estimated coefficients are unknown in both known input fields.
    sdPlmt = putcoeffs(sdPlmt, 0, coeffLocs);
    sdStdPlmt = putcoeffs(sdStdPlmt, 0, coeffLocs);
end

%% Subfunctions
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
    addParameter(ip, 'RotationFeedback', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x) && isreal(x) && ismember(x, [0, 1]));
    addParameter(ip, 'OceanKernel', [], @isnumeric);
    addParameter(ip, 'EarthResponseKernel', [], @isnumeric);
    addParameter(ip, 'ForcingKernel', [], @isnumeric);
    addParameter(ip, 'RelativeTolerance', 1e-12, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
    addParameter(ip, 'Tolerance', 1e-6, valdTol);
    addParameter(ip, 'SleMaxIter', 100, valdLimit);
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

    if ip.Results.Tolerance == 0
        error('ULMO:solvedegree1:InvalidTolerance', 'Tolerance must be positive.');
    end

    solveOpts = struct('RelativeTolerance', ip.Results.RelativeTolerance, ...
        'Tolerance', ip.Results.Tolerance, 'SleMaxIter', ip.Results.SleMaxIter, ...
        'RotationFeedback', logical(ip.Results.RotationFeedback), ...
        'OceanKernel', ip.Results.OceanKernel, ...
        'EarthResponseKernel', ip.Results.EarthResponseKernel, ...
        'ForcingKernel', ip.Results.ForcingKernel);
    solveOpts.HasKernels = ~isempty(solveOpts.OceanKernel) || ...
        ~isempty(solveOpts.EarthResponseKernel) || ~isempty(solveOpts.ForcingKernel);
    validatekernels(solveOpts, Lsle);

    if ~isempty(regexp(lower(giaMdl), '^ice[-]?6g[_-]?d$', 'once'))
        giaMdl = 'ICE6GD';
    elseif contains(lower(giaMdl), 'caron')
        giaMdl = replace(giaMdl, 'caron', 'Caron');
    elseif contains(lower(giaMdl), 'paulson')
        giaMdl = replace(giaMdl, 'paulson', 'Paulson');
    end

    varargout = ...
        {{pcenter, rlevel, Ldata}, Lused, Lsle, giaMdl, ocnDom, ...
          inclC20, inclC30, inclC21S21, mthd, timelim, unit, solveOpts, ...
          forceNew, saveData, beQuiet, callChain};
end

% Validate caller-supplied matrices before reading any coefficient cache.
function validatekernels(options, L)
    names = {'OceanKernel', 'EarthResponseKernel', 'ForcingKernel'};
    nCoeffs = (L + 1) ^ 2;

    for k = 1:numel(names)
        matrix = options.(names{k});

        if isempty(matrix)
            continue
        end

        if ~isreal(matrix) || ~isequal(size(matrix), [nCoeffs, nCoeffs])
            error('ULMO:solvedegree1:InvalidKernel', ...
                '%s must be a real (Lsle+1)^2 by (Lsle+1)^2 matrix.', names{k});
        end

        % Bound validation memory for dense high-degree ocean kernels.
        for first = 1:256:nCoeffs

            if any(~isfinite(matrix(:, first:min(first + 255, nCoeffs))), 'all')
                error('ULMO:solvedegree1:InvalidKernel', ...
                    '%s must contain only finite values.', names{k});
            end

        end

    end

    if ~isempty(options.OceanKernel) && options.OceanKernel(1, 1) <= 0
        error('ULMO:solvedegree1:InvalidKernel', ...
        'OceanKernel must have positive ocean area in (1,1).');
    end

    if xor(isempty(options.EarthResponseKernel), isempty(options.ForcingKernel))
        error('ULMO:solvedegree1:IncompleteSleKernels', ...
        'Supply EarthResponseKernel and ForcingKernel together.');
    end

end

% Get the input and output file names
function outPath = getoutputfile( ...
        prod, giaMdl, inclC20, inclC30, inclC21S21, Lused, Lsle, ocnDom, mthd, doRotation)
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

            if ~doRotation
                mthdFlag = '-NRot';
            end

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
function out = formatoutput( ...
        coeffs, coeffStds, dates, timelim, unit, nOut, includeC2030, product, beQuiet)
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
