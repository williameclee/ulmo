%% SOLVESLE - Solves the elastic sea-level equation for a fixed ocean domain.
% Computes the relative sea-level load and barystatic sea level from
% spherical harmonic forcing, with optional formal uncertainty propagation.
%
% Syntax
%   rslSd = solvesle(forcing)
%   [rslSd, rslSdStd, gmsl, gmslStd] = solvesle(forcing, forcingStd, L, ocean)
%   [rslSd, rslSdStd, gmsl, gmslStd] = solvesle(forcing, [], L, ocean)
%   rslSd = solvesle(__, "Name", value)
%
% Input arguments
%   forcing - Surface mass density spherical harmonic coefficients.
%       Forcing is not additionally localised; any ocean signal is retained.
%       M = addmup(L). A single epoch may be 2-D.
%       Size: (M, 4, N) in lmcosi order or (M, 2, N) in cosi order.
%       Unit: kg/m^2, equivalent to mm of fresh water.
%   forcingStd (optional) - Formal standard deviations of forcing.
%       An empty value disables uncertainty propagation.
%       Size: same as forcing.
%       Unit: same as forcing for coefficient columns.
%       Default field: [].
%   L (optional) - SLE bandwidth.
%       Default bandwidth: inferred from forcing.
%   ocean (optional) - Ocean domain.
%       Accepts a GeoDomain or a longitude/latitude polygon.
%       Size: (n, 2) for a polygon.
%       Unit: degrees for polygon coordinates.
%       Default domain: GeoDomain('alloceans', "Buffer", 0.5).
%   Frame (name-value) - Reference frame.
%       'CF': centre of figure; 'CM': centre of mass.
%       Default frame: 'CF'.
%   LoveNumSource (name-value) - Source of the Love numbers.
%       Options: 'ISSM' or 'ALMA3 <model>'. The 'Wahr' set lacks the
%       required Love numbers.
%       Default source: 'ISSM'.
%   RotationFeedback (name-value) - Whether to include rotation feedback.
%       Default option: true.
%   Tolerance (name-value) - Absolute tolerance.
%       Checks returned and ocean-localised coefficients, and standard
%       deviations when supplied, at every epoch. Set 0 for exactly
%       MaxIter steps.
%       Unit: m.
%       Default tolerance: 1e-9.
%   MaxIter (name-value) - Maximum iterations.
%       Default iteration limit: 100.
%   OceanKernel (name-value) - Precomputed ocean-localisation matrix.
%       The first column gives the ocean function; entry (1,1) gives
%       the positive ocean area fraction used for mass conservation.
%       An empty value computes or loads the kernel.
%       Size: ((L+1)^2, (L+1)^2).
%       Default kernel: [].
%   OceanVarianceKernel (name-value) - Precomputed mass-constrained variance kernel.
%       An empty value computes the kernel.
%       Size: ((L+1)^2, (L+1)^2).
%       Default kernel: [].
%   OceanFunction (name-value) - Legacy ocean function; accepted but unused.
%       OceanKernel supplies the geometry.
%       Default function: [].
%   KernelOrder (name-value) - Indices mapping cosi to kernel order.
%       An empty value uses KERNELORDER(L).
%       Size: vector with (L+1)^2 elements.
%       Default order: [].
%   InitialRSL (name-value) - Initial relative sea-level coefficients.
%       Uses the same layout as forcing. An empty value uses an internal
%       initial guess.
%       Size: same as forcing.
%       Unit: m for coefficient columns.
%       Default field: [].
%   BeQuiet (name-value) - Whether to suppress the progress waitbar.
%       Default option: false.
%
% Output arguments
%   rslSd - Relative sea-level load coefficients, physical over the ocean.
%       Matches the input column layout. A single epoch is a 2-D matrix.
%       Size: (addmup(L), 4, N) or (addmup(L), 2, N).
%       Unit: kg/m^2 for coefficient columns, equivalent to mm of fresh water.
%   rslSdStd - Formal standard deviations of rslSd.
%       Empty when forcingStd is not supplied.
%       Size: same as rslSd.
%       Unit: same as rslSd for coefficient columns.
%   gmsl - Barystatic sea level over the specified ocean domain.
%       Size: (N, 1).
%       Unit: mm.
%   gmslStd - Formal standard deviation of gmsl.
%       Empty when forcingStd is not supplied.
%       Size: (N, 1).
%       Unit: mm.
%
% Notes
%   The ocean function is fixed. Its domain should include the polar oceans
%   for global mass balance.
%   Uncertainties use diagonal propagation. Correlations between forcing
%   and iterated ocean load, between provisional response harmonics, and
%   between target mass and that response are neglected.
%   Positive tolerance requires convergence; reaching MaxIter without it
%   raises an error. Tolerance=0 explicitly requests a fixed iteration count.
%
% References
%   Adhikari, S et al. (2016). ISSM-SESAW v1.0: mesh-based computation
%       of gravitationally consistent sea-level and geodetic signatures caused
%       by cryosphere and climate driven mass change. Geoscientific Model
%       Development, 9(3):1087-1109. doi: 10.5194/gmd-9-1087-2016.
%   Adhikari, S et al. (2019). Sea-level fingerprints emergent from GRACE
%       mission data. Earth System Science Data, 11(2):629-646.
%       doi: 10.5194/essd-11-629-2019.
%
% See also
%   GRACE2FINGERPRINT, SOLVEDEGREE1, KERNELCP_NEW, KERNELORDER, LOVENUMBER,
%   PREPAREOCEANVARIANCEKERNEL
%
% Created by
%   2024/11/20, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function varargout = solvesle(varargin)
    %% Initialisation
    % Load all missing variables and convert to correct formats

    % Parse input arguments
    [forcSdPlm, forcSdStdPlm, inclDO, doError, ...
         L, ocean, frame, loveNumSrc, doRotation, maxIter, ...
         ocnKernel, ocnVarKer, kernelOrder, initRslPlm, beQuiet, tol] = ...
        parseinputs(varargin{:});

    % Create progress bar
    if ~beQuiet
        wbar = waitbar(0, 'Initialising', ...
            "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');
        cleanup = onCleanup(@() deleteWaitbar(wbar));
    end

    % Load the ocean kernel if not provided
    if isempty(ocnKernel)

        if ~beQuiet
            waitbar(0, wbar, 'Computing ocean kernel');
        end

        ocnKernel = kernelcp_new(L, ocean, "BeQuiet", true);
    end

    % Multiplying a constant field by the ocean mask gives this column.
    % Use the same geometry for localization, mass balance and ocean area.
    ocnFunPlms = ocnKernel(:, 1);
    ocnArea = ocnFunPlms(1);

    if ~isfinite(ocnArea) || ocnArea <= 0
        error('ULMO:solvesle:InvalidOceanArea', ...
        'OceanKernel must have a finite, positive ocean area in (1,1).');
    end

    if doError
        % The mass adjustment correlates the constant term with the other
        % harmonics. The variance kernel is not simply squaring the ocean kernel
        ocnWghts = ocnFunPlms / ocnArea;

        if isempty(ocnVarKer)
            ocnVarKer = prepareoceanvariancekernel(ocnKernel);
        end

    end

    % Load the permutation order from cosi -> flat
    if isempty(kernelOrder)
        kernelOrder = kernelorder(L);
    end

    kernelOrder = kernelOrder(:);

    % Flatten the SH arrays to 2D
    forcSdPlm = ensureplmdegree(forcSdPlm, L, 'traditional');
    forcSdStdPlm = ensureplmdegree(forcSdStdPlm, L, 'traditional');

    forcSdPlm = reshape(forcSdPlm, [], size(forcSdPlm, 3));
    forcSdStdPlm = reshape(forcSdStdPlm, [], size(forcSdStdPlm, 3));
    forcSdPlms = forcSdPlm(kernelOrder, :);
    forcSdStdPlms = forcSdStdPlm(kernelOrder, :);

    if ~isempty(initRslPlm)
        initRslPlm = ensureplmdegree(initRslPlm, L, 'traditional');
        initRslPlm = reshape(initRslPlm, [], size(initRslPlm, 3));
        initRslPlms = initRslPlm(kernelOrder, :);
    end

    %% Constants
    % Define constants and precompute SLE kernel

    if ~beQuiet
        waitbar(0, wbar, 'Constructing SLE kernel');
    end

    [~, degree] = addmon(L);
    degrees = [degree, degree];
    degrees = degrees(kernelOrder);
    llnGeoids = lovenumber(degrees, 'LLN geoid', frame, "Source", loveNumSrc);
    llnVlms = lovenumber(degrees, 'LLN VLM', frame, "Source", loveNumSrc);
    llnFactors = 1 + llnGeoids - llnVlms;
    tlnGeoids = lovenumber(degrees, 'TLN geoid', frame, "Source", loveNumSrc);
    tlnVlms = lovenumber(degrees, 'TLN VLM', frame, "Source", loveNumSrc);
    tlnFactors = 1 + tlnGeoids - tlnVlms;
    llnGeoid2 = lovenumber(2, 'LLN geoid', frame, "Source", loveNumSrc);

    % Many values from Adhikari et al. (2016)
    WATER_DENSITY = 1000; % to be consistent with other functions(mainly MASS2WEQ) and Adhikari et al. (2016)
    EARTH_DENSITY = 5517; % to be consistent with PLM2POT
    EARTH_RADIUS = 6371e3;
    EARTH_ANGULAR_VELOCITY = 7.2921e-5;
    EQUATORIAL_INERTIA = 8.0077e37;
    POLAR_INERTIA = 8.0345e37;
    SECULAR_LOVE_NUMBER = 0.942; % 'Secular (fluid) Love number '
    CHANDLER_WOBBLE_FREQUENCY = EARTH_ANGULAR_VELOCITY ...
        * (POLAR_INERTIA - EQUATORIAL_INERTIA) / EQUATORIAL_INERTIA ...
        * (1 - tlnGeoids(5) / SECULAR_LOVE_NUMBER); % Should use the elastic Earth value
    GRAVITY = 9.81;

    % Construct the SLE kernel
    % Gravitational potential
    sleKernel = spdiags(3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1), 0, ...
        numel(degrees), numel(degrees));

    if doRotation
        cmnFac = (4 * pi / 3) * (1 + llnGeoid2) * EARTH_RADIUS ^ 6 * EARTH_ANGULAR_VELOCITY ^ 2 / GRAVITY;
        c00c20Fac = 4/3 / POLAR_INERTIA;
        c21s21Fac = EARTH_ANGULAR_VELOCITY / 5 / EQUATORIAL_INERTIA / CHANDLER_WOBBLE_FREQUENCY;

        % C00
        sleKernel(1, 1) = sleKernel(1, 1) - cmnFac * c00c20Fac;
        sleKernel(1, 5) = sleKernel(1, 5) + cmnFac * c00c20Fac / sqrt(5);
        % C20
        sleKernel(5, 1) = sleKernel(5, 1) + cmnFac * c00c20Fac / sqrt(5) * tlnFactors(5);
        sleKernel(5, 5) = sleKernel(5, 5) - cmnFac * c00c20Fac / 5 * tlnFactors(5);
        % C21 and S21
        sleKernel(6, 6) = sleKernel(6, 6) + cmnFac * c21s21Fac * tlnFactors(6);
        sleKernel(7, 7) = sleKernel(7, 7) + cmnFac * c21s21Fac * tlnFactors(7);
    end

    %% Iteration
    % Solve the SLE iteratively

    % Initial guess
    gmsl =- forcSdPlms(1, :) / ocnFunPlms(1, :) / WATER_DENSITY;

    if ~isempty(initRslPlm)
        rslPlms = initRslPlms;
    else
        rslPlms = ocnFunPlms * gmsl;
    end

    rslOcnPlms = rslPlms;

    if doError
        % Carry variances through the iteration (to avoid taking square roots repeatedly); convert to STD only at output.
        forcVarPlms = forcSdStdPlms .^ 2;
        targetOcnVar = forcVarPlms(1, :) / WATER_DENSITY ^ 2;
        sleKernelSq = sleKernel .^ 2;
        ocnWghtsSq = ocnWghts .^ 2;
        ocnMeanWghtsSq = ocnKernel(1, 2:end) .^ 2;
        rslOcnVarPlms = ocnWghtsSq * targetOcnVar;
    end

    converged = false;
    maxChange = Inf;
    prevStd = [];

    for iIter = 1:maxIter
        prevRslPlms = rslPlms;
        prevRslOcnPlms = rslOcnPlms;

        if ~beQuiet
            waitbar(iIter / maxIter, wbar, ...
                sprintf('Solving SLE iteratively (%d/%d)', iIter, maxIter));

            if getappdata(wbar, 'canceling')
                error('ULMO:ProcessCancelledByUser', 'SLE computation cancelled');
            end

        end

        sdPlms = forcSdPlms + rslOcnPlms * WATER_DENSITY;

        rslPlms = sleKernel * sdPlms;
        rslOcnPlms = ocnKernel * rslPlms;

        % Add a constant sea-level offset and its complete ocean footprint.
        % Keep the iterated field equal to oceanKernel * the returned field.
        massOffset = (-forcSdPlms(1, :) / WATER_DENSITY ...
            - rslOcnPlms(1, :)) / ocnArea;
        rslPlms(1, :) = rslPlms(1, :) + massOffset;
        rslOcnPlms = rslOcnPlms + ocnFunPlms * massOffset;

        % Require convergence of both returned and localized coefficients.
        maxChange = max(absoluteiterationchange(rslPlms, prevRslPlms), ...
            absoluteiterationchange(rslOcnPlms, prevRslOcnPlms));

        if doError
            sdVarPlms = forcVarPlms + WATER_DENSITY ^ 2 * rslOcnVarPlms;
            rslVarPlms = sleKernelSq * sdVarPlms;
            rslOcnVarPlms = ocnVarKer * rslVarPlms + ocnWghtsSq * targetOcnVar;
            rslOcnVarPlms(1, :) = targetOcnVar;
            rslVarPlms(1, :) = (targetOcnVar + ocnMeanWghtsSq * rslVarPlms(2:end, :)) / ocnArea ^ 2;
            crntStd = sqrt([rslVarPlms; rslOcnVarPlms]);

            if isempty(prevStd)
                maxChange = Inf;
            else
                maxChange = max(maxChange, absoluteiterationchange(crntStd, prevStd));
            end

            prevStd = crntStd;

            if any(~isfinite(crntStd), 'all')
                error('ULMO:solvesle:NonFiniteIterate', 'SLE uncertainties became nonfinite.');
            end

        end

        if any(~isfinite(rslPlms), 'all') || any(~isfinite(rslOcnPlms), 'all')
            error('ULMO:solvesle:NonFiniteIterate', 'SLE coefficients became nonfinite.');
        elseif tol > 0 && maxChange <= tol
            converged = true;
            break
        end

    end

    if tol > 0 && ~converged
        error('ULMO:solvesle:NonConvergence', ...
            'SLE did not converge in %d iterations (maximum absolute change %.6g; tolerance %.6g).', ...
            maxIter, maxChange, tol);
    end

    %% Post-processing
    if ~beQuiet
        waitbar(1, wbar, 'Post-processing results');
    end

    rslSdPlms = rslPlms * WATER_DENSITY;
    gmsl = rslOcnPlms(1, :) ./ ocnFunPlms(1, :) * 1e3; % m -> mm (= kg/m^2)

    if isvector(rslSdPlms)
        rslSdPlm = zeros([addmup(L), 2]);
        rslSdPlm(kernelOrder) = rslSdPlms;
    else
        rslSdPlm = zeros([addmup(L) * 2, size(rslSdPlms, 2)]);
        rslSdPlm(kernelOrder, :) = rslSdPlms;
        rslSdPlm = ...
            reshape(rslSdPlm, [addmup(L), 2, size(rslSdPlms, 2)]);
    end

    if inclDO
        [order, degree] = addmon(L);

        if isvector(rslSdPlms)
            rslSdPlm = [degree, order, rslSdPlm];
        else
            degree = repmat(degree, [1, 1, size(rslSdPlm, 3)]);
            order = repmat(order, [1, 1, size(rslSdPlm, 3)]);
            rslSdPlm = cat(2, degree, order, rslSdPlm);
        end

    end

    if ~doError
        varargout = {rslSdPlm, [], gmsl(:), []};
        return
    end

    % Post-process the standard deviation
    rslSdStdPlms = sqrt(rslVarPlms) * WATER_DENSITY;
    gmslStd = sqrt(targetOcnVar) / ocnArea * 1e3; % m -> mm

    if isvector(rslSdStdPlms)
        rslSdStdPlm = zeros([addmup(L), 2]);
        rslSdStdPlm(kernelOrder) = rslSdStdPlms;
    else
        rslSdStdPlm = zeros([addmup(L) * 2, size(rslSdStdPlms, 2)]);
        rslSdStdPlm(kernelOrder, :) = rslSdStdPlms;
        rslSdStdPlm = ...
            reshape(rslSdStdPlm, [addmup(L), 2, size(rslSdStdPlms, 2)]);
    end

    if inclDO

        if isvector(rslSdStdPlms)
            rslSdStdPlm = [degree, order, rslSdStdPlm];
        else
            rslSdStdPlm = cat(2, degree, order, rslSdStdPlm);
        end

    end

    varargout = {rslSdPlm, rslSdStdPlm, gmsl(:), gmslStd(:)};
end

%% Subfunctions
% Parse input arguments
function varargout = parseinputs(varargin)
    ip = inputParser;
    addRequired(ip, 'ForcingPlm', ...
        @(x) isnumeric(x));
    addOptional(ip, 'ForcingStdPlm', [], ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'L', [], ...
        @(x) (isnumeric(x) && isscalar(x) && x > 0) || isempty(x));
    addOptional(ip, 'Ocean', GeoDomain('alloceans', "Buffer", 0.5), ...
        @(x) isa(x, 'GeoDomain') || (isnumeric(x) && ismatrix(x) && size(x, 2) == 2));
    addOptional(ip, 'frame', 'CF', ...
        @(x) ischar(validatestring(upper(x), {'CM', 'CF'})));
    addOptional(ip, 'RotationFeedback', true, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'maxIter', 100, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == fix(x));
    addParameter(ip, 'Tolerance', 1e-9, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(ip, 'LoveNumSource', 'ISSM', ...
        @(x) ischar(x) || isstring(x));
    addParameter(ip, 'OceanKernel', [], @(x) isnumeric(x));
    addParameter(ip, 'OceanVarianceKernel', [], ...
        @(x) isa(x, 'double') && isreal(x) && ismatrix(x));
    addParameter(ip, 'OceanFunction', [], @(x) isnumeric(x));
    addParameter(ip, 'KernelOrder', [], @(x) isnumeric(x));
    addParameter(ip, 'InitialRSL', [], @(x) isnumeric(x));
    addParameter(ip, 'BeQuiet', false, @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    parse(ip, varargin{:});
    forcPlm = squeeze(ip.Results.ForcingPlm);
    forcStdPlm = squeeze(ip.Results.ForcingStdPlm);
    L = ip.Results.L;
    ocnDom = ip.Results.Ocean;
    frame = upper(ip.Results.frame);
    loveNumSrc = char(ip.Results.LoveNumSource);
    doRotation = ip.Results.RotationFeedback;
    maxIter = ip.Results.maxIter;
    ocnKernel = ip.Results.OceanKernel;
    ocnVarKer = ip.Results.OceanVarianceKernel;
    kernelOrder = ip.Results.KernelOrder;
    initRslPlm = ip.Results.InitialRSL;
    beQuiet = logical(ip.Results.BeQuiet);

    if strcmpi(loveNumSrc, "Wahr")
        error(['Wahr''s (and slepian''s) love numbers are only provided for loading gravitational potential, which would not be sufficient for solving the SLE. ', ...
               'Consider using ISSM (preferred) or ALMA3 love numbers instead.']);
    end

    % Check input sizes
    inclDO = size(forcPlm, 2) == 4;

    if isempty(forcStdPlm)
        forcStdPlm = zeros(size(forcPlm));
        doError = false;
    else

        if ~isequal(size(forcStdPlm), size(forcPlm))
            error('ULMO:solvesle:InvalidInput:InputSizeNotMatch', ...
                'Input size of FORCINGSTDPLM ([%s]) does not match FORCINGPLM ([%s])', ...
                strjoin(string(size(forcStdPlm)), ', '), ...
                strjoin(string(size(forcPlm)), ', '));
        end

        doError = true;

    end

    if ~isempty(initRslPlm) && ~isequal(size(initRslPlm), size(forcPlm))
        error( ...
            sprintf('%s:InvalidInput:InputSizeNotMatch', upper(mfilename)), ...
            'Input size of INITIAL RSL ([%s]) does not match FORCINGPLM ([%s])', ...
            strjoin(string(size(initRslPlm)), ', '), ...
            strjoin(string(size(forcPlm)), ', '));
    end

    % Find the degree of the spherical harmonic coefficients
    if inclDO

        if isempty(L)
            L = max(forcPlm(:, 1, 1));
        end

        forcPlm = forcPlm(:, 3:4, :);
        forcStdPlm = forcStdPlm(:, 3:4, :);

        if ~isempty(initRslPlm)
            initRslPlm = initRslPlm(:, 3:4, :);
        end

    elseif isempty(L)
        L = finddegree(forcPlm(:, :, 1));
    end

    if ~isempty(ocnVarKer)

        if isempty(ocnKernel)
            error('ULMO:solvesle:MissingOceanKernel', ...
            'OceanVarianceKernel requires an explicit matching OceanKernel.');
        end

        expdSize = [(L + 1) ^ 2, (L + 1) ^ 2];

        if ~isequal(size(ocnVarKer), size(ocnKernel), expdSize)
            error('ULMO:solvesle:InvalidOceanVarianceKernelSize', ...
            'OceanVarianceKernel and OceanKernel must both have size (L+1)^2 by (L+1)^2.');
        end

    end

    varargout = ...
        {forcPlm, forcStdPlm, inclDO, doError, ...
         L, ocnDom, frame, loveNumSrc, doRotation, maxIter, ...
         ocnKernel, ocnVarKer, kernelOrder, initRslPlm, beQuiet, ip.Results.Tolerance};
end

% Find the degree of the spherical harmonic coefficients based on the
% number of terms
function L = finddegree(Plm)
    terms = size(Plm, 1);
    syms x y
    eq = x ^ 2 + 3 * x + (2 - 2 * y) == 0;
    eq = subs(eq, y, terms);
    Ls = solve(eq, x);
    L = Ls(Ls > 0);
end
