%% SOLVESLE
% Solves the (elastic) sea level equation (SLE) for a given forcing and
% ocean domain.
% The ocean function is assumed constant.
%
% Syntax
%   [rslSd, [], gmsl, []] = SOLVESLE(forcing, [], L, ocean)
%   [rslSd, rslSdStd, gmsl, gmslStd] = SOLVESLE(forcing, forcingStd, L, ocean)
%   rslSd = SOLVESLE(__, "Name", Value)
%
% Input arguments
%   forcing - SH coefficients of the forcing
%       The forcing will not be localised, i.e. signal in the ocean is
%       retained even when enforcing mass conservation (could be useful
%       for, e.g. GAD).
%       Unit: kg/m^2 (SD, equivalent to mm of fresh water)
%       Data type: DOUBLE
%       Dimension: [addmup(L), 4] (lmcosi) | [addmup(L) x 2] (cosi) |
%           [addmup(L) x 4 x N] (lmcosi, t) | [addmup(L) x 2 x N] (cosi, t)
%           (detected automatically, the GRACE2PLMT output format should be
%           set to 'traditional' (lmcosi, t))
%   forcingStd - Standard deviation of the forcing coefficients
%       The data type and dimension are the same as for forcing.
%   L - Degree to solve the SLE
%       The default degree is the degree of forcing (recommended to be
%       >= 96).
%       Data type: DOUBLE
%   ocean - Ocean domain to contruct the ocean function
%       Note this ocean domain should include, e.g. the Arctics and is not
%       the traditional sense of 'global' ocean in altimetry studies.
%       The default domain is all oceans (ALLOCEANS) with a buffer of 0.5°.
%       Data type: GEODOMAIN
%           or whatever is accepted by GEOBOXCAP and KERNELCP(_NEW)
%   Frame - Reference frame
%       Centre of mass (CM) frame or centre of figure (CF) frame.
%       The default frame is the CF frame, which is the frame used in
%       SLEPIAN_DELTA and most gravimetry and altimetry studies.
%       Data type: CHAR
%   LoveNumSource - Source of load Love numbers
%       - 'ISSM': Love numbers from the Ice Sheet System Model (ISSM)
%       - 'ALMA3 <model>': Love numbers from the ALMA3 models (currently
%           has to be added manually)
%       Note: The 'Wahr' Love number set used in the SLEPIAN_DELTA package
%       is not supported for SOLVESLE and will result in an error if used.
%       The default source is 'ISSM'.
%   RotationFeedback - Logical flag to include rotation feedback
%       The default option is true.
%       Data type: LOGICAL
%   Tolerance (name-value) - Absolute tolerance.
%       Estimated coefficient and standard-deviation vectors are checked
%       separately; every coefficient in every epoch must pass.
%       Set 0 to use exactly MaxIter steps.
%       Default absolute tolerance is 1e-6 mm or 1e-9 m.
%   MaxIter (name-value) - Maximum number of SLE iterations.
%       Default number of iterations is 100.
%   OceanKernel - Precomputed ocean-localization matrix.
%       Providing it avoids loading or computing the kernel on each call.
%       Its first column defines the ocean function and its (1,1) entry
%       defines the ocean area fraction used for mass conservation.
%   OceanFunction - Legacy precomputed ocean function, accepted but unused.
%       The kernel supplies consistent geometry for the mass constraint.
%   InitialRSL - Initial guess of the relative sea level.
%       Supplies a starting point for the SLE iteration so that an
%       available prior estimate can accelerate convergence.
%       Unit: m
%
% Output arguments
%   rslSd - SH coefficients of the mass load from the forcing
%       Note this value is only physical in the ocean domain.
%       Unit: kg/m^2 (SD, equivalent to mm of fresh water)
%       The dimension should be the same as the input forcing.
%   rslSdStd - Standard deviation of the mass load coefficients
%       When the standard deviation of the forcing is not provided, this
%       output is empty.
%       The data type and dimension should be the same as rslSd.
%   gmsl - Barystatic sea level from the forcing
%       Since the ocean domain is the true global ocean, the GMSL is not
%       comparible with values reported in latitude-bounded oceans.
%       Unit: mm
%       Dimension: scalar | [N x 1]
%           (depending on the input forcing)
%   gmslStd - Standard deviation of the barystatic sea level
%       When the standard deviation of the forcing is not provided, this
%       output is empty.
%       The data type and dimension is the same as gmsl.
%
% Notes
%   In the implementation, all SH coefficient variables with a 's' suffix
%   means they are flattened to a column vector.
%   Uncertainties use diagonal propagation: correlations between the forcing
%   and iterated ocean load, and between provisional response harmonics, are
%   neglected, as is covariance of the target mass with that response.
%
% References
%   Adhikari, S. et al. (2019). Sea-level fingerprints emergent from GRACE
%       mission data. Earth System Science Data, 11(2):629–646.
%       doi: 10.5194/essd-11-629-2019
%   Adhikari, S. et al. (2016). ISSM-SESAW v1.0: mesh-based computation of
%       gravitationally consistent sea-level and geodetic signatures caused
%       by cryosphere and climate driven mass change. Geoscientific Model
%       Development, 9(3):1087–1109. doi: 10.5194/gmd-9-1087-2016
%
% Author
%   2024/11/20, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified by
%   2026/10/26, En-Chi Lee (williameclee@arizona.edu)

function varargout = solvesle(varargin)
    %% Initialisation
    % Load all missing variables and convert to correct formats

    % Parse input arguments
    [forcSdPlm, forcSdStdPlm, inclDO, doError, ...
         L, ocean, frame, loveNumSrc, doRotation, maxIter, ...
         ocnKernel, ~, kernelOrder, initRslPlm, beQuiet, tol] = ...
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
        % harmonics. Apply that constraint before squaring the operator.
        ocnWghts = ocnFunPlms / ocnArea;
        cstredOcnKernel = ocnKernel - ocnWghts * ocnKernel(1, :);
        cstredOcnKernel(:, 1) = 0;
        cstredOceanKernelSq = cstredOcnKernel .^ 2;
        clear cstredOceanKernel
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

    % Many values from Adhikari et al. (2016)
    WATER_DENSITY = 1000; % to be consistent with other functions, mainly MASS2WEQ
    EARTH_DENSITY = 5517; % to be consistent with PLM2POT
    EARTH_RADIUS = 6371e3;
    EARTH_ANGULAR_VELOCITY = 7.2921e-5;
    EQUATORIAL_INERTIA = 8.0077e37;
    POLAR_INERTIA = 8.0345e37;
    CHANDLER_WOBBLE_FREQUENCY = 2.4405e-7;
    GRAVITY = 9.81;

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

    % Construct the SLE kernel
    % Gravitational potential
    sleKernel = spdiags(3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1), 0, ...
        numel(degrees), numel(degrees));

    if doRotation
        % Rotation potential
        c00c20Factor = 2/3 * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 ...
            * (- (1 + llnGeoid2) / POLAR_INERTIA) * 8 * pi / 3 * EARTH_RADIUS ^ 4;
        c21s21Factor = -1 / sqrt(15) * EARTH_RADIUS ^ 2 * EARTH_ANGULAR_VELOCITY ^ 2 ...
            * EARTH_ANGULAR_VELOCITY * (1 + llnGeoid2) ...
            / (EQUATORIAL_INERTIA * CHANDLER_WOBBLE_FREQUENCY) ...
            * (-4 * pi / sqrt(15)) * EARTH_RADIUS ^ 4;
        c21s21Factor = c21s21Factor / GRAVITY;
        % C00
        sleKernel(1, 1) = sleKernel(1, 1) + c00c20Factor;
        sleKernel(1, 5) = sleKernel(1, 5) + (-1 / sqrt(5)) * c00c20Factor;
        % C20
        sleKernel(5, 1) = sleKernel(5, 1) + (-1 / sqrt(5)) * c00c20Factor;
        sleKernel(5, 5) = sleKernel(5, 5) + (1/5) * c00c20Factor;
        % C21
        sleKernel(6, 6) = sleKernel(6, 6) + c21s21Factor * tlnFactors(6);
        sleKernel(7, 7) = sleKernel(7, 7) + c21s21Factor * tlnFactors(7);
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
            rslOcnVarPlms = cstredOceanKernelSq * rslVarPlms ...
                + ocnWghtsSq * targetOcnVar;
            rslOcnVarPlms(1, :) = targetOcnVar;
            rslVarPlms(1, :) = (targetOcnVar ...
                + ocnMeanWghtsSq * rslVarPlms(2:end, :)) / ocnArea ^ 2;
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
    ocnFunPlm = ip.Results.OceanFunction;
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
            error( ...
                sprintf('%s:InvalidInput:InputSizeNotMatch', upper(mfilename)), ...
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

    varargout = ...
        {forcPlm, forcStdPlm, inclDO, doError, ...
         L, ocnDom, frame, loveNumSrc, doRotation, maxIter, ...
         ocnKernel, ocnFunPlm, kernelOrder, initRslPlm, beQuiet, ip.Results.Tolerance};
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
