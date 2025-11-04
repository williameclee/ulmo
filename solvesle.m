%% SOLVESLE
% Solves the (elastic) sea level equation (SLE) for a given forcing and
% ocean domain.
% The ocean function is assumed constant.
%
% Syntax
%   [rslLoad, [], gmsl, []] = SOLVESLE(forcing, [], L, ocean)
%   [rslLoad, rslLoadStd, gmsl, gmslStd] = SOLVESLE(forcing, forcingStd, L, ocean)
%   rslLoad = SOLVESLE(__, "Name", Value)
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
%   RotationFeedback - Logical flag to include rotation feedback
%       The default option is true.
%       Data type: LOGICAL
%   MaxIter - Maximum number of iterations
%       The default value is 10.
%       Data type: [NUMERIC]
%   OceanKernel/OceanFunction - Precomputed ocean kernel or ocean function
%       If this function is to be called recursively, providing these
%       variables prevents having to read these from the disk.
%
% Output arguments
%   rslLoad - SH coefficients of the mass load from the forcing
%       Note this value is only physical in the ocean domain.
%       Unit: kg/m^2 (SD, equivalent to mm of fresh water)
%       The dimension should be the same as the input forcing.
%   rslLoadStd - Standard deviation of the mass load coefficients
%       When the standard deviation of the forcing is not provided, this
%       output is empty.
%       The data type and dimension should be the same as rslLoad.
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
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/11/02, williameclee@arizona.edu (@williameclee)

function varargout = solvesle(varargin)
    %% Initialisation
    % Load all missing variables and convert to correct formats

    % Parse input arguments
    [forcingLoadPlm, forcingLoadStdPlm, includesDO, computeError, ...
         L, ocean, frame, doRotationFeedback, maxIter, ...
         oceanKernel, oceanFunPlm, kernelOrder, initialPlm, beQuiet] = ...
        parseinputs(varargin{:});

    % Create progress bar
    if ~beQuiet
        wbar = waitbar(0, 'Initialising', ...
            "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');
    end

    % Load the ocean function if not provided
    if isempty(oceanFunPlm)

        if ~beQuiet
            waitbar(0, wbar, 'Computing ocean function');
        end

        [~, ~, ~, ~, ~, oceanFunPlm] = geoboxcap(L, ocean, "BeQuiet", true);
    end

    % Remove the d/o columns if present
    if size(oceanFunPlm, 2) == 4
        oceanFunPlm = oceanFunPlm(:, 3:4);
    end

    % Load the ocean kernel if not provided
    if isempty(oceanKernel)

        if ~beQuiet
            waitbar(0, wbar, 'Computing ocean kernel');
        end

        oceanKernel = kernelcp_new(L, ocean, "BeQuiet", true);
    end

    % Precompute the squared ocean kernel for error computation
    oceanKernelSquare = oceanKernel .^ 2;

    % Load the permutation order from cosi -> flat
    if isempty(kernelOrder)
        kernelOrder = kernelorder(L);
    end

    kernelOrder = kernelOrder(:);

    % Truncate SH input to correct d/o
    if ismatrix(forcingLoadPlm)

        if size(forcingLoadPlm, 1) < addmup(L)
            forcingLoadPlm(addmup(L), 2) = 0;
            forcingLoadStdPlm(addmup(L), 2) = 0;
            initialPlm(addmup(L), 2) = 0;
        elseif size(forcingLoadPlm, 1) > addmup(L)
            forcingLoadPlm = forcingLoadPlm(1:addmup(L), :);
            forcingLoadStdPlm = forcingLoadStdPlm(1:addmup(L), :);
            initialPlm = initialPlm(1:addmup(L), :);
        end

        forcingLoadPlms = forcingLoadPlm(kernelOrder);
        forcingLoadStdPlms = forcingLoadStdPlm(kernelOrder);
        initialPlms = initialPlm(kernelOrder);
    else

        if size(forcingLoadPlm, 1) < addmup(L)
            forcingLoadPlm(addmup(L), 2, 1) = 0;
            forcingLoadStdPlm(addmup(L), 2, 1) = 0;
            initialPlm(addmup(L), 2, 1) = 0;
        elseif size(forcingLoadPlm, 1) > addmup(L)
            forcingLoadPlm = forcingLoadPlm(1:addmup(L), :, 1);
            forcingLoadStdPlm = forcingLoadStdPlm(1:addmup(L), :, 1);
            initialPlm = initialPlm(1:addmup(L), :, 1);
        end

        forcingLoadPlm = reshape(forcingLoadPlm, ...
            [prod(size(forcingLoadPlm, 1:2)), size(forcingLoadPlm, 3)]);
        forcingLoadStdPlm = reshape(forcingLoadStdPlm, ...
            [prod(size(forcingLoadStdPlm, 1:2)), size(forcingLoadStdPlm, 3)]);
        initialPlm = reshape(initialPlm, ...
            [prod(size(initialPlm, 1:2)), size(initialPlm, 3)]);
        forcingLoadPlms = forcingLoadPlm(kernelOrder, :);
        forcingLoadStdPlms = forcingLoadStdPlm(kernelOrder, :);
        initialPlms = initialPlm(kernelOrder, :);
    end

    oceanFunPlms = oceanFunPlm(kernelOrder);

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
    llnGeoids = lovenumber(degrees, 'LLN geoid', frame);
    llnVlms = lovenumber(degrees, 'LLN VLM', frame);
    llnFactors = 1 + llnGeoids - llnVlms;
    tlnGeoids = lovenumber(degrees, 'TLN geoid', frame);
    tlnVlms = lovenumber(degrees, 'TLN VLM', frame);
    tlnFactors = 1 + tlnGeoids - tlnVlms;
    llnGeoid2 = lovenumber(2, 'LLN geoid', frame);

    % Construct the SLE kernel
    % Gravitational potential
    sleKernel = sparse(diag(3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1)));

    if doRotationFeedback
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
    gmsl =- forcingLoadPlms(1, :) / oceanFunPlms(1, :) / WATER_DENSITY;

    if ~isempty(initialPlm)
        rslPlms = initialPlms;
    else
        rslPlms = oceanFunPlms * gmsl;
    end

    rslOceanPlms = rslPlms;

    if computeError
        gmslStd = ...
            -forcingLoadStdPlms(1, :) / oceanFunPlms(1, :) / WATER_DENSITY;
        rslOceanStdPlms = oceanFunPlms * gmslStd;
    end

    for iIter = 1:maxIter

        if ~beQuiet
            waitbar(iIter / maxIter, wbar, ...
                sprintf('Solving SLE iteratively (%d/%d)', iIter, maxIter));

            if getappdata(wbar, 'canceling')
                delete(wbar);
                error(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
                'Processing cancelled');
            end

        end

        loadPlms = forcingLoadPlms + rslOceanPlms * WATER_DENSITY;

        rslPlms = sleKernel * loadPlms;
        rslOceanPlms = oceanKernel * rslPlms;

        % Enforce mass conservation
        rslOceanPlms(1, :) = -forcingLoadPlms(1, :) / WATER_DENSITY;
        rslPlms(1, :) = rslOceanPlms(1, :) / oceanFunPlms(1, :);

        % Compute error if required
        if ~computeError
            continue
        end

        loadStdPlms = sqrt( ...
            forcingLoadStdPlms .^ 2 + (rslOceanStdPlms * WATER_DENSITY) .^ 2);
        rslStdPlms = sqrt((sleKernel .^ 2) * (loadStdPlms .^ 2));
        rslOceanStdPlms = sqrt(oceanKernelSquare * (rslStdPlms .^ 2));
        rslOceanStdPlms(1, :) = -forcingLoadStdPlms(1, :) / WATER_DENSITY;
        rslStdPlms(1, :) = rslOceanStdPlms(1, :) ./ oceanFunPlms(1, :);

    end

    %% Post-processing
    if ~beQuiet
        waitbar(1, wbar, 'Post-processing results');
    end

    rslLoadPlms = rslPlms * WATER_DENSITY;
    gmsl = rslOceanPlms(1, :) ./ oceanFunPlms(1, :) * 1e3; % m -> mm (= kg/m^2)

    if isvector(rslLoadPlms)
        rslLoadPlm = zeros([addmup(L), 2]);
        rslLoadPlm(kernelOrder) = rslLoadPlms;
    else
        rslLoadPlm = zeros([addmup(L) * 2, size(rslLoadPlms, 2)]);
        rslLoadPlm(kernelOrder, :) = rslLoadPlms;
        rslLoadPlm = ...
            reshape(rslLoadPlm, [addmup(L), 2, size(rslLoadPlms, 2)]);
    end

    if includesDO
        [order, degree] = addmon(L);

        if isvector(rslLoadPlms)
            rslLoadPlm = [degree, order, rslLoadPlm];
        else
            degree = repmat(degree, [1, 1, size(rslLoadPlm, 3)]);
            order = repmat(order, [1, 1, size(rslLoadPlm, 3)]);
            rslLoadPlm = cat(2, degree, order, rslLoadPlm);
        end

    end

    if ~computeError
        varargout = {rslLoadPlm, [], gmsl(:), []};

        if ~beQuiet
            delete(wbar);
        end

        return
    end

    % Post-process the standard deviation
    rslLoadStdPlms = rslStdPlms * WATER_DENSITY;
    gmslStd = rslOceanStdPlms(1, :) ./ oceanFunPlms(1, :) * 1e3; % m -> mm

    if isvector(rslLoadStdPlms)
        rslLoadStdPlm = zeros([addmup(L), 2]);
        rslLoadStdPlm(kernelOrder) = rslLoadStdPlms;
    else
        rslLoadStdPlm = zeros([addmup(L) * 2, size(rslLoadStdPlms, 2)]);
        rslLoadStdPlm(kernelOrder, :) = rslLoadStdPlms;
        rslLoadStdPlm = ...
            reshape(rslLoadStdPlm, [addmup(L), 2, size(rslLoadStdPlms, 2)]);
    end

    if includesDO

        if isvector(rslLoadStdPlms)
            rslLoadStdPlm = [degree, order, rslLoadStdPlm];
        else
            rslLoadStdPlm = cat(2, degree, order, rslLoadStdPlm);
        end

    end

    varargout = {rslLoadPlm, rslLoadStdPlm, gmsl(:), gmslStd(:)};

    if ~beQuiet
        delete(wbar);
    end

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
    addOptional(ip, 'maxIter', 10, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'OceanKernel', [], @(x) isnumeric(x));
    addParameter(ip, 'OceanFunction', [], @(x) isnumeric(x));
    addParameter(ip, 'KernelOrder', [], @(x) isnumeric(x));
    addParameter(ip, 'InitialCondition', [], @(x) isnumeric(x));
    addParameter(ip, 'BeQuiet', false, @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    parse(ip, varargin{:});
    forcingPlm = squeeze(ip.Results.ForcingPlm);
    forcingStdPlm = squeeze(ip.Results.ForcingStdPlm);
    L = ip.Results.L;
    ocean = ip.Results.Ocean;
    frame = upper(ip.Results.frame);
    doRotationFeedback = ip.Results.RotationFeedback;
    maxIter = ip.Results.maxIter;
    oceanKernel = ip.Results.OceanKernel;
    oceanFunSph = ip.Results.OceanFunction;
    kernelOrder = ip.Results.KernelOrder;
    initialPlm = ip.Results.InitialCondition;
    beQuiet = logical(ip.Results.BeQuiet);

    % Check input sizes
    includesDO = size(forcingPlm, 2) == 4;

    if isempty(forcingStdPlm) || ...
            all(forcingStdPlm(:, end - 1:end) == 0, "all")
        forcingStdPlm = zeros(size(forcingPlm));
        computeError = false;
    else

        if ~isequal(size(forcingStdPlm), size(forcingPlm))
            error( ...
                sprintf('%s:InvalidInput:InputSizeNotMatch', upper(mfilename)), ...
                'Input size of FORCINGSTDPLM ([%s]) does not match FORCINGPLM ([%s])', ...
                strjoin(string(size(forcingStdPlm)), ', '), ...
                strjoin(string(size(forcingPlm)), ', '));
        end

        computeError = true;

    end

    if isempty(initialPlm) || ...
            all(initialPlm(:, end - 1:end) == 0, "all")
        initialPlm = zeros(size(forcingPlm));
    else

        if ~isequal(size(initialPlm), size(forcingPlm))
            error( ...
                sprintf('%s:InvalidInput:InputSizeNotMatch', upper(mfilename)), ...
                'Input size of INITIAL CONDITION ([%s]) does not match FORCINGPLM ([%s])', ...
                strjoin(string(size(initialPlm)), ', '), ...
                strjoin(string(size(forcingPlm)), ', '));
        end

    end

    % Find the degree of the spherical harmonic coefficients
    if includesDO

        if isempty(L)
            L = max(forcingPlm(:, 1, 1));
        end

        forcingPlm = forcingPlm(:, 3:4, :);
        forcingStdPlm = forcingStdPlm(:, 3:4, :);

        if ~isempty(initialPlm)
            initialPlm = initialPlm(:, 3:4, :);
        end

    elseif isempty(L)
        L = finddegree(forcingPlm(:, :, 1));
    end

    varargout = ...
        {forcingPlm, forcingStdPlm, includesDO, computeError, ...
         L, ocean, frame, doRotationFeedback, maxIter, ...
         oceanKernel, oceanFunSph, kernelOrder, initialPlm, beQuiet};
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
