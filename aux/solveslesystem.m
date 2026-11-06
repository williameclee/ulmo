%% SOLVESLESYSTEM - Solves the SLE system or its transpose from numeric kernels.
% The forward system gives the RSL load induced by a given forcing.
% The transposed system gives an auxiliary response used to prepare the
% geocentre and other low-degree calculations. Both use GMRES.
%
% Syntax
%   [rslLoads, rltvResids, iterCnts] = solveslesystem(ocnKer, earthRespKer, forcIndLoads)
%   adjointResp = solveslesystem(ocnKer, earthRespKer, unknOcnKer.', Transpose=true)
%
% Input arguments
%   ocnKer - Ocean localisation matrix in kernel order.
%       Size: (M, M), M = (L + 1)^2.
%       Unit: dimensionless.
%   earthRespKer - Sparse Earth response matrix G+R.
%       Maps a load to its globally extended RSL load before mass adjustment.
%       Size: same as ocnKer.
%       Unit: dimensionless.
%   forcIndLoads - Right-hand sides, with coefficients in rows and fields in columns.
%       In the forward system, supply KF * F.
%       With Transpose=true, supply unknOcnKer.' = O^T * EU.
%       Size: (M, N); N = U for the selected ocean response.
%   Transpose (name-value) - Whether to solve the transposed SLE system.
%       Default option: false.
%   RelativeTolerance (name-value) - Relative Euclidean residual tolerance.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum GMRES iterations per right-hand side.
%       Default limit: 100.
%
% Output arguments
%   rslLoads - Solution coefficients S in the forward system.
%       With Transpose=true, returns H^T (adjointResp in the syntax above).
%       Size: same as the supplied right-hand sides.
%       Unit: same as the supplied right-hand sides.
%   rltvResids - Relative Euclidean residuals, one per right-hand side.
%       Zero for a zero right-hand side.
%       Size: (1, N).
%       Unit: dimensionless.
%   iterCnts - GMRES iteration counts, one per right-hand side.
%       Zero for a zero right-hand side.
%       Size: (1, N).
%
% Notes
%   The forward system is
%       A * S = KF * F,
%       A = I - KS,
%       KS = (I - e00 * e00^T * O / O00) * (G + R) * O.
%   where:
%     - O: Ocean localisation kernel (ocnKer)
%     - G+R: Earth's response kernel to load (earthRespKer)
%     - KS: RSL load's response kernel to itself
%     - S: RSL load coefficients (rslLoads)
%     - KF * F: Forcing-induced RSL load, supplied as the RHS (forcIndLoads)
%     - O00: Ocean area fraction (ocnKer(1,1))
%     - e00: Unit vector with its first entry equal to 1
%
%   For geocentre and other unknown coefficients, the required response is
%       D = EU^T * O * A^(-1) * KF * (I - O).
%   Define H = EU^T * O * A^(-1). Then
%       H * A = EU^T * O,
%       A^T * H^T = O^T * EU.
%   The second form lets GMRES solve for each column of H^T. Afterwards
%   (see PREPARELOWDEGREEMAP),
%       D = H * KF * (I - O).
%
%   Using GMRES gives the adventage that neither A = (I - KS) nor its 
%   inverse is constructed.
%
% See also
%   PREPARESLEKERNELS, SOLVESLE, PREPARELOWDEGREEMAP
%
% Created by
%   2026/10/28, En-Chi Lee (williameclee@arizona.edu)

function [rslLoads, rltvResids, iterCnts] = solveslesystem( ...
        ocnKer, earthRespKer, forcIndLoads, options)

    arguments (Input)
        ocnKer (:, :) double {mustBeReal}
        earthRespKer (:, :) double {mustBeReal}
        forcIndLoads (:, :) double {mustBeReal, mustBeFinite}
        options.Transpose (1, 1) logical = false
        options.RelativeTolerance (1, 1) double ...
            {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end

    %% Initialisation
    nCoeffs = size(ocnKer, 1);

    if nCoeffs == 0 || size(ocnKer, 2) ~= nCoeffs || ...
            ~isequal(size(earthRespKer), size(ocnKer)) || ...
            ~isfinite(ocnKer(1, 1)) || ocnKer(1, 1) <= 0
        error('ULMO:solveslesystem:InvalidKernels', ...
        'Kernels must be matching square matrices with positive finite ocean area.');
    end

    if size(forcIndLoads, 1) ~= size(ocnKer, 1)
        error('ULMO:solveslesystem:InvalidSize', 'rightHandSides must have one row per harmonic.');
    end

    rslLoads = zeros(size(forcIndLoads));
    rltvResids = zeros(1, size(forcIndLoads, 2));
    iterCnts = zeros(1, size(forcIndLoads, 2));
    applySystem = @(rslLoads) systemproduct( ...
        ocnKer, earthRespKer, rslLoads, options.Transpose);

    %% Linear solve
    for iField = 1:size(forcIndLoads, 2)
        forcIndLoad = full(forcIndLoads(:, iField));
        forcIndLoadNorm = norm(forcIndLoad);

        if forcIndLoadNorm == 0
            continue
        end

        % Normalise the RHS and restore its scale after solving.
        % applySystem multiplies by A or A^T without constructing the matrix.
        [normRslLoad, convFlag, ~, solverIterCnts] = ...
            gmres(applySystem, forcIndLoad / forcIndLoadNorm, [], ...
            options.RelativeTolerance, min(options.MaxIter, size(forcIndLoads, 1)));
        rltvResid = norm(applySystem(normRslLoad) - forcIndLoad / forcIndLoadNorm);

        if convFlag ~= 0 || ~isfinite(rltvResid) || rltvResid > options.RelativeTolerance
            error('ULMO:solveslesystem:NonConvergence', ...
                'SLE system failed at right-hand side %d (relative residual %.6g).', ...
                iField, rltvResid);
        end

        rslLoads(:, iField) = normRslLoad * forcIndLoadNorm;
        rltvResids(iField) = rltvResid;
        iterCnts(iField) = solverIterCnts(end);
    end

end

%% SYSTEMPRODUCT - Applies the SLE system or its transpose.
%
% Notes
%   The forward product applies I - KS:
%       SS' = (G + R) * O * S,
%       SS = KS * S = SS' - e00 * (e00^T * O * SS') / O00,
%       (I - KS) * S = S - SS.
%   where:
%     - O: Ocean localisation kernel (ocnKer)
%     - G+R: Earth's response kernel to load (earthRespKer)
%     - S: Input RSL load coefficients (rslLoads)
%     - SS': Self-induced RSL load before mass adjustment (selfIndLoads)
%     - SS: Self-induced RSL load after mass adjustment (selfIndLoads)
%     - (I - KS) * S: System response (forcIndLoads); equals KF * F at the solution
%     - O00: Ocean area fraction (ocnKer(1,1))
%     - e00: Unit vector with its first entry equal to 1
%   The prime in SS' labels the response before mass adjustment, not a transpose.
%
%   Transposition reverses the order of the factors in KS:
%       P = I - e00 * e00^T * O / O00,
%       KS^T = O^T * (G + R)^T * P^T.
%   The transposed product therefore applies P^T first:
%       Z = H^T - O^T * e00 * (e00^T * H^T) / O00,
%       (I - KS)^T * H^T = H^T - O^T * (G + R)^T * Z.
%   where:
%     - P: Ocean mass adjustment, applied without constructing its matrix
%     - H^T: Input adjoint response coefficients (adjResp)
%     - Z: Adjoint coefficients after mass adjustment (prjAdjResp)
%     - (I - KS)^T * H^T: System response (adjSysResp); equals O^T * EU at the solution
%
%   Formulas use unnormalised coefficients. During GMRES, the inputs can
%   be search vectors; these linear products apply at any scale.
function varargout = systemproduct(varargin)
    ocnKer = varargin{1};
    earthRespKer = varargin{2};
    doTranspose = varargin{4};

    if doTranspose
        adjResp = varargin{3};
        prjAdjResp = adjResp - ocnKer(1, :).' * (adjResp(1, :) / ocnKer(1, 1));
        adjSysResp = adjResp - ocnKer.' * (earthRespKer.' * prjAdjResp);
        varargout = {adjSysResp};
        return
    end

    rslLoads = varargin{3};
    selfIndLoads = earthRespKer * (ocnKer * rslLoads);
    selfIndLoads(1, :) = selfIndLoads(1, :) - ocnKer(1, :) * selfIndLoads / ocnKer(1, 1);
    forcIndLoads = rslLoads - selfIndLoads;
    varargout = {forcIndLoads};
end
