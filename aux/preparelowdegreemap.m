%% PREPARELOWDEGREEMAP - Builds the linear map from known to unknown loads.
%
% Syntax
%   [known2unknownMap, knownIds] = preparelowdegreemap(ocnKer, earthRespKer, forcKer, unknownIds)
%   [known2unknownMap, knownIds, ocnRespMap, unknownConstraint, ...
%       knownConstraint, reciprocalCondition, relativeResiduals, iterCnts] = ...
%       preparelowdegreemap(ocnKer, earthRespKer, forcKer, unknownIds)
%   known2unknownMap = preparelowdegreemap(ocnKer, [], [], unknownIds, Method="uniform")
%
% Input arguments
%   ocnKer - Ocean localisation matrix in kernel order.
%       Size: (M, M), M = (L + 1)^2.
%       Unit: dimensionless.
%   earthRespKer - Sparse G+R matrix that maps a load to its induced global-
%       extended RSL load.
%       Size: Same as ocnKer.
%       Unit: dimensionless.
%   forcKer - Mass-constrained forcing matrix from PREPARESLEKERNELS.
%       For Method="uniform", pass [] for earthRespKer and forcKer;
%       only ocnKer is used, with no Earth response or SLE solve.
%       Size: same as ocnKer.
%       Unit: dimensionless.
%   unknownIds - Distinct kernel-order indices of the estimated coefficients.
%       Order is retained in the output.
%       Size: (1, U).
%   Method (name-value) - Ocean constraint, fingerprint or uniform.
%       Default method: fingerprint.
%   RelativeTolerance (name-value) - Relative residual tolerance for SLE solves.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum iterations in each SLE system solve.
%       Default limit: 100.
%
% Output arguments
%   known2unknownMap - Maps known load coefficients to estimated coefficients.
%       Rows follow unknownIds; columns follow knownIds.
%       Size: (U, M - U).
%       Unit: dimensionless.
%   knownIds - Kernel-order indices of known coefficients, in ascending order.
%       Size: (1, M - U).
%   ocnRespMap - Maps a full load to its selected ocean response.
%       Rows follow unknownIds; columns use full kernel order.
%       Size: (U, M).
%       Unit: dimensionless.
%   unknownConstraint - Unknown-coefficient block of the ocean constraint.
%       Size: (U, U).
%       Unit: dimensionless.
%   knownConstraint - Known-coefficient block on the right-hand side.
%       Size: (U, M-U).
%       Unit: dimensionless.
%   reciprocalCondition - Reciprocal condition number of unknownConstraint.
%       Size: scalar.
%       Unit: dimensionless.
%   relativeResiduals - Relative residual of each transposed SLE solve.
%       Zero for uniform loading, which needs no SLE solve.
%       Size: (1, U).
%   iterCnts - GMRES iterations for each selected output direction.
%       Zero for uniform loading.
%       Size: (1, U).
%
% Notes
%   The selected ocean response to a full load is
%       D = EU^T * O * A^(-1) * KF * (I - O),
%       A = I - KS.
%   This is evaluated by solving for H^T, then forming D:
%       A^T * H^T = O^T * EU,
%       D = H * KF * (I - O).
%   where:
%     - O: Ocean localisation kernel (ocnKer)
%     - G+R: Earth's response kernel to load (earthRespKer)
%     - KS: RSL load's response kernel to itself, applied by SOLVESLESYSTEM
%     - KF: RSL load's response kernel to forcing (forcKer)
%     - EU: Unit vectors selecting the unknown coefficients (unknownIds)
%     - EU^T * O: Selected rows of the ocean kernel (selectedOcnKer)
%     - H^T: Adjoint response from the transposed solve (adjointResp)
%     - H * KF: Selected response to forcing (forcResp)
%     - D: Selected ocean response to a full load (ocnRespMap)
%
%   Matching the modelled and observed ocean coefficients gives
%       DU * LU = DnotU * LnotU,
%       DU = D * EU - OU,
%       DnotU = OnotU - D * EnotU,
%       T = DU \ DnotU,
%       LU = T * LnotU.
%   where:
%     - EnotU: Unit vectors selecting the known coefficients (knownIds)
%     - OU = EU^T * O * EU: selectedOcnKer(:,unknownIds)
%     - OnotU = EU^T * O * EnotU: selectedOcnKer(:,knownIds)
%     - DU: Unknown-coefficient constraint matrix (unknownConstraint)
%     - DnotU: Known-coefficient constraint matrix (knownConstraint)
%     - T: Map from known to unknown load coefficients (known2unknownMap)
%     - LU, LnotU: Unknown and known load coefficients, supplied by the caller
%   U and notU denote the unknown set and its complement. The selection
%   matrices EU and EnotU are applied by indexing, without constructing them.
%
%   For uniform loading, the selected ocean response is instead
%       D = -(EU^T * O * e00 / O00) * e00^T * (I - O).
%   where:
%     - e00^T * (I - O): Land C00 projection (landMassRow)
%     - O00: Ocean area fraction (ocnKer(1,1))
%     - e00: Unit vector with its first entry equal to 1
%   The same DU, DnotU and T construction applies to this response.
%
%   All epochs use the same T. Independent input variances propagate as
%       Var(LU) = (T.^2) * Var(LnotU),
%   where .^2 squares each entry. Reuse requires the same ocean matrix,
%   Earth response, bandwidth, coefficient selection and numerical accuracy.
%   Finite-bandwidth O is not assumed idempotent.
%
% See also
%   PREPARESLEKERNELS, SOLVESLESYSTEM, SOLVEDEGREE1
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/28, En-Chi Lee (williameclee@arizona.edu)

function [known2unknownMap, knownIds, ocnRespMap, unknownConstraint, knownConstraint, reciprocalCondition, relativeResiduals, iterCnts] = ...
        preparelowdegreemap(ocnKer, earthRespKer, forcKer, unknownIds, options)

    arguments (Input)
        ocnKer (:, :) double {mustBeReal}
        earthRespKer (:, :) double {mustBeReal}
        forcKer (:, :) double {mustBeReal}
        unknownIds (1, :) double {mustBeInteger, mustBePositive}
        options.Method {mustBeTextScalar, mustBeMember(options.Method, ["fingerprint", "uniform"])} = "fingerprint"
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end

    %% Initialisation
    numCoeffs = size(ocnKer, 1);

    if isempty(unknownIds) || numel(unique(unknownIds)) ~= numel(unknownIds) || any(unknownIds > numCoeffs)
        error('ULMO:preparelowdegreemap:InvalidIndices', ...
        'Unknown indices must be distinct and within the bandwidth.');
    end

    if numCoeffs == 0 || size(ocnKer, 2) ~= numCoeffs || ~isfinite(ocnKer(1, 1)) || ocnKer(1, 1) <= 0
        error('ULMO:preparelowdegreemap:InvalidocnKer', ...
        'Ocean kernel must be square with positive finite area.');
    end

    % Select the known coefficients and the ocean rows for the unknowns.
    knownIds = setdiff(1:numCoeffs, unknownIds);
    selectedOcnKer = ocnKer(unknownIds, :);
    relativeResiduals = zeros(1, numel(unknownIds));
    iterCnts = zeros(1, numel(unknownIds));

    %% Ocean response map
    if strcmpi(options.Method, "uniform")
        % Uniform loading: land C00 supplies the compensating ocean mass.
        landMassRow = -ocnKer(1, :);
        landMassRow(1) = landMassRow(1) + 1;
        ocnRespMap =- (ocnKer(unknownIds, 1) / ocnKer(1, 1)) * landMassRow;
    else

        if ~isequal(size(forcKer), size(ocnKer))
            error('ULMO:preparelowdegreemap:InvalidForcingKernel', ...
            'forcKer must match the ocean matrix size.');
        end

        % Solve A^T * H^T = O^T * EU for the selected output directions.
        [adjointResp, relativeResiduals, iterCnts] = ...
            solveslesystem(ocnKer, earthRespKer, selectedOcnKer.', Transpose = true, ...
            RelativeTolerance = options.RelativeTolerance, MaxIter = options.MaxIter);
        % D = H * KF * (I - O).
        forcResp = adjointResp.' * forcKer;
        ocnRespMap = forcResp - forcResp * ocnKer;
    end

    %% Ocean constraint
    % DU = D * EU - OU; DnotU = OnotU - D * EnotU.
    unknownConstraint = ocnRespMap(:, unknownIds) - selectedOcnKer(:, unknownIds);
    knownConstraint = selectedOcnKer(:, knownIds) - ocnRespMap(:, knownIds);
    reciprocalCondition = rcond(unknownConstraint);

    if ~isfinite(reciprocalCondition) || reciprocalCondition <= eps
        error('ULMO:preparelowdegreemap:SingularSystem', ...
        'The selected ocean constraint is singular to working precision.');
    elseif reciprocalCondition < 1e-10
        warning('ULMO:preparelowdegreemap:IllConditionedSystem', ...
            'The selected ocean constraint has reciprocal condition %.6g.', reciprocalCondition);
    end

    %% Load response map
    % T = DU \ DnotU, so LU = T * LnotU for each epoch.
    known2unknownMap = unknownConstraint \ knownConstraint;

    if any(~isfinite(known2unknownMap), 'all')
        error('ULMO:preparelowdegreemap:NonFiniteMap', 'The response map is nonfinite.');
    end

end
