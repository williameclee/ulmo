%% PREPARELOWDEGREEMAP - Builds the linear map from known to unknown loads.
%
% Syntax
%   [knownToUnknownMap, knownIds] = preparelowdegreemap(oceanKernel, earthResponseKernel, forcingKernel, unknownIds)
%   [knownToUnknownMap, knownIds, oceanResponseMap, unknownConstraint, ...
%       knownConstraint, reciprocalCondition, relativeResiduals, iterationCounts] = ...
%       preparelowdegreemap(oceanKernel, earthResponseKernel, forcingKernel, unknownIds)
%   knownToUnknownMap = preparelowdegreemap(oceanKernel, [], [], unknownIds, Method="uniform")
%
% Input arguments
%   oceanKernel - Ocean localization matrix in kernel order.
%   earthResponseKernel - G+R matrix from PREPARESLEKERNELS.
%   forcingKernel - Mass-constrained forcing matrix from PREPARESLEKERNELS.
%       All kernels have size (M,M), M=(L+1)^2, and are dimensionless.
%       For Method="uniform", pass [] for earthResponseKernel and forcingKernel;
%       only oceanKernel is used, with no Earth response or SLE solve.
%   unknownIds - Distinct kernel-order indices of the estimated coefficients.
%       Size: (1, U). Order is retained in the output.
%   Method (name-value) - Ocean constraint, fingerprint or uniform.
%       Default method: fingerprint.
%   RelativeTolerance (name-value) - Relative residual tolerance for SLE solves.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum iterations in each SLE system solve.
%       Default limit: 100.
%
% Output arguments
%   knownToUnknownMap - Maps known load coefficients to estimated coefficients.
%       Rows follow unknownIds; columns follow knownIds.
%       Size: (U, M-U). Unit: dimensionless.
%   knownIds - Kernel-order indices of known coefficients, in ascending order.
%       Size: (1, M-U).
%   oceanResponseMap - Maps a full load to its selected ocean response.
%       Rows follow unknownIds; columns use full kernel order.
%       Size: (U, M). Unit: dimensionless.
%   unknownConstraint - Unknown-coefficient block of the ocean constraint.
%       Size: (U, U). Unit: dimensionless.
%   knownConstraint - Known-coefficient block on the right-hand side.
%       Size: (U, M-U). Unit: dimensionless.
%   reciprocalCondition - Reciprocal condition number of unknownConstraint.
%       Size: scalar. Unit: dimensionless.
%   relativeResiduals - Relative residual of each transposed SLE solve.
%       Size: (1, U). Zero for uniform loading, which needs no SLE solve.
%   iterationCounts - GMRES iterations for each selected output direction.
%       Size: (1, U). Zero for uniform loading.
%
% Notes
%   For fingerprint loading, D=(O(unknownIds,:)/A)*KF*(I-O).
%   oceanResponseMap is D in this expression. Only U transposed systems are
%   solved. knownToUnknownMap = unknownConstraint \ knownConstraint.
%   All epochs use knownToUnknownMap*knownLoads. Independent input variances
%   propagate as (knownToUnknownMap.^2)*knownVariances.
%   Reuse this map only with the same ocean matrix, Earth response, bandwidth,
%   coefficient selection and numerical accuracy. It is independent of epochs.
%   Finite-bandwidth O is not assumed idempotent.
%
% See also
%   PREPARESLEKERNELS, SOLVESLESYSTEM, SOLVEDEGREE1
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function [knownToUnknownMap, knownIds, oceanResponseMap, unknownConstraint, ...
          knownConstraint, reciprocalCondition, relativeResiduals, iterationCounts] = ...
        preparelowdegreemap(oceanKernel, earthResponseKernel, forcingKernel, unknownIds, options)
    arguments (Input)
        oceanKernel (:, :) double {mustBeReal}
        earthResponseKernel (:, :) double {mustBeReal}
        forcingKernel (:, :) double {mustBeReal}
        unknownIds (1, :) double {mustBeInteger, mustBePositive}
        options.Method (1, 1) string {mustBeMember(options.Method, ["fingerprint", "uniform"])} = "fingerprint"
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end
    numCoeffs = size(oceanKernel, 1);
    if isempty(unknownIds) || numel(unique(unknownIds)) ~= numel(unknownIds) || any(unknownIds > numCoeffs)
        error('ULMO:preparelowdegreemap:InvalidIndices', 'Unknown indices must be distinct and within the bandwidth.');
    end
    if numCoeffs == 0 || size(oceanKernel, 2) ~= numCoeffs || ~isfinite(oceanKernel(1, 1)) || oceanKernel(1, 1) <= 0
        error('ULMO:preparelowdegreemap:InvalidOceanKernel', 'oceanKernel must be square with positive finite area.');
    end
    knownIds = setdiff(1:numCoeffs, unknownIds);
    selectedOceanKernel = oceanKernel(unknownIds, :);
    relativeResiduals = zeros(1, numel(unknownIds));
    iterationCounts = zeros(1, numel(unknownIds));
    if options.Method == "uniform"
        landMassRow = -oceanKernel(1, :);
        landMassRow(1) = landMassRow(1) + 1;
        oceanResponseMap = -(oceanKernel(unknownIds, 1) / oceanKernel(1, 1)) * landMassRow;
    else
        if ~isequal(size(forcingKernel), size(oceanKernel))
            error('ULMO:preparelowdegreemap:InvalidForcingKernel', ...
                'forcingKernel must match the ocean matrix size.');
        end
        [adjointResponse, solverDiagnostics] = solveslesystem(oceanKernel, earthResponseKernel, selectedOceanKernel.', Transpose=true, ...
            RelativeTolerance=options.RelativeTolerance, MaxIter=options.MaxIter);
        relativeResiduals = solverDiagnostics.RelativeResidual;
        iterationCounts = solverDiagnostics.Iterations;
        forcingResponse = adjointResponse.' * forcingKernel;
        oceanResponseMap = forcingResponse - forcingResponse * oceanKernel;
    end
    unknownConstraint = oceanResponseMap(:, unknownIds) - selectedOceanKernel(:, unknownIds);
    knownConstraint = selectedOceanKernel(:, knownIds) - oceanResponseMap(:, knownIds);
    reciprocalCondition = rcond(unknownConstraint);
    if ~isfinite(reciprocalCondition) || reciprocalCondition <= eps
        error('ULMO:preparelowdegreemap:SingularSystem', 'The selected ocean constraint is singular to working precision.');
    elseif reciprocalCondition < 1e-10
        warning('ULMO:preparelowdegreemap:IllConditionedSystem', ...
            'The selected ocean constraint has reciprocal condition %.6g.', reciprocalCondition);
    end
    knownToUnknownMap = unknownConstraint \ knownConstraint;
    if any(~isfinite(knownToUnknownMap), 'all')
        error('ULMO:preparelowdegreemap:NonFiniteMap', 'The response map is nonfinite.');
    end
end
