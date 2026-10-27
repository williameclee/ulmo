%% PREPARELOWDEGREEMAP - Builds the linear map from known to unknown loads.
%
% Syntax
%   [known2unknownMap, knownIds] = preparelowdegreemap(ocnKer, earthRespKer, forcKer, unknownIds)
%   [known2unknownMap, knownIds, oceanResponseMap, unknownConstraint, ...
%       knownConstraint, reciprocalCondition, relativeResiduals, iterationCounts] = ...
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
%       All kernels have size (M,M), M=(L+1)^2, and are dimensionless.
%       For Method="uniform", pass [] for earthRespKer and forcKer;
%       only ocnKer is used, with no Earth response or SLE solve.
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
%   known2unknownMap - Maps known load coefficients to estimated coefficients.
%       Rows follow unknownIds; columns follow knownIds.
%       Size: (U, M - U).
%       Unit: dimensionless.
%   knownIds - Kernel-order indices of known coefficients, in ascending order.
%       Size: (1, M - U).
%   oceanResponseMap - Maps a full load to its selected ocean response.
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
%   For fingerprint loading, D=(O(unknownIds,:)/A)*KF*(I-O).
%   oceanResponseMap is D in this expression. Only U transposed systems are
%   solved. known2unknownMap = unknownConstraint \ knownConstraint.
%   All epochs use known2unknownMap*knownLoads. Independent input variances
%   propagate as (known2unknownMap.^2)*knownVariances.
%   Reuse this map only with the same ocean matrix, Earth response, bandwidth,
%   coefficient selection and numerical accuracy. It is independent of epochs.
%   Finite-bandwidth O is not assumed idempotent.
%
% See also
%   PREPARESLEKERNELS, SOLVESLESYSTEM, SOLVEDEGREE1
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

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

    numCoeffs = size(ocnKer, 1);

    if isempty(unknownIds) || numel(unique(unknownIds)) ~= numel(unknownIds) || any(unknownIds > numCoeffs)
        error('ULMO:preparelowdegreemap:InvalidIndices', ...
        'Unknown indices must be distinct and within the bandwidth.');
    end

    if numCoeffs == 0 || size(ocnKer, 2) ~= numCoeffs || ~isfinite(ocnKer(1, 1)) || ocnKer(1, 1) <= 0
        error('ULMO:preparelowdegreemap:InvalidocnKer', ...
        'Ocean kernel must be square with positive finite area.');
    end

    knownIds = setdiff(1:numCoeffs, unknownIds);
    selectedOcnKer = ocnKer(unknownIds, :);
    relativeResiduals = zeros(1, numel(unknownIds));
    iterCnts = zeros(1, numel(unknownIds));

    if strcmpi(options.Method, "uniform")
        landMassRow = -ocnKer(1, :);
        landMassRow(1) = landMassRow(1) + 1;
        ocnRespMap =- (ocnKer(unknownIds, 1) / ocnKer(1, 1)) * landMassRow;
    else

        if ~isequal(size(forcKer), size(ocnKer))
            error('ULMO:preparelowdegreemap:InvalidforcKer', ...
            'forcKer must match the ocean matrix size.');
        end

        [adjointResp, relativeResiduals, iterCnts] = ...
            solveslesystem(ocnKer, earthRespKer, selectedOcnKer.', Transpose = true, ...
            RelativeTolerance = options.RelativeTolerance, MaxIter = options.MaxIter);
        forcResp = adjointResp.' * forcKer;
        ocnRespMap = forcResp - forcResp * ocnKer;
    end

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

    known2unknownMap = unknownConstraint \ knownConstraint;

    if any(~isfinite(known2unknownMap), 'all')
        error('ULMO:preparelowdegreemap:NonFiniteMap', 'The response map is nonfinite.');
    end

end
