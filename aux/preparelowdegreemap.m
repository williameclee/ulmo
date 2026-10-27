%% PREPARELOWDEGREEMAP - Builds the linear map from known to unknown loads.
%
% Syntax
%   map = preparelowdegreemap(kernels, unknownIds)
%   map = preparelowdegreemap(kernels, unknownIds, Method="uniform")
%
% Input arguments
%   kernels - Output of PREPARESLEKERNELS. For Method="uniform", only O
%       is required; no Earth response or SLE solve is used.
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
%   map - Numeric response matrices and index metadata.
%       D: selected ocean response to a full load, size (U,M).
%       Dunknown, Dknown: the small-system blocks, sizes (U,U), (U,M-U).
%       T: maps known loads to estimated loads, size (U,M-U).
%       unknownIds, knownIds: row and column index order for T.
%       reciprocalCondition: rcond(Dunknown); diagnostics: SLE solve accuracy.
%
% Notes
%   For fingerprint loading, D=(O(unknownIds,:)/A)*KF*(I-O).
%   Only U transposed systems are solved. All epochs then use T*knownLoads.
%   Independent input variances propagate as (T.^2)*knownVariances.
%   Reuse this map only with the same ocean matrix, Earth response, bandwidth,
%   coefficient selection and numerical accuracy. It is independent of epochs.
%   Finite-bandwidth O is not assumed idempotent.
%
% See also
%   PREPARESLEKERNELS, SOLVESLESYSTEM, SOLVEDEGREE1
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function map = preparelowdegreemap(kernels, unknownIds, options)
    arguments (Input)
        kernels (1, 1) struct
        unknownIds (1, :) double {mustBeInteger, mustBePositive}
        options.Method (1, 1) string {mustBeMember(options.Method, ["fingerprint", "uniform"])} = "fingerprint"
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end
    O = kernels.O;
    M = size(O, 1);
    if isempty(unknownIds) || numel(unique(unknownIds)) ~= numel(unknownIds) || any(unknownIds > M)
        error('ULMO:preparelowdegreemap:InvalidIndices', 'Unknown indices must be distinct and within the bandwidth.');
    end
    if size(O, 2) ~= M || ~isfinite(O(1, 1)) || O(1, 1) <= 0
        error('ULMO:preparelowdegreemap:InvalidOceanKernel', 'O must be square with positive finite area.');
    end
    knownIds = setdiff(1:M, unknownIds);
    Orows = O(unknownIds, :);
    diagnostics = struct;
    if options.Method == "uniform"
        landC00 = -O(1, :);
        landC00(1) = landC00(1) + 1;
        D = -(O(unknownIds, 1) / O(1, 1)) * landC00;
    else
        [Ht, diagnostics] = solveslesystem(kernels, Orows.', Transpose=true, ...
            RelativeTolerance=options.RelativeTolerance, MaxIter=options.MaxIter);
        H = Ht.' * kernels.KF;
        D = H - H * O;
    end
    Dunknown = D(:, unknownIds) - Orows(:, unknownIds);
    Dknown = Orows(:, knownIds) - D(:, knownIds);
    reciprocalCondition = rcond(Dunknown);
    if ~isfinite(reciprocalCondition) || reciprocalCondition <= eps
        error('ULMO:preparelowdegreemap:SingularSystem', 'The selected ocean constraint is singular to working precision.');
    elseif reciprocalCondition < 1e-10
        warning('ULMO:preparelowdegreemap:IllConditionedSystem', ...
            'The selected ocean constraint has reciprocal condition %.6g.', reciprocalCondition);
    end
    T = Dunknown \ Dknown;
    if any(~isfinite(T), 'all')
        error('ULMO:preparelowdegreemap:NonFiniteMap', 'The response map is nonfinite.');
    end
    map = struct('D', D, 'Dunknown', Dunknown, 'Dknown', Dknown, 'T', T, ...
        'unknownIds', unknownIds, 'knownIds', knownIds, ...
        'reciprocalCondition', reciprocalCondition, 'diagnostics', diagnostics);
end
