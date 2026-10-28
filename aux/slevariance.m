%% SLEVARIANCE - Propagates independent forcing errors through the complete SLE.
%
% Syntax
%   var = slevariance(ocnKer, earthRespKer, forcKer, forcVar)
%
% Input arguments
%   ocnKer - Ocean localization matrix in kernel order.
%       Size: (M, M), M = (L + 1)^2.
%       Unit: dimensionless.
%   earthRespKer - G+R matrix from PREPARESLEKERNELS.
%       Size: same as ocnKer.
%       Unit: dimensionless.
%   forcKer - Mass-constrained forcing matrix from PREPARESLEKERNELS.
%       Size: same as ocnKer.
%       Unit: dimensionless.
%   forcVar - Independent input variances, with epochs in columns.
%       Size: (M, N).
%       Unit: (kg/m^2)^2.
%   RelativeTolerance (name-value) - SLE relative residual tolerance.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum GMRES iterations per right-hand side.
%       Default limit: 100.
%   BlockSize (name-value) - Number of response columns held at once.
%       Default size: 32.
%
% Output arguments
%   var - Variances of the globally extended RSL load coefficients.
%       Rows use kernel order; columns follow the input epochs.
%       Size: same as forcVar.
%       Unit: (kg/m^2)^2.
%
% Notes
%   For each epoch, the forcing and RSL load perturbations satisfy
%       A * deltaS = KF * deltaF,
%       A = I - KS,
%       T = A \ KF,
%       deltaS = T * deltaF.
%   The covariance and variance propagation are
%       C[S] = T * C[F] * T^T,
%       V[S] = (T.^2) * V[F],
%   where:
%     - O: Ocean localisation kernel (ocnKer)
%     - G+R: Earth's response kernel to load (earthRespKer)
%     - KS: RSL load's response kernel to itself, applied by SOLVESLESYSTEM
%     - KF: RSL load's response kernel to forcing (forcKer)
%     - deltaF, deltaS: Perturbations in the forcing and RSL load coefficients
%     - T: Full response from forcing to RSL load; not constructed in full
%     - C[F], C[S]: Input and output covariance matrices for one epoch
%     - V[F]: Input coefficient variances (forcVar)
%     - V[S]: Output coefficient variances (var)
%   The variance formula assumes C[F] is diagonal; .^2 squares each entry.
%   For multiple epochs, V[F] and V[S] have one column per epoch.
%
%   Compute the response in blocks of input coefficients J:
%       A * T(:,J) = KF(:,J),
%       V[S] = sum_J (T(:,J).^2) * V[F](J,:).
%   where:
%     - J: Input coefficient indices in the current block (ids)
%     - T(:,J): Response columns from the current solve (forcResp)
%     - KF(:,J): Right-hand sides for the current solve (forcKer(:,ids))
%     - V[F](J,:): Variances for those input coefficients (forcVar(ids,:))
%   actv contains coefficients with nonzero variance in at least one epoch.
%   The sum covers disjoint blocks of actv and accumulates into var.
%
%   Only input columns with nonzero variance are solved. Full-field uncertainty
%   can require M large solves; blocking limits memory, not total work.
%   Correlations induced by the SLE and mass constraint are retained.
%
% See also
%   SOLVESLE, SOLVESLESYSTEM, PREPARESLEKERNELS
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/28, En-Chi Lee (williameclee@arizona.edu)

function var = slevariance(ocnKer, earthRespKer, forcKer, forcVar, options)

    arguments (Input)
        ocnKer (:, :) double {mustBeReal}
        earthRespKer (:, :) double {mustBeReal}
        forcKer (:, :) double {mustBeReal}
        forcVar (:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.RelativeTolerance (1, 1) double ...
            {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBePositive, mustBeInteger} = 100
        options.BlockSize (1, 1) double {mustBePositive, mustBeInteger} = 32
    end

    arguments (Output)
        var (:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    %% Initialisation
    if ~isequal(size(forcKer), size(ocnKer), size(earthRespKer))
        error('ULMO:slevariance:InvalidSize', 'Kernels must have matching sizes.');
    end

    if size(forcVar, 1) ~= size(ocnKer, 1)
        error('ULMO:slevariance:InvalidSize', 'Variance must have one row per harmonic.');
    end

    var = zeros(size(forcVar));
    % Omit only coefficients with zero variance at every epoch.
    actv = find(any(forcVar ~= 0, 2));

    %% Variance propagation
    for first = 1:options.BlockSize:numel(actv)
        ids = actv(first:min(first + options.BlockSize - 1, numel(actv)));
        % Solve A * T(:,J) = KF(:,J) for the current input block.
        forcResp = solveslesystem(ocnKer, earthRespKer, forcKer(:, ids), ...
            RelativeTolerance = options.RelativeTolerance, MaxIter = options.MaxIter);
        % Add this block's contribution to V[S] at every epoch.
        var = var + (forcResp .^ 2) * forcVar(ids, :);
    end

end
