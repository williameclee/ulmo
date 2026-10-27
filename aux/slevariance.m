%% SLEVARIANCE - Propagates independent forcing errors through the complete SLE.
%
% Syntax
%   variance = slevariance(oceanKernel, earthResponseKernel, forcingKernel, forcingVariance)
%
% Input arguments
%   oceanKernel - Ocean localization matrix in kernel order.
%   earthResponseKernel - G+R matrix from PREPARESLEKERNELS.
%   forcingKernel - Mass-constrained forcing matrix from PREPARESLEKERNELS.
%       All kernels have size ((L+1)^2,(L+1)^2) and are dimensionless.
%   forcingVariance - Independent input variances, size ((L+1)^2,N).
%       Unit: (kg/m^2)^2.
%   RelativeTolerance (name-value) - SLE relative residual tolerance.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum GMRES iterations per right-hand side.
%       Default limit: 100.
%   BlockSize (name-value) - Number of response columns held at once.
%       Default size: 32.
%
% Output arguments
%   variance - Output diagonal of forcingResponse*diag(forcingVariance)*forcingResponse.', where
%       forcingResponse=A\KF. Same size and units as forcingVariance.
%
% Notes
%   Only input columns with nonzero variance are solved. Full-field uncertainty
%   can require M large solves; blocking limits memory, not total work.
%   Correlations induced by the SLE and mass constraint are retained.
%
% See also
%   SOLVESLE, SOLVESLESYSTEM, PREPARESLEKERNELS
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function variance = slevariance(oceanKernel, earthResponseKernel, forcingKernel, forcingVariance, options)
    arguments (Input)
        oceanKernel (:, :) double {mustBeReal}
        earthResponseKernel (:, :) double {mustBeReal}
        forcingKernel (:, :) double {mustBeReal}
        forcingVariance (:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBePositive, mustBeInteger} = 100
        options.BlockSize (1, 1) double {mustBePositive, mustBeInteger} = 32
    end
    if ~isequal(size(forcingKernel), size(oceanKernel), size(earthResponseKernel))
        error('ULMO:slevariance:InvalidSize', 'Kernels must have matching sizes.');
    end
    if size(forcingVariance, 1) ~= size(oceanKernel, 1)
        error('ULMO:slevariance:InvalidSize', 'Variance must have one row per harmonic.');
    end
    variance = zeros(size(forcingVariance));
    active = find(any(forcingVariance ~= 0, 2));
    for first = 1:options.BlockSize:numel(active)
        ids = active(first:min(first+options.BlockSize-1, numel(active)));
        forcingResponse = solveslesystem(oceanKernel, earthResponseKernel, forcingKernel(:, ids), ...
            RelativeTolerance=options.RelativeTolerance, MaxIter=options.MaxIter);
        variance = variance + (forcingResponse.^2) * forcingVariance(ids, :);
    end
end
