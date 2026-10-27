%% SLEVARIANCE - Propagates independent forcing errors through the complete SLE.
%
% Syntax
%   variance = slevariance(kernels, forcingVariance)
%
% Input arguments
%   kernels - Output of PREPARESLEKERNELS.
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
%   variance - Output diagonal of T*diag(forcingVariance)*T.', where
%       T=A\KF. Same size and units as forcingVariance.
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

function variance = slevariance(kernels, forcingVariance, options)
    arguments (Input)
        kernels (1, 1) struct
        forcingVariance (:, :) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBePositive, mustBeInteger} = 100
        options.BlockSize (1, 1) double {mustBePositive, mustBeInteger} = 32
    end
    if size(forcingVariance, 1) ~= size(kernels.O, 1)
        error('ULMO:slevariance:InvalidSize', 'Variance must have one row per harmonic.');
    end
    variance = zeros(size(forcingVariance));
    active = find(any(forcingVariance ~= 0, 2));
    for first = 1:options.BlockSize:numel(active)
        ids = active(first:min(first+options.BlockSize-1, numel(active)));
        T = solveslesystem(kernels, kernels.KF(:, ids), ...
            RelativeTolerance=options.RelativeTolerance, MaxIter=options.MaxIter);
        variance = variance + (T.^2) * forcingVariance(ids, :);
    end
end
