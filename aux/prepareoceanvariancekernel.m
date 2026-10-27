%% PREPAREOCEANVARIANCEKERNEL - Squares the mass-constrained ocean operator.
% Constructs the operator in column blocks without storing the complete
% unsquared constrained matrix.
%
% Syntax
%   Q = prepareoceanvariancekernel(K)
%
% Input arguments
%   K - Real double ocean-localisation kernel in kernel order.
%       Size: (n, n), where n = (L+1)^2.
%
% Output arguments
%   Q - Elementwise square of K - K(:,1)/K(1,1)*K(1,:), with Q(:,1)=0.
%       Applies diagonal variance propagation after the mass constraint.
%       Size: (n, n).
%
% See also
%   SOLVESLE, SOLVEDEGREE1, KERNELCP_NEW
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@gmail.com)

function Q = prepareoceanvariancekernel(K)

    arguments (Input)
        K (:, :) double {mustBeReal}
    end

    arguments (Output)
        Q (:, :) double {mustBeReal}
    end

    N = size(K, 1);

    if N == 0 || size(K, 2) ~= N
        error('ULMO:prepareoceanvariancekernel:InvalidSize', ...
        'Ocean kernel must be a nonempty square matrix.');
    end

    if ~isfinite(K(1, 1)) || K(1, 1) <= 0
        error('ULMO:prepareoceanvariancekernel:InvalidOceanArea', ...
        'Ocean kernel must have a finite, positive area in (1, 1).');
    end

    % Actual formula:
    % Qij = (Kij - Ki0 * K0j / K00) ^ 2 for j /= 0
    % Qi0 = 0 because of the mass conservation constraint
    % The second term comes from the spatially constant offset to conserve ocean mass
    % Ki0 / K00 is constant for every column, so precompute and store it in w
    % The strightforward formula is
    %   C = K - w * K(1, :);
    %   Q = C .^ 2;
    % but here we do the calculation in blocks to consume less memory for temporary arrays

    Q = zeros(N, N);
    w = K(:, 1) / K(1, 1);
    blockSize = 256;

    for jstart = 1:blockSize:N
        js = jstart:min(jstart + blockSize - 1, N);
        block = K(:, js) - w * K(1, js);

        if any(~isfinite(block), 'all')
            error('ULMO:prepareoceanvariancekernel:NonFiniteKernel', ...
            'The constrained ocean kernel must be finite.');
        end

        Q(:, js) = block .^ 2;
    end

    Q(:, 1) = 0;

end
