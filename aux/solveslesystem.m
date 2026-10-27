%% SOLVESLESYSTEM - Solves the SLE system or its transpose from numeric kernels.
%
% Syntax
%   [X, diagnostics] = solveslesystem(kernels, B)
%   X = solveslesystem(kernels, B, Transpose=true)
%
% Input arguments
%   kernels - Output of PREPARESLEKERNELS.
%   B - Right-hand sides, with coefficients in rows and fields in columns.
%       Size: ((L+1)^2, N). Unit: inherited by X.
%   Transpose (name-value) - Whether to solve A.'*X=B instead of A*X=B.
%       Default option: false.
%   RelativeTolerance (name-value) - Relative Euclidean residual tolerance.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum GMRES iterations per right-hand side.
%       Default limit: 100.
%
% Output arguments
%   X - Solutions with the same size and units as B.
%   diagnostics - RelativeResidual and Iterations, one per right-hand side.
%
% Notes
%   A=I-P*(G+R)*O. GMRES uses products of the supplied ordinary matrices;
%   no additional dense A, factorization, or inverse is constructed.
%   Each right-hand side is checked independently. Failure raises an error.
%   The forward and transpose products preserve the full rotation matrix.
%
% See also
%   PREPARESLEKERNELS, SOLVESLE, PREPARELOWDEGREEMAP
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function [X, diagnostics] = solveslesystem(kernels, B, options)
    arguments (Input)
        kernels (1, 1) struct
        B (:, :) double {mustBeReal, mustBeFinite}
        options.Transpose (1, 1) logical = false
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end
    if size(B, 1) ~= size(kernels.O, 1)
        error('ULMO:solveslesystem:InvalidSize', 'B must have one row per harmonic.');
    end
    X = zeros(size(B));
    diagnostics.RelativeResidual = zeros(1, size(B, 2));
    diagnostics.Iterations = zeros(1, size(B, 2));
    apply = @(v) systemproduct(kernels, v, options.Transpose);
    for j = 1:size(B, 2)
        b = full(B(:, j));
        scale = norm(b);
        if scale == 0
            continue
        end
        % Normalize before solving to make accuracy independent of load units.
        [x, flag, ~, iterations] = gmres(apply, b / scale, [], ...
            options.RelativeTolerance, min(options.MaxIter, size(B, 1)));
        residual = norm(apply(x) - b / scale);
        if flag ~= 0 || ~isfinite(residual) || residual > options.RelativeTolerance
            error('ULMO:solveslesystem:NonConvergence', ...
                'SLE system failed at right-hand side %d (relative residual %.6g).', j, residual);
        end
        X(:, j) = x * scale;
        diagnostics.RelativeResidual(j) = residual;
        diagnostics.Iterations(j) = iterations(end);
    end
end

function y = systemproduct(kernels, x, transpose)
    if transpose
        p = x - kernels.O(1, :).' * (x(1, :) / kernels.oceanArea);
        y = x - kernels.O.' * (kernels.response.' * p);
    else
        r = kernels.response * (kernels.O * x);
        r(1, :) = r(1, :) - kernels.O(1, :) * r / kernels.oceanArea;
        y = x - r;
    end
end
