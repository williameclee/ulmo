%% SOLVESLESYSTEM - Solves the SLE system or its transpose from numeric kernels.
%
% Syntax
%   [solutionCoefficients, relativeResiduals, iterationCounts] = ...
%       solveslesystem(oceanKernel, earthResponseKernel, rightHandSides)
%   solutionCoefficients = solveslesystem(oceanKernel, earthResponseKernel, ...
%       rightHandSides, Transpose=true)
%
% Input arguments
%   oceanKernel - Ocean localization matrix in kernel order.
%       Size: (M,M), M=(L+1)^2. Unit: dimensionless.
%   earthResponseKernel - Sparse G+R matrix from PREPARESLEKERNELS.
%       Size: (M,M). Unit: dimensionless.
%   rightHandSides - Right-hand sides, with coefficients in rows and fields in columns.
%       Size: ((L+1)^2, N). Unit: inherited by solutionCoefficients.
%   Transpose (name-value) - Whether to solve the transposed SLE system.
%       Default option: false.
%   RelativeTolerance (name-value) - Relative Euclidean residual tolerance.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum GMRES iterations per right-hand side.
%       Default limit: 100.
%
% Output arguments
%   solutionCoefficients - Solutions with the same size and units as rightHandSides.
%   relativeResiduals - Relative Euclidean residuals, one per right-hand side.
%       Size: (1,N). Zero for a zero right-hand side.
%   iterationCounts - GMRES iteration counts, one per right-hand side.
%       Size: (1,N). Zero for a zero right-hand side.
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

function [solutionCoefficients, relativeResiduals, iterationCounts] = ...
        solveslesystem(oceanKernel, earthResponseKernel, rightHandSides, options)
    arguments (Input)
        oceanKernel (:, :) double {mustBeReal}
        earthResponseKernel (:, :) double {mustBeReal}
        rightHandSides (:, :) double {mustBeReal, mustBeFinite}
        options.Transpose (1, 1) logical = false
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end
    numCoeffs = size(oceanKernel, 1);
    if numCoeffs == 0 || size(oceanKernel, 2) ~= numCoeffs || ...
            ~isequal(size(earthResponseKernel), size(oceanKernel)) || ...
            ~isfinite(oceanKernel(1, 1)) || oceanKernel(1, 1) <= 0
        error('ULMO:solveslesystem:InvalidKernels', ...
            'Kernels must be matching square matrices with positive finite ocean area.');
    end
    if size(rightHandSides, 1) ~= size(oceanKernel, 1)
        error('ULMO:solveslesystem:InvalidSize', 'rightHandSides must have one row per harmonic.');
    end
    solutionCoefficients = zeros(size(rightHandSides));
    relativeResiduals = zeros(1, size(rightHandSides, 2));
    iterationCounts = zeros(1, size(rightHandSides, 2));
    applySystem = @(inputCoefficients) systemproduct( ...
        oceanKernel, earthResponseKernel, inputCoefficients, options.Transpose);
    for rightHandSideIndex = 1:size(rightHandSides, 2)
        rightHandSide = full(rightHandSides(:, rightHandSideIndex));
        rightHandSideNorm = norm(rightHandSide);
        if rightHandSideNorm == 0
            continue
        end
        % Normalize before solving to make accuracy independent of load units.
        [normalizedSolution, convergenceFlag, ~, solverIterations] = ...
            gmres(applySystem, rightHandSide / rightHandSideNorm, [], ...
            options.RelativeTolerance, min(options.MaxIter, size(rightHandSides, 1)));
        relativeResidual = norm(applySystem(normalizedSolution) - rightHandSide / rightHandSideNorm);
        if convergenceFlag ~= 0 || ~isfinite(relativeResidual) || relativeResidual > options.RelativeTolerance
            error('ULMO:solveslesystem:NonConvergence', ...
                'SLE system failed at right-hand side %d (relative residual %.6g).', rightHandSideIndex, relativeResidual);
        end
        solutionCoefficients(:, rightHandSideIndex) = normalizedSolution * rightHandSideNorm;
        relativeResiduals(rightHandSideIndex) = relativeResidual;
        iterationCounts(rightHandSideIndex) = solverIterations(end);
    end
end

function systemResponse = systemproduct(oceanKernel, earthResponseKernel, inputCoefficients, doTranspose)
    if doTranspose
        projectedCoefficients = inputCoefficients - oceanKernel(1, :).' * (inputCoefficients(1, :) / oceanKernel(1, 1));
        systemResponse = inputCoefficients - oceanKernel.' * (earthResponseKernel.' * projectedCoefficients);
    else
        loadResponse = earthResponseKernel * (oceanKernel * inputCoefficients);
        loadResponse(1, :) = loadResponse(1, :) - oceanKernel(1, :) * loadResponse / oceanKernel(1, 1);
        systemResponse = inputCoefficients - loadResponse;
    end
end
