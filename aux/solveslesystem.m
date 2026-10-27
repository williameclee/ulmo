%% SOLVESLESYSTEM - Solves the SLE system or its transpose from numeric kernels.
%
% Syntax
%   [solutionCoefficients, relativeResiduals, iterationCounts] = ...
%       solveslesystem(ocnKer, earthRespKer, rightHandSides)
%   solutionCoefficients = solveslesystem(ocnKer, earthRespKer, ...
%       rightHandSides, Transpose=true)
%
% Input arguments
%   ocnKer - Ocean localisation matrix in kernel order.
%       Size: (M, M), M = (L + 1)^2. 
%       Unit: dimensionless.
%   earthRespKer - Sparse G+R matrix that maps a load to its induced global-
%       extended RSL load.
%       Size: Same as ocnKer. 
%       Unit: dimensionless.
%   rhss - Right-hand sides, with coefficients in rows and fields in columns.
%       Size: ((L+1)^2, N). 
%   Transpose (name-value) - Whether to solve the transposed SLE system.
%       Default option: false.
%   RelativeTolerance (name-value) - Relative Euclidean residual tolerance.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum GMRES iterations per right-hand side.
%       Default limit: 100.
%
% Output arguments
%   coeffs - Solutions with the same size and units as rhss.
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

function [coeffs, relativeResiduals, iterCnts] = solveslesystem( ...
        ocnKer, earthRespKer, rhss, options)

    arguments (Input)
        ocnKer (:, :) double {mustBeReal}
        earthRespKer (:, :) double {mustBeReal}
        rhss (:, :) double {mustBeReal, mustBeFinite}
        options.Transpose (1, 1) logical = false
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
    end

    numCoeffs = size(ocnKer, 1);

    if numCoeffs == 0 || size(ocnKer, 2) ~= numCoeffs || ...
            ~isequal(size(earthRespKer), size(ocnKer)) || ...
            ~isfinite(ocnKer(1, 1)) || ocnKer(1, 1) <= 0
        error('ULMO:solveslesystem:InvalidKernels', ...
        'Kernels must be matching square matrices with positive finite ocean area.');
    end

    if size(rhss, 1) ~= size(ocnKer, 1)
        error('ULMO:solveslesystem:InvalidSize', 'rightHandSides must have one row per harmonic.');
    end

    coeffs = zeros(size(rhss));
    relativeResiduals = zeros(1, size(rhss, 2));
    iterCnts = zeros(1, size(rhss, 2));
    applySystem = @(inputCoefficients) systemproduct( ...
        ocnKer, earthRespKer, inputCoefficients, options.Transpose);

    for rhsIdx = 1:size(rhss, 2)
        rhs = full(rhss(:, rhsIdx));
        rhsNorm = norm(rhs);

        if rhsNorm == 0
            continue
        end

        % Normalise before solving to make accuracy independent of load units.
        [normalizedSolution, convergenceFlag, ~, solverIterations] = ...
            gmres(applySystem, rhs / rhsNorm, [], ...
            options.RelativeTolerance, min(options.MaxIter, size(rhss, 1)));
        relativeResidual = norm(applySystem(normalizedSolution) - rhs / rhsNorm);

        if convergenceFlag ~= 0 || ~isfinite(relativeResidual) || relativeResidual > options.RelativeTolerance
            error('ULMO:solveslesystem:NonConvergence', ...
                'SLE system failed at right-hand side %d (relative residual %.6g).', rhsIdx, relativeResidual);
        end

        coeffs(:, rhsIdx) = normalizedSolution * rhsNorm;
        relativeResiduals(rhsIdx) = relativeResidual;
        iterCnts(rhsIdx) = solverIterations(end);
    end

end

function sysResp = systemproduct(ocnKer, earthRespKernel, inputCoeffs, doTranspose)

    if doTranspose
        projectedCoefficients = inputCoeffs - ocnKer(1, :).' * (inputCoeffs(1, :) / ocnKer(1, 1));
        sysResp = inputCoeffs - ocnKer.' * (earthRespKernel.' * projectedCoefficients);
    else
        loadResponse = earthRespKernel * (ocnKer * inputCoeffs);
        loadResponse(1, :) = loadResponse(1, :) - ocnKer(1, :) * loadResponse / ocnKer(1, 1);
        sysResp = inputCoeffs - loadResponse;
    end

end
