%% SOLVESLEMASSCONSERVATIONTEST - Checks returned SLE mass and constrained errors.
% Uses a hemispheric ocean kernel and synthetic Love numbers in an isolated
% temporary data root. No production masks, data or caches are modified.
% Add slepian_alpha to the path before running this class.
%
% Syntax
%   results = runtests('unittests/solvesleMassConservationTest.m')
%
% Created by
%   2026/10/24, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

classdef solvesleMassConservationTest < matlab.unittest.TestCase

    properties (TestParameter)
        epochCount = {1, 3}
        iterations = {1, 5, 20}
        rotation = {false, true}
        sourceDegree = {1, 3}
        columnCount = {2, 4}
    end

    methods (TestClassSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            testCase.applyFixture(PathFixture(fullfile(repo, 'legacy')));
            testCase.assertNotEmpty(which('addmon'), 'Add slepian_alpha to the path.');
            testCase.assertEqual(which('solvesle'), fullfile(repo, 'solvesle.m'));
            folder = testCase.applyFixture(TemporaryFolderFixture);
            previous = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', previous));
            setenv('IFILES', folder.Folder);
            mkdir(fullfile(folder.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; 0.021; -0.3];
            loadingverticaldisplacement = [0; -0.269; -0.6];
            tidalgravitationalpotential = [0; 0; 0.3];
            tidalverticaldisplacement = [0; 0; 0.6];
            save(fullfile(folder.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential', 'loadingverticaldisplacement', ...
                'tidalgravitationalpotential', 'tidalverticaldisplacement');
        end

    end

    methods (Test)

        function returnedLoadConservesMass(testCase, epochCount, iterations, rotation)
            K = hemisphereKernel();
            f = (1:9)' * (1:epochCount);
            initial = flipud(f) / 1000;
            [r, empty, meanLevel] = solve(f, [], K, iterations, rotation, initial);
            ocean = K * r;
            testCase.verifyEqual(ocean(1, :), -f(1, :), 'AbsTol', 1e-11);
            testCase.verifyEqual(meanLevel, (-f(1, :) / K(1, 1))', 'AbsTol', 1e-11);
            testCase.verifyEmpty(empty);
        end

        function resizingPreservesEveryEpoch(testCase, sourceDegree, columnCount, epochCount)
            % Compare implicit padding/truncation with explicit target-degree
            % inputs for forcing, uncertainties and the initial RSL together.
            K = hemisphereKernel();
            n = (sourceDegree + 1) * (sourceDegree + 2) / 2;
            values = reshape(1:n * 2 * epochCount, n, 2, epochCount);
            input = {values, values / 10, values / 1000};
            expected = cell(1, 3);
            [order, degree] = addmon(sourceDegree);

            for j = 1:3
                expected{j} = zeros(6, 2, epochCount);
                count = min(n, 6);
                expected{j}(1:count, :, :) = input{j}(1:count, :, :);

                if columnCount == 4
                    input{j} = cat(2, repmat([degree, order], 1, 1, epochCount), input{j});
                end

            end

            opts = {"OceanKernel", K, "Tolerance", 0, "MaxIter", 5, "BeQuiet", true};
            actual = cell(1, 4); reference = cell(1, 4);
            [actual{:}] = solvesle_iteration(input{1}, input{2}, 2, GeoDomain('alloceans'), ...
                opts{:}, "InitialRSL", input{3});
            [reference{:}] = solvesle_iteration(expected{1}, expected{2}, 2, GeoDomain('alloceans'), ...
                opts{:}, "InitialRSL", expected{3});

            if columnCount == 4
                [targetOrder, targetDegree] = addmon(2);
                metadata = repmat([targetDegree, targetOrder], 1, 1, epochCount);
                testCase.verifyEqual(actual{1}(:, 1:2, :), metadata);
                testCase.verifyEqual(actual{2}(:, 1:2, :), metadata);
                actual{1} = actual{1}(:, 3:4, :);
                actual{2} = actual{2}(:, 3:4, :);
            end

            for j = 1:4
                testCase.verifyEqual(actual{j}, reference{j}, 'AbsTol', 1e-12);
            end

        end

        function convergesToIndependentLinearSolution(testCase)
            K = hemisphereKernel();
            f = reshape(sin(1:27), 9, 3);
            [~, ~, ~, ~, ~, ~, ~, degrees] = addmon(2);
            k = lovenumber(degrees, 'LLN geoid', 'CF');
            h = lovenumber(degrees, 'LLN VLM', 'CF');
            S = diag(3/5517 * (1 + k - h) ./ (2 * degrees + 1));
            C = eye(9); C(1, :) = C(1, :) - K(1, :) / K(1, 1);
            target = zeros(size(f)); target(1, :) = -f(1, :) / 1000 / K(1, 1);
            expected = (eye(9) - 1000 * C * S * K) \ (C * S * f + target) * 1000;
            actual = solve(f, [], K, 80, false, zeros(size(f)));
            otherStart = solve(f, [], K, 80, false, ones(size(f)));
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
            testCase.verifyEqual(otherStart, expected, 'AbsTol', 1e-10);
        end

        function independentOceanExpansionCannotChangeMass(testCase)
            K = hemisphereKernel(); f = (1:9)';
            reference = solve(f, [], K, 5, false, zeros(size(f)));
            p = unflatten(f); wrongOcean = ones(6, 2) * 0.9;
            % Hold initialisation fixed while varying the ocean expansion.
            r = solvesle_iteration(p, [], 2, GeoDomain('alloceans'), ...
                "OceanKernel", K, "OceanFunction", wrongOcean, ...
                "InitialRSL", zeros(size(p)), ...
                "RotationFeedback", false, "Tolerance", 0, "MaxIter", 5, "BeQuiet", true);
            testCase.verifyEqual(flatten(r), reference, 'AbsTol', 1e-12);
        end

        function marginalErrorsIncludeMassAdjustment(testCase)
            % With zero degree-zero forcing error and one nonrotating step,
            % independent input perturbations give an exact variance oracle.
            K = hemisphereKernel(); sigma = (0:8)' / 10;
            f = zeros(9, 1);
            [~, stdLoad, ~, stdMean] = solve(f, sigma, K, 1, false, f);
            responses = solve(diag(sigma), [], K, 1, false, zeros(9));
            expected = sqrt(sum(responses .^ 2, 2));
            testCase.verifyGreaterThan(expected(1), 0);
            testCase.verifyEqual(stdLoad, expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(stdMean, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(K(1, :) * responses, zeros(1, 9), 'AbsTol', 1e-12);
        end

        function standardDeviationsRemainNonnegative(testCase, epochCount)
            K = hemisphereKernel(); f = (1:9)' * (1:epochCount);
            sigma = abs(f) / 10;
            [r, stdLoad, ~, stdMean] = solve(f, sigma, K, 10, true, zeros(size(f)));
            testCase.verifyTrue(all(isfinite(stdLoad), 'all'));
            testCase.verifyGreaterThanOrEqual(stdLoad, zeros(size(stdLoad)));
            testCase.verifyEqual(stdMean, (sigma(1, :) / K(1, 1))', 'AbsTol', 1e-12);
            testCase.verifyEqual(K(1, :) * r, -f(1, :), 'AbsTol', 1e-11);
        end

        function invalidOceanAreaIsRejected(testCase)
            K = hemisphereKernel(); K(1, 1) = 0;
            testCase.verifyError(@() solve(zeros(9, 1), [], K, 1, false, zeros(9, 1)), ...
            'ULMO:solvesle:InvalidOceanArea');
        end

        function degreeOrderColumnsArePreserved(testCase)
            K = hemisphereKernel(); f = reshape(1:27, 9, 3); p = unflatten(f);
            [order, degree] = addmon(2);
            p = cat(2, repmat([degree, order], 1, 1, 3), p);
            r = solvesle_iteration(p, [], 2, GeoDomain('alloceans'), ...
                "OceanKernel", K, "OceanFunction", unflatten(K(:, 1)), "Tolerance", 0, "MaxIter", 3, "BeQuiet", true);
            testCase.verifyEqual(r(:, 1:2, :), p(:, 1:2, :));
            testCase.verifyEqual(K(1, :) * flatten(r(:, 3:4, :)), -f(1, :), 'AbsTol', 1e-11);
        end

    end

end

function K = hemisphereKernel()
    % Gauss-Legendre quadrature over x=cos(theta) in [0,1], all longitudes.
    n = 8; j = (1:n - 1)'; beta = j ./ sqrt(4 * j .^ 2 - 1);
    [V, X] = eig(diag(beta, 1) + diag(beta, -1));
    x = (diag(X) + 1) / 2; w = V(1, :)' .^ 2;
    phi = (0:15) * (2 * pi / 16);
    [~, ~, ~, ~, ~, ~, m, l] = addmon(2);
    Y = zeros(n * numel(phi), 9);

    for i = 1:9
        P = legendre(l(i), x', 'sch') * sqrt(2 * l(i) + 1);
        if m(i) > 0, trig = sin(m(i) * phi); else, trig = cos(m(i) * phi); end
        sample = P(abs(m(i)) + 1, :)' * trig;
        Y(:, i) = sample(:);
    end

    weights = repmat(w, numel(phi), 1) / (2 * numel(phi));
    K = Y' * (Y .* weights);
end

function [r, sigma, g, gStd] = solve(f, s, K, count, rotation, initial)
    if ~isempty(s), s = unflatten(s); end
    [p, stdP, g, gStd] = solvesle_iteration(unflatten(f), s, 2, GeoDomain('alloceans'), ...
        "OceanKernel", K, "OceanFunction", unflatten(K(:, 1)), ...
        "RotationFeedback", rotation, "Tolerance", 0, "MaxIter", count, ...
        "InitialRSL", unflatten(initial), "BeQuiet", true);
    r = flatten(p); sigma = [];
    if ~isempty(stdP), sigma = flatten(stdP); end
end

function p = unflatten(v)
    ko = kernelorder(2); p = zeros(12, size(v, 2)); p(ko, :) = v;
    p = reshape(p, 6, 2, size(v, 2));
end

function v = flatten(p)
    p = reshape(p, 12, []); v = p(kernelorder(2), :);
end
