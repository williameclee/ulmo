%% solvesleConvergenceTest - Compares adaptive stopping with independent equations and fixed-step limits.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

classdef solvesleConvergenceTest < matlab.unittest.TestCase

    methods (TestClassSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
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

        function adaptiveAgreesWithIndependentSolution(testCase)
            K = hemisphereKernel();
            f = reshape(sin(1:27), 9, 3);
            [~, ~, ~, ~, ~, ~, ~, degree] = addmon(2);
            k = lovenumber(degree, 'LLN geoid', 'CF');
            h = lovenumber(degree, 'LLN VLM', 'CF');
            S = diag(3/5517 * (1 + k - h) ./ (2 * degree + 1));
            C = eye(9);
            C(1, :) = C(1, :) - K(1, :) / K(1, 1);
            b = zeros(size(f)); b(1, :) = -f(1, :) / 1000 / K(1, 1);
            exact = (eye(9) - 1000 * C * S * K) \ (C * S * f + b) * 1000;
            p = solvesle(unflatten(f), [], 2, GeoDomain('alloceans'), ...
                "OceanKernel", K, "RotationFeedback", false, "Tolerance", 1e-12, "BeQuiet", true);
            testCase.verifyEqual(flatten(p), exact, 'AbsTol', 1e-8);
        end

        function stoppingIsInvariantToEpochScaling(testCase)
            f = reshape(sin(1:27), 9, 3);
            sigma = abs(f) / 10;
            scales = [1e-9, 1, 1e9];
            opts = {"OceanKernel", hemisphereKernel(), "RotationFeedback", false, "BeQuiet", true};
            [a, sa] = solvesle(unflatten(f), unflatten(sigma), 2, GeoDomain('alloceans'), opts{:});
            [b, sb] = solvesle(unflatten(f .* scales), unflatten(sigma .* scales), ...
                2, GeoDomain('alloceans'), opts{:});
            testCase.verifyEqual(flatten(b) ./ scales, flatten(a), 'AbsTol', 1e-11);
            testCase.verifyEqual(flatten(sb) ./ scales, flatten(sa), 'AbsTol', 1e-11);
        end

        function checksEveryEpochAndUncertainties(testCase)
            K = hemisphereKernel(); f = zeros(9, 3); sigma = zeros(9, 3); sigma(:, 3) = (1:9)';
            opts = {"OceanKernel", K, "RotationFeedback", false, "BeQuiet", true};
            [p, s] = solvesle(unflatten(f), unflatten(sigma), 2, GeoDomain('alloceans'), ...
                opts{:}, "Tolerance", 1e-12);
            [reference, referenceStd] = solvesle(unflatten(f), unflatten(sigma), 2, GeoDomain('alloceans'), ...
                opts{:}, "Tolerance", 0, "MaxIter", 100);
            testCase.verifyEqual(p, reference, 'AbsTol', 1e-8);
            testCase.verifyEqual(s, referenceStd, 'AbsTol', 1e-8);
        end

        function capFailureAndFixedModeAreDistinct(testCase)
            p = unflatten((1:9)');
            opts = {"OceanKernel", hemisphereKernel(), "MaxIter", 1, "BeQuiet", true};
            testCase.verifyError(@() solvesle(p, [], 2, GeoDomain('alloceans'), opts{:}, "Tolerance", 1e-15), ...
            'ULMO:solvesle:NonConvergence');
        end

        function nonfiniteInputCannotConverge(testCase)
            p = zeros(6, 2);
            p(2, 1) = NaN;
            testCase.verifyError(@() solvesle(p, [], 2, GeoDomain('alloceans'), ...
                "OceanKernel", hemisphereKernel(), "BeQuiet", true), 'ULMO:solvesle:NonFiniteIterate');
        end

    end

end

function K = hemisphereKernel()
    % Gauss-Legendre quadrature over x=cos(theta) in [0,1], all longitudes.
    n = 8;
    j = (1:n - 1)';
    beta = j ./ sqrt(4 * j .^ 2 - 1);
    [V, X] = eig(diag(beta, 1) + diag(beta, -1));
    x = (diag(X) + 1) / 2;
    w = V(1, :)' .^ 2;
    phi = (0:15) * (2 * pi / 16);
    [~, ~, ~, ~, ~, ~, m, l] = addmon(2);
    Y = zeros(n * numel(phi), 9);

    for i = 1:9
        P = legendre(l(i), x', 'sch') * sqrt(2 * l(i) + 1);

        if m(i) > 0
            trig = sin(m(i) * phi);
        else
            trig = cos(m(i) * phi);
        end

        sample = P(abs(m(i)) + 1, :)' * trig;
        Y(:, i) = sample(:);
    end

    weights = repmat(w, numel(phi), 1) / (2 * numel(phi));
    K = Y' * (Y .* weights);
end

function p = unflatten(v)
    ko = kernelorder(2);
    p = zeros(12, size(v, 2));
    p(ko, :) = v;
    p = reshape(p, 6, 2, size(v, 2));
end

function v = flatten(p)
    p = reshape(p, 12, []);
    v = p(kernelorder(2), :);
end
