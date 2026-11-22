%% GLMALPHAPOLARCAPTEST - Checks polar-cap concentration and cache handling.
%
% Created by
%   2026/11/22, En-Chi Lee (williameclee@arizona.edu)

classdef glmalphaPolarCapTest < matlab.unittest.TestCase

    properties
        Folder
    end

    methods (TestMethodSetup)

        function setup(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root, IncludingSubfolders = true));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', testCase.Folder));
        end

    end

    methods (Test)

        function lowpassDispatchUsesScalarRadius(testCase)
            [G, V, ~, ~, N] = glmalpha_new(60, 4, SaveData = false, BeQuiet = true);
            testCase.verifySize(G, [25, 25]);
            testCase.verifyEqual(G.' * G, eye(25), AbsTol = 1e-12);
            testCase.verifyNumElements(V, 25);
            testCase.verifyEqual(N, 25 * (1 - cosd(60)) / 2, AbsTol = 1e-12);
            testCase.verifyFalse(isfolder(fullfile(testCase.Folder, 'GLMALPHA')));
        end

        function continuationUsesItsOwnCache(testCase)
            baseline = glmalpha_new(60, 4, SaveData = false, BeQuiet = true);
            continued = glmalpha_new(60, 4, 1, 0, 0.01, BeQuiet = true);
            cached = glmalpha_new(60, 4, 1, 0, 0.01, SaveData = false, BeQuiet = true);
            testCase.verifyNotEqual(continued, baseline);
            testCase.verifyEqual(cached, continued);
            testCase.verifyNumElements(dir(fullfile(testCase.Folder, 'GLMALPHAPTO', '*.mat')), 1);
        end

        function fractionalRadiiHaveDistinctCaches(testCase)
            [~, ~, ~, ~, first] = glmalpha_new(20.1, 4, BeQuiet = true);
            [~, ~, ~, ~, second] = glmalpha_new(20.2, 4, BeQuiet = true);
            testCase.verifyNotEqual(first, second);
            folder = fullfile(testCase.Folder, 'GLMALPHA');
            testCase.verifyTrue(isfile(fullfile(folder, 'glmalpha-20p1-4-1-0.mat')));
            testCase.verifyTrue(isfile(fullfile(folder, 'glmalpha-20p2-4-1-0.mat')));
        end

        function integerRadiusReusesExistingCacheName(testCase)
            folder = fullfile(testCase.Folder, 'GLMALPHA');
            mkdir(folder);
            G = eye(25);
            V = zeros(1, 25);
            [EM, EL] = addmout(4);
            N = 0;
            save(fullfile(folder, 'glmalpha-20-4-1-0.mat'), 'G', 'V', 'EL', 'EM', 'N');
            [cachedG, cachedV] = glmalpha_new(20, 4, SaveData = false, BeQuiet = true);
            testCase.verifyEqual(cachedG, G);
            testCase.verifyEqual(cachedV, V);
        end

        function numbersRoundTripThroughFilenameTokens(testCase)
            testCase.verifyEqual(num2filename(121.48), '121p48');
            testCase.verifyEqual(num2filename(-30.72), 'n30p72');
            testCase.verifyEqual(num2filename(-0), '0');
            values = [pi, -1e-20, eps, realmin, realmax, 1.000000000000001, 20, 180];
            values = [values, 0.1 - eps(0.1), 1 - eps(1), 10 - eps(10)];

            for value = values
                token = num2filename(value);
                testCase.verifyEqual(str2double(strrep(strrep(token, 'p', '.'), 'n', '-')), value);
                testCase.verifyEmpty(regexp(token, '[.+-]', 'once'));
            end

            testCase.verifyNotEqual(num2filename(121.480001), num2filename(121.480002));
        end

        function bandpassWithZeroMinimumMatchesLowpass(testCase)
            [G, V, ~, ~, N] = glmalpha_new(60, [0, 4], SaveData = false, BeQuiet = true);
            [lowpass, lowpassV, ~, ~, lowpassN] = glmalpha_new(60, 4, SaveData = false, BeQuiet = true);
            testCase.verifyEqual(G.' * G, eye(25), AbsTol = 1e-12);
            testCase.verifyEqual(sort(V), sort(lowpassV), AbsTol = 1e-12);
            testCase.verifyEqual(N, lowpassN, AbsTol = 1e-12);
            testCase.verifyEqual(G * diag(V) * G.', ...
                lowpass * diag(lowpassV) * lowpass.', AbsTol = 1e-12);
            testCase.verifyEmpty(dir(fullfile(testCase.Folder, '**', '*.mat')));
        end

        function positiveMinimumMatchesIndependentKernel(testCase)
            minDegree = 2;
            maxDegree = 4;
            radius = 60;
            [G, V, ~, ~, N] = glmalpha_new(radius, [minDegree, maxDegree], BeQuiet = true);
            firstRetainedRow = minDegree ^ 2 + 1;
            nModes = (maxDegree + 1) ^ 2 - minDegree ^ 2;
            testCase.verifySize(G, [25, nModes]);
            testCase.verifyEqual(G(1:firstRetainedRow - 1, :), zeros(minDegree ^ 2, nModes));
            testCase.verifyEqual(G.' * G, eye(nModes), AbsTol = 1e-12);
            testCase.verifyEqual(N, nModes * (1 - cosd(radius)) / 2, AbsTol = 1e-12);

            [weights, nodes] = gausslegendrecof(20);
            z = (1 + cosd(radius)) / 2 + (1 - cosd(radius)) / 2 * nodes(:);
            weights = weights(:) * (1 - cosd(radius)) / 2;
            longitude = (0:19).' * 2 * pi / 20;
            [theta, phi] = ndgrid(acos(z), longitude);
            orders = addmout(maxDegree);
            harmonics = ylm([0, maxDegree], [], theta(:), phi(:), [], [], 0, 1);
            harmonics = harmonics .* ((-1) .^ orders(:));
            harmonics = harmonics(firstRetainedRow:end, :);
            kernel = (harmonics .* repmat(weights, 20, 1).') * harmonics.' * 2 * pi / 20;
            retained = G(firstRetainedRow:end, :);
            testCase.verifyEqual(kernel * retained, retained .* V, AbsTol = 1e-11);

            [cachedG, cachedV] = glmalpha_new(radius, [minDegree, maxDegree], SaveData = false, BeQuiet = true);
            testCase.verifyEqual(cachedG, G);
            testCase.verifyEqual(cachedV, V);
            testCase.verifyNumElements(dir(fullfile(testCase.Folder, '**', '*.mat')), 1);
            testCase.verifyFalse(isfolder(fullfile(testCase.Folder, 'SDWCAP')));
        end

    end

end
