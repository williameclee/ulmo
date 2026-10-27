%% degree1ConvergenceTest - Exercises the actual outer loop with isolated readers and a small kernel.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/10/26, En-Chi Lee (williameclee@gmail.com)

classdef degree1ConvergenceTest < matlab.unittest.TestCase

    properties
        Folder
    end

    methods (TestMethodSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            temp = testCase.applyFixture(TemporaryFolderFixture); testCase.Folder = temp.Folder;
            fixtures = fullfile(repo, 'unittests', 'fixtures', 'degree1convergence');
            files = dir(fullfile(fixtures, '*.txt'));

            for k = 1:numel(files)
                copyfile(fullfile(fixtures, files(k).name), fullfile(temp.Folder, erase(files(k).name, '.txt')));
            end

            testCase.applyFixture(PathFixture(temp.Folder));
            old = getenv('GRACEDATA'); setenv('GRACEDATA', temp.Folder);
            testCase.addTeardown(@() setenv('GRACEDATA', old));
            testCase.addTeardown(@() clearFixture());
            previousData = getenv('IFILES'); setenv('IFILES', temp.Folder);
            testCase.addTeardown(@() setenv('IFILES', previousData));
            mkdir(fullfile(temp.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; .021; - .3];
            loadingverticaldisplacement = [0; - .269; - .6];
            tidalgravitationalpotential = [0; 0; .3];
            tidalverticaldisplacement = [0; 0; .6];
            save(fullfile(temp.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential', 'loadingverticaldisplacement', ...
                'tidalgravitationalpotential', 'tidalverticaldisplacement');
            global DEGREE1_CONVERGENCE_FIXTURE
            DEGREE1_CONVERGENCE_FIXTURE = struct('x', reshape(sin(1:27), 9, 3), 'sigma', ones(9, 3) * .1);
            testCase.assertEqual(which('solvedegree1'), fullfile(repo, 'solvedegree1.m'));
        end

    end

    methods (Test)

        function firstUpdateIgnoresImportedUnknowns(testCase)
            global DEGREE1_CONVERGENCE_FIXTURE
            [a, sa] = runSolver("Tolerance", 0, "MaxIter", 1);
            DEGREE1_CONVERGENCE_FIXTURE.x(2:4, :) = 1e8;
            DEGREE1_CONVERGENCE_FIXTURE.sigma(2:4, :) = 1e7;
            [b, sb] = runSolver("Tolerance", 0, "MaxIter", 1);
            testCase.verifyEqual(a, b, 'AbsTol', 1e-14); testCase.verifyEqual(sa, sb, 'AbsTol', 1e-14);
            testCase.verifyGreaterThan(max(abs(a), [], 'all'), 0);
        end

        function absoluteToleranceAllowsSmallUpdates(testCase)
            global DEGREE1_CONVERGENCE_FIXTURE
            DEGREE1_CONVERGENCE_FIXTURE.x = DEGREE1_CONVERGENCE_FIXTURE.x * 1e-12;
            DEGREE1_CONVERGENCE_FIXTURE.sigma = DEGREE1_CONVERGENCE_FIXTURE.sigma * 1e-12;
            [a, sa] = runSolver("Tolerance", 1e-6, "MaxIter", 1);
            [b, sb] = runSolver("Tolerance", 0, "MaxIter", 1);
            testCase.verifyEqual(a, b); testCase.verifyEqual(sa, sb);
        end

        function adaptiveMatchesLongRun(testCase)
            [a, sa] = runSolver("Tolerance", 1e-10);
            [b, sb] = runSolver("Tolerance", 0, "MaxIter", 100);
            testCase.verifyEqual(a, b, 'AbsTol', 1e-9); testCase.verifyEqual(sa, sb, 'AbsTol', 1e-9);
        end

        function nestedFingerprintMatchesLongRun(testCase)
            [a, sa] = runSolver("Method", 'fingerprint', "Tolerance", 1e-8, "SleTolerance", 1e-12);
            [b, sb] = runSolver("Method", 'fingerprint', "Tolerance", 0, "MaxIter", 50, ...
                "SleTolerance", 0, "SleMaxIter", 80);
            testCase.verifyEqual(a, b, 'AbsTol', 1e-7); testCase.verifyEqual(sa, sb, 'AbsTol', 1e-7);
        end

        function jointFingerprintMatchesLongRun(testCase)
            % Exercise merged joint-coefficient support through both solver loops.
            opts = {"IncludeC21S21", true, "Method", 'fingerprint'};
            [a, sa] = runSolver(opts{:});
            [b, sb] = runSolver(opts{:}, "Tolerance", 0, "MaxIter", 100, ...
                "SleTolerance", 0, "SleMaxIter", 100);
            testCase.verifySize(a, [3, 5]);
            testCase.verifyEqual(a, b, 'AbsTol', 1e-5);
            testCase.verifyEqual(sa, sb, 'AbsTol', 1e-5);
        end

        function jointFirstUpdateIgnoresImportedUnknowns(testCase)
            global DEGREE1_CONVERGENCE_FIXTURE
            opts = {"IncludeC21S21", true, "Tolerance", 0, "MaxIter", 1};
            [a, sa] = runSolver(opts{:});
            DEGREE1_CONVERGENCE_FIXTURE.x([2:4, 6:7], :) = 1e8;
            DEGREE1_CONVERGENCE_FIXTURE.sigma([2:4, 6:7], :) = 1e7;
            [b, sb] = runSolver(opts{:});
            testCase.verifyEqual(a, b, 'AbsTol', 1e-14);
            testCase.verifyEqual(sa, sb, 'AbsTol', 1e-14);
        end

        function innerCapFailureDoesNotSave(testCase)
            testCase.verifyError(@() ...
                runSolver("Method", 'fingerprint', "SleMaxIter", 1, ...
                "SleTolerance", 1e-15, "SaveData", true), 'ULMO:solvesle:NonConvergence');
            testCase.verifyEmpty(dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat')));
        end

        function capFailureDoesNotSave(testCase)
            testCase.verifyError(@() ...
                runSolver("Tolerance", 1e-15, "MaxIter", 1, "SaveData", true), ...
            'ULMO:solvedegree1:NonConvergence');
            testCase.verifyEmpty(dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat')));
        end

        function iterationSettingsShareCache(testCase)
            % Cached solutions are assumed converged at sufficiently small tolerances.
            settings = {{"Tolerance", 1e-8}, {"SleTolerance", 1e-10}, ...
                            {"MaxIter", 50}, {"SleMaxIter", 80}};
            cachePattern = fullfile(testCase.Folder, 'Degree1', 'new', '*.mat');

            for joint = [false, true]
                opts = {"IncludeC21S21", joint, "SaveData", true};
                runSolver(opts{:});
                files = dir(cachePattern);
                % Degree-1-only and joint solutions must have distinct caches.
                testCase.assertNumElements(files, 1 + double(joint));

                for k = 1:numel(settings)
                    runSolver(opts{:}, "ForceNew", true, settings{k}{:});
                    actual = dir(cachePattern);
                    testCase.verifyEqual({actual.name}, {files.name});
                end

            end

        end

        function forceNewAppliesChangedIterationSettings(testCase)
            % Unusual settings take effect only when explicitly bypassing the cache.
            for joint = [false, true]
                opts = {"IncludeC21S21", joint};
                [coeffs, coeffStds, dates] = runSolver(opts{:}, "Tolerance", 1e-10, "SaveData", true);
                limited = {"Tolerance", 1e-15, "MaxIter", 1};
                [actual, actualStd, actualDates] = runSolver(opts{:}, limited{:});
                testCase.verifyEqual(actual, coeffs);
                testCase.verifyEqual(actualStd, coeffStds);
                testCase.verifyEqual(actualDates, dates);
                testCase.verifyError(@() runSolver(opts{:}, limited{:}, ...
                    "ForceNew", true, "SaveData", true), 'ULMO:solvedegree1:NonConvergence');
                % A failed fresh run must leave the converged cache intact.
                [actual, actualStd] = runSolver(opts{:});
                testCase.verifyEqual(actual, coeffs);
                testCase.verifyEqual(actualStd, coeffStds);
            end

        end

        function ignoresLegacyGadCache(testCase)
            [expected, expectedStd, dates] = runSolver("SaveData", true);
            files = dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat'));
            testCase.assertNumElements(files, 1);
            testCase.verifyFalse(contains(files.name, '_GAD-'));
            cachePath = fullfile(files.folder, files.name);
            legacyPath = fullfile(files.folder, strrep(files.name, 'CSRRL062_', 'CSRRL062_GAD-'));
            movefile(cachePath, legacyPath);
            coeffs = ones(size(expected)) * 1e6;
            coeffStds = coeffs;
            save(legacyPath, 'coeffs', 'coeffStds', 'dates');
            [actual, actualStd] = runSolver("SaveData", true);
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(actualStd, expectedStd);
            testCase.verifyTrue(isfile(legacyPath));
            testCase.verifyTrue(isfile(cachePath));
        end

        function reusesGsmCache(testCase)
            [coeffs, coeffStds, dates] = runSolver("SaveData", true);
            files = dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat'));
            testCase.assertNumElements(files, 1);
            coeffs(:) = 42;
            cachePath = fullfile(files.folder, files.name);
            save(cachePath, 'coeffs', 'coeffStds', 'dates');
            testCase.verifyEqual(runSolver(), coeffs);
        end

        function checksLastEpoch(testCase)
            global DEGREE1_CONVERGENCE_FIXTURE
            DEGREE1_CONVERGENCE_FIXTURE.x(:, 1:2) = 0;
            DEGREE1_CONVERGENCE_FIXTURE.sigma(:) = 0;
            testCase.verifyError(@() runSolver("Tolerance", 1e-15, "MaxIter", 1), ...
            'ULMO:solvedegree1:NonConvergence');
        end

    end

end

function varargout = runSolver(varargin)
    [varargout{1:nargout}] = legacysolvedegree1('CSR', 'RL06', 2, 2, 2, 'ICE6GD', GeoDomain('alloceans'), ...
        "Method", 'uniform', "BeQuiet", true, "SaveData", false, varargin{:});
end

function clearFixture()
    clear global DEGREE1_CONVERGENCE_FIXTURE
end

% Keep the historical iteration contract explicit in regression tests.
function varargout = legacysolvedegree1(varargin)
    [varargout{1:nargout}] = solvedegree1(varargin{:}, "Solver", "iteration");
end
