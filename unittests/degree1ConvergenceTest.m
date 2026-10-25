%% degree1ConvergenceTest - Exercises the actual outer loop with isolated readers and a small kernel.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

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

        function relativeStoppingPreservesSmallEpochs(testCase)
            global DEGREE1_CONVERGENCE_FIXTURE
            [a, sa] = runSolver();
            scales = [1e-9, 1, 1e9];
            DEGREE1_CONVERGENCE_FIXTURE.x = DEGREE1_CONVERGENCE_FIXTURE.x .* scales;
            DEGREE1_CONVERGENCE_FIXTURE.sigma = DEGREE1_CONVERGENCE_FIXTURE.sigma .* scales;
            [b, sb] = runSolver();
            testCase.verifyEqual(b ./ scales', a, 'AbsTol', 1e-12);
            testCase.verifyEqual(sb ./ scales', sa, 'AbsTol', 1e-12);
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

        function toleranceSettingsHaveSeparateCaches(testCase)
            a = runSolver("OnlyId", true); b = runSolver("OnlyId", true, "Tolerance", 1e-8);
            c = runSolver("OnlyId", true, "SleTolerance", 1e-10);
            d = runSolver("OnlyId", true, "MaxIter", 50);
            testCase.verifyEqual(numel(unique({a, b, c, d})), 4);
            testCase.verifyTrue(contains(a, 'zero-relative-v2'));
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
    [varargout{1:nargout}] = solvedegree1('CSR', 'RL06', 2, 2, 2, 'ICE6GD', GeoDomain('alloceans'), ...
        "ReplaceWithGAD", false, "Method", 'uniform', "BeQuiet", true, "SaveData", false, varargin{:});
end

function clearFixture()
    clear global DEGREE1_CONVERGENCE_FIXTURE
end
