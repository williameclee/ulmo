%% fingerprintGsmTest - Checks GSM forcing.
% Uses synthetic readers and an SLE boundary observer; it does not test SLE physics.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

classdef fingerprintGsmTest < matlab.unittest.TestCase

    properties
        Folder
    end

    properties (TestParameter)
        redoDegree1 = {false, true}
        regional = {false, true}
    end

    methods (TestMethodSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            temp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.Folder = temp.Folder;

            for name = ["degree1convergence", "fingerprintgsm"]
                source = fullfile(repo, 'unittests', 'fixtures', name);
                files = dir(fullfile(source, '*.txt'));

                for k = 1:numel(files)
                    copyfile(fullfile(source, files(k).name), ...
                        fullfile(temp.Folder, erase(files(k).name, '.txt')));
                end

            end

            testCase.applyFixture(PathFixture(temp.Folder));
            previous = getenv('IFILES');
            setenv('IFILES', temp.Folder);
            testCase.addTeardown(@() setenv('IFILES', previous));
            testCase.addTeardown(@() clearFixture());
            global DEGREE1_CONVERGENCE_FIXTURE
            DEGREE1_CONVERGENCE_FIXTURE = struct( ...
                'x', reshape(sin(1:27), 9, 3), 'sigma', ones(9, 3) * .1);
            testCase.assertEqual(which('grace2fingerprint'), fullfile(repo, 'grace2fingerprint.m'));
        end

    end

    methods (Test)

        function preservesGsmForcing(testCase, regional)
            args = {};

            if regional
                args = {"LandDomain", GeoDomain('alloceans')};
            end

            [actual, actualStd] = runFingerprint(args{:}, "SaveData", false);
            [p, s] = grace2plmt_new();
            K = kernelcp_new();

            if ~regional
                K = eye(size(K)) - K;
            end

            expected = localise(permute(p, [2, 3, 1]), "L", 2, "K", K);
            expectedStd = localise(permute(s, [2, 3, 1]), "L", 2, "K", K, "IsError", true);
            testCase.verifyEqual(actual(:, 3:4, :), expected(:, 3:4, :), 'AbsTol', 1e-14);
            testCase.verifyEqual(actualStd(:, 3:4, :), expectedStd(:, 3:4, :), 'AbsTol', 1e-14);
        end

    end

end

function varargout = runFingerprint(varargin)
    [varargout{1:nargout}] = grace2fingerprint({'CSR', 'RL06', 2}, 2, ...
        GeoDomain('alloceans'), "GIA", false, "RotationFeedback", false, ...
        "OutputFormat", 'traditional', "TimeFormat", 'datetime', ...
        "BeQuiet", true, varargin{:});
end

function clearFixture()
    clear global DEGREE1_CONVERGENCE_FIXTURE
end
