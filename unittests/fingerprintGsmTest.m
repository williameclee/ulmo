%% fingerprintGsmTest - Checks GIA-corrected GSM forcing.
% Uses synthetic readers and an SLE boundary observer; it does not test SLE physics.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/10/26, En-Chi Lee (williameclee@gmail.com)

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
            % Include nonzero low-degree GIA to catch accidental correction of degree 1.
            DEGREE1_CONVERGENCE_FIXTURE.gia = zeros(3, 6, 4);
            DEGREE1_CONVERGENCE_FIXTURE.gia(:, :, 3:4) = reshape(1:36, 3, 6, 2) / 100;
            testCase.assertEqual(which('grace2fingerprint'), fullfile(repo, 'grace2fingerprint.m'));
        end

    end

    methods (Test)

        function appliesGiaToHigherDegreeForcing(testCase, regional)
            global DEGREE1_CONVERGENCE_FIXTURE
            args = {};

            if regional
                args = {"LandDomain", GeoDomain('alloceans')};
            end

            [actual, actualStd] = runFingerprint(args{:}, "SaveData", false);
            [p, s, dates] = grace2plmt_new();
            testCase.verifyEqual(DEGREE1_CONVERGENCE_FIXTURE.giaModel, 'synthetic');
            testCase.verifyEqual(DEGREE1_CONVERGENCE_FIXTURE.giaDates, dates);
            % Only degrees >= 2 are corrected; formal GRACE errors are unchanged.
            p(:, 4:end, 3:4) = p(:, 4:end, 3:4) - ...
                DEGREE1_CONVERGENCE_FIXTURE.gia(:, 4:end, 3:4);
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
        GeoDomain('alloceans'), "GIA", 'synthetic', "RotationFeedback", false, ...
        "OutputFormat", 'traditional', "TimeFormat", 'datetime', ...
        "BeQuiet", true, varargin{:});
end

function clearFixture()
    clear global DEGREE1_CONVERGENCE_FIXTURE
end
