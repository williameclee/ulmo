%% GRACEEPOCHSTEST - Checks product discovery, gap flags and optional reference data.
% Small GRCOF2 fixtures exercise actual parsing without production data.
% Existing one/two-output calls must not require reference header metadata.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

classdef graceepochsTest < matlab.unittest.TestCase
    properties
        Folder
        Windows = [datetime(2003, 1, 1), datetime(2003, 2, 1); ...
            datetime(2003, 3, 1), datetime(2003, 4, 1)]
    end
    properties (TestParameter)
        product = {{'CSR', 'RL06', 60}, {'GFZ', 'RL06', 96}, {'JPL', 'RL06', 60}, ...
            {'CSR', 'RL05', 60}, {'GFZ', 'RL05', 90}, {'JPL', 'RL05', 90}, ...
            {'CSR mascon', 'RL06', 720}}
    end
    methods (TestMethodSetup)
        function setup(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root, IncludingSubfolders=true));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', testCase.Folder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GRACEDATA', fullfile(testCase.Folder, 'GRACE')));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'ORIGINALGRACEDATA', fullfile(testCase.Folder, 'source')));
        end
    end
    methods (Test)
        function discoversProductsAndPreservesActualWindows(testCase, product)
            testCase.writeFiles(product, testCase.Windows, [6378136.3; 6378136.3], ...
                getenv('ORIGINALGRACEDATA'), true);
            [epochs, synthetic, radius] = graceepochs(product);
            [expected, expectedSynthetic] = fillgraceepochs(testCase.Windows);
            testCase.verifyEqual(epochs, expected);
            testCase.verifyEqual(synthetic, expectedSynthetic);
            testCase.verifyEqual(epochs(~synthetic, :), testCase.Windows);
            testCase.verifyEqual(sum(synthetic), 1);
            testCase.verifyEqual(radius, 6378136.3);
        end
        function oldOutputFormsDoNotReadReferenceHeaders(testCase)
            spec = {'CSR', 'RL06', 60};
            testCase.writeFiles(spec, testCase.Windows, [0; 0], getenv('ORIGINALGRACEDATA'), false);
            expected = fillgraceepochs(testCase.Windows);
            testCase.verifyEqual(graceepochs(spec), expected);
            [epochs, synthetic] = graceepochs(spec);
            testCase.verifyEqual(epochs, expected);
            testCase.verifyEqual(epochs(~synthetic, :), testCase.Windows);
            testCase.verifyError(@() readreference(spec), 'ULMO:parsegracesourcefile:InvalidHeader');
        end
        function radiusCheckIsOptionalAndIncludesTheFullRecord(testCase)
            spec = {'CSR', 'RL06', 60};
            testCase.writeFiles(spec, testCase.Windows, [6378136.3; 6378137.3], ...
                getenv('ORIGINALGRACEDATA'), true);
            expected = fillgraceepochs(testCase.Windows);
            [epochs, ~] = graceepochs(spec);
            testCase.verifyEqual(epochs, expected);
            testCase.verifyError(@() readreference(spec), 'ULMO:Epochs:Reference');
            testCase.verifyError(@() readreference(spec, ...
                [datetime(2003, 1, 1), datetime(2003, 2, 1)]), 'ULMO:Epochs:Reference');
        end
        function rangeKeepsSyntheticWindowsAndRadius(testCase)
            spec = {'CSR', 'RL06', 60};
            testCase.writeFiles(spec, testCase.Windows, [6378136.3; 6378136.3], ...
                getenv('ORIGINALGRACEDATA'), true);
            range = [datetime(2003, 2, 1), datetime(2003, 2, 28)];
            [epochs, synthetic, radius] = graceepochs(spec, datenum(range)); %#ok<DATNM>
            testCase.verifyEqual(epochs, [datetime(2003, 2, 1), datetime(2003, 3, 1)]);
            testCase.verifyTrue(synthetic);
            testCase.verifyEqual(radius, 6378136.3);
        end
        function aliasesAndFixedProductDegreesFollowTheLoader(testCase)
            testCase.writeFiles({'CSR', 'RL05', 60}, testCase.Windows, [1; 1], ...
                getenv('ORIGINALGRACEDATA'), true);
            testCase.verifyEqual(graceepochs({'csr', 5, 96}), fillgraceepochs(testCase.Windows));
            testCase.writeFiles({'CSR mascon', 'RL06', 720}, testCase.Windows, [1; 1], ...
                getenv('ORIGINALGRACEDATA'), true);
            testCase.verifyEqual(graceepochs({'csr mascon', 'rl06', 60}), fillgraceepochs(testCase.Windows));
        end
        function rootsUseTheExistingPriority(testCase)
            spec = {'CSR', 'RL06', 60};
            roots = {getenv('ORIGINALGRACEDATA'), fullfile(getenv('GRACEDATA'), 'raw'), ...
                fullfile(getenv('IFILES'), 'GRACE', 'raw')};
            % Give GRACEDATA a separate root to distinguish both fallbacks.
            setenv('GRACEDATA', fullfile(testCase.Folder, 'alternate'));
            roots{2} = fullfile(getenv('GRACEDATA'), 'raw');
            for iRoot = 1:3
                testCase.writeFiles(spec, testCase.Windows, [iRoot; iRoot], roots{iRoot}, true);
            end
            [~, ~, radius] = graceepochs(spec);
            testCase.verifyEqual(radius, 1);
            setenv('ORIGINALGRACEDATA', '');
            [~, ~, radius] = graceepochs(spec);
            testCase.verifyEqual(radius, 2);
            setenv('GRACEDATA', '');
            [~, ~, radius] = graceepochs(spec);
            testCase.verifyEqual(radius, 3);
        end
        function missingAndUnsupportedProductsFailClearly(testCase)
            testCase.verifyError(@() graceepochs({'CSR', 'RL06', 60}), 'ULMO:Epochs:MissingReference');
            for spec = {{'CSR', 'RL04', 60}, {'CSR', 'RL06', 90}, {'CSR mascon', 'RL05', 720}}
                testCase.verifyError(@() graceepochs(spec{1}), 'ULMO:Epochs:UnsupportedProduct');
            end
        end
        function multipleSolutionsInOneMonthArePreserved(testCase)
            spec = {'CSR', 'RL06', 60};
            windows = [datetime(2003, 1, 1), datetime(2003, 1, 20); ...
                datetime(2003, 1, 15), datetime(2003, 2, 1); testCase.Windows(2, :)];
            testCase.writeFiles(spec, windows, [1; 1; 1], getenv('ORIGINALGRACEDATA'), true);
            [epochs, synthetic] = graceepochs(spec);
            testCase.verifyEqual(epochs(~synthetic, :), windows);
        end
    end
    methods
        function writeFiles(~, spec, windows, radii, root, header)
            folder = fullfile(root, spec{2}, spec{1});
            mkdir(folder);
            for iDate = 1:size(windows, 1)
                dateTag = sprintf('2003%03d-2003%03d', iDate * 30, iDate * 30 + 29);
                if strcmp(spec{1}, 'CSR mascon')
                    name = ['GSU-2_', dateTag, '_B---_0603'];
                elseif strcmp(spec{2}, 'RL06')
                    if spec{3} == 60, bandwidth = 'BA01'; else, bandwidth = 'BB01'; end
                    name = ['GSM-2_', dateTag, '_', bandwidth, '_0600'];
                else
                    switch spec{1}
                        case 'CSR', suffix = 'UTCSR_0060_0005';
                        case 'GFZ', suffix = 'EIGEN_G---_005a';
                        case 'JPL', suffix = 'JPLEM_0000_0005';
                    end
                    name = ['GSM-2_', dateTag, '_TEST_', suffix];
                end
                fid = fopen(fullfile(folder, name), 'w');
                cleanup = onCleanup(@() fclose(fid));
                if header, fprintf(fid, 'EARTH 3.986004415e14 %.9f\n', radii(iDate)); end
                fprintf(fid, 'GRCOF2 0 0 0 0 0 0 %s %s 0\n', ...
                    string(windows(iDate, 1), 'yyyyMMdd.HHmm'), string(windows(iDate, 2), 'yyyyMMdd.HHmm'));
                clear cleanup
            end
        end
    end
end

function radius = readreference(spec, range)
    arguments (Input)
        spec
        range = []
    end
    [~, ~, radius] = graceepochs(spec, range);
end
