%% epochAveragingTest - Unit tests for epoch averaging, temporal interpolation, and GRACE epoch utilities.
%
% Last modified
%	2026/11/01, En-Chi Lee (williameclee@arizona.edu)

classdef epochAveragingTest < matlab.unittest.TestCase

    methods (Test)

        function testAverageEpochsDatetime(testCase)
            % Synthetic data checks independent of external product files using datetime.
            dates = datetime(2020, 1, 1) + days(0:59)';
            mesh = reshape(1:60, 1, 1, []);
            epochs = [datetime(2020, 1, 1), datetime(2020, 2, 1); ...
                          datetime(2020, 2, 1), datetime(2020, 3, 1); ...
                          datetime(2020, 4, 1), datetime(2020, 5, 1)];
            [vals, labels] = averageepochs(dates, mesh, epochs);

            testCase.verifyEqual(squeeze(vals(1, 1, 1:2)), [16; 46]);
            testCase.verifyTrue(isnan(vals(1, 1, 3)));
            testCase.verifyEqual(labels, mean(epochs, 2));
        end

        function testAverageEpochsDatenum(testCase)
            % Synthetic data checks using numeric datenum and handling NaNs.
            dates = datetime(2020, 1, 1) + days(0:59)';
            mesh = reshape(1:60, 1, 1, []);
            epochs = [datetime(2020, 1, 1), datetime(2020, 2, 1); ...
                          datetime(2020, 2, 1), datetime(2020, 3, 1); ...
                          datetime(2020, 4, 1), datetime(2020, 5, 1)];
            mesh(1, 1, 1) = NaN;
            vals = averageepochs(datenum(dates), mesh, datenum(epochs)); %#ok<DATNM>

            testCase.verifyEqual(vals(1, 1, 1), 16.5);
        end

        function testMonthlyMeanWindows(testCase)
            % Selecting a midpoint keeps the entire leap-year February window.
            dates = datetime(2020, 1, 1) + days(0:90);
            limits = datetime(2020, 2, [14, 16]);
            epochs = temporalmeanepochs(dates, 'midmonth', limits);
            testCase.verifyEqual(epochs, [datetime(2020, 2, 1), datetime(2020, 3, 1)]);
        end

        function testFixedStepMeanWindows(testCase)
            % Boundary samples belong to the window on their right.
            dates = datetime(2020, 1, 1) + days(0:8);
            epochs = temporalmeanepochs(dates, days(4), []);
            [values, labels] = averageepochs(dates, reshape(0:8, 1, 1, []), epochs);
            testCase.verifyEqual(labels, dates([1, 5, 9])');
            testCase.verifyEqual(values(:), [0.5; 3.5; 7]);
        end

        function testSSHGraceCache(testCase)
            % A subset request caches the full spatial grid and all epoch means.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', folder.Folder));
            sshFolder = fullfile(folder.Folder, 'SSH');
            refFolder = fullfile(folder.Folder, 'RL06', 'CSR');
            mkdir(sshFolder);
            mkdir(refFolder);
            refFile = fullfile(refFolder, 'GSM-2_test_BA01_0600');
            writelines(["# End of YAML header"; ...
                        "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], refFile);

            dates = (datetime(2020, 1, 1):datetime(2020, 2, 29))';
            lon = [0; 180];
            lat = [-90; 90];
            sshs = lon + reshape(1:numel(dates), 1, 1, []) + 0 * lat';
            sshErrors = 2 * sshs;
            rawPath = fullfile(sshFolder, 'MEaSUREs.mat');
            save(rawPath, 'dates', 'lon', 'lat', 'sshs', 'sshErrors');

            % A spatial-only cache must be ignored on a miss for the GRACE grid.
            sshs(:) = NaN;
            save(fullfile(sshFolder, 'MEaSUREs-H90-linear.mat'), 'dates', 'lon', 'lat', 'sshs', 'sshErrors');

            limits = datetime(2020, 1, [16, 17]);
            [subset, ~, subsetDates] = ssh2lonlatt('MEaSUREs', 'GRACE', 90, limits, BeQuiet = true);
            testCase.verifyEqual(subsetDates, datetime(2020, 1, 16, 12, 0, 0));
            testCase.verifyEqual(subset(2, 2), single(90 + 16));

            testCase.assertTrue(isfile(fullfile(sshFolder, 'MEaSUREs-Tgrace-H90-linear.mat')));

            % A cache hit must work without the raw grid or reference epochs.
            delete(rawPath);
            delete(refFile);
            [ssh, sigma, labels, x, y] = ssh2lonlatt('MEaSUREs', "GRACE", 90, BeQuiet = true);
            testCase.verifyEqual(labels, [subsetDates; datetime(2020, 2, 15, 12, 0, 0)]);
            testCase.verifyEqual(x, 0:90:360);
            testCase.verifyEqual(y, (-90:90:90)');
            testCase.verifyEqual(ssh(2, 2, 2), single(90 + mean(32:60)));
            testCase.verifyEqual(sigma(2, 2, 2), single(sqrt(mean((2 * (90 + (32:60))) .^ 2))), ...
                'RelTol', single(1e-6));
        end

        function testStericGraceCache(testCase)
            % Derived components use paired native samples before epoch averaging.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', folder.Folder));
            stericFolder = fullfile(folder.Folder, 'STERIC', 'EN4c14');
            refFolder = fullfile(folder.Folder, 'RL06', 'CSR');
            mkdir(stericFolder);
            mkdir(refFolder);
            refFile = fullfile(refFolder, 'GSM-2_test_BA01_0600');
            writelines(["# End of YAML header"; ...
                        "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], refFile);

            dates = datetime(2020, [1; 1; 2; 2], [1; 20; 1; 29]);
            lon = [0, 180];
            lat = [-90; 0; 90];
            stericSl = zeros(size(lat)) + lon + reshape([10, 100, 20, 40], 1, 1, []);
            halostericSl = repmat(reshape([3, NaN, 4, 6], 1, 1, []), 3, 2);
            rawPath = fullfile(stericFolder, 'EN4c14-StericSeaLevel.mat');
            save(rawPath, 'dates', 'lon', 'lat', 'stericSl', 'halostericSl');

            % A spatial-only cache with different values must be ignored.
            stericSl(:) = NaN;
            save(fullfile(stericFolder, 'EN4c14-StericSeaLevel-H90-linear.mat'), ...
                'dates', 'lon', 'lat', 'stericSl', 'halostericSl');
            limits = datetime(2020, 1, [16, 17]);
            args = {'Type', 'non-halosteric', 'BeQuiet', true};
            subset = steric2lonlatt('EN4c14', 'GRACE', 90, limits, args{:});
            testCase.verifyEqual(subset(2, 2), 97);
            total = steric2lonlatt('EN4c14', 'GRACE', 90, limits, Type = 'total', BeQuiet = true);
            testCase.verifyEqual(total(2, 2), 145);

            % The wider request must succeed using only the finished component cache.
            delete(rawPath);
            delete(refFile);
            [steric, sigma, labels, x, y] = steric2lonlatt('EN4c14', "GRACE", 90, [], args{:});
            testCase.verifyEqual(labels, datetime(2020, [1; 2], [16; 15], 12, 0, 0));
            testCase.verifyEqual(squeeze(steric(2, 2, :)), [97; 115]);
            testCase.verifyEqual(x, 0:90:360);
            testCase.verifyEqual(y, (-90:90:90)');
            % Unavailable steric uncertainty is represented by a zero grid.
            testCase.verifyEqual(sigma, zeros(size(steric), 'like', steric));
        end

        function testStericGraceTruncation(testCase)
            % Truncation selects the same midpoints and uses the same native samples as a full run.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', folder.Folder));
            stericFolder = fullfile(folder.Folder, 'STERIC', 'EN4c14');
            refFolder = fullfile(folder.Folder, 'RL06', 'CSR');
            mkdir(stericFolder);
            mkdir(refFolder);
            writelines(["# End of YAML header"; ...
                        "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], ...
                fullfile(refFolder, 'GSM-2_test_BA01_0600'));
            dates = [datetime(2020, 1, 16, 12, 0, 0); datetime(2020, 2, 1); datetime(2020, 2, 15, 12, 0, 0)];
            lon = [0, 180];
            lat = [-90; 0; 90];
            stericSl = repmat(reshape([4, 8, 12], 1, 1, []), 3, 2);
            save(fullfile(stericFolder, 'EN4c14-StericSeaLevel.mat'), 'dates', 'lon', 'lat', 'stericSl');

            args = {'Type', 'total', 'BeQuiet', true};
            [fullGrid, ~, labels] = steric2lonlatt('EN4c14', 'GRACE', [], [], args{:});
            testCase.verifyEqual(labels, dates([1, 3]));
            testCase.verifyEqual(squeeze(fullGrid(1, 1, :)), [4; 10]);
            ranges = [datetime(2022, 1, 1), datetime(2022, 3, 1); ...
                      datetime(2020, 2, 10), datetime(2020, 4, 1); ...
                      datetime(2019, 12, 1), dates(1)];
            for iRange = 1:size(ranges, 1)
                limits = ranges(iRange, :);
                keep = labels >= limits(1) & labels <= limits(2);
                [grid, ~, selected] = steric2lonlatt('EN4c14', 'GRACE', [], limits, args{:}, ...
                    TruncationBeforeInterpolation = true);
                testCase.verifyEqual(selected, labels(keep));
                testCase.verifyEqual(grid, fullGrid(:, :, keep));
            end
        end

        function testInterpTemporal(testCase)
            dates = datetime(2020, 1, 1) + days(0:99)';
            mesh = cat(1, reshape(1:100, 1, 1, []), reshape(101:200, 1, 1, []));
            [vals, labels] = interptemporal(dates, mesh, days(10), 'linear', BeQuiet = true);
            [~, group] = min(abs(dates(:)' - labels(:)), [], 1);

            for k = 1:numel(labels)
                testCase.verifyEqual(vals(:, :, k), mean(mesh(:, :, group == k), 3));
            end

        end

        function testFillGraceEpochs(testCase)
            refs = [datetime(2017, 6, 3), datetime(2017, 6, 28, 18, 0, 0); ...
                        datetime(2018, 6, 5, 6, 0, 0), datetime(2018, 6, 29)];
            [filled, synthetic] = fillgraceepochs(refs, [datetime(2017, 6, 1), datetime(2018, 7, 1)]);
            gap = filled(synthetic, :);

            testCase.verifyEqual(gap(1, 1), refs(1, 2));
            testCase.verifyEqual(gap(end, 2), refs(2, 1));
            testCase.verifyEqual(gap(1:end - 1, 2), gap(2:end, 1));
            testCase.verifyEqual(filled(~synthetic, :), refs);
        end

        function testFillGraceEpochsSubset(testCase)
            refs = [datetime(2017, 6, 3), datetime(2017, 6, 28, 18, 0, 0); ...
                        datetime(2018, 6, 5, 6, 0, 0), datetime(2018, 6, 29)];
            [filled, ~] = fillgraceepochs(refs, [datetime(2017, 6, 1), datetime(2018, 7, 1)]);
            [subset, ~] = fillgraceepochs(refs, datetime(2017, [9 11], 16));
            labels = mean(filled, 2);
            keep = labels >= datetime(2017, 9, 16) & labels <= datetime(2017, 11, 16);

            testCase.verifyEqual(subset, filled(keep, :));
        end

        function testMustBeTimeStep(testCase)
            testCase.verifyWarningFree(@() mustBeTimeStep([]));
            testCase.verifyWarningFree(@() mustBeTimeStep('midmonth'));
            testCase.verifyWarningFree(@() mustBeTimeStep('GRACE'));
            testCase.verifyWarningFree(@() mustBeTimeStep(days(10)));

            testCase.verifyError(@() mustBeTimeStep('invalid'), 'ULMO:notValidTimeStep');
            testCase.verifyError(@() mustBeTimeStep(-days(5)), 'ULMO:notValidTimeStep');
            testCase.verifyError(@() mustBeTimeStep([days(1), days(2)]), 'ULMO:notValidTimeStep');
        end

    end

end
