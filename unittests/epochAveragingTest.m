%% epochAveragingTest - Unit tests for epoch averaging, temporal interpolation, and GRACE epoch utilities.
%
% Last modified by
%	2026/09/26, En-Chi Lee (williameclee@arizona.edu)

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
