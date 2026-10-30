%% EPOCHSTATISTICSTEST - Tests generic epoch means, completeness and uncertainty.
% Uses multiple sampling cadences and independent analytic error propagation.
%
% Created by
%   2026/10/30

classdef epochStatisticsTest < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function suppressExpectedWarnings(testCase)
            state = warning('off', 'ULMO:averageepochs:IncompleteEpoch');
            testCase.addTeardown(@() warning(state));
        end
    end
    methods (Test)
        function calendarMonthsSpanLeapYearAndGaps(testCase)
            dates = [datetime(2004, 3, 1), datetime(2004, 1, 31)];
            epochs = calendarmonthepochs(dates);
            starts = datetime(2004, (1:3)', 1);
            testCase.verifyEqual(epochs, [starts, starts + calmonths(1)]);
            testCase.verifyEqual(calendarmonthepochs(datenum(dates)), epochs);
            testCase.verifyEqual(epochs(2, 2) - epochs(2, 1), days(29));
        end
        function calendarMonthsRejectInvalidDates(testCase)
            testCase.verifyError(@() calendarmonthepochs(NaT), 'ULMO:calendarmonthepochs:InvalidDates');
            testCase.verifyError(@() calendarmonthepochs(Inf), 'ULMO:calendarmonthepochs:InvalidDates');
        end
        function irregularOverlappingWindowsUseSuppliedBounds(testCase)
            dates = datetime(2003, 1, 31, 18, 0, 0) + hours((0:4)' * 6);
            data = reshape(1:5, 1, 1, []);
            epochs = dates(1) + hours([7, 25; 0, 12; 6, 18]);
            [p, d, s] = averageepochs(dates, data, epochs, hours(6), 'CheckCompleteness', true);
            testCase.verifyEqual(p(:), [4; 1.5; 2.5]);
            testCase.verifyEqual(d, mean(epochs, 2));
            testCase.verifyEqual(s(:), [1 / sqrt(3); 0.5; 0.5], 'AbsTol', 1e-14);
        end
        function cadenceAndPhaseAreCallerSupplied(testCase)
            origin = datetime(2020, 1, 1, 0, 15, 0);
            dates = origin + hours((0:7)' * 7);
            epochs = [origin, origin + hours(56)];
            data = reshape(1:8, 1, 1, []);
            [p, ~] = averageepochs(dates, data, epochs, hours(7), ...
                'CheckCompleteness', true, 'SamplingOrigin', origin);
            testCase.verifyEqual(p, 4.5);
            missing = averageepochs(dates(2:end), data(:, :, 2:end), epochs, hours(7), ...
                'CheckCompleteness', true, 'SamplingOrigin', origin);
            testCase.verifyTrue(isnan(missing));
        end
        function checkingCanBeDisabledForIrregularSamples(testCase)
            origin = datetime(2020, 1, 1);
            dates = origin + hours([0; 1; 4]);
            data = reshape([2, 4, 9], 1, 1, []);
            epochs = [origin, origin + hours(5)];
            unchecked = averageepochs(dates, data, epochs, hours(1), 'CheckCompleteness', false);
            checked = averageepochs(dates, data, epochs, hours(1), 'CheckCompleteness', true);
            testCase.verifyEqual(unchecked, 5);
            testCase.verifyTrue(isnan(checked));
        end
        function independentFormalErrorsPropagateToMean(testCase)
            dates = datetime(2020, 1, 1) + days((0:1)');
            data = reshape([100, -100], 1, 1, []);
            sigma = reshape([3, 4], 1, 1, []);
            epochs = [dates(1), dates(2) + days(1)];
            [p, ~, s] = averageepochs(dates, data, epochs, 1, ...
                'CheckCompleteness', true, 'SampleStd', sigma);
            testCase.verifyEqual(p, 0);
            testCase.verifyEqual(s, 2.5); % sqrt(3^2+4^2)/2; not temporal spread.
        end
        function singleSampleRetainsFormalSigma(testCase)
            date = datetime(2020, 1, 1);
            [p, ~, s] = averageepochs(date, 7, [date, date + days(1)], ...
                'SampleStd', 3);
            testCase.verifyEqual(p, 7);
            testCase.verifyEqual(s, 3);
        end
        function unknownSigmaUsesEmpiricalStandardError(testCase)
            dates = datetime(2020, 1, 1) + days((0:1)');
            [p, ~, s] = averageepochs(dates, reshape([2, 4], 1, 1, []), ...
                [dates(1), dates(2) + days(1)], 'SampleStd', reshape([3, NaN], 1, 1, []));
            testCase.verifyEqual(p, 3);
            testCase.verifyEqual(s, 1);
        end
        function empiricalEstimateUsesSampleVarianceAndValidCount(testCase)
            dates = datetime(2020, 1, 1) + days((0:3)');
            data = reshape([2, NaN, 6, 10], 1, 1, []);
            epoch = [dates(1), dates(end) + days(1)];
            [p, ~, s] = averageepochs(dates, data, epoch);
            testCase.verifyEqual(p, 6);
            testCase.verifyEqual(s, 4 / sqrt(3), 'RelTol', 1e-14);
            [~, ~, s] = averageepochs(dates(1), 2, epoch);
            testCase.verifyTrue(isnan(s));
            [~, ~, s] = averageepochs(dates, ones(size(data)), epoch);
            testCase.verifyEqual(s, 0);
        end
        function missingValuesUseSameMaskForMeanAndSigma(testCase)
            dates = datetime(2020, 1, 1) + days((0:2)');
            data = reshape([2, NaN, 6], 1, 1, []);
            sigma = reshape([3, NaN, 4], 1, 1, []);
            epoch = [dates(1), dates(3) + days(1)];
            [p, ~, s] = averageepochs(dates, data, epoch, 'SampleStd', sigma);
            testCase.verifyEqual(p, 4);
            testCase.verifyEqual(s, 2.5);
            [p, ~, s] = averageepochs(dates, data, epoch, days(1), ...
                'SampleStd', sigma, 'CheckCompleteness', true);
            testCase.verifyTrue(isnan(p) && isnan(s));
        end
        function arbitraryGridShapeAndSinglePrecision(testCase)
            dates = datetime(2020, 1, 1) + hours((0:2)');
            data = single(reshape(1:18, [2, 3, 3]));
            epoch = [dates(1), dates(3) + hours(1)];
            [p, ~, s] = averageepochs(dates, data, epoch, 'SampleStd', ones(size(data), 'single'));
            testCase.verifyClass(p, 'single');
            testCase.verifySize(p, [2, 3]);
            testCase.verifyEqual(p, mean(data, 3));
            testCase.verifyEqual(s, ones(2, 3, 'single') / sqrt(single(3)), 'RelTol', single(1e-6));
        end
        function emptyEpochsAndUncoveredWindows(testCase)
            date = datetime(2020, 1, 1);
            [p, d, s] = averageepochs(date, 1, NaT(0, 2));
            testCase.verifySize(p, [1, 1, 0]);
            testCase.verifySize(s, [1, 1, 0]);
            testCase.verifySize(d, [0, 1]);
            [p, ~, s] = averageepochs(date, 1, [date + days(1), date + days(2)], 'SampleStd', 2);
            testCase.verifyTrue(isnan(p) && isnan(s));
        end
        function numericDatesAndStepAreSupported(testCase)
            dates = datetime(2020, 1, 1) + days((0:2)');
            data = reshape(1:3, 1, 1, []);
            epoch = [dates(1), dates(3) + days(1)];
            [p, d] = averageepochs(datenum(dates), data, datenum(epoch), 1, ...
                'CheckCompleteness', true, 'SamplingOrigin', datenum(dates(1)));
            testCase.verifyEqual(p, 2);
            testCase.verifyEqual(d, mean(epoch, 2));
        end
        function invalidInputsAreRejected(testCase)
            date = datetime(2020, 1, 1);
            epoch = [date, date + days(1)];
            testCase.verifyError(@() averageepochs(date, 1, fliplr(epoch)), 'ULMO:averageepochs:InvalidEpochs');
            testCase.verifyError(@() averageepochs(date, 1, [date, NaT]), 'ULMO:averageepochs:InvalidDates');
            testCase.verifyError(@() averageepochs([date;date], ones(1,1,2), epoch), 'ULMO:averageepochs:InvalidDates');
            testCase.verifyError(@() averageepochs(date, ones(1,1,2), epoch), 'ULMO:averageepochs:DimensionMismatch');
            testCase.verifyError(@() averageepochs(date, 1, epoch, 'CheckCompleteness', true), 'ULMO:averageepochs:InvalidTimeStep');
            testCase.verifyError(@() averageepochs(date, 1, epoch, 'SampleStd', -1), 'ULMO:averageepochs:InvalidUncertainty');
            testCase.verifyError(@() averageepochs(date, 1, epoch, 'SampleStd', [1,2]), 'ULMO:averageepochs:DimensionMismatch');
        end
    end
end
