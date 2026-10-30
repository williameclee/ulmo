%% ESMEPOCHSTEST - Checks epoch generation and calendar-independent averaging.
% Uses explicit irregular windows to exercise the interface needed by a
% future GRACE epoch provider without enabling that time-step option.
%
% Created by
%   2026/10/30

classdef esmEpochsTest < matlab.unittest.TestCase
    methods (Test)
        function calendarEpochsSpanLeapYearAndGaps(testCase)
            nativeDates = [datetime(2004, 1, 31); datetime(2004, 3, 1)];
            epochs = esmepochs('midmonth', nativeDates);
            starts = datetime(2004, (1:3)', 1);
            testCase.verifyEqual(epochs, [starts, starts + calmonths(1)]);
            testCase.verifyEqual(epochs(2, 2) - epochs(2, 1), days(29));
        end
        function irregularOverlappingWindowsUseSuppliedBounds(testCase)
            [dates, raw] = testCase.fields();
            start = dates(1);
            % Windows cross a month boundary, overlap and have off-grid bounds.
            epochs = [start + hours(7), start + hours(25); ...
                start, start + hours(12); start + hours(6), start + hours(18)];
            [p, s, d] = averageesmepochs(dates, raw, epochs);
            testCase.verifyEqual(p(:, 1, 3), [4; 1.5; 2.5]);
            testCase.verifyEqual(p(:, 1, 4), -[4; 1.5; 2.5]);
            testCase.verifyEqual(d, epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2);
            testCase.verifyEqual(s(:, :, 1:2), p(:, :, 1:2));
            testCase.verifyTrue(all(isnan(s(:, :, 3:4)), 'all'));
        end
        function incompleteWindowsDoNotAffectOtherEpochs(testCase)
            [dates, raw] = testCase.fields();
            epochs = [dates(1), dates(3); dates(1), dates(4); ...
                dates(1) + hours(1), dates(1) + hours(2)];
            keep = [1, 2, 4, 5];
            call = @() averageesmepochs(dates(keep), raw(keep, :, :), epochs);
            testCase.verifyWarning(call, 'ULMO:esm2plmt:IncompleteEpoch');
            warningState = warning('off', 'ULMO:esm2plmt:IncompleteEpoch');
            cleanup = onCleanup(@() warning(warningState));
            p = call();
            testCase.verifyEqual(p(1, 1, 3), 1.5);
            testCase.verifyTrue(all(isnan(p(2:3, :, 3:4)), 'all'));
            testCase.verifyEqual(p(:, :, 1:2), zeros(3, 1, 2));
        end
        function emptyEpochListReturnsEmptyArrays(testCase)
            [dates, raw] = testCase.fields();
            [p, s, d] = averageesmepochs(dates, raw, NaT(0, 2));
            testCase.verifySize(p, [0, 1, 4]);
            testCase.verifySize(s, [0, 1, 4]);
            testCase.verifySize(d, [0, 1]);
        end
        function invalidEpochsAndDatesFail(testCase)
            [dates, raw] = testCase.fields();
            testCase.verifyError(@() averageesmepochs(dates, raw, [dates(2), dates(1)]), ...
                'ULMO:averageesmepochs:InvalidEpochs');
            testCase.verifyError(@() averageesmepochs(dates, raw, [dates(1), NaT]), ...
                'ULMO:averageesmepochs:InvalidEpochs');
            testCase.verifyError(@() averageesmepochs(dates(1:2), raw, [dates(1), dates(3)]), ...
                'ULMO:averageesmepochs:DimensionMismatch');
            dates(2) = dates(1);
            testCase.verifyError(@() averageesmepochs(dates, raw, [dates(1), dates(3)]), ...
                'ULMO:averageesmepochs:InvalidDates');
        end
        function graceOptionRemainsUnsupported(testCase)
            testCase.verifyError(@() esmepochs('GRACE', datetime(2003, 1, 1)), ...
                'ULMO:esm2plmt:InvalidTimeStep');
            testCase.verifyError(@() esm2plmt('O', 'GRACE'), 'ULMO:esm2plmt:InvalidTimeStep');
        end
    end
    methods (Static, Access = private)
        function [dates, raw] = fields()
            dates = datetime(2003, 1, 31, 18, 0, 0) + hours((0:4)' * 6);
            raw = zeros(5, 1, 4);
            raw(:, 1, 3) = (1:5)';
            raw(:, 1, 4) = -(1:5)';
        end
    end
end
