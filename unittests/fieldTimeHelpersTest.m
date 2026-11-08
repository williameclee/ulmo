%% FIELDTIMEHELPERSTEST - Tests temporal field processing and projection cache tokens.
% Covers native date handling, independent monthly means, selection by
% window centre, fixed-step interpolation and averaging, and cache tokens.
%
% Syntax
%   results = runtests('unittests/fieldTimeHelpersTest.m')
%
% Created by
%   2026/11/07, En-Chi Lee (williameclee@gmail.com)

classdef fieldTimeHelpersTest < matlab.unittest.TestCase

    methods (Test)

        function normalisationKeepsProjectionTokens(t)
            [step, method, useMean, token] = normalisetimeconfig(0.5, 'LINEAR');
            t.verifyEqual(step, hours(12)); t.verifyEqual(method, 'linear');
            t.verifyFalse(useMean); t.verifyEqual(token, 'T43200seconds-linear');
            [step, ~, useMean, token] = normalisetimeconfig('GRACE', 'pchip');
            t.verifyEqual(step, "grace"); t.verifyTrue(useMean); t.verifyEqual(token, 'grace');
            [~, ~, useMean, token] = normalisetimeconfig([], 'mean');
            t.verifyFalse(useMean); t.verifyEqual(token, 'native');
            [~, ~, useMean, token] = normalisetimeconfig('MIDMONTH', 'linear');
            t.verifyTrue(useMean); t.verifyEqual(token, 'midmonth');
            [~, ~, ~, a] = normalisetimeconfig(days(1), 'linear');
            [~, ~, ~, b] = normalisetimeconfig(seconds(86400 +1e-7), 'linear');
            t.verifyNotEqual(a, b);
            [~, ~, ~, b] = normalisetimeconfig(days(1), 'mean');
            t.verifyNotEqual(a, b);
        end

        function nativeKeepsDataAndIgnoresLimits(t)
            dates = datetime(2020, 1, 1:5);
            input.a = reshape([1 NaN 3 4 5], 1, 1, []);
            [actual, when] = processfieldtime(input, dates, [], 'unused', dates([2 3]));
            t.verifyEqual(actual, input); t.verifyEqual(when, dates(:));
        end

        function monthlyFieldsUseIndependentSamples(t)
            dates = datetime(2020, 1, (1:60).');
            input.a = reshape(1:60, 1, 1, []); input.b = 2 * input.a;
            input.a(:, :, 1:10) = NaN; input.b(:, :, 20:31) = NaN;
            [actual, when] = processfieldtime(input, dates, 'MIDMONTH', 'linear', []);
            t.verifyEqual(when, [datetime(2020, 1, 16, 12, 0, 0); datetime(2020, 2, 15, 12, 0, 0)]);
            t.verifyEqual(actual.a(:), [mean(11:31); mean(32:60)], AbsTol = 1e-13);
            t.verifyEqual(actual.b(:), 2 * [mean(1:19); mean(32:60)], AbsTol = 1e-13);
        end

        function limitsSelectCentresWithoutClippingSamples(t)
            dates = datetime(2020, 1, (1:31).');
            input.a = reshape(1:31, 1, 1, []);
            limits = datetime(2020, 1, [15 17]);
            [actual, when] = processfieldtime(input, dates, 'midmonth', 'mean', limits);
            t.verifyEqual(actual.a, 16); t.verifyEqual(when, datetime(2020, 1, 16, 12, 0, 0));
            [actual, when] = processfieldtime(input, dates, 'midmonth', 'mean', datetime(2020, 1, [1 2]));
            t.verifyEmpty(when); t.verifySize(actual.a, [1 1 0]);
        end

        function fixedStepsMatchKnownLinearSignal(t)
            dates = datetime(2020, 1, (1:2:11).');
            input.a = reshape(2 * (0:2:10) + 3, 1, 1, []);
            [actual, when] = processfieldtime(input, dates, 1, 'linear', []);
            t.verifyEqual(when, datetime(2020, 1, (1:11).'));
            t.verifyEqual(actual.a(:), (2 * (0:10) + 3).', AbsTol = 1e-12);
        end

        function fixedMeanUsesHalfOpenWindows(t)
            dates = datetime(2020, 1, (1:5).');
            input.a = reshape(1:5, 1, 1, []);
            [actual, when] = processfieldtime(input, dates, days(2), 'mean', []);
            t.verifyEqual(when, dates([1 3 5]));
            t.verifyEqual(actual.a(:), [1; 2.5; 4.5]);
        end

    end

end
