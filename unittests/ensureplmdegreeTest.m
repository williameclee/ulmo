%% ENSUREPLMDEGREETEST - Verifies harmonic resizing and array-format preservation.
% Uses synthetic coefficients without loading scientific data or caches.
% Covers single/multiple epochs, degree zero, padding, truncation, no-op
% resizing, default format, format recovery and malformed input shapes.
%
% Syntax
%   results = runtests('unittests/ensureplmdegreeTest.m');
%
% Notes
%   Add slepian_alpha to the MATLAB path when it is not beside this checkout.
%
% Created by
%   2026/10/24, En-Chi Lee (williameclee@gmail.com)

classdef ensureplmdegreeTest < matlab.unittest.TestCase

    properties (TestParameter)
        format = {'traditional', 'timefirst'}
        epochCount = {1, 3}
        targetDegree = {0, 1, 2, 4}
    end

    methods (TestClassSetup)

        function configurePaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            dependency = fullfile(fileparts(repo), 'slepian_alpha');

            if isfolder(dependency)
                testCase.applyFixture(PathFixture(dependency));
            end

            testCase.assertNotEmpty(which('addmon'), ...
                'Add slepian_alpha to the MATLAB path before running these tests.');
            testCase.assertNotEmpty(which('addmup'), ...
                'Add slepian_alpha to the MATLAB path before running these tests.');
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            % Guard against the experiment-local helper shadowing production.
            testCase.assertEqual(which('ensureplmdegree'), ...
                fullfile(repo, 'aux', 'ensureplmdegree.m'));
        end

    end

    methods (Test)

        function resizePreservesValuesAndFormat(testCase, format, epochCount, targetDegree)
            input = makeCoefficients(2, epochCount);
            expected = makeCoefficients(targetDegree, epochCount);

            if targetDegree > 2
                expected(7:end, 3:4, :) = 0;
            end

            input = formatCoefficients(input, format);
            expected = formatCoefficients(expected, format);

            actual = ensureplmdegree(input, targetDegree, format);

            testCase.verifySize(actual, size(expected));
            testCase.verifyEqual(actual, expected);
        end

        function truncateGraceDegree60To45(testCase, format)
            input = formatCoefficients(makeCoefficients(60, 3), format);
            expected = formatCoefficients(makeCoefficients(45, 3), format);
            testCase.verifyEqual(ensureplmdegree(input, 45, format), expected);
        end

        function defaultFormatPreservesSingleEpoch(testCase)
            input = formatCoefficients(makeCoefficients(2, 1), 'timefirst');
            expected = formatCoefficients(makeCoefficients(1, 1), 'timefirst');
            testCase.verifyEqual(ensureplmdegree(input, 1), expected);
        end

        function recoverMislabeledFormat(testCase, format)
            input = formatCoefficients(makeCoefficients(2, 3), format);
            expected = formatCoefficients(makeCoefficients(1, 3), format);

            if strcmp(format, 'traditional')
                requestedFormat = 'timefirst';
            else
                requestedFormat = 'traditional';
            end

            actual = testCase.verifyWarning( ...
                @() ensureplmdegree(input, 1, requestedFormat), ...
            'ULMO:InvalidInputSize');
            testCase.verifyEqual(actual, expected);
        end

        function rejectMalformedShape(testCase, format)
            testCase.verifyError(@() ensureplmdegree(zeros(2, 6, 3), 1, format), ...
            'ULMO:InvalidInputSize');
        end

    end

end

function plm = makeCoefficients(L, epochCount)
    % Construct expected harmonic ordering independently of ADDMON/ADDMUP.
    count = (L + 1) * (L + 2) / 2;
    plm = zeros(count, 4, epochCount);
    row = 0;

    for degree = 0:L

        for order = 0:degree
            row = row + 1;

            for epoch = 1:epochCount
                plm(row, :, epoch) = ...
                    [degree, order, 10000 * epoch + row, -20000 * epoch - row];
            end

        end

    end

end

function plm = formatCoefficients(plm, format)

    if strcmp(format, 'timefirst')
        plm = permute(plm, [3, 1, 2]);
    end

end
