%% absoluteIterationChangeTest - Checks absolute updates across coefficients and epochs.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

classdef absoluteIterationChangeTest < matlab.unittest.TestCase

    methods (Test)

        function zeroAndFirstStep(testCase)
            testCase.verifyEqual(absoluteiterationchange(zeros(3, 2), zeros(3, 2)), 0);
            testCase.verifyEqual(absoluteiterationchange([0; 1e-200], zeros(2, 1)), 1e-200);
            testCase.verifyEqual(absoluteiterationchange(zeros(2, 1), [0; 1e-200]), 1e-200);
        end

        function checksEveryElement(testCase)
            previous = [1e9, 0; 0, 0]; current = previous; current(2, 2) = 2e-9;
            testCase.verifyEqual(absoluteiterationchange(current, previous), 2e-9);
        end

        function nonfiniteCannotConverge(testCase)
            testCase.verifyEqual(absoluteiterationchange([NaN; 0], [0; 0]), Inf);
        end

    end

end
