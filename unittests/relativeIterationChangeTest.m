classdef relativeIterationChangeTest < matlab.unittest.TestCase
    % Cover zero vectors and independent epoch normalization without a unit floor.
    methods (Test)

        function zeroAndFirstStep(testCase)
            testCase.verifyEqual(relativeiterationchange(zeros(3, 2), zeros(3, 2)), 0);
            testCase.verifyEqual(relativeiterationchange([0; 1e-200], zeros(2, 1)), 1);
            testCase.verifyEqual(relativeiterationchange(zeros(2, 1), [0; 1e-200]), 1);
        end

        function smallEpochCannotBeHidden(testCase)
            previous = [1e9, 1e-9; 0, 0]; current = previous; current(1, 2) = 2e-9;
            testCase.verifyEqual(relativeiterationchange(current, previous), .5, 'AbsTol', eps);
        end

        function nonfiniteCannotConverge(testCase)
            testCase.verifyEqual(relativeiterationchange([NaN; 0], [0; 0]), Inf);
        end

    end

end
