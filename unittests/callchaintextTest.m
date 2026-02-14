classdef callchaintextTest < matlab.unittest.TestCase
    % Unit tests for the callchaintext function

    methods (Test)

        function testEmptyInput(testCase)
            % Test with empty input
            actualOutput = callchaintext('');
            testCase.verifyEmpty(actualOutput);
            actualOutput = callchaintext("");
            testCase.verifyEmpty(actualOutput);
            actualOutput = callchaintext({});
            testCase.verifyEmpty(actualOutput);
        end

        function testSingleFunctionInput(testCase)
            % Test with a single function name
            fun = "sum";
            expectedOutput = sprintf('<a href="matlab: open(''%s'')">%s</a>', which(fun), fun);
            actualOutput = callchaintext(fun);
            testCase.verifyEqual(actualOutput, expectedOutput);
        end

        function testMultipleFunctionsInput(testCase)
            % Test with multiple function names
            funs = {"sum", "numel", "ndims"};
            expectedOutput = sprintf('<a href="matlab: open(''%s'')">%s</a>>', which(funs{1}), funs{1});

            for i = 2:length(funs)
                expectedOutput = strcat(expectedOutput, sprintf('<a href="matlab: open(''%s'')">%s</a>>', which(funs{i}), funs{i}));
            end

            expectedOutput = expectedOutput(1:end - 1); % Remove trailing '>'
            actualOutput = callchaintext(funs);
            testCase.verifyEqual(actualOutput, expectedOutput);
        end

        function testColumnCellVectorInput(testCase)
            % Test with a column cell vector input
            funs = {"sum"; "numel"};
            expectedOutput = sprintf('<a href="matlab: open(''%s'')">%s</a>>', which(funs{1}), funs{1});
            expectedOutput = strcat(expectedOutput, sprintf('<a href="matlab: open(''%s'')">%s</a>', which(funs{2}), funs{2}));
            actualOutput = callchaintext(funs);
            testCase.verifyEqual(actualOutput, expectedOutput);
        end

        function testConsecutiveIdenticalFunctions(testCase)
            % Test with consecutive identical function names
            funs = {"sum", "sum", "numel", "numel", "ndims"};
            expectedOutput = sprintf('<a href="matlab: open(''%s'')">%s</a>>', which(funs{1}), funs{1});
            expectedOutput = strcat(expectedOutput, sprintf('<a href="matlab: open(''%s'')">%s</a>>', which(funs{3}), funs{3}));
            expectedOutput = strcat(expectedOutput, sprintf('<a href="matlab: open(''%s'')">%s</a>', which(funs{5}), funs{5}));
            actualOutput = callchaintext(funs);
            testCase.verifyEqual(actualOutput, expectedOutput);
        end

        function testInvalidInput(testCase)
            % Test with invalid input
            testCase.verifyError(@() callchaintext(123), 'MATLAB:validators:mustBeA');
        end

    end

end
