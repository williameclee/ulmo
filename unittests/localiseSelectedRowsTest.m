%% LOCALISESELECTEDROWSTEST - Compare selected coefficients with full localization.
% Covers signal and diagonal STD propagation, padding/truncation, metadata,
% single/multiple epochs, inverse kernels and requested row ordering.
% Uses synthetic arrays without reading or writing scientific data.
%
% Created by
%   2026/10/26, OpenAI Codex

classdef localiseSelectedRowsTest < matlab.unittest.TestCase

    properties (TestParameter)
        sourceDegree = {2, 5, 8}
        epochCount = {1, 4}
        metadata = {false, true}
        isError = {false, true}
        inverse = {false, true}
    end

    methods (TestClassSetup)
        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            testCase.assertEqual(which('localise'), fullfile(repo, 'aux', 'localise.m'));
        end
    end

    methods (Test)
        function matchesFullLocalization(testCase, sourceDegree, epochCount, metadata, isError, inverse)
            L = 5; n = (L + 1)^2;
            K = reshape(sin(1:n^2), n, n)/n + 0.6*eye(n);
            rows = (sourceDegree + 1)*(sourceDegree + 2)/2;
            input = reshape(cos(1:rows*2*epochCount), rows, 2, epochCount);
            if metadata
                [order, degree] = addmon(sourceDegree);
                input = cat(2, repmat([degree, order], 1, 1, epochCount), input);
            end
            opts = {'L', L, 'K', K, 'IsError', isError, 'Inverse', inverse};
            [full, expectedKernel] = localise(input, opts{:});
            if metadata
                full = full(:, 3:4, :);
            end
            full = reshape(full, (L + 1)*(L + 2), epochCount);
            order = kernelorder(L);
            full = full(order, :);
            % Include repeated/nonascending indices and a single selected row.
            for selected = {[n, 4, 1, 3, 4], 2, (n:-1:1)'}
                ids = selected{1};
                [actual, actualKernel] = localise(input, opts{:}, 'OutputIndices', ids);
                testCase.verifySize(actual, [numel(ids), epochCount]);
                testCase.verifyEqual(actual, full(ids, :), 'AbsTol', 1e-12);
                testCase.verifyEqual(actualKernel, expectedKernel);
            end
        end

        function emptySelectionPreservesDefaultFormat(testCase)
            input = reshape(1:48, 6, 4, 2);
            K = reshape(sin(1:81), 9, 9);
            expected = localise(input, 'L', 2, 'K', K);
            [actual, actualKernel] = localise(input, 'L', 2, 'K', K, 'OutputIndices', []);
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(actualKernel, K);
        end

        function rejectsInvalidIndices(testCase)
            for selected = {0, -1, 1.5, NaN, Inf, 10, [1, 2; 3, 4], 1 + 1i}
                testCase.verifyError(@() localise(ones(6, 2), 'L', 2, 'K', eye(9), ...
                    'OutputIndices', selected{1}), 'ULMO:localise:InvalidOutputIndices');
            end
        end
    end
end
