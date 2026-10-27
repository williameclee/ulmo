%% OCEANVARIANCEKERNELTEST - Checks blocked construction and reusable SLE kernels.
% Compares with the original dense constrained expression and actual SLE
% outputs. Synthetic Love numbers and temporary fixtures isolate all writes.
% A counting fixture checks one construction per fingerprint degree-1 run
% and no construction for uniform loading.
%
% Syntax
%   results = runtests('unittests/oceanVarianceKernelTest.m')
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@gmail.com)

classdef oceanVarianceKernelTest < matlab.unittest.TestCase

    properties (TestParameter)
        kernelSize = {1, 9, 256, 289, 576}
        epochCnt = {1, 3}
        rotation = {false, true}
        tol = {0, 1e-12}
        columns = {2, 4}
        mthd = {'fingerprint', 'uniform'}
    end

    methods (TestClassSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            testCase.applyFixture(PathFixture(fullfile(repo, 'domains')));
            testCase.assertEqual(which('solvesle'), fullfile(repo, 'solvesle.m'));
            folder = testCase.applyFixture(TemporaryFolderFixture);
            previous = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', previous));
            setenv('IFILES', folder.Folder);
            mkdir(fullfile(folder.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; .021; - .3];
            loadingverticaldisplacement = [0; - .269; - .6];
            tidalgravitationalpotential = [0; 0; .3];
            tidalverticaldisplacement = [0; 0; .6];
            save(fullfile(folder.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential', 'loadingverticaldisplacement', ...
                'tidalgravitationalpotential', 'tidalverticaldisplacement');
        end

    end

    methods (Test)

        function degree1BuildsOnce(testCase, mthd)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            folder = testCase.applyFixture(TemporaryFolderFixture);
            fixtures = fullfile(repo, 'unittests', 'fixtures', 'degree1convergence');
            files = dir(fullfile(fixtures, '*.txt'));

            for k = 1:numel(files)
                copyfile(fullfile(fixtures, files(k).name), ...
                    fullfile(folder.Folder, erase(files(k).name, '.txt')));
            end

            fid = fopen(fullfile(folder.Folder, 'prepareoceanvariancekernel.m'), 'w');
            fprintf(fid, ['function Q = prepareoceanvariancekernel(K)\n' ...
                              'global DEGREE1_CONVERGENCE_FIXTURE\n' ...
                              'DEGREE1_CONVERGENCE_FIXTURE.builds = DEGREE1_CONVERGENCE_FIXTURE.builds + 1;\n' ...
                          'C = K - K(:,1)/K(1,1)*K(1,:); C(:,1)=0; Q=C.^2;\nend\n']);
            fclose(fid);
            testCase.applyFixture(PathFixture(folder.Folder));
            previous = getenv('GRACEDATA');
            testCase.addTeardown(@() setenv('GRACEDATA', previous));
            setenv('GRACEDATA', folder.Folder);
            testCase.addTeardown(@() clearDegree1Fixture());
            % The existing isolated readers share this fixture state.
            global DEGREE1_CONVERGENCE_FIXTURE %#ok<GVMIS>
            DEGREE1_CONVERGENCE_FIXTURE = struct('x', reshape(sin(1:27), 9, 3), ...
                'sigma', ones(9, 3) * .1, 'builds', 0);
            [coeffs, sigma] = solvedegree1('CSR', 'RL06', 2, 2, 2, 'ICE6GD', ...
                GeoDomain('alloceans'), 'Method', mthd, 'Tolerance', 0, 'MaxIter', 3, ...
                'SleTolerance', 0, 'SleMaxIter', 3, 'ForceNew', true, ...
                'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(DEGREE1_CONVERGENCE_FIXTURE.builds, double(strcmp(mthd, 'fingerprint')));
            testCase.verifyTrue(all(isfinite([coeffs, sigma]), 'all'));
        end

        function blockedMatchesDenseExpression(testCase, kernelSize)
            ids = 1:kernelSize;
            K = 0.6 * eye(kernelSize) + sin(ids' * ids) / (10 * kernelSize);
            expd = denseVarianceKernel(K);
            actl = prepareoceanvariancekernel(K);
            testCase.verifyEqual(actl, expd, 'AbsTol', 1e-15);
            testCase.verifyEqual(actl(:, 1), zeros(kernelSize, 1));
            var = reshape(abs(cos(1:3 * kernelSize)), kernelSize, 3);
            testCase.verifyEqual(actl * var, expd * var, 'AbsTol', 1e-13);
        end

        function suppliedMatchesFallback(testCase, epochCnt, rotation, tol, columns)
            K = syntheticKernel();
            Q = denseVarianceKernel(K);
            f = reshape(sin(1:9 * epochCnt), 9, epochCnt);
            sigma = abs(f) / 10;

            if epochCnt > 1
                sigma(:, end) = 0;
            end

            forc = unflatten(f, columns);
            unc = unflatten(sigma, columns);
            opts = {'OceanKernel', K, 'RotationFeedback', rotation, ...
                        'Tolerance', tol, 'MaxIter', 100, 'BeQuiet', true};

            if tol == 0
                opts = [opts, {'MaxIter', 5}];
            end

            actl = cell(1, 4);
            ref = cell(1, 4);
            [ref{:}] = solvesle(forc, unc, 2, GeoDomain('alloceans'), opts{:});
            [actl{:}] = solvesle(forc, unc, 2, GeoDomain('alloceans'), ...
                opts{:}, 'OceanVarianceKernel', Q);

            for k = 1:4
                testCase.verifyEqual(actl{k}, ref{k}, 'AbsTol', 1e-12);
            end

            ocn = K * flatten(actl{1});
            testCase.verifyEqual(ocn(1, :), -f(1, :), 'AbsTol', 1e-11);
            testCase.verifyEqual(actl{4}, (sigma(1, :) / K(1, 1))', 'AbsTol', 1e-12);
            testCase.verifyEqual(Q, denseVarianceKernel(K));
        end

        function preparedOrNoErrorsSkipsConstruction(testCase)
            % A throwing builder detects unnecessary fallback allocations.
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            folder = testCase.applyFixture(TemporaryFolderFixture);
            fid = fopen(fullfile(folder.Folder, 'prepareoceanvariancekernel.m'), 'w');
            fprintf(fid, ['function Q = prepareoceanvariancekernel(varargin)\n' ...
                          'error(''ULMO:UnexpectedVarianceConstruction'', ''Unexpected builder call'');\nend\n']);
            fclose(fid);
            testCase.applyFixture(PathFixture(folder.Folder));
            K = syntheticKernel();
            forcing = unflatten((1:9)', 2);
            opts = {'OceanKernel', K, 'Tolerance', 0, 'MaxIter', 3, 'BeQuiet', true};
            [~, sigma] = solvesle(forcing, [], 2, GeoDomain('alloceans'), opts{:});
            testCase.verifyEmpty(sigma);
            [~, sigma] = solvesle(forcing, forcing / 10, 2, GeoDomain('alloceans'), ...
                opts{:}, 'OceanVarianceKernel', denseVarianceKernel(K));
            testCase.verifyTrue(all(isfinite(sigma), 'all'));
            testCase.verifyError(@() solvesle(forcing, forcing / 10, 2, ...
                GeoDomain('alloceans'), opts{:}), 'ULMO:UnexpectedVarianceConstruction');
        end

        function suppliedRequiresMatchingDimensions(testCase)
            p = zeros(6, 2);
            common = {p, p, 2, GeoDomain('alloceans'), 'BeQuiet', true};
            testCase.verifyError(@() solvesle(common{:}, 'OceanVarianceKernel', eye(9)), ...
            'ULMO:solvesle:MissingOceanKernel');
            testCase.verifyError(@() solvesle(common{:}, 'OceanKernel', eye(9), ...
                'OceanVarianceKernel', eye(4)), 'ULMO:solvesle:InvalidOceanVarianceKernelSize');
            testCase.verifyError(@() solvesle(common{:}, 'OceanKernel', eye(4), ...
                'OceanVarianceKernel', eye(4)), 'ULMO:solvesle:InvalidOceanVarianceKernelSize');
        end

        function invalidConstructionInputsAreRejected(testCase)
            testCase.verifyError(@() prepareoceanvariancekernel([]), ...
            'ULMO:prepareoceanvariancekernel:InvalidSize');
            testCase.verifyError(@() prepareoceanvariancekernel(ones(3, 2)), ...
            'ULMO:prepareoceanvariancekernel:InvalidSize');

            for area = [0, -1, NaN, Inf]
                K = eye(9); K(1, 1) = area;
                testCase.verifyError(@() prepareoceanvariancekernel(K), ...
                'ULMO:prepareoceanvariancekernel:InvalidOceanArea');
            end

            K = eye(9); K(8, 7) = Inf;
            testCase.verifyError(@() prepareoceanvariancekernel(K), ...
            'ULMO:prepareoceanvariancekernel:NonFiniteKernel');
        end

    end

end

function K = syntheticKernel()
    ids = 1:9;
    K = 0.6 * eye(9) + 0.002 * (sin(ids' * ids) + cos(ids' + ids));
end

function clearDegree1Fixture()
    clear global DEGREE1_CONVERGENCE_FIXTURE
end

function Q = denseVarianceKernel(K)
    constrained = K - K(:, 1) / K(1, 1) * K(1, :);
    constrained(:, 1) = 0;
    Q = constrained .^ 2;
end

function p = unflatten(v, columns)
    ko = kernelorder(2);
    p = zeros(12, size(v, 2)); p(ko, :) = v;
    p = reshape(p, 6, 2, size(v, 2));

    if columns == 4
        [order, degree] = addmon(2);
        p = cat(2, repmat([degree, order], 1, 1, size(v, 2)), p);
    end

end

function v = flatten(p)

    if size(p, 2) == 4
        p = p(:, 3:4, :);
    end

    p = reshape(p, 12, []);
    v = p(kernelorder(2), :);
end
