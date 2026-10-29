%% LINEARSLESYSTEMTEST - Checks full and selected responses against dense equations.
% Uses isolated synthetic Love numbers, ocean kernels and GRACE/GIA readers.
% Tests induced covariance, transposed solves, joint coefficient sets, cache
% reuse, physical mass balance and failure handling without production data.
%
% Created by
%   2026/10/28, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/29, En-Chi Lee (williameclee@arizona.edu)

classdef linearSleSystemTest < matlab.unittest.TestCase

    properties
        Folder
    end

    properties (TestParameter)
        kernelSource = {"ocean", "sle", "all"}
        rotation = {false, true}
        method = {"fingerprint", "uniform"}
        c20 = {false, true}
        c30 = {false, true}
        cs21 = {false, true}
    end

    methods (TestMethodSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            testCase.applyFixture(PathFixture(fullfile(repo, 'legacy')));
            testCase.applyFixture(PathFixture(fullfile(repo, 'domains')));
            temp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.Folder = temp.Folder;
            fixtures = fullfile(repo, 'unittests', 'fixtures', 'jointlowdegreesolver');
            files = dir(fullfile(fixtures, '*.txt'));

            for j = 1:numel(files)
                copyfile(fullfile(fixtures, files(j).name), ...
                    fullfile(temp.Folder, erase(files(j).name, '.txt')));
            end

            testCase.applyFixture(PathFixture(temp.Folder));

            for name = {'IFILES', 'GRACEDATA'}
                previous = getenv(name{1});
                testCase.addTeardown(@() setenv(name{1}, previous));
                setenv(name{1}, temp.Folder);
            end

            mkdir(fullfile(temp.Folder, 'EARTHMODELS', 'CONSTANTS'));
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(temp.Folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
            mkdir(fullfile(temp.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; .021; - .3; - .2];
            loadingverticaldisplacement = [0; - .269; - .6; - .4];
            tidalgravitationalpotential = [0; 0; .3; .1];
            tidalverticaldisplacement = [0; 0; .6; .2];
            save(fullfile(temp.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential', 'loadingverticaldisplacement', ...
                'tidalgravitationalpotential', 'tidalverticaldisplacement');
            testCase.addTeardown(@() clearFixture());
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, 3);
            p = cat(2, meta, unpack(reshape(sin(1:48), 16, 3)));
            s = cat(2, meta, unpack(reshape(abs(cos(1:48)) / 10, 16, 3)));
            g = cat(2, meta, unpack(reshape(cos(1:48) / 20, 16, 3)));
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(g, [3, 1, 2]), ...
                'dates', datetime(2003, (1:3)', 15), 'kernel', ocean());
        end

    end

    methods (Test)

        function forwardAndTransposeMatchDense(testCase, rotation)
            ocnKer = ocean();
            [earthRespKer, forcKer, waterDen] = ...
                prepareslekernels(3, ocnKer, RotationFeedback = rotation);
            [A, KF] = explicitSystem(ocnKer, earthRespKer);
            rhss = reshape(sin(1:48), 16, 3);
            rhss(:, 2) = 0;
            [solCoeffs, relResids, iterCnt] = solveslesystem(ocnKer, earthRespKer, rhss);
            [transSolCoeffs, transRelResids, transIterCnt] = ...
                solveslesystem(ocnKer, earthRespKer, rhss, Transpose = true);
            testCase.verifyEqual(solCoeffs, A \ rhss, 'AbsTol', 2e-11);
            testCase.verifyEqual(transSolCoeffs, A.' \ rhss, 'AbsTol', 2e-11);
            testCase.verifyEqual(full(forcKer), KF, 'AbsTol', 1e-14);
            testCase.verifyClass(relResids, 'double');
            testCase.verifyClass(iterCnt, 'double');
            testCase.verifySize(relResids, [1, 3]);
            testCase.verifySize(iterCnt, [1, 3]);
            testCase.verifySize(transRelResids, [1, 3]);
            testCase.verifySize(transIterCnt, [1, 3]);
            testCase.verifyEqual(relResids(2), 0);
            testCase.verifyEqual(iterCnt(2), 0);
            testCase.verifyEqual(transRelResids(2), 0);
            testCase.verifyEqual(transIterCnt(2), 0);
            testCase.verifyGreaterThan(iterCnt([1, 3]), zeros(1, 2));
            testCase.verifyGreaterThan(transIterCnt([1, 3]), zeros(1, 2));
            testCase.verifyLessThanOrEqual(relResids, 1e-12);
            testCase.verifyLessThanOrEqual(transRelResids, 1e-12);

            for column = [1, 3]
                rhs = rhss(:, column);
                testCase.verifyEqual(relResids(column), ...
                    norm(A * solCoeffs(:, column) - rhs) / norm(rhs), 'AbsTol', 2e-15);
                testCase.verifyEqual(transRelResids(column), ...
                    norm(A.' * transSolCoeffs(:, column) - rhs) / norm(rhs), 'AbsTol', 2e-15);
            end

            testCase.verifyTrue(issparse(earthRespKer) && issparse(forcKer));
            testCase.verifyClass(earthRespKer, 'double');
            testCase.verifyClass(forcKer, 'double');
            testCase.verifySize(earthRespKer, [16, 16]);
            testCase.verifySize(forcKer, [16, 16]);
            testCase.verifyEqual(waterDen, 1000);
            testCase.verifyGreaterThan(norm(ocnKer * ocnKer - ocnKer, 'fro'), 0.1);
        end

        function massAndFullCovariance(testCase, rotation)
            ocnKer = ocean();
            [earthRespKer, ~] = ...
                prepareslekernels(3, ocnKer, RotationFeedback = rotation);
            [A, KF] = explicitSystem(ocnKer, earthRespKer);
            fullSleResp = A \ KF;
            f = reshape(sin(1:48), 16, 3);
            sigma = reshape(abs(cos(1:48)), 16, 3) / 10;
            sigma(:, 1) = 0;
            [s, ss, meanLevel, meanStd] = solvesle(unpack(f), unpack(sigma), ...
                3, GeoDomain('alloceans'), 'OceanKernel', ocnKer, ...
                'RotationFeedback', rotation, 'BeQuiet', true);
            testCase.verifyEqual(pack(s), fullSleResp * f, 'AbsTol', 2e-11);
            testCase.verifyEqual(pack(ss) .^ 2, (fullSleResp .^ 2) * (sigma .^ 2), 'AbsTol', 2e-13);
            local = ocnKer * pack(s);
            testCase.verifyEqual(local(1, :), -f(1, :), 'AbsTol', 2e-11);
            testCase.verifyEqual(meanLevel, (-f(1, :) / ocnKer(1, 1)).', 'AbsTol', 2e-11);
            testCase.verifyEqual(meanStd, (sigma(1, :) / ocnKer(1, 1)).', 'AbsTol', 2e-13);

            for epoch = 1:3
                C = fullSleResp * diag(sigma(:, epoch) .^ 2) * fullSleResp.';
                testCase.verifyEqual(pack(ss(:, :, epoch)) .^ 2, diag(C), 'AbsTol', 2e-13);
            end

        end

        function selectedMapsAndPublicSolver(testCase, method, c20, c30, cs21, rotation)
            global JOINT_LOW_DEGREE_FIXTURE
            f = JOINT_LOW_DEGREE_FIXTURE;
            ocnKer = f.kernel;
            [earthRespKer, forcKer] = prepareslekernels(3, ocnKer, RotationFeedback = rotation);
            [A, KF] = explicitSystem(ocnKer, earthRespKer);
            [locs, ids] = lowdegreecoefficients(c20, c30, cs21);

            if method == "fingerprint"
                expOcnResp = f.kernel(ids, :) * (A \ KF) * (eye(16) - f.kernel);
            else
                % The uniform response consumes only the land C00 row.
                expOcnResp = -f.kernel(ids, 1) / f.kernel(1, 1) * ([1, zeros(1, 15)] - f.kernel(1, :));
            end

            [kn2unknMap, knIds, ocnRespMsp, unknCstr, knCtsr, rCond, relResid, iterCnt] = ...
                preparelowdegreemap(ocnKer, earthRespKer, forcKer, ids, Method = method);
            known = setdiff(1:16, ids);
            expLoadResp = (expOcnResp(:, ids) - f.kernel(ids, ids)) \ (f.kernel(ids, known) - expOcnResp(:, known));
            testCase.verifyEqual(ocnRespMsp, expOcnResp, 'AbsTol', 2e-11);
            testCase.verifyEqual(kn2unknMap, expLoadResp, 'AbsTol', 2e-11);
            testCase.verifyEqual(knIds, known);
            testCase.verifyEqual(unknCstr, ...
                expOcnResp(:, ids) - f.kernel(ids, ids), 'AbsTol', 2e-11);
            testCase.verifyEqual(knCtsr, ...
                f.kernel(ids, known) - expOcnResp(:, known), 'AbsTol', 2e-11);
            testCase.verifyEqual(rCond, rcond(unknCstr));
            testCase.verifySize(relResid, [1, numel(ids)]);
            testCase.verifySize(iterCnt, [1, numel(ids)]);
            testCase.verifyLessThanOrEqual(max(relResid), 1e-12);

            if method == "uniform"
                testCase.verifyEqual(iterCnt, zeros(1, numel(ids)));
            else
                testCase.verifyGreaterThan(iterCnt, zeros(1, numel(ids)));
            end

            p = pack(permute(f.p(:, :, 3:4), [2, 3, 1]));
            sigma = pack(permute(f.s(:, :, 3:4), [2, 3, 1]));
            gia = pack(permute(f.gia(:, :, 3:4), [2, 3, 1]));
            p(5:end, :) = p(5:end, :) - gia(5:end, :);
            expLoad = expLoadResp * p(known, :);
            expLoadStd = sqrt((expLoadResp .^ 2) * (sigma(known, :) .^ 2));

            if numel(ids) > 3
                expLoad(4:end, :) = expLoad(4:end, :) + gia(ids(4:end), :);
            end

            opts = {'IncludeC20', c20, 'IncludeC30', c30, 'IncludeC21S21', cs21, ...
                        'Method', char(method), 'ForceNew', true, 'SaveData', false, 'BeQuiet', true};
            [actlCoeffs, actlCoeffStds, dates, ~, actlLocs] = solvedegree1( ...
                'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), opts{:}, ...
                'RotationFeedback', rotation);
            testCase.verifyEqual(actlCoeffs, expLoad.', 'AbsTol', 2e-11);
            testCase.verifyEqual(actlCoeffStds, expLoadStd.', 'AbsTol', 2e-11);
            testCase.verifyEqual(dates, f.dates);
            testCase.verifyEqual(actlLocs, locs);
            % Independently check the central coefficients against converged iteration.
            if rotation || method == "uniform"
                legacy = solvedegree1_iteration('CSR', 'RL06', 3, 3, 3, 'ICE6GD', ...
                    GeoDomain('alloceans'), opts{:}, ...
                    'Tolerance', 1e-12, 'SleTolerance', 1e-14);
                testCase.verifyEqual(actlCoeffs, legacy, 'AbsTol', 2e-9);
            end
        end

        function suppliedKernelsSkipConstructionAndCache(testCase, rotation, kernelSource)
            global JOINT_LOW_DEGREE_FIXTURE
            args = solverArgs();
            [expected, expectedStd, dates] = solvedegree1(args{:}, ...
                'RotationFeedback', rotation, 'SaveData', true);
            files = dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat'));
            testCase.assertNumElements(files, 1);
            % A sentinel cache must never hide explicitly supplied matrices.
            coeffs = expected + 100;
            coeffStds = expectedStd + 100;
            save(fullfile(files.folder, files.name), 'coeffs', 'coeffStds', 'dates');
            ocnKer = ocean();
            [earthRespKer, forcKer] = prepareslekernels(3, ocnKer, RotationFeedback = rotation);
            supplied = {};

            if kernelSource == "ocean" || kernelSource == "all"
                supplied = {'OceanKernel', ocnKer};
                % Invalid fallback geometry proves the supplied ocean matrix is used.
                JOINT_LOW_DEGREE_FIXTURE.kernel = NaN(16);
            end

            if kernelSource == "sle" || kernelSource == "all"
                supplied = [supplied, {'EarthResponseKernel', earthRespKer, 'ForcingKernel', forcKer}];
                % Missing Love numbers prove the SLE matrices are not rebuilt.
                movefile(fullfile(testCase.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                    fullfile(testCase.Folder, 'saved-lovenumbers.mat'));
            end

            [actual, actualStd, actualDates] = solvedegree1(args{:}, supplied{:}, ...
                'RotationFeedback', rotation, 'SaveData', false);
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-11);
            testCase.verifyEqual(actualStd, expectedStd, 'AbsTol', 2e-11);
            testCase.verifyEqual(actualDates, dates);
            stored = load(fullfile(files.folder, files.name));
            testCase.verifyEqual(stored.coeffs, coeffs);
        end

        function rotationCachesAreSeparate(testCase)
            global JOINT_LOW_DEGREE_FIXTURE
            args = solverArgs();
            [withRotation, withStd] = solvedegree1(args{:}, 'SaveData', true);
            [withoutRotation, withoutStd] = solvedegree1(args{:}, ...
                'RotationFeedback', false, 'SaveData', true);
            testCase.verifyGreaterThan(norm(withRotation - withoutRotation), 1e-8);
            files = dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat'));
            dom = GeoDomain('alloceans');
            stem = sprintf('CSRRL063_ICE6GD-Ld3_Ls3-%s', dom.Id);
            testCase.verifyEqual(sort({files.name}), sort({[stem, '-SD.mat'], [stem, '-NRot-SD.mat']}));
            % Remove valid source data to verify each filename is read independently.
            JOINT_LOW_DEGREE_FIXTURE.p = [];
            [actual, actualStd] = solvedegree1(args{:}, 'RotationFeedback', true);
            testCase.verifyEqual(actual, withRotation);
            testCase.verifyEqual(actualStd, withStd);
            [actual, actualStd] = solvedegree1(args{:}, 'RotationFeedback', false);
            testCase.verifyEqual(actual, withoutRotation);
            testCase.verifyEqual(actualStd, withoutStd);
        end

        function uniformIgnoresRotationAndSleKernels(testCase)
            args = solverArgs();
            [expected, expectedStd] = solvedegree1(args{:}, 'Method', 'uniform', 'SaveData', true);
            [actual, actualStd] = solvedegree1(args{:}, 'Method', 'uniform', ...
                'RotationFeedback', false, 'OceanKernel', ocean(), ...
                'EarthResponseKernel', zeros(16), 'ForcingKernel', zeros(16), 'SaveData', true);
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(actualStd, expectedStd);
            files = dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat'));
            testCase.verifyNumElements(files, 1);
            testCase.verifyTrue(endsWith(files.name, '-uniform-SD.mat'));
        end

        function invalidKernelsFailBeforeCacheRead(testCase)
            args = solverArgs();
            solvedegree1(args{:}, 'SaveData', true);
            names = {'OceanKernel', 'EarthResponseKernel', 'ForcingKernel'};
            bad = {eye(9), NaN(16), complex(eye(16), eye(16)), Inf(16), ones(16, 16, 2)};

            for k = 1:numel(names)
                for j = 1:numel(bad)
                    testCase.verifyError(@() solvedegree1(args{:}, names{k}, bad{j}), ...
                        'ULMO:solvedegree1:InvalidKernel');
                end
            end

            testCase.verifyError(@() solvedegree1(args{:}, 'OceanKernel', zeros(16)), ...
                'ULMO:solvedegree1:InvalidKernel');
            testCase.verifyError(@() solvedegree1(args{:}, 'EarthResponseKernel', eye(16)), ...
                'ULMO:solvedegree1:IncompleteSleKernels');
            testCase.verifyError(@() solvedegree1(args{:}, 'ForcingKernel', eye(16)), ...
                'ULMO:solvedegree1:IncompleteSleKernels');
        end

        function paddedInputsIgnoreImportedUnknowns(testCase, method)
            global JOINT_LOW_DEGREE_FIXTURE
            f = JOINT_LOW_DEGREE_FIXTURE;
            [~, ids] = lowdegreecoefficients(false, true, true);
            ocnKer = f.kernel;
            [earthRespKer, forcKer] = prepareslekernels(3, ocnKer);
            [kn2unknMap, knIds] = preparelowdegreemap(ocnKer, earthRespKer, forcKer, ids, Method = method);
            x = pack(permute(f.p(:, :, 3:4), [2, 3, 1]));
            sigma = pack(permute(f.s(:, :, 3:4), [2, 3, 1]));
            gia = pack(permute(f.gia(:, :, 3:4), [2, 3, 1]));
            x(10:end, :) = 0; sigma(10:end, :) = 0;
            x(5:end, :) = x(5:end, :) - gia(5:end, :);
            expd = kn2unknMap * x(knIds, :);
            expd(4:end, :) = expd(4:end, :) + gia(ids(4:end), :);
            args = {'CSR', 'RL06', 3, 2, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                        'IncludeC30', true, 'IncludeC21S21', true, 'Method', char(method), ...
                        'SaveData', false, 'BeQuiet', true};
            [actl, actlStd] = solvedegree1(args{:});
            testCase.verifyEqual(actl, expd.', 'AbsTol', 2e-11);
            testCase.verifyEqual(actlStd, sqrt((kn2unknMap .^ 2) * (sigma(knIds, :) .^ 2)).', 'AbsTol', 2e-11);
            [locs, ~] = lowdegreecoefficients(false, true, true);

            for j = 1:size(locs, 1)
                JOINT_LOW_DEGREE_FIXTURE.p(:, locs(j, 1), locs(j, 2) + 2) = 1e9;
                JOINT_LOW_DEGREE_FIXTURE.s(:, locs(j, 1), locs(j, 2) + 2) = 1e8;
            end

            [again, againStd] = solvedegree1(args{:});
            testCase.verifyEqual(again, actl);
            testCase.verifyEqual(againStd, actlStd);
        end

        function varianceBlocksAndZeroRhs(testCase)
            ocnKer = ocean();
            [earthRespKer, forcKer] = prepareslekernels(3, ocnKer);
            v = zeros(16, 3); v([1, 7, 16], 3) = [1; 1e-20; 2];
            a = slevariance(ocnKer, earthRespKer, forcKer, v, BlockSize = 1);
            b = slevariance(ocnKer, earthRespKer, forcKer, v, BlockSize = 16);
            [A, KF] = explicitSystem(ocnKer, earthRespKer);
            testCase.verifyEqual(a, ((A \ KF) .^ 2) * v, 'AbsTol', 2e-12);
            testCase.verifyEqual(a, b, 'AbsTol', 2e-14);
            testCase.verifyEqual(solveslesystem(ocnKer, earthRespKer, zeros(16, 2)), zeros(16, 2));
        end

        function linearFailureContracts(testCase)
            ocnKer = ocean();
            [earthRespKer, forcKer] = prepareslekernels(3, ocnKer);
            testCase.verifyError(@() solveslesystem(ocnKer, earthRespKer, (1:16)', MaxIter = 1), ...
            'ULMO:solveslesystem:NonConvergence');
            testCase.verifyError(@() solveslesystem(ocnKer, eye(3), ones(16, 1)), ...
            'ULMO:solveslesystem:InvalidKernels');
            testCase.verifyError(@() slevariance(ocnKer, earthRespKer, eye(3), zeros(16, 1)), ...
            'ULMO:slevariance:InvalidSize');
            testCase.verifyError(@() preparelowdegreemap(ocnKer, earthRespKer, eye(3), 2:4), ...
            'ULMO:preparelowdegreemap:InvalidForcingKernel');
            testCase.verifyError(@() preparelowdegreemap(ocnKer, earthRespKer, forcKer, [2, 2]), ...
            'ULMO:preparelowdegreemap:InvalidIndices');
            badOcnKer = zeros(16); badOcnKer(1, 1) = 1;
            testCase.verifyError(@() preparelowdegreemap(badOcnKer, [], [], 2:4, Method = "uniform"), ...
            'ULMO:preparelowdegreemap:SingularSystem');
            f = unpack((1:16)');
            args = {f, [], 3, GeoDomain('alloceans'), 'OceanKernel', ocnKer, ...
                        'Tolerance', 0, 'MaxIter', 3, 'BeQuiet', true};
            testCase.verifyError(@() solvesle(args{:}), ...
            'ULMO:solvesle:InvalidTolerance');
            testCase.verifyError(@() solvedegree1('Tolerance', 0, 'BeQuiet', true), ...
            'ULMO:solvedegree1:InvalidTolerance');
        end

        function mapFailureIsReported(testCase)
            ocnKer = ocean();
            [earthRespKer, forcKer] = prepareslekernels(3, ocnKer);
            testCase.verifyError(@() preparelowdegreemap(ocnKer, earthRespKer, forcKer, 2:4, MaxIter = 1), ...
            'ULMO:solveslesystem:NonConvergence');
        end

        function originalCacheNamesAndExplicitRecomputation(testCase)
            args = {'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                        'Method', 'uniform', 'SaveData', true, 'BeQuiet', true};
            [legacy, legacyStd] = solvedegree1_iteration(args{:});
            files = dir(fullfile(testCase.Folder, 'Degree1', 'new', '*.mat'));
            testCase.verifyNumElements(files, 1);
            dom = GeoDomain('alloceans');
            testCase.verifyEqual(files.name, ...
                sprintf('CSRRL063_ICE6GD-Ld3_Ls3-%s-uniform-SD.mat', dom.Id));
            [loaded, loadedStd] = solvedegree1(args{:});
            testCase.verifyEqual(loaded, legacy);
            testCase.verifyEqual(loadedStd, legacyStd);
            [fresh, freshStd] = solvedegree1(args{:}, 'ForceNew', true, 'SaveData', false);
            stored = load(fullfile(files.folder, files.name));
            testCase.verifyEqual(stored.coeffStds, legacyStd);
            testCase.verifyGreaterThan(norm(freshStd - legacyStd), 1e-8);
            solvedegree1(args{:}, 'ForceNew', true);
            [loaded, loadedStd] = solvedegree1(args{:}, 'RelativeTolerance', 1e-10, 'Tolerance', 1e-5);
            testCase.verifyEqual(loaded, fresh);
            testCase.verifyEqual(loadedStd, freshStd);
            testCase.verifyNumElements(dir(fullfile(files.folder, '*.mat')), 1);
            testCase.verifyFalse(isfolder(fullfile(testCase.Folder, 'Degree1', 'maps')));
        end

    end

end

function O = ocean()
    i = 1:16;
    O = .6 * eye(16) + .003 * (sin(i' * i) + cos(i' + i));
end

function [A, KF] = explicitSystem(oceanKernel, earthResponseKernel)
    e = eye(size(oceanKernel, 1), 1);
    P = eye(size(oceanKernel)) - e * oceanKernel(1, :) / oceanKernel(1, 1);
    A = eye(size(oceanKernel)) - P * full(earthResponseKernel) * oceanKernel;
    KF = P * full(earthResponseKernel) - e * e.' / oceanKernel(1, 1);
end

function p = unpack(v)
    p = zeros(20, size(v, 2)); p(kernelorder(3), :) = v;
    p = reshape(p, 10, 2, []);
end

function v = pack(p)
    v = reshape(p, 20, []); v = v(kernelorder(3), :);
end

function clearFixture()
    clear global JOINT_LOW_DEGREE_FIXTURE
end

% Default fingerprint call with synthetic degree-3 input and output in SD.
function args = solverArgs()
    args = {'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), 'BeQuiet', true};
end
