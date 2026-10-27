%% LINEARSLESYSTEMTEST - Checks full and selected responses against dense equations.
% Uses isolated synthetic Love numbers, ocean kernels and GRACE/GIA readers.
% Tests induced covariance, transposed solves, joint coefficient sets, cache
% separation, physical mass balance and failure handling without production data.
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

classdef linearSleSystemTest < matlab.unittest.TestCase
    properties
        Folder
    end
    properties (TestParameter)
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
            loadinggravitationalpotential = [0; .021; -.3; -.2];
            loadingverticaldisplacement = [0; -.269; -.6; -.4];
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
            s = cat(2, meta, unpack(reshape(abs(cos(1:48))/10, 16, 3)));
            g = cat(2, meta, unpack(reshape(cos(1:48)/20, 16, 3)));
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3,1,2]), ...
                's', permute(s, [3,1,2]), 'gia', permute(g, [3,1,2]), ...
                'dates', datetime(2003, (1:3)', 15), 'kernel', ocean());
        end
    end
    methods (Test)
        function forwardAndTransposeMatchDense(testCase, rotation)
            k = prepareslekernels(3, ocean(), RotationFeedback=rotation);
            [A, KF] = explicitSystem(k);
            b = reshape(sin(1:48), 16, 3);
            [x, info] = solveslesystem(k, b);
            xt = solveslesystem(k, b, Transpose=true);
            testCase.verifyEqual(x, A\b, 'AbsTol', 2e-11);
            testCase.verifyEqual(xt, A.'\b, 'AbsTol', 2e-11);
            testCase.verifyEqual(full(k.KF), KF, 'AbsTol', 1e-14);
            testCase.verifyLessThanOrEqual(max(info.RelativeResidual), 1e-12);
            testCase.verifyTrue(issparse(k.response) && issparse(k.KF));
            testCase.verifyFalse(any(structfun(@(x) isa(x, 'function_handle'), k)));
            testCase.verifyGreaterThan(norm(k.O*k.O-k.O, 'fro'), 0.1);
        end
        function massAndFullCovariance(testCase, rotation)
            k = prepareslekernels(3, ocean(), RotationFeedback=rotation);
            [A, KF] = explicitSystem(k);
            fullSleResponse = A\KF;
            f = reshape(sin(1:48), 16, 3);
            sigma = reshape(abs(cos(1:48)), 16, 3)/10;
            sigma(:, 1) = 0;
            [s, ss, meanLevel, meanStd] = solvesle(unpack(f), unpack(sigma), ...
                3, GeoDomain('alloceans'), 'OceanKernel', k.O, ...
                'RotationFeedback', rotation, 'BeQuiet', true);
            testCase.verifyEqual(pack(s), fullSleResponse*f, 'AbsTol', 2e-11);
            testCase.verifyEqual(pack(ss).^2, (fullSleResponse.^2)*(sigma.^2), 'AbsTol', 2e-13);
            local = k.O*pack(s);
            testCase.verifyEqual(local(1, :), -f(1, :), 'AbsTol', 2e-11);
            testCase.verifyEqual(meanLevel, (-f(1,:)/k.oceanArea).', 'AbsTol', 2e-11);
            testCase.verifyEqual(meanStd, (sigma(1,:)/k.oceanArea).', 'AbsTol', 2e-13);
            for epoch = 1:3
                C = fullSleResponse*diag(sigma(:,epoch).^2)*fullSleResponse.';
                testCase.verifyEqual(pack(ss(:, :, epoch)).^2, diag(C), 'AbsTol', 2e-13);
            end
        end
        function selectedMapsAndPublicSolver(testCase, method, c20, c30, cs21)
            global JOINT_LOW_DEGREE_FIXTURE
            f = JOINT_LOW_DEGREE_FIXTURE;
            k = prepareslekernels(3, f.kernel);
            [A, KF] = explicitSystem(k);
            [locs, ids] = lowdegreecoefficients(c20, c30, cs21);
            if method == "fingerprint"
                expectedOceanResponse = f.kernel(ids,:)*(A\KF)*(eye(16)-f.kernel);
            else
                % The uniform response consumes only the land C00 row.
                expectedOceanResponse = -f.kernel(ids,1)/f.kernel(1,1)*( [1, zeros(1,15)]-f.kernel(1,:) );
            end
            [knownToUnknownMap, knownIds, oceanResponseMap, unknownConstraint, ...
                knownConstraint, reciprocalCondition, relativeResiduals, iterationCounts] = ...
                preparelowdegreemap(k, ids, Method=method);
            known = setdiff(1:16, ids);
            expectedLoadResponse = (expectedOceanResponse(:,ids)-f.kernel(ids,ids)) \ (f.kernel(ids,known)-expectedOceanResponse(:,known));
            testCase.verifyEqual(oceanResponseMap, expectedOceanResponse, 'AbsTol', 2e-11);
            testCase.verifyEqual(knownToUnknownMap, expectedLoadResponse, 'AbsTol', 2e-11);
            testCase.verifyEqual(knownIds, known);
            testCase.verifyEqual(unknownConstraint, ...
                expectedOceanResponse(:,ids)-f.kernel(ids,ids), 'AbsTol', 2e-11);
            testCase.verifyEqual(knownConstraint, ...
                f.kernel(ids,known)-expectedOceanResponse(:,known), 'AbsTol', 2e-11);
            testCase.verifyEqual(reciprocalCondition, rcond(unknownConstraint));
            testCase.verifySize(relativeResiduals, [1, numel(ids)]);
            testCase.verifySize(iterationCounts, [1, numel(ids)]);
            testCase.verifyLessThanOrEqual(max(relativeResiduals), 1e-12);
            if method == "uniform"
                testCase.verifyEqual(iterationCounts, zeros(1, numel(ids)));
            else
                testCase.verifyGreaterThan(iterationCounts, zeros(1, numel(ids)));
            end
            p = pack(permute(f.p(:,:,3:4), [2,3,1]));
            sigma = pack(permute(f.s(:,:,3:4), [2,3,1]));
            gia = pack(permute(f.gia(:,:,3:4), [2,3,1]));
            p(5:end,:) = p(5:end,:)-gia(5:end,:);
            expected = expectedLoadResponse*p(known,:);
            expectedStd = sqrt((expectedLoadResponse.^2)*(sigma(known,:).^2));
            if numel(ids) > 3
                expected(4:end,:) = expected(4:end,:)+gia(ids(4:end),:);
            end
            opts = {'IncludeC20', c20, 'IncludeC30', c30, 'IncludeC21S21', cs21, ...
                'Method', char(method), 'ForceNew', true, 'SaveData', false, 'BeQuiet', true};
            [actual, actualStd, dates, ~, actualLocs] = solvedegree1( ...
                'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), opts{:});
            testCase.verifyEqual(actual, expected.', 'AbsTol', 2e-11);
            testCase.verifyEqual(actualStd, expectedStd.', 'AbsTol', 2e-11);
            testCase.verifyEqual(dates, f.dates);
            testCase.verifyEqual(actualLocs, locs);
            % Independently check the central coefficients against converged iteration.
            legacy = solvedegree1('CSR', 'RL06', 3, 3, 3, 'ICE6GD', ...
                GeoDomain('alloceans'), opts{:}, 'Solver', 'iteration', ...
                'Tolerance', 1e-12, 'SleTolerance', 1e-14);
            testCase.verifyEqual(actual, legacy, 'AbsTol', 2e-9);
        end
        function paddedInputsIgnoreImportedUnknowns(testCase, method)
            global JOINT_LOW_DEGREE_FIXTURE
            f = JOINT_LOW_DEGREE_FIXTURE;
            [~, ids] = lowdegreecoefficients(false, true, true);
            k = prepareslekernels(3, f.kernel);
            [knownToUnknownMap, knownIds] = preparelowdegreemap(k, ids, Method=method);
            x = pack(permute(f.p(:,:,3:4), [2,3,1]));
            sigma = pack(permute(f.s(:,:,3:4), [2,3,1]));
            gia = pack(permute(f.gia(:,:,3:4), [2,3,1]));
            x(10:end,:) = 0; sigma(10:end,:) = 0;
            x(5:end,:) = x(5:end,:)-gia(5:end,:);
            expected = knownToUnknownMap*x(knownIds,:);
            expected(4:end,:) = expected(4:end,:)+gia(ids(4:end),:);
            args = {'CSR','RL06',3,2,3,'ICE6GD',GeoDomain('alloceans'), ...
                'IncludeC30',true,'IncludeC21S21',true,'Method',char(method), ...
                'SaveData',false,'BeQuiet',true};
            [actual, actualStd] = solvedegree1(args{:});
            testCase.verifyEqual(actual,expected.','AbsTol',2e-11);
            testCase.verifyEqual(actualStd,sqrt((knownToUnknownMap.^2)*(sigma(knownIds,:).^2)).','AbsTol',2e-11);
            [locs, ~] = lowdegreecoefficients(false,true,true);
            for j = 1:size(locs,1)
                JOINT_LOW_DEGREE_FIXTURE.p(:,locs(j,1),locs(j,2)+2) = 1e9;
                JOINT_LOW_DEGREE_FIXTURE.s(:,locs(j,1),locs(j,2)+2) = 1e8;
            end
            [again, againStd] = solvedegree1(args{:});
            testCase.verifyEqual(again,actual);
            testCase.verifyEqual(againStd,actualStd);
        end
        function varianceBlocksAndZeroRhs(testCase)
            k = prepareslekernels(3, ocean());
            v = zeros(16,3); v([1,7,16],3) = [1; 1e-20; 2];
            a = slevariance(k, v, BlockSize=1);
            b = slevariance(k, v, BlockSize=16);
            [A, KF] = explicitSystem(k);
            testCase.verifyEqual(a, ((A\KF).^2)*v, 'AbsTol', 2e-12);
            testCase.verifyEqual(a,b,'AbsTol',2e-14);
            testCase.verifyEqual(solveslesystem(k, zeros(16,2)), zeros(16,2));
        end
        function failureAndFixedStepContracts(testCase)
            k = prepareslekernels(3, ocean());
            testCase.verifyError(@() solveslesystem(k, (1:16)', MaxIter=1), ...
                'ULMO:solveslesystem:NonConvergence');
            testCase.verifyError(@() preparelowdegreemap(k, [2,2]), ...
                'ULMO:preparelowdegreemap:InvalidIndices');
            bad = k; bad.O(:) = 0; bad.O(1,1) = 1;
            testCase.verifyError(@() preparelowdegreemap(bad, 2:4, Method="uniform"), ...
                'ULMO:preparelowdegreemap:SingularSystem');
            f = unpack((1:16)');
            args = {f, [], 3, GeoDomain('alloceans'), 'OceanKernel', k.O, ...
                'Tolerance', 0, 'MaxIter', 3, 'BeQuiet', true};
            testCase.verifyEqual(solvesle(args{:}), solvesle(args{:}, 'Solver', 'iteration'));
            testCase.verifyError(@() solvesle(args{:}, 'Solver', 'linear'), ...
                'ULMO:solvesle:FixedStepsRequireIteration');
        end
        function mapFailureIsReported(testCase)
            k = prepareslekernels(3,ocean());
            testCase.verifyError(@() preparelowdegreemap(k,2:4,MaxIter=1), ...
                'ULMO:solveslesystem:NonConvergence');
        end
        function linearAndLegacyCachesAreSeparate(testCase)
            args = {'CSR','RL06',3,3,3,'ICE6GD',GeoDomain('alloceans'), ...
                'Method','uniform','SaveData',true,'BeQuiet',true};
            [a, sa] = solvedegree1(args{:});
            solvedegree1(args{:}, 'Solver','iteration');
            files = dir(fullfile(testCase.Folder,'Degree1','new','*.mat'));
            testCase.verifyNumElements(files,2);
            testCase.verifyEqual(sum(contains({files.name},'-linear-v1-')),1);
            [b,sb] = solvedegree1(args{:});
            testCase.verifyEqual(a,b); testCase.verifyEqual(sa,sb);
            testCase.verifyFalse(isfolder(fullfile(testCase.Folder,'Degree1','maps')));
        end
    end
end
function O = ocean()
    i = 1:16;
    O = .6*eye(16) + .003*(sin(i'*i)+cos(i'+i));
end
function [A,KF] = explicitSystem(k)
    e = eye(size(k.O,1),1);
    P = eye(size(k.O))-e*k.O(1,:)/k.O(1,1);
    A = eye(size(k.O))-P*full(k.response)*k.O;
    KF = P*full(k.response)-e*e.'/k.O(1,1);
end
function p = unpack(v)
    p = zeros(20,size(v,2)); p(kernelorder(3),:) = v;
    p = reshape(p,10,2,[]);
end
function v = pack(p)
    v = reshape(p,20,[]); v = v(kernelorder(3),:);
end
function clearFixture()
    clear global JOINT_LOW_DEGREE_FIXTURE
end
