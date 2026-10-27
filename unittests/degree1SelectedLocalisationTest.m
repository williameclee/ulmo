%% DEGREE1SELECTEDLOCALISATIONTEST - Compare solver outputs with full localization.
% Uses the actual SLE and isolated reader/kernel fixtures. The oracle
% retains full-field localisation before selecting coefficients. All optional
% coefficient combinations, both methods, and fixed/adaptive outer loops are
% checked, including uncertainties, GIA restoration and output ordering.
% Dense signed kernels and zero uncertainties also exercise reusable land
% forcing against full localisation at every reference iteration.
% Uniform cases retain unestimated high-degree contributions to land C00,
% including GIA above the input truncation degree and nonzero input C00.
% Fingerprint cases at degrees 15, 16 and 23 check rectangular variance
% products and both signs of the diagonal correction against full fields.
% Padded-input cases check exact variance support across epochs, empty
% support, tiny nonzero uncertainties, and higher-degree GIA signals.
% Independent signal/variance supports exercise empty, half-full and dense
% projections, including support present only after the first epoch.
%
% Created by
%   2026/10/26, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

classdef degree1SelectedLocalisationTest < matlab.unittest.TestCase

    properties (TestParameter)
        inclC20 = {false, true}
        inclC30 = {false, true}
        inclCS21 = {false, true}
        mthd = {'fingerprint', 'uniform'}
        epochCnt = {1, 3}
        adaptive = {false, true}
        kernelCase = {'original', 'dense', 'denseZeroStd'}
        inputDegree = {2, 3}
        fingerprintDegree = {15, 16, 23}
        oceanDiagonal = {0.25, 0.6, 0.95}
        varianceSupport = {'empty', 'mixed', 'tiny'}
        signalCount = {0, 8, 9}
        varianceCount = {0, 8, 9}
    end

    methods (TestClassSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(repo));
            testCase.applyFixture(PathFixture(fullfile(repo, 'aux')));
            testCase.applyFixture(PathFixture(fullfile(repo, 'domains')));
            folder = testCase.applyFixture(TemporaryFolderFixture);
            fixtures = fullfile(repo, 'unittests', 'fixtures', 'jointlowdegreesolver');
            files = dir(fullfile(fixtures, '*.txt'));

            for k = 1:numel(files)
                copyfile(fullfile(fixtures, files(k).name), ...
                    fullfile(folder.Folder, erase(files(k).name, '.txt')));
            end

            testCase.applyFixture(PathFixture(folder.Folder));

            for name = {'IFILES', 'GRACEDATA'}
                previous = getenv(name{1});
                testCase.addTeardown(@() setenv(name{1}, previous));
                setenv(name{1}, folder.Folder);
            end

            testCase.addTeardown(@() clearFixture());
            mkdir(fullfile(folder.Folder, 'EARTHMODELS', 'CONSTANTS'));
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(folder.Folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
            mkdir(fullfile(folder.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; .021; - .3; repmat(- .2, 21, 1)];
            loadingverticaldisplacement = [0; - .269; - .6; repmat(- .4, 21, 1)];
            tidalgravitationalpotential = [0; 0; .3; repmat(.1, 21, 1)];
            tidalverticaldisplacement = [0; 0; .6; repmat(.2, 21, 1)];
            save(fullfile(folder.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential', 'loadingverticaldisplacement', ...
                'tidalgravitationalpotential', 'tidalverticaldisplacement');
            testCase.assertEqual(which('solvedegree1'), fullfile(repo, 'solvedegree1.m'));
            testCase.assertEqual(which('solvesle'), fullfile(repo, 'solvesle.m'));
        end

    end

    methods (Test)

        function projectionSupportThresholdMatchesReference(testCase, mthd, epochCnt, signalCount, varianceCount)
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, epochCnt);
            ko = kernelorder(3);
            available = [16, 1, 9, 5, 12, 6, 15, 7, 13, 8, 10, 11, 14];
            epochs = 1:epochCnt;

            if epochCnt > 1
                epochs(1) = 0;
            end

            values = zeros(20, epochCnt);
            sigma = zeros(20, epochCnt);
            values(ko(available(1:signalCount)), :) = (1:signalCount)' * epochs;
            sigma(ko(available(1:varianceCount)), :) = (1:varianceCount)' * epochs / 10;
            giaValues = zeros(10, 2, epochCnt);
            giaValues(10, 2, :) = 0.125;
            p = cat(2, meta, reshape(values, 10, 2, epochCnt) + giaValues);
            s = cat(2, meta, reshape(sigma, 10, 2, epochCnt));
            gia = cat(2, meta, giaValues);
            dates = datetime(2003, (1:epochCnt)', 15);
            indices = 1:16;
            K = 0.6 * eye(16) + 0.002 * (sin(indices' * indices) + cos(indices' + indices));
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), ...
                'dates', dates, 'kernel', K);
            [expc, expcStd] = fullFieldReference(p, s, gia, ...
                [2, 1; 3, 1; 3, 2], [2, 3, 4], mthd, 0, 3);
            [actl, actlStd, actlDates] = solvedegree1legacy( ...
                'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                'Method', mthd, 'Tolerance', 0, 'MaxIter', 3, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actl, expc, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlStd, expcStd, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlDates, dates);

            if epochCnt > 1
                testCase.verifyEqual(actl(1, :), zeros(1, 3));
                testCase.verifyEqual(actlStd(1, :), zeros(1, 3));
            end

        end

        function paddedVariancesPreserveFullSignal(testCase, mthd, epochCnt, varianceSupport)
            global JOINT_LOW_DEGREE_FIXTURE
            L = 16;
            [order, degree] = addmon(L);
            meta = repmat([degree, order], 1, 1, epochCnt);
            values = zeros(numel(degree), 2, epochCnt);
            values(1:10, :, :) = reshape(sin(1:20 * epochCnt), 10, 2, epochCnt);
            sigma = zeros(size(values));

            if ~strcmp(varianceSupport, 'empty')
                scale = 0.2;

                if strcmp(varianceSupport, 'tiny')
                    scale = 1e-100;
                end

                sigma(1, 1, :) = scale;
                sigma(6, 2, :) = 2 * scale; % S22 is not jointly estimated.

                if epochCnt > 1
                    % The first epoch alone must not determine retained columns.
                    sigma(:, :, 1) = 0;
                    sigma(6, 2, 2) = 0;
                    sigma(1, 1, 3) = 0;
                end

            end

            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            giaValues = zeros(size(values));
            giaValues(end, 2, :) = reshape(1:epochCnt, 1, 1, []);
            gia = cat(2, meta, giaValues);
            dates = datetime(2003, (1:epochCnt)', 15);
            indices = 1:(L + 1) ^ 2;
            K = 0.6 * eye(numel(indices)) ...
                + 0.001 * (sin(indices' * indices) + cos(indices' + indices));
            K(1, end) = -0.04; K(end, 1) = K(1, end);
            % The reader supplies degree 3; estimation truncates it to 2,
            % then pads it to 16 before subtracting the degree-16 GIA field.
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p(1:10, :, :), [3, 1, 2]), ...
                's', permute(s(1:10, :, :), [3, 1, 2]), ...
                'gia', permute(gia, [3, 1, 2]), 'dates', dates, 'kernel', K);
            p(degree > 2, 3:4, :) = 0;
            s(degree > 2, 3:4, :) = 0;
            locs = [2, 1; 3, 1; 3, 2; 4, 1; 7, 1; 5, 1; 5, 2];
            ids = [2, 3, 4, 5, 10, 6, 7];
            [expc, expcStd] = fullFieldReference(p, s, gia, ...
                locs, ids, mthd, 0, 4);
            [actl, actlStd, actlDates] = solvedegree1legacy( ...
                'CSR', 'RL06', 3, 2, L, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', true, 'IncludeC30', true, 'IncludeC21S21', true, ...
                'Method', mthd, 'Tolerance', 0, 'MaxIter', 4, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actl, expc, 'AbsTol', 1e-10);
            % Relative tolerance detects accidental removal of tiny variances.
            testCase.verifyEqual(actlStd, expcStd, 'RelTol', 1e-11);
            testCase.verifyEqual(actlDates, dates);

            if strcmp(varianceSupport, 'empty')
                testCase.verifyEqual(actlStd, zeros(epochCnt, 7));
            else
                testCase.verifyGreaterThan(max(actlStd, [], 'all'), 0);
            end

            if epochCnt > 1
                testCase.verifyEqual(actlStd(1, :), zeros(1, 7));
            end

        end

        function fingerprintMatchesAcrossBandwidths(testCase, fingerprintDegree, epochCnt, adaptive, oceanDiagonal)
            global JOINT_LOW_DEGREE_FIXTURE
            L = fingerprintDegree;
            [order, degree] = addmon(L);
            meta = repmat([degree, order], 1, 1, epochCnt);
            values = reshape(sin(1:2 * numel(degree) * epochCnt), ...
                numel(degree), 2, epochCnt);
            sigma = 0.1 + abs(values) / 5;

            if epochCnt > 1
                sigma(:, :, end) = 0;
            end

            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            gia = cat(2, meta, values / 7);
            dates = datetime(2003, (1:epochCnt)', 15);
            idxs = 1:(L + 1) ^ 2;
            K = oceanDiagonal * eye(numel(idxs)) ...
                + 0.001 * (sin(idxs' * idxs) + cos(idxs' + idxs));
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), ...
                'dates', dates, 'kernel', K);
            locs = [2, 1; 3, 1; 3, 2; 4, 1; 7, 1; 5, 1; 5, 2];
            ids = [2, 3, 4, 5, 10, 6, 7];

            if adaptive
                tol = 1e-6;
                lim = 100;
            else
                tol = 0;
                lim = 12;
            end

            [expc, expcStd] = fullFieldReference(p, s, gia, ...
                locs, ids, 'fingerprint', tol, lim);
            [actl, actlStd, actlDates] = solvedegree1legacy( ...
                'CSR', 'RL06', L, L, L, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', true, 'IncludeC30', true, 'IncludeC21S21', true, ...
                'Method', 'fingerprint', 'Tolerance', tol, 'MaxIter', lim, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actl, expc, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlStd, expcStd, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlDates, dates);
            testCase.verifyGreaterThanOrEqual(actlStd, zeros(size(actlStd)));

            if epochCnt > 1
                testCase.verifyEqual(actlStd(end, :), zeros(1, 7));
            end

        end

        function uniformRetainsUnestimatedHarmonics(testCase, epochCnt, adaptive, inputDegree)
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, epochCnt);
            values = zeros(10, 2, epochCnt);
            sigma = values;
            giaValues = values;
            % C00 and S33 are not among the seven jointly estimated terms.
            values(1, 1, :) = reshape(1:epochCnt, 1, 1, []);
            values(10, 2, :) = reshape(-3 * (1:epochCnt), 1, 1, []);
            sigma(1, 1, :) = 0.2;
            sigma(10, 2, :) = 0.7;
            giaValues(10, 2, :) = reshape(0.4 * (1:epochCnt), 1, 1, []);

            if epochCnt > 1
                sigma(:, :, end) = 0;
            end

            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            gia = cat(2, meta, giaValues);
            dates = datetime(2003, (1:epochCnt)', 15);
            K = 0.6 * eye(16);
            K(1, 2:16) = 0.008 * (-1) .^ (2:16);
            K(1, 16) = -0.08;
            K(2:16, 1) = K(1, 2:16)';
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), ...
                'dates', dates, 'kernel', K);
            locs = [2, 1; 3, 1; 3, 2; 4, 1; 7, 1; 5, 1; 5, 2];
            ids = [2, 3, 4, 5, 10, 6, 7];
            refP = p;
            refS = s;
            refP(degree > inputDegree, 3:4, :) = 0;
            refS(degree > inputDegree, 3:4, :) = 0;

            if adaptive
                tol = 1e-6;
                lim = 100;
            else
                tol = 0;
                lim = 12;
            end

            [expc, expcStd] = fullFieldReference(refP, refS, ...
                gia, locs, ids, 'uniform', tol, lim);
            [actl, actlStd, actlDates] = solvedegree1legacy( ...
                'CSR', 'RL06', 3, inputDegree, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', true, 'IncludeC30', true, 'IncludeC21S21', true, ...
                'Method', 'uniform', 'Tolerance', tol, 'MaxIter', lim, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actl, expc, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlStd, expcStd, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlDates, dates);
            testCase.verifyGreaterThan(max(abs(actl), [], 'all'), 0);

            if epochCnt > 1
                testCase.verifyEqual(actlStd(end, :), zeros(1, 7));
            end

        end

        function matchesFullFieldReference(testCase, inclC20, inclC30, inclCS21, mthd, epochCnt, adaptive, kernelCase)
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, epochCnt);
            values = reshape(sin(1:20 * epochCnt), 10, 2, epochCnt);
            sigma = reshape(1:20 * epochCnt, 10, 2, epochCnt) / 100;

            if strcmp(kernelCase, 'denseZeroStd')
                sigma(:) = 0;
            elseif strcmp(kernelCase, 'dense') && epochCnt > 1
                % A zero-variance epoch must not inherit another epoch's errors.
                sigma(:, :, 1) = 0;
            end

            gia = cat(2, meta, values / 7);
            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            dates = datetime(2003, (1:epochCnt)', 15);
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), 'dates', dates);

            if ~strcmp(kernelCase, 'original')
                idxs = 1:16;
                JOINT_LOW_DEGREE_FIXTURE.kernel = .6 * eye(16) ...
                    + .002 * (sin(idxs' * idxs) + cos(idxs' + idxs));
            end

            locs = [2, 1; 3, 1; 3, 2];
            ids = [2, 3, 4];

            if inclC20
                locs(end + 1, :) = [4, 1];
                ids(end + 1) = 5;
            end

            if inclC30
                locs(end + 1, :) = [7, 1];
                ids(end + 1) = 10;
            end

            if inclCS21
                locs = [locs; 5, 1; 5, 2];
                ids = [ids, 6, 7];
            end

            if adaptive
                tol = 1e-6;
                lim = 100;
            else
                tol = 0;
                lim = 3;

                if ~strcmp(kernelCase, 'original')
                    lim = 12;
                end

            end

            [expc, expcStd] = fullFieldReference(p, s, gia, locs, ids, mthd, tol, lim);
            [actl, actlStd, actlDates, field, actlLocs] = solvedegree1legacy( ...
                'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', inclC20, 'IncludeC30', inclC30, ...
                'IncludeC21S21', inclCS21, 'Method', mthd, ...
                'Tolerance', tol, 'MaxIter', lim, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actl, expc, 'AbsTol', 1e-10);
            testCase.verifyEqual(actlStd, expcStd, 'AbsTol', 1e-10);
            testCase.verifyTrue(all(isfinite(actl), 'all'));
            testCase.verifyTrue(all(isfinite(actlStd), 'all'));
            testCase.verifyGreaterThanOrEqual(actlStd, zeros(size(actlStd)));
            testCase.verifyEqual(actlDates, dates);
            testCase.verifyEqual(actlLocs, locs);

            for k = 1:size(locs, 1)
                testCase.verifyEqual(field(:, locs(k, 1), 2 + locs(k, 2)), expc(:, k), 'AbsTol', 1e-10);
            end

        end

    end

end

function [coeffs, stds] = fullFieldReference(p, s, gia, locs, ids, mthd, tol, lim)
    % Keep the original full-localisation arithmetic as the comparison oracle.
    L = max(p(:, 1, 1));
    p(4:end, 3:4, :) = p(4:end, 3:4, :) - gia(4:end, 3:4, :);

    for k = 1:size(locs, 1)
        p(locs(k, 1), 2 + locs(k, 2), :) = 0;
        s(locs(k, 1), 2 + locs(k, 2), :) = 0;
    end

    K = kernelcp_new(L); land = eye((L + 1) ^ 2) - K; block = K(ids, ids);
    known = select(localise(p, 'L', L, 'K', K), locs);
    knownStd = select(localise(s, 'L', L, 'K', K, 'IsError', true), locs);
    coeffs = zeros(size(known)); stds = coeffs;

    for iter = 1:lim
        old = coeffs;
        oldStd = stds;
        f = localise(p, 'L', L, 'K', land);
        fs = localise(s, 'L', L, 'K', land, 'IsError', true);

        if strcmp(mthd, 'fingerprint')
            [r, rs] = solveslelegacy(f, fs, L, GeoDomain('alloceans'), ...
                'OceanKernel', K, 'KernelOrder', kernelorder(L), ...
                'InitialRSL', p / 1000, 'RotationFeedback', true, ...
                'Tolerance', 1e-9, 'MaxIter', 100, 'BeQuiet', true);
            ocn = select(localise(r, 'L', L, 'K', K), locs);
            ocnStd = select(localise(rs, 'L', L, 'K', K, 'IsError', true), locs);
        else
            ocn = K(ids, 1) * reshape(-f(1, 3, :) / K(1, 1), 1, []);
            ocnStd = K(ids, 1) * reshape(-fs(1, 3, :) / K(1, 1), 1, []);
        end

        coeffs = block \ (ocn - known);
        stds = sqrt((block \ eye(numel(ids))) .^ 2 * (ocnStd .^ 2 + knownStd .^ 2));

        for k = 1:size(locs, 1)
            p(locs(k, 1), 2 + locs(k, 2), :) = coeffs(k, :);
            s(locs(k, 1), 2 + locs(k, 2), :) = stds(k, :);
        end

        change = max([abs(coeffs - old); abs(stds - oldStd)], [], 'all');

        if tol > 0 && change <= tol
            break
        end

    end

    assert(tol == 0 || change <= tol, 'Full-field reference failed to converge.');

    for k = 4:size(locs, 1)
        coeffs(k, :) = coeffs(k, :) + reshape(gia(locs(k, 1), 2 + locs(k, 2), :), 1, []);
    end

    coeffs = coeffs'; stds = stds';
end

function values = select(p, locs)
    values = zeros(size(locs, 1), size(p, 3));

    for k = 1:size(locs, 1)
        values(k, :) = reshape(p(locs(k, 1), 2 + locs(k, 2), :), 1, []);
    end

end

function clearFixture()
    clear global JOINT_LOW_DEGREE_FIXTURE
end
