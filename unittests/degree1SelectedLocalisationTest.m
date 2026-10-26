%% DEGREE1SELECTEDLOCALISATIONTEST - Compare solver outputs with full localization.
% Uses the actual SLE and isolated reader/kernel fixtures. The oracle
% retains full-field localisation before selecting coefficients. All optional
% coefficient combinations, both methods, and fixed/adaptive outer loops are
% checked, including uncertainties, GIA restoration and output ordering.
% Dense signed kernels and zero uncertainties also exercise reusable land
% forcing against full localisation at every reference iteration.
% Uniform cases retain unestimated high-degree contributions to land C00,
% including GIA above the input truncation degree and nonzero input C00.
% Fingerprint cases at degrees 15, 16 and 23 exercise full and partial
% 256-row variance blocks against the original full-field calculation.
%
% Created by
%   2026/10/26, OpenAI Codex

classdef degree1SelectedLocalisationTest < matlab.unittest.TestCase

    properties (TestParameter)
        includeC20 = {false, true}
        includeC30 = {false, true}
        includeC21S21 = {false, true}
        method = {'fingerprint', 'uniform'}
        epochCount = {1, 3}
        adaptive = {false, true}
        kernelCase = {'original', 'dense', 'denseZeroStd'}
        inputDegree = {2, 3}
        fingerprintDegree = {15, 16, 23}
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
            loadinggravitationalpotential = [0; .021; - .3; repmat(-.2, 21, 1)];
            loadingverticaldisplacement = [0; - .269; - .6; repmat(-.4, 21, 1)];
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

        function fingerprintMatchesAcrossVarianceBlocks(testCase, fingerprintDegree, epochCount, adaptive)
            global JOINT_LOW_DEGREE_FIXTURE
            L = fingerprintDegree;
            [order, degree] = addmon(L);
            meta = repmat([degree, order], 1, 1, epochCount);
            values = reshape(sin(1:2 * numel(degree) * epochCount), ...
                numel(degree), 2, epochCount);
            sigma = 0.1 + abs(values) / 5;

            if epochCount > 1
                sigma(:, :, end) = 0;
            end

            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            gia = cat(2, meta, values / 7);
            dates = datetime(2003, (1:epochCount)', 15);
            indices = 1:(L + 1)^2;
            K = 0.6 * eye(numel(indices)) ...
                + 0.001 * (sin(indices' * indices) + cos(indices' + indices));
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), ...
                'dates', dates, 'kernel', K);
            locs = [2, 1; 3, 1; 3, 2; 4, 1; 7, 1; 5, 1; 5, 2];
            ids = [2, 3, 4, 5, 10, 6, 7];

            if adaptive
                tolerance = 1e-6; limit = 100;
            else
                tolerance = 0; limit = 12;
            end

            [expected, expectedStd] = fullFieldReference(p, s, gia, ...
                locs, ids, 'fingerprint', tolerance, limit);
            [actual, actualStd, actualDates] = solvedegree1( ...
                'CSR', 'RL06', L, L, L, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', true, 'IncludeC30', true, 'IncludeC21S21', true, ...
                'Method', 'fingerprint', 'Tolerance', tolerance, 'MaxIter', limit, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
            testCase.verifyEqual(actualStd, expectedStd, 'AbsTol', 1e-10);
            testCase.verifyEqual(actualDates, dates);
            testCase.verifyGreaterThanOrEqual(actualStd, zeros(size(actualStd)));

            if epochCount > 1
                testCase.verifyEqual(actualStd(end, :), zeros(1, 7));
            end

        end

        function uniformRetainsUnestimatedHarmonics(testCase, epochCount, adaptive, inputDegree)
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, epochCount);
            values = zeros(10, 2, epochCount);
            sigma = values;
            giaValues = values;
            % C00 and S33 are not among the seven jointly estimated terms.
            values(1, 1, :) = reshape(1:epochCount, 1, 1, []);
            values(10, 2, :) = reshape(-3 * (1:epochCount), 1, 1, []);
            sigma(1, 1, :) = 0.2;
            sigma(10, 2, :) = 0.7;
            giaValues(10, 2, :) = reshape(0.4 * (1:epochCount), 1, 1, []);

            if epochCount > 1
                sigma(:, :, end) = 0;
            end

            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            gia = cat(2, meta, giaValues);
            dates = datetime(2003, (1:epochCount)', 15);
            K = 0.6 * eye(16);
            K(1, 2:16) = 0.008 * (-1) .^ (2:16);
            K(1, 16) = -0.08;
            K(2:16, 1) = K(1, 2:16)';
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), ...
                'dates', dates, 'kernel', K);
            locs = [2, 1; 3, 1; 3, 2; 4, 1; 7, 1; 5, 1; 5, 2];
            ids = [2, 3, 4, 5, 10, 6, 7];
            referenceP = p;
            referenceS = s;
            referenceP(degree > inputDegree, 3:4, :) = 0;
            referenceS(degree > inputDegree, 3:4, :) = 0;

            if adaptive
                tolerance = 1e-6; limit = 100;
            else
                tolerance = 0; limit = 12;
            end

            [expected, expectedStd] = fullFieldReference(referenceP, referenceS, ...
                gia, locs, ids, 'uniform', tolerance, limit);
            [actual, actualStd, actualDates] = solvedegree1( ...
                'CSR', 'RL06', 3, inputDegree, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', true, 'IncludeC30', true, 'IncludeC21S21', true, ...
                'Method', 'uniform', 'Tolerance', tolerance, 'MaxIter', limit, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
            testCase.verifyEqual(actualStd, expectedStd, 'AbsTol', 1e-10);
            testCase.verifyEqual(actualDates, dates);
            testCase.verifyGreaterThan(max(abs(actual), [], 'all'), 0);

            if epochCount > 1
                testCase.verifyEqual(actualStd(end, :), zeros(1, 7));
            end

        end

        function matchesFullFieldReference(testCase, includeC20, includeC30, includeC21S21, method, epochCount, adaptive, kernelCase)
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, epochCount);
            values = reshape(sin(1:20 * epochCount), 10, 2, epochCount);
            sigma = reshape(1:20 * epochCount, 10, 2, epochCount) / 100;

            if strcmp(kernelCase, 'denseZeroStd')
                sigma(:) = 0;
            elseif strcmp(kernelCase, 'dense') && epochCount > 1
                % A zero-variance epoch must not inherit another epoch's errors.
                sigma(:, :, 1) = 0;
            end

            gia = cat(2, meta, values / 7);
            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            dates = datetime(2003, (1:epochCount)', 15);
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), 'dates', dates);

            if ~strcmp(kernelCase, 'original')
                indices = 1:16;
                JOINT_LOW_DEGREE_FIXTURE.kernel = .6 * eye(16) ...
                    + .002 * (sin(indices' * indices) + cos(indices' + indices));
            end

            locs = [2, 1; 3, 1; 3, 2]; ids = [2, 3, 4];

            if includeC20
                locs(end + 1, :) = [4, 1]; ids(end + 1) = 5;
            end

            if includeC30
                locs(end + 1, :) = [7, 1]; ids(end + 1) = 10;
            end

            if includeC21S21
                locs = [locs; 5, 1; 5, 2]; ids = [ids, 6, 7];
            end

            if adaptive
                tolerance = 1e-6; limit = 100;
            else
                tolerance = 0; limit = 3;

                if ~strcmp(kernelCase, 'original')
                    limit = 12;
                end

            end

            [expected, expectedStd] = fullFieldReference(p, s, gia, locs, ids, method, tolerance, limit);
            [actual, actualStd, actualDates, field, actualLocs] = solvedegree1( ...
                'CSR', 'RL06', 3, 3, 3, 'ICE6GD', GeoDomain('alloceans'), ...
                'IncludeC20', includeC20, 'IncludeC30', includeC30, ...
                'IncludeC21S21', includeC21S21, 'Method', method, ...
                'Tolerance', tolerance, 'MaxIter', limit, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
            testCase.verifyEqual(actualStd, expectedStd, 'AbsTol', 1e-10);
            testCase.verifyTrue(all(isfinite(actual), 'all'));
            testCase.verifyTrue(all(isfinite(actualStd), 'all'));
            testCase.verifyGreaterThanOrEqual(actualStd, zeros(size(actualStd)));
            testCase.verifyEqual(actualDates, dates);
            testCase.verifyEqual(actualLocs, locs);

            for k = 1:size(locs, 1)
                testCase.verifyEqual(field(:, locs(k, 1), 2 + locs(k, 2)), expected(:, k), 'AbsTol', 1e-10);
            end

        end

    end

end

function [coeffs, stds] = fullFieldReference(p, s, gia, locs, ids, method, tolerance, limit)
    % Keep the original full-localisation arithmetic as the comparison oracle.
    L = max(p(:, 1, 1));
    p(4:end, 3:4, :) = p(4:end, 3:4, :) - gia(4:end, 3:4, :);

    for k = 1:size(locs, 1)
        p(locs(k, 1), 2 + locs(k, 2), :) = 0;
        s(locs(k, 1), 2 + locs(k, 2), :) = 0;
    end

    K = kernelcp_new(L); land = eye((L + 1)^2) - K; block = K(ids, ids);
    known = select(localise(p, 'L', L, 'K', K), locs);
    knownStd = select(localise(s, 'L', L, 'K', K, 'IsError', true), locs);
    coeffs = zeros(size(known)); stds = coeffs;

    for iteration = 1:limit
        old = coeffs; oldStd = stds;
        f = localise(p, 'L', L, 'K', land);
        fs = localise(s, 'L', L, 'K', land, 'IsError', true);

        if strcmp(method, 'fingerprint')
            [r, rs] = solvesle(f, fs, L, GeoDomain('alloceans'), ...
                'OceanKernel', K, 'KernelOrder', kernelorder(L), ...
                'InitialRSL', p / 1000, 'RotationFeedback', true, ...
                'Tolerance', 1e-9, 'MaxIter', 100, 'BeQuiet', true);
            ocean = select(localise(r, 'L', L, 'K', K), locs);
            oceanStd = select(localise(rs, 'L', L, 'K', K, 'IsError', true), locs);
        else
            ocean = K(ids, 1) * reshape(-f(1, 3, :) / K(1, 1), 1, []);
            oceanStd = K(ids, 1) * reshape(-fs(1, 3, :) / K(1, 1), 1, []);
        end

        coeffs = block \ (ocean - known);
        stds = sqrt((block \ eye(numel(ids))) .^ 2 * (oceanStd .^ 2 + knownStd .^ 2));

        for k = 1:size(locs, 1)
            p(locs(k, 1), 2 + locs(k, 2), :) = coeffs(k, :);
            s(locs(k, 1), 2 + locs(k, 2), :) = stds(k, :);
        end

        change = max([abs(coeffs - old); abs(stds - oldStd)], [], 'all');

        if tolerance > 0 && change <= tolerance
            break
        end

    end

    assert(tolerance == 0 || change <= tolerance, 'Full-field reference failed to converge.');

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
