%% DEGREE1SELECTEDLOCALIZATIONTEST - Compare solver outputs with full localization.
% Uses the actual SLE and isolated degree-3 reader/kernel fixtures. The oracle
% retains full-field localization before selecting coefficients. All optional
% coefficient combinations, both methods, and fixed/adaptive outer loops are
% checked, including uncertainties, GIA restoration and output ordering.
%
% Created by
%   2026/10/26, OpenAI Codex

classdef degree1SelectedLocalizationTest < matlab.unittest.TestCase
    properties (TestParameter)
        includeC20 = {false, true}
        includeC30 = {false, true}
        includeC21S21 = {false, true}
        method = {'fingerprint', 'uniform'}
        epochCount = {1, 3}
        adaptive = {false, true}
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
            loadinggravitationalpotential = [0; .021; -.3; -.2];
            loadingverticaldisplacement = [0; -.269; -.6; -.4];
            tidalgravitationalpotential = [0; 0; .3; .1];
            tidalverticaldisplacement = [0; 0; .6; .2];
            save(fullfile(folder.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential', 'loadingverticaldisplacement', ...
                'tidalgravitationalpotential', 'tidalverticaldisplacement');
            testCase.assertEqual(which('solvedegree1'), fullfile(repo, 'solvedegree1.m'));
            testCase.assertEqual(which('solvesle'), fullfile(repo, 'solvesle.m'));
        end
    end

    methods (Test)
        function matchesFullFieldReference(testCase, includeC20, includeC30, includeC21S21, method, epochCount, adaptive)
            global JOINT_LOW_DEGREE_FIXTURE
            [order, degree] = addmon(3);
            meta = repmat([degree, order], 1, 1, epochCount);
            values = reshape(sin(1:20*epochCount), 10, 2, epochCount);
            sigma = reshape(1:20*epochCount, 10, 2, epochCount)/100;
            gia = cat(2, meta, values/7);
            p = cat(2, meta, values);
            s = cat(2, meta, sigma);
            dates = datetime(2003, (1:epochCount)', 15);
            JOINT_LOW_DEGREE_FIXTURE = struct('p', permute(p, [3, 1, 2]), ...
                's', permute(s, [3, 1, 2]), 'gia', permute(gia, [3, 1, 2]), 'dates', dates);
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
            testCase.verifyEqual(actualDates, dates);
            testCase.verifyEqual(actualLocs, locs);
            for k = 1:size(locs, 1)
                testCase.verifyEqual(field(:, locs(k, 1), 2 + locs(k, 2)), expected(:, k), 'AbsTol', 1e-10);
            end
        end
    end
end

function [coeffs, stds] = fullFieldReference(p, s, gia, locs, ids, method, tolerance, limit)
    % Keep the original full-localization arithmetic as the comparison oracle.
    p(4:end, 3:4, :) = p(4:end, 3:4, :) - gia(4:end, 3:4, :);
    for k = 1:size(locs, 1)
        p(locs(k, 1), 2 + locs(k, 2), :) = 0;
        s(locs(k, 1), 2 + locs(k, 2), :) = 0;
    end
    K = kernelcp_new(3); land = eye(16) - K; block = K(ids, ids);
    known = select(localise(p, 'L', 3, 'K', K), locs);
    knownStd = select(localise(s, 'L', 3, 'K', K, 'IsError', true), locs);
    coeffs = zeros(size(known)); stds = coeffs;
    for iteration = 1:limit
        old = coeffs; oldStd = stds;
        f = localise(p, 'L', 3, 'K', land);
        fs = localise(s, 'L', 3, 'K', land, 'IsError', true);
        if strcmp(method, 'fingerprint')
            [r, rs] = solvesle(f, fs, 3, GeoDomain('alloceans'), ...
                'OceanKernel', K, 'KernelOrder', kernelorder(3), ...
                'InitialRSL', p/1000, 'RotationFeedback', true, ...
                'Tolerance', 1e-9, 'MaxIter', 100, 'BeQuiet', true);
            ocean = select(localise(r, 'L', 3, 'K', K), locs);
            oceanStd = select(localise(rs, 'L', 3, 'K', K, 'IsError', true), locs);
        else
            ocean = K(ids, 1) * reshape(-f(1, 3, :)/K(1, 1), 1, []);
            oceanStd = K(ids, 1) * reshape(-fs(1, 3, :)/K(1, 1), 1, []);
        end
        coeffs = block \ (ocean - known);
        stds = sqrt((block \ eye(numel(ids))).^2 * (oceanStd.^2 + knownStd.^2));
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
