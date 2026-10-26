%% jointLowDegreeIntegrationTest - Checks joint coefficients across reader and solver boundaries.
% Uses isolated caches and synthetic degree-3 fields; no production data are read.
%
% Created by
%   2026/10/26, En-Chi Lee (williameclee@gmail.com)

classdef jointLowDegreeIntegrationTest < matlab.unittest.TestCase

    properties
        Folder
        Repo
    end

    properties (TestParameter)
        includeC20 = {false, true}
        includeC30 = {false, true}
        unit = {'SD', 'POT', 'GRAV'}
        filtered = {false, true}
        dateCase = {'reordered', 'extra'}
    end

    methods (TestMethodSetup)

        function configure(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.Repo = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(PathFixture(testCase.Repo));
            testCase.applyFixture(PathFixture(fullfile(testCase.Repo, 'aux')));
            temp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.Folder = temp.Folder;

            for name = {'GRACEDATA', 'ORIGINALGRACEDATA', 'IFILES'}
                original = getenv(name{1});
                testCase.addTeardown(@() setenv(name{1}, original));
                setenv(name{1}, temp.Folder);
            end

            testCase.addTeardown(@() clearFixture());

            % The real reader's default freshness check must prefer the seeded cache.
            inputFolder = fullfile(temp.Folder, 'RL06', 'CSR');
            mkdir(inputFolder);
            inputDirectory = java.io.File(inputFolder);
            testCase.assertTrue(inputDirectory.setLastModified(946684800000));

            mkdir(fullfile(temp.Folder, 'EARTHMODELS', 'CONSTANTS'));
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(temp.Folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
            mkdir(fullfile(temp.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; .021; - .2; - .3];
            save(fullfile(temp.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
            'loadinggravitationalpotential');
            mkdir(fullfile(temp.Folder, 'Degree1', 'new'));

            global JOINT_LOW_DEGREE_FIXTURE
            p = zeros(3, 10, 4);
            p(:, :, 1) = repmat([0, 1, 1, 2, 2, 2, 3, 3, 3, 3], 3, 1);
            p(:, :, 2) = repmat([0, 0, 1, 0, 1, 2, 0, 1, 2, 3], 3, 1);
            p(:, :, 3:4) = reshape(1:60, 3, 10, 2) / 10;
            s = p;
            s(:, :, 3:4) = reshape(1:60, 3, 10, 2) / 1000;
            JOINT_LOW_DEGREE_FIXTURE = struct('p', p, 's', s, ...
                'dates', datetime(2003, (1:3)', 15), 'gia', zeros(size(p)), ...
                'coeffs', reshape(1:21, 3, 7), 'coeffStds', reshape(1:21, 3, 7) / 100, ...
                'solverDates', datetime(2003, (1:3)', 15), ...
                'locs', [2, 1; 3, 1; 3, 2; 4, 1; 7, 1; 5, 1; 5, 2]);
        end

    end

    methods (Test)

        function insertsCoefficientsByDate(testCase, dateCase)
            global JOINT_LOW_DEGREE_FIXTURE
            testCase.installFixtures('jointlowdegreereader');
            testCase.assertEqual(which('grace2plmt_new'), fullfile(testCase.Repo, 'grace2plmt_new.m'));
            testCase.seedGraceCache('SD');
            f = JOINT_LOW_DEGREE_FIXTURE;

            if strcmp(dateCase, 'reordered')
                order = [3, 1, 2];
            else
                order = [4, 3, 1, 2];
                f.coeffs(4, :) = 100 + (1:7);
                f.coeffStds(4, :) = 1 + (1:7) / 100;
            end

            allDates = datetime(2003, (1:4)', 15);
            JOINT_LOW_DEGREE_FIXTURE.coeffs = f.coeffs(order, :);
            JOINT_LOW_DEGREE_FIXTURE.coeffStds = f.coeffStds(order, :);
            JOINT_LOW_DEGREE_FIXTURE.solverDates = allDates(order);
            [actual, actualStd, dates] = readJoint();
            expected = f.p;
            expectedStd = f.s;

            for k = 1:7
                expected(:, f.locs(k, 1), 2 + f.locs(k, 2)) = f.coeffs(1:3, k);
                expectedStd(:, f.locs(k, 1), 2 + f.locs(k, 2)) = f.coeffStds(1:3, k);
            end

            % Whole-array comparison also checks that unrelated fields are untouched.
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(actualStd, expectedStd);
            testCase.verifyEqual(dates, f.dates);
        end

        function rejectsMissingSolverDates(testCase)
            global JOINT_LOW_DEGREE_FIXTURE
            testCase.installFixtures('jointlowdegreereader');
            testCase.seedGraceCache('SD');
            JOINT_LOW_DEGREE_FIXTURE.solverDates(2) = [];
            JOINT_LOW_DEGREE_FIXTURE.coeffs(2, :) = [];
            JOINT_LOW_DEGREE_FIXTURE.coeffStds(2, :) = [];
            testCase.verifyError(@() readJoint(), 'ULMO:MissingLowDegreeDates');
        end

        function restoresNonzeroGia(testCase, includeC20, includeC30)
            global JOINT_LOW_DEGREE_FIXTURE
            testCase.installFixtures('jointlowdegreesolver');
            testCase.assertEqual(which('solvedegree1'), fullfile(testCase.Repo, 'solvedegree1.m'));
            opts = solverOptions(includeC20, includeC30);
            [baseline, baselineStd, dates] = solvedegree1({'CSR', 'RL06', 3}, ...
                opts{:}, 'ForceNew', true, 'SaveData', false);
            gia = zeros(size(JOINT_LOW_DEGREE_FIXTURE.p));
            gia(:, :, 3:4) = reshape(1:60, 3, 10, 2) / 20;
            JOINT_LOW_DEGREE_FIXTURE.gia = gia;
            % Keep the corrected forcing identical to baseline. Nonzero GIA at
            % degrees 0/1 must be ignored, while estimated higher degrees return.
            JOINT_LOW_DEGREE_FIXTURE.p(:, 4:end, 3:4) = ...
                JOINT_LOW_DEGREE_FIXTURE.p(:, 4:end, 3:4) + gia(:, 4:end, 3:4);
            [actual, actualStd, actualDates] = solvedegree1({'CSR', 'RL06', 3}, ...
                opts{:}, 'ForceNew', true, 'SaveData', false);
            locs = expectedLocations(includeC20, includeC30);
            expected = baseline;

            for k = 4:size(locs, 1)
                expected(:, k) = expected(:, k) + gia(:, locs(k, 1), 2 + locs(k, 2));
            end

            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(actualStd, baselineStd, 'AbsTol', 1e-12);
            testCase.verifyEqual(actualDates, dates);
        end

        function formatsMixedDegreeOutputs(testCase, includeC20, includeC30, unit, filtered)
            global JOINT_LOW_DEGREE_FIXTURE
            % Use the real reader and solver together, backed only by synthetic caches.
            testCase.assertEqual(which('solvedegree1'), fullfile(testCase.Repo, 'solvedegree1.m'));
            testCase.assertEqual(which('grace2plmt_new'), fullfile(testCase.Repo, 'grace2plmt_new.m'));
            testCase.seedGraceCache(unit);
            [locs, degrees] = expectedLocations(includeC20, includeC30);
            n = size(locs, 1);
            coeffs = reshape(1:3 * n, 3, n) .* (-1) .^ (1:n);
            coeffStds = reshape(1:3 * n, 3, n) / 100;
            dates = JOINT_LOW_DEGREE_FIXTURE.dates;
            domain = GeoDomain('alloceans');
            flag = 'C21S21';

            if includeC30
                flag = ['C30', flag];
            end

            if includeC20
                flag = ['C20', flag];
            end

            filename = sprintf('CSRRL063_ICE6GD-W%s-Ld3_Ls3-%s-uniform-SD.mat', flag, domain.Id);
            save(fullfile(testCase.Folder, 'Degree1', 'new', filename), 'coeffs', 'coeffStds', 'dates');
            opts = solverOptions(includeC20, includeC30);
            keep = 1:3;

            if filtered
                keep = 2:3;
                opts = [opts, {'TimeRange', dates([2, 3])}];
            end

            [actual, actualStd, actualDates, actualPlmt, actualLocs] = ...
                solvedegree1({'CSR', 'RL06', 3}, opts{:}, 'Unit', unit);
            factors = densityFactors(degrees, unit);
            expected = coeffs(keep, :) .* factors;
            expectedStd = coeffStds(keep, :) .* abs(factors);
            expectedPlmt = convertField(JOINT_LOW_DEGREE_FIXTURE.p, unit);
            expectedPlmt = expectedPlmt(keep, :, :);

            for k = 1:n
                expectedPlmt(:, locs(k, 1), 2 + locs(k, 2)) = expected(:, k);
            end

            testCase.verifyEqual(actual, expected, 'RelTol', 1e-13);
            testCase.verifyEqual(actualStd, expectedStd, 'RelTol', 1e-13);
            testCase.verifyEqual(actualDates, dates(keep));
            testCase.verifyEqual(actualLocs, locs);
            testCase.verifyEqual(actualPlmt, expectedPlmt, 'RelTol', 1e-13);
            [fourC, fourS, fourDates, fourPlmt] = ...
                solvedegree1({'CSR', 'RL06', 3}, opts{:}, 'Unit', unit);
            testCase.verifyEqual(fourC, actual);
            testCase.verifyEqual(fourS, actualStd);
            testCase.verifyEqual(fourDates, actualDates);
            testCase.verifyEqual(fourPlmt, actualPlmt);
        end

    end

    methods (Access = private)

        function installFixtures(testCase, group)
            import matlab.unittest.fixtures.PathFixture
            folder = fullfile(testCase.Folder, group);
            mkdir(folder);
            source = fullfile(testCase.Repo, 'unittests', 'fixtures', group);
            files = dir(fullfile(source, '*.txt'));

            for k = 1:numel(files)
                copyfile(fullfile(source, files(k).name), fullfile(folder, erase(files(k).name, '.txt')));
            end

            testCase.applyFixture(PathFixture(folder));
        end

        function seedGraceCache(testCase, unit)
            global JOINT_LOW_DEGREE_FIXTURE
            gracePlmt = convertField(JOINT_LOW_DEGREE_FIXTURE.p, unit);
            graceStdPlmt = convertField(JOINT_LOW_DEGREE_FIXTURE.s, unit);
            dates = JOINT_LOW_DEGREE_FIXTURE.dates;

            switch unit
                case 'SD'
                    filename = 'CSR_RL06_3_alldata_SD_ISSM.mat';
                case 'GRAV'
                    filename = 'CSR_RL06_alldata_3_GRAV.mat';
                case 'POT'
                    filename = 'CSR_RL06_alldata_3.mat';
            end

            save(fullfile(testCase.Folder, filename), 'gracePlmt', 'graceStdPlmt', 'dates');
        end

    end

end

function opts = solverOptions(c20, c30)
    opts = {'Ltruncation', 3, 'Lsle', 3, 'GIAModel', 'ICE6GD', ...
                'OceanDomain', GeoDomain('alloceans'), 'Method', 'uniform', ...
                'IncludeC20', c20, 'IncludeC30', c30, 'IncludeC21S21', true, 'BeQuiet', true};
end

function varargout = readJoint()
    opts = solverOptions(true, true);
    [varargout{1:nargout}] = grace2plmt_new('CSR', 'RL06', 3, 'Unit', 'SD', ...
        'RecomputeDegree1', opts, 'TimeFormat', 'datetime', ...
        'ForceNew', false, 'SaveData', false, 'BeQuiet', true);
end

function [locs, degrees] = expectedLocations(c20, c30)
    % Explicit harmonic locations, independent of lowdegreecoefficients.
    locs = [2, 1; 3, 1; 3, 2];
    degrees = [1, 1, 1];

    if c20
        locs = [locs; 4, 1];
        degrees = [degrees, 2];
    end

    if c30
        locs = [locs; 7, 1];
        degrees = [degrees, 3];
    end

    locs = [locs; 5, 1; 5, 2];
    degrees = [degrees, 2, 2];
end

function factors = densityFactors(degrees, unit)
    % Analytic SD conversion for the fixture's Love numbers and Earth constants.
    love = [0, .021, - .2, - .3];
    factors = ones(size(degrees));

    if ~strcmp(unit, 'SD')
        factors = 3/5517 .* (1 + love(degrees + 1)) ./ (2 * degrees + 1);

        if strcmp(unit, 'GRAV')
            factors = factors / 6378136.3;
        end

    end

end

function p = convertField(p, unit)
    factors = densityFactors([0, 1, 1, 2, 2, 2, 3, 3, 3, 3], unit);
    p(:, :, 3:4) = p(:, :, 3:4) .* factors;
end

function clearFixture()
    clear global JOINT_LOW_DEGREE_FIXTURE
end
