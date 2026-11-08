%% CONVERTSTERICSOURCETSTEST - Tests native T/S conversion and forced reprocessing.
% Compares stored CT, SA, pressure, depth, and bottom values with direct
% GSW calculations using synthetic monthly MAT files.
%
% Syntax
%   results = runtests(fullfile('unittests', 'convertStericSourceTSTest.m'));
%   assertSuccess(results);
%
% Output arguments
%   results - Array of TestResult objects containing the outcome of each test method.
%
% See also
%   convertStericSourceTS
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

classdef convertStericSourceTSTest < matlab.unittest.TestCase

    properties (Access = private)
        MonthPath
    end

    methods (TestClassSetup)
        % Make the local converter available and restore the path after the suite.
        function addSourcePath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'aux')));
        end

    end

    methods (TestMethodSetup)
        % Allocate an isolated monthly MAT file that is removed after each test.
        function createFixture(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            tc.MonthPath = fullfile(folder.Folder, 'month.mat');
        end

    end

    methods (Test)
        % Check pressure-to-depth conversion for shared and latitude-dependent levels.
        function testPressureVectorAndMatrix(tc)
            tc.checkConversion([0; 1000], [-60; 20], 'p');
            tc.checkConversion([0, 1000; 10, 1200], [-60; 20], 'p');
        end

        % Check depth-to-pressure conversion without flattening coordinate matrices.
        function testDepthVectorAndMatrix(tc)
            tc.checkConversion([0, 1000], [-60; 20], 'z');
            tc.checkConversion([0, 1000; 10, 1200], [-60; 20], 'z');
        end

        % Compare every supported temperature/salinity interpretation with GSW.
        function testTemperatureAndSalinityTypes(tc)

            for ttype = {'T', 'PT', 'CT'}

                for stype = {'SP', 'SA'}
                    tc.checkConversion([0; 1200; 2400], [-70; 0; 70], ...
                        'z', 3, ttype{1}, stype{1});
                end

            end

        end

        % Distinguish a single level per latitude from a shared multi-level vector.
        function testSingletonDimensions(tc)
            tc.checkConversion([0, 1000], 20, 'p');
            tc.checkConversion([100; 200], [-60; 20], 'p', 1);
            tc.checkConversion(100, [-60; 20], 'z', 1);
        end

        % Reject coordinates whose shape cannot represent the source field levels.
        function testInvalidVerticalShape(tc)
            data = convertStericSourceTSTest.nativeData(ones(3, 2), [-60; 20], 'p', 2);
            save(tc.MonthPath, '-struct', 'data');
            tc.verifyError(@() convertStericSourceTS(tc.MonthPath), ...
            'ULMO:convertStericSourceTS:InvalidInputSize');
        end

        % Reject descending levels before computing pressure or seawater properties.
        function testNonIncreasingLevels(tc)
            data = convertStericSourceTSTest.nativeData([1000, 0], [-60; 20], 'p', 2);
            save(tc.MonthPath, '-struct', 'data');
            tc.verifyError(@() convertStericSourceTS(tc.MonthPath), ...
            'ULMO:convertStericSourceTS:InvalidVerticalCoordinate');
        end

    end

    methods (Access = private)
        % Compare the first conversion with GSW and verify stable forced reruns.
        function checkConversion(tc, z, lat, ztype, nLevels, ttype, stype)

            arguments (Input)
                tc
                z
                lat
                ztype
                nLevels = 2
                ttype = 'T'
                stype = 'SP'
            end

            data = convertStericSourceTSTest.nativeData(z, lat, ztype, nLevels);
            data.ttype = ttype;
            data.stype = stype;
            save(tc.MonthPath, '-struct', 'data');

            if ~isequal(size(z), [numel(lat), nLevels])
                z = repmat(z(:)', numel(lat), 1);
            end

            if strcmp(ztype, 'p')
                expectedP = z;
                expectedZ = -gsw_z_from_p(z, lat(:));
            else
                expectedP = gsw_p_from_z(-z, lat(:));
                expectedZ = z;
            end

            expectedS = nan(size(data.S));
            expectedT = nan(size(data.T));

            for k = 1:nLevels

                if strcmp(stype, 'SP')
                    expectedS(:, :, k) = gsw_SA_from_SP(data.S(:, :, k), ...
                        expectedP(:, k), data.lon, lat(:));
                else
                    expectedS(:, :, k) = data.S(:, :, k);
                end

                switch ttype
                    case 'T'
                        expectedT(:, :, k) = gsw_CT_from_t(expectedS(:, :, k), ...
                            data.T(:, :, k), expectedP(:, k));
                    case 'PT'
                        expectedT(:, :, k) = gsw_CT_from_pt(expectedS(:, :, k), data.T(:, :, k));
                    case 'CT'
                        expectedT(:, :, k) = data.T(:, :, k);
                end

            end

            convertStericSourceTS(tc.MonthPath);
            first = load(tc.MonthPath);
            tc.verifyEqual(first.T, single(expectedT));
            tc.verifyEqual(first.S, single(expectedS));
            tc.verifyEqual(first.p, expectedP, AbsTol = 1e-8);
            tc.verifyEqual(first.z, expectedZ, AbsTol = 1e-8);
            tc.verifyEqual(first.bottom, expectedZ(:, end), AbsTol = 1e-8);
            tc.verifyEqual({first.ttype, first.stype, first.ztype}, {'CT', 'SA', 'z'});

            for k = 1:2
                convertStericSourceTS(tc.MonthPath, ForceNew = true);
                repeated = load(tc.MonthPath);
                tc.verifyEqual(repeated.T, first.T);
                tc.verifyEqual(repeated.S, first.S);
                tc.verifyEqual(repeated.p, first.p, AbsTol = 1e-8);
                tc.verifyEqual(repeated.z, first.z);
                tc.verifyEqual(repeated.bottom, first.bottom);
            end

        end

    end

    methods (Static, Access = private)
        % Build synthetic native fields with a missing column to test NaN preservation.
        function data = nativeData(z, lat, ztype, nLevels)
            data = struct('T', 10 * ones(numel(lat), 3, nLevels), ...
                'S', 35 * ones(numel(lat), 3, nLevels), 'z', z, ...
                'lon', [0, 120, 240], 'lat', lat(:), ...
                'date', datetime(2005, 1, 16), 'ttype', 'T', 'stype', 'SP', 'ztype', ztype);
            data.T(:, 3, :) = NaN;
            data.S(:, 3, :) = NaN;
        end

    end

end
