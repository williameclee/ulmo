%% STERIC2GMSLTEST - Native mean, schema and cache tests using temporary files.
%
% Syntax
%   results = runtests('unittests/steric2gmslTest.m')
%
% Notes
%   Exercises component/layer means on uniform and nonuniform latitude grids,
%   finite support, native schemas, output formats and cache reuse/rebuilding.
%   Latitude cases include descending order and clipped polar boundaries.
%   Each test uses a temporary IFILES folder; production files are not read
%   or changed.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

classdef steric2gmslTest < matlab.unittest.TestCase
    properties
        nativePath
        data
    end

    methods (TestMethodSetup)
        function prepareNativeFile(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.nativePath = stericfilepath('EN4c14');
            mkdir(fileparts(testCase.nativePath));
            data.lon = [0, 120, 240];
            data.lat = [0; 60];
            data.dates = datetime(2020, 1, (1:5)');
            base = [1, 3, 5; 2, 4, 6] + reshape(0:4, 1, 1, []);
            data.shallowHalostericSl = base;
            data.deepHalostericSl = 0.5 * base - 0.5;
            data.halostericSl = data.shallowHalostericSl + data.deepHalostericSl;
            data.shallowThermostericSl = 2 * base;
            data.deepThermostericSl = base + 1;
            data.thermostericSl = data.shallowThermostericSl + data.deepThermostericSl;
            data.shallowStericSl = 10 * base;
            data.deepStericSl = base + 2;
            data.stericSl = data.shallowStericSl + data.deepStericSl;
            data.preservedNote = 'An unrelated native variable';
            testCase.data = data;
            testCase.writeNative;
        end
    end

    methods (Test)
        function weightedComponentsAndResiduals(testCase)
            % Independent latitude-weighted base mean: (9 + 0.5*12)/4.5.
            baseMean = 10/3 + (0:4)';
            expected = [11*baseMean+2, 3*baseMean+1, 1.5*baseMean-0.5, ...
                8*baseMean+1, 9.5*baseMean+2.5];
            types = ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"];
            for iType = 1:numel(types)
                [actual, dates] = steric2gmsl('EN4c14', Type=types(iType), BeQuiet=true);
                testCase.verifyEqual(actual, expected(:, iType), AbsTol=1e-12);
                testCase.verifyEqual(dates, testCase.data.dates);
            end
            stored = load(testCase.nativePath);
            testCase.verifyEqual(stored.preservedNote, testCase.data.preservedNote);
            original = fieldnames(testCase.data);
            for iName = 1:numel(original)
                testCase.verifyEqual(stored.(original{iName}), testCase.data.(original{iName}));
            end
            cached = who('-file', testCase.nativePath);
            testCase.verifyEqual(nnz(endsWith(cached, 'Gmsl')), 15);
            testCase.verifyEqual(stored.nonThermostericGmsl, stored.stericGmsl-stored.thermostericGmsl, AbsTol=1e-12);
            testCase.verifyNotEqual(stored.nonThermostericGmsl, stored.halostericGmsl);
        end

        function directLayersRetainZeroContributions(testCase)
            shallow = steric2gmsl('EN4c14', Type="halosteric", Depth="shallow", BeQuiet=true);
            deep = steric2gmsl('EN4c14', Type="halosteric", Depth="deep", BeQuiet=true);
            full = steric2gmsl('EN4c14', Type="halosteric", BeQuiet=true);
            testCase.verifyEqual(shallow, 10/3+(0:4)', AbsTol=1e-12);
            testCase.verifyEqual(deep, 7/6+0.5*(0:4)', AbsTol=1e-12);
            testCase.verifyEqual(full, shallow+deep);
        end

        function commonSupportAndFiniteNonzeroCounts(testCase)
            % Reference cell (1,1) has only three finite nonzero samples;
            % a zero and an Inf must exclude it for the entire record.
            testCase.data.halostericSl(1,1,1) = Inf;
            testCase.data.halostericSl(1,1,2) = 0;
            % A different component removes one more cell at date 3 only.
            testCase.data.thermostericSl(2,3,3) = NaN;
            testCase.writeNative;
            actual = steric2gmsl('EN4c14', Type="halosteric", Depth="shallow", SaveData=false, BeQuiet=true);
            expected = 4+(0:4)'; % Remaining weights sum to 3.5.
            expected(3) = (5+7+0.5*(4+6))/3;
            testCase.verifyEqual(actual, expected, AbsTol=1e-12);
        end

        function noObservedSupportReturnsNaN(testCase)
            testCase.data.thermostericSl(:,:,2) = NaN;
            testCase.writeNative;
            actual = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            testCase.verifyTrue(isnan(actual(2)));
            testCase.verifyTrue(all(isfinite(actual([1,3:5]))));
        end

        function completeCacheDoesNotReadGrids(testCase)
            expected = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            writelines(["function out = matfile(varargin)"; ...
                "error('test:GridRead', 'A complete cache must bypass grid reads.');"; "end"], ...
                fullfile(folder.Folder, 'matfile.m'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder.Folder));
            actual = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            testCase.verifyEqual(actual, expected);
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total", ForceNew=true, BeQuiet=true), 'test:GridRead');
        end

        function incompleteAndMalformedCachesRebuildTogether(testCase)
            halostericGmsl = 999*ones(5,1);
            save(testCase.nativePath, 'halostericGmsl', '-append');
            actual = steric2gmsl('EN4c14', Type="halosteric", BeQuiet=true);
            testCase.verifyEqual(actual, 4.5+1.5*(0:4)', AbsTol=1e-12);
            stericGmsl = zeros(1,5);
            halostericGmsl = 999*ones(5,1);
            save(testCase.nativePath, 'stericGmsl', 'halostericGmsl', '-append');
            rebuilt = steric2gmsl('EN4c14', Type="halosteric", BeQuiet=true);
            testCase.verifyEqual(rebuilt, actual);
            stored = load(testCase.nativePath, 'stericGmsl');
            testCase.verifySize(stored.stericGmsl, [5,1]);
        end

        function saveFalsePreservesBytesAndForceNewRefreshes(testCase)
            before = testCase.fileBytes;
            expected = steric2gmsl('EN4c14', Type="total", SaveData=false, BeQuiet=true);
            testCase.verifyEqual(testCase.fileBytes, before);
            steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            % No automatic source invalidation: a changed grid keeps the cache
            % until an explicit ForceNew call; all refreshed means stay coherent.
            source = matfile(testCase.nativePath, 'Writable', true);
            source.stericSl = testCase.data.stericSl + 3;
            source.shallowStericSl = testCase.data.shallowStericSl + 3;
            before = testCase.fileBytes;
            cached = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            testCase.verifyEqual(cached, expected);
            fresh = steric2gmsl('EN4c14', Type="total", ForceNew=true, SaveData=false, BeQuiet=true);
            testCase.verifyEqual(fresh, expected+3, AbsTol=1e-12);
            testCase.verifyEqual(testCase.fileBytes, before);
            steric2gmsl('EN4c14', Type="total", ForceNew=true, BeQuiet=true);
            saved = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            testCase.verifyEqual(saved, fresh);
        end

        function formatsChangeOnlyReturnedVectors(testCase)
            [actual, dates] = steric2gmsl('EN4c14', Type="halosteric", Unit="mm", TimeFormat="datenum", BeQuiet=true);
            testCase.verifySize(actual, [5,1]);
            testCase.verifyEqual(actual, 1000*(4.5+1.5*(0:4)'), AbsTol=1e-9);
            testCase.verifyEqual(dates, datenum(testCase.data.dates)); %#ok<DATNM>
            [charValues, charDates] = steric2gmsl('EN4c14', Type="halosteric", Unit='mm', TimeFormat='datenum', BeQuiet=true);
            testCase.verifyEqual(charValues, actual);
            testCase.verifyEqual(charDates, dates);
            [metres, nativeDates] = steric2gmsl('EN4c14', Type="halosteric", Unit='m', TimeFormat='datetime', BeQuiet=true);
            testCase.verifyEqual(metres, actual/1000, AbsTol=1e-12);
            testCase.verifyEqual(nativeDates, testCase.data.dates);
            stored = load(testCase.nativePath, 'halostericGmsl', 'dates');
            testCase.verifyEqual(stored.halostericGmsl, actual/1000, AbsTol=1e-12);
            testCase.verifyClass(stored.dates, 'datetime');
        end

        function thermalOnlyProductAndUnavailableDepth(testCase)
            data = testCase.data;
            testCase.data = struct('lon', data.lon, 'lat', data.lat, 'dates', data.dates, ...
                'thermostericSl', data.thermostericSl, 'deepThermostericSl', []);
            testCase.writeNative;
            actual = steric2gmsl('EN4c14', Type="thermosteric", BeQuiet=true);
            testCase.verifyEqual(actual, 11+3*(0:4)', AbsTol=1e-12);
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="thermosteric", Depth="deep"), 'ULMO:steric2gmsl:UnsupportedField');
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="halosteric"), 'ULMO:steric2gmsl:UnsupportedField');
            stored = who('-file', testCase.nativePath);
            testCase.verifyEqual(nnz(endsWith(stored,'Gmsl')), 1);
        end

        function partialSplitFailsBeforeWriting(testCase)
            testCase.data = rmfield(testCase.data, 'deepHalostericSl');
            testCase.writeNative;
            before = testCase.fileBytes;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:InvalidSchema');
            testCase.verifyEqual(testCase.fileBytes, before);
        end

        function inconsistentDepthFailsBeforeWriting(testCase)
            testCase.data.stericSl = testCase.data.stericSl + 1;
            testCase.writeNative;
            before = testCase.fileBytes;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:DepthMismatch');
            testCase.verifyEqual(testCase.fileBytes, before);
        end

        function cyclicEndpointAndTransposedGrids(testCase)
            expected = steric2gmsl('EN4c14', Type="total", SaveData=false, BeQuiet=true);
            testCase.data.lon = [0,120,240,360];
            names = fieldnames(testCase.data);
            for iName = 1:numel(names)
                if endsWith(names{iName}, 'Sl')
                    field = testCase.data.(names{iName});
                    field(:,4,:) = field(:,1,:);
                    testCase.data.(names{iName}) = permute(field,[2,1,3]);
                end
            end
            testCase.writeNative;
            actual = steric2gmsl('EN4c14', Type="total", SaveData=false, BeQuiet=true);
            testCase.verifyEqual(actual, expected);
        end

        function invalidGeometryAndDimensions(testCase)
            testCase.data.lon = [0,110,240];
            testCase.writeNative;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:InvalidGrid');
            testCase.data.lon = [0,120,240];
            testCase.data.stericSl = testCase.data.stericSl(:,:,1:4);
            testCase.writeNative;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:InvalidSchema');
        end

        function fractionalSingleLongitudeCoordinates(testCase)
            % Single-precision rounding on a regular 1/12-degree grid must
            % not reject equivalent negative, zero or positive origins.
            dates = testCase.data.dates;
            for lonOrigin = [-180, 0, 20]
                testCase.data = struct('lon', single(lonOrigin + (0:4319)/12), ...
                    'lat', [-1; 1], 'dates', dates, ...
                    'thermostericSl', ones(2, 4320, numel(dates), 'single'));
                testCase.writeNative;
                actual = steric2gmsl('EN4c14', Type="thermosteric", SaveData=false, BeQuiet=true);
                testCase.verifyEqual(actual, ones(numel(dates), 1), AbsTol=1e-12);
            end
        end

        function nonuniformLatitudeMidpointAreas(testCase)
            % Explicit boundaries cover polar clipping and limited coverage;
            % neither irregular spacing nor descending order needs regridding.
            latitudes = {[-80; -20; 30; 85], [-40; -10; 10; 35]};
            boundaries = {[-90; -50; 5; 57.5; 90], [-55; -25; 0; 22.5; 47.5]};
            rowValues = [1; 4; 9; 16];
            dates = testCase.data.dates;
            for iGrid = 1:numel(latitudes)
                area = diff(sind(boundaries{iGrid}));
                expected = sum(area .* rowValues) / sum(area) + (0:4)';
                field = repmat(rowValues, 1, 3, 5) + reshape(0:4, 1, 1, []);
                for descending = [false, true]
                    latitude = latitudes{iGrid};
                    grid = field;
                    if descending
                        latitude = flip(latitude);
                        grid = flip(grid, 1);
                    end
                    testCase.data = struct('lon', [0, 120, 240], 'lat', latitude, ...
                        'dates', dates, 'stericSl', grid, ...
                        'shallowStericSl', 0.25 * grid, 'deepStericSl', 0.75 * grid);
                    testCase.writeNative;
                    full = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
                    shallow = steric2gmsl('EN4c14', Type="total", Depth="shallow", BeQuiet=true);
                    deep = steric2gmsl('EN4c14', Type="total", Depth="deep", BeQuiet=true);
                    testCase.verifyEqual(full, expected, AbsTol=1e-12);
                    testCase.verifyEqual(shallow, 0.25 * expected, AbsTol=1e-12);
                    testCase.verifyEqual(deep, 0.75 * expected, AbsTol=1e-12);
                    testCase.verifyEqual(full, shallow + deep);
                end
            end
        end

        function nearlyUniformLatitudeRetainsCosineWeights(testCase)
            % Coordinate rounding below the spacing tolerance retains the
            % existing uniform-grid convention, including at the outer rows.
            latitude = [-60; 0; 60 + 1e-5];
            rowValues = [1; 5; 13];
            testCase.data = struct('lon', [0, 90, 180, 270], 'lat', latitude, ...
                'dates', testCase.data.dates, 'stericSl', repmat(rowValues, 1, 4, 5));
            testCase.writeNative;
            actual = steric2gmsl('EN4c14', Type="total", SaveData=false, BeQuiet=true);
            expected = sum(cosd(latitude) .* rowValues) / sum(cosd(latitude));
            testCase.verifyEqual(actual, repmat(expected, 5, 1), AbsTol=1e-12);
        end

        function latitudeMustBeStrictlyMonotone(testCase)
            testCase.data = struct('lon', [0, 90, 180, 270], 'lat', [0; 60; 30], ...
                'dates', testCase.data.dates, 'stericSl', ones(3, 4, 5));
            testCase.writeNative;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:InvalidGrid');
            testCase.data.lat = [0; 30; 30];
            testCase.writeNative;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:InvalidGrid');
        end

        function numericAndInvalidDates(testCase)
            expected = testCase.data.dates;
            testCase.data.dates = datenum(expected); %#ok<DATNM>
            testCase.writeNative;
            [~, dates] = steric2gmsl('EN4c14', Type="total", SaveData=false, BeQuiet=true);
            testCase.verifyEqual(dates, expected);
            testCase.data.dates(2) = testCase.data.dates(1);
            testCase.writeNative;
            testCase.verifyError(@() steric2gmsl('EN4c14', Type="total"), 'ULMO:steric2gmsl:InvalidDates');
        end

        function singleSampleSinglePrecisionAndZeroDeep(testCase)
            % Native single-precision storage permits rounding in depth closure;
            % a singleton time axis and an entirely zero layer remain supported.
            names = fieldnames(testCase.data);
            for iName = 1:numel(names)
                if endsWith(names{iName}, 'Sl')
                    testCase.data.(names{iName}) = single(testCase.data.(names{iName})(:,:,1)/7);
                end
            end
            testCase.data.dates = testCase.data.dates(1);
            testCase.data.deepHalostericSl(:) = 0;
            testCase.data.halostericSl = testCase.data.shallowHalostericSl;
            testCase.writeNative;
            actual = steric2gmsl('EN4c14', Type="halosteric", Depth="deep", BeQuiet=true);
            testCase.verifyEqual(actual, 0);
            total = steric2gmsl('EN4c14', Type="total", BeQuiet=true);
            testCase.verifyEqual(total, (11*10/3+2)/7, AbsTol=1e-6);
        end

        function explicitTypeAndMissingFile(testCase)
            testCase.verifyError(@() steric2gmsl('EN4c14'), 'ULMO:steric2gmsl:TypeRequired');
            testCase.verifyError(@() steric2gmsl('NOAA', Type="total"), 'ULMO:steric2gmsl:FileNotFound');
        end
    end

    methods (Access=private)
        function writeNative(testCase)
            data = testCase.data;
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
        end

        function bytes = fileBytes(testCase)
            fid = fopen(testCase.nativePath, 'rb');
            cleanup = onCleanup(@() fclose(fid));
            bytes = fread(fid, Inf, '*uint8');
        end
    end
end
