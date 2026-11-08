%% STERICRAWFIELDSTEST - Raw loader, temporal GMSL and cache contract tests.
% Uses temporary native grids and isolated IFILES. No production data is read.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/08, OpenAI Codex
classdef stericRawFieldsTest < matlab.unittest.TestCase
    properties
        fixtureData
        nativePath
    end
    methods (TestMethodSetup)
        function prepare(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.nativePath = stericfilepath('EN4c14');
            mkdir(fileparts(testCase.nativePath));
            data.lon = [0, 120, 240];
            data.lat = [-60; -20; 20; 60];
            data.dates = datetime(2020, 1, (1:60)');
            base = [1, 2, 3; 4, 5, 6; 7, 8, 9; 10, 11, 12] + reshape((0:59)/60, 1, 1, []);
            data.shallowStericSl = 10 * base;
            data.deepStericSl = 2 * base;
            data.stericSl = data.shallowStericSl + data.deepStericSl;
            data.shallowHalostericSl = base;
            data.deepHalostericSl = 0.5 * base;
            data.halostericSl = data.shallowHalostericSl + data.deepHalostericSl;
            data.thermostericSl = 3 * base;
            data.shallowThermostericSl = 2 * base;
            data.deepThermostericSl = base;
            data.note = 'preserve native metadata';
            testCase.fixtureData = data;
            testCase.writeNative;
        end
    end
    methods (Test)
        function rawDefaultAndFormats(testCase)
            [raw, sigma, dates, lon, lat, halo] = steric2lonlatt('EN4c14', ...
                Unit='mm', TimeFormat='datenum', OutputFormat='ndgrid', BeQuiet=true);
            testCase.verifyEqual(raw, 1000 * permute(testCase.fixtureData.stericSl, [2, 1, 3]));
            testCase.verifyEqual(sigma, zeros(size(raw)));
            testCase.verifyEqual(dates, datenum(testCase.fixtureData.dates)); %#ok<DATNM>
            testCase.verifyEqual(lon, testCase.fixtureData.lon');
            testCase.verifyEqual(lat, testCase.fixtureData.lat');
            testCase.verifyEqual(halo, 1000 * (9.75 + (0:59)'/40), AbsTol=1e-11);
            testCase.verifyEqual(steric2lonlatt('EN4c14', Type='total', BeQuiet=true), testCase.fixtureData.stericSl);
        end
        function removedTypes(testCase)
            for type = ["default", "full", "no halosteric drift"]
                testCase.verifyError(@() steric2lonlatt('EN4c14', Type=type), 'MATLAB:validators:mustBeMember');
            end
        end
        function rawOnlyWithoutHalo(testCase)
            data = rmfield(testCase.fixtureData, {'halostericSl', 'shallowHalostericSl', 'deepHalostericSl'});
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
            raw = steric2lonlatt('EN4c14', BeQuiet=true);
            testCase.verifyEqual(raw, data.stericSl);
            testCase.verifyFalse(any(endsWith(who('-file', testCase.nativePath), 'Gmsl')));
            testCase.verifyError(@() sixth('EN4c14'), 'ULMO:steric2lonlatt:UnsupportedField');
        end
        function independentMonthlyMeans(testCase)
            % One missing cell changes its effective correction-before-average.
            testCase.fixtureData.stericSl(1, 1, 1:10) = NaN;
            testCase.fixtureData.shallowStericSl(1, 1, 1:10) = NaN;
            testCase.fixtureData.deepStericSl(1, 1, 1:10) = NaN;
            testCase.writeNative;
            [raw, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', 'midmonth', BeQuiet=true);
            nativeHalo = steric2gmsl('EN4c14', 'halosteric', BeQuiet=true);
            testCase.verifyEqual(dates, datetime(2020, [1; 2], [16; 15], 12, 0, 0));
            testCase.verifyEqual(raw(:, :, 1), mean(testCase.fixtureData.stericSl(:, :, 1:31), 3, 'omitmissing'));
            testCase.verifyEqual(halo(1), mean(nativeHalo(1:31)), AbsTol=1e-12);
            before = mean(testCase.fixtureData.stericSl(1, 1, 1:31) - reshape(nativeHalo(1:31), 1, 1, []), 3, 'omitmissing');
            testCase.verifyGreaterThan(abs(raw(1, 1, 1) - halo(1) - before), 0.01);
        end
        function allMeansAndFieldsInEpochCache(testCase)
            sixth('EN4c14', 'midmonth');
            path = stericfilepath('EN4c14', 'midmonth', [], [], 'linear', false, [], '-full-total');
            stored = load(path);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Gmsl')), 15);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Sl')), 9);
            testCase.verifySize(stored.halostericGmsl, [2, 1]);
            testCase.verifyEqual(stored.halostericGmsl, stored.shallowHalostericGmsl + stored.deepHalostericGmsl);
        end
        function nonlinearDepthAdditivity(testCase)
            dates = testCase.fixtureData.dates([1, 16, 31, 46, 60]);
            testCase.fixtureData.dates = dates;
            names = fieldnames(testCase.fixtureData);
            for iName = 1:numel(names)
                if endsWith(names{iName}, 'Sl'), testCase.fixtureData.(names{iName}) = testCase.fixtureData.(names{iName})(:, :, [1,16,31,46,60]); end
            end
            testCase.fixtureData.shallowHalostericSl = ones(4, 3) .* reshape([1, 4, 2, 8, 3], 1, 1, []);
            testCase.fixtureData.deepHalostericSl = ones(4, 3) .* reshape([5, 2, 7, 1, 9], 1, 1, []);
            testCase.fixtureData.halostericSl = testCase.fixtureData.shallowHalostericSl + testCase.fixtureData.deepHalostericSl;
            testCase.writeNative;
            args = {'EN4c14', days(5), [], [], 'Interpolation', 'pchip', 'SaveData', false};
            full = sixth(args{:});
            shallow = sixth(args{:}, 'Depth', 'shallow');
            deep = sixth(args{:}, 'Depth', 'deep');
            testCase.verifyEqual(full, shallow + deep);
            expected = interp1(dates, [1;4;2;8;3], dates(1):days(5):dates(end), 'pchip') + ...
                interp1(dates, [5;2;7;1;9], dates(1):days(5):dates(end), 'pchip');
            testCase.verifyEqual(full, expected(:), AbsTol=1e-12);
        end
        function spatialGridDoesNotChangeMean(testCase)
            original = sixth('EN4c14');
            spatial = sixth('EN4c14', [], 60);
            testCase.verifyEqual(spatial, original);
        end
        function truncatedAndFilteredDates(testCase)
            limits = datetime(2020, 1, [10, 30]);
            [~, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', days(7), [], limits, ...
                TruncationBeforeInterpolation=true, Unit='mm', BeQuiet=true);
            testCase.verifyEqual(dates, datetime(2020, 1, [15;22;29]));
            testCase.verifyEqual(halo, 1000 * (9.75 + [14;21;28]/40), AbsTol=1e-11);
        end
        function augmentRawCacheWithoutGridRebuild(testCase)
            raw = steric2lonlatt('EN4c14', days(10), BeQuiet=true); %#ok<NASGU>
            path = stericfilepath('EN4c14', days(10), [], [], 'linear');
            cached = load(path);
            % A sentinel proves cached raw grids remain in use during augmentation.
            stericSl = cached.stericSl + 123;
            save(path, 'stericSl', '-append');
            [raw, ~, ~, ~, ~, halo] = steric2lonlatt('EN4c14', days(10), BeQuiet=true);
            testCase.verifyEqual(raw, stericSl);
            testCase.verifyEqual(halo, 9.75 + (0:10:59)'/40, AbsTol=1e-12);
            stored = load(path);
            testCase.verifyEqual(stored.stericSl, stericSl);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Gmsl')), 15);
        end
        function forcedFiveOutputsRebuildEverything(testCase)
            sixth('EN4c14', 'midmonth');
            old = load(testCase.nativePath);
            names = fieldnames(testCase.fixtureData);
            for iName = 1:numel(names)
                if endsWith(names{iName}, 'Sl'), testCase.fixtureData.(names{iName}) = 2 * testCase.fixtureData.(names{iName}); end
            end
            data = testCase.fixtureData;
            save(testCase.nativePath, '-struct', 'data', '-append');
            [raw, ~, ~, ~, ~] = steric2lonlatt('EN4c14', 'midmonth', ForceNew=true, BeQuiet=true);
            stored = load(stericfilepath('EN4c14', 'midmonth', [], [], 'linear', false, [], '-full-total'));
            native = load(testCase.nativePath);
            testCase.verifyEqual(native.halostericGmsl, 2 * old.halostericGmsl, AbsTol=1e-12);
            testCase.verifyEqual(native.deepNonThermostericGmsl, 2 * old.deepNonThermostericGmsl, AbsTol=1e-12);
            testCase.verifyEqual(stored.deepThermostericSl(:, :, 1), mean(data.deepThermostericSl(:, :, 1:31), 3));
            testCase.verifyEqual(raw, stored.stericSl);
        end
        function forcedUnsavedRebuildPreservesFiles(testCase)
            sixth('EN4c14', 'midmonth');
            data = testCase.fixtureData;
            data.shallowHalostericSl = 2 * data.shallowHalostericSl;
            data.deepHalostericSl = 2 * data.deepHalostericSl;
            data.halostericSl = data.shallowHalostericSl + data.deepHalostericSl;
            save(testCase.nativePath, '-struct', 'data', '-append');
            cache = stericfilepath('EN4c14', 'midmonth', [], [], 'linear', false, [], '-full-total');
            beforeNative = bytes(testCase.nativePath); beforeCache = bytes(cache);
            actual = sixth('EN4c14', 'midmonth', [], [], 'ForceNew', true, 'SaveData', false);
            expected = 2 * [mean(9.75+(0:30)/40); mean(9.75+(31:59)/40)];
            testCase.verifyEqual(actual, expected, AbsTol=1e-12);
            testCase.verifyEqual(bytes(testCase.nativePath), beforeNative);
            testCase.verifyEqual(bytes(cache), beforeCache);
        end
        function correctedLegacyCacheIsNotRaw(testCase)
            data = testCase.fixtureData;
            data.steric = -999 * ones(4, 3, 2);
            data.dates = datetime(2020, [1;2], [16;15], 12, 0, 0);
            legacy = stericfilepath('EN4c14', 'midmonth', [], [], 'linear', false, [], '-full-no-halosteric-drift');
            save(legacy, '-struct', 'data');
            raw = steric2lonlatt('EN4c14', 'midmonth', SaveData=false, BeQuiet=true);
            testCase.verifyEqual(raw(:, :, 1), mean(testCase.fixtureData.stericSl(:, :, 1:31), 3));
        end
        function rawComponentsAndRequestedDepth(testCase)
            types = ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"];
            expected = {testCase.fixtureData.deepStericSl, testCase.fixtureData.deepThermostericSl, ...
                testCase.fixtureData.deepHalostericSl, ...
                testCase.fixtureData.deepStericSl-testCase.fixtureData.deepThermostericSl, ...
                testCase.fixtureData.deepStericSl-testCase.fixtureData.deepHalostericSl};
            for iType = 1:numel(types)
                [raw, ~, ~, ~, ~, halo] = steric2lonlatt('EN4c14', ...
                    Type=types(iType), Depth='deep', BeQuiet=true);
                testCase.verifyEqual(raw, expected{iType});
                testCase.verifyEqual(halo, 3.25+(0:59)'/120, AbsTol=1e-12);
            end
        end
        function nativeAndDerivedMeansAreReused(testCase)
            expected = sixth('EN4c14', days(10));
            % Keep grid metadata intact but destroy values. Neither cached
            % native means nor a complete derived cache may recalculate them.
            data = testCase.fixtureData;
            fields = fieldnames(data);
            for iField = 1:numel(fields)
                if endsWith(fields{iField}, 'Sl'), data.(fields{iField})(:) = NaN; end
            end
            save(testCase.nativePath, '-struct', 'data', '-append');
            actual = sixth('EN4c14', days(10));
            testCase.verifyEqual(actual, expected);
            native = sixth('EN4c14');
            testCase.verifyEqual(native, 9.75+(0:59)'/40, AbsTol=1e-12);
        end
        function forceWithoutHaloRebuildsAvailableMeans(testCase)
            fields = string(fieldnames(testCase.fixtureData));
            data = rmfield(testCase.fixtureData, cellstr(fields(contains(lower(fields), 'halosteric'))));
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
            raw = steric2lonlatt('EN4c14', days(10), ForceNew=true, BeQuiet=true);
            stored = load(testCase.nativePath);
            testCase.verifyEqual(raw, data.stericSl(:, :, 1:10:end));
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Gmsl')), 9);
            testCase.verifyFalse(isfield(stored, 'halostericGmsl'));
        end
        function graceAndFixedStepMeanDiagnostics(testCase)
            folder = fileparts(fileparts(fileparts(testCase.nativePath)));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', folder));
            reference = fullfile(folder, 'RL06', 'CSR'); mkdir(reference);
            writelines(["# End of YAML header"; ...
                "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], ...
                fullfile(reference, 'GSM-2_test_BA01_0600'));
            [~, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', 'GRACE', BeQuiet=true);
            testCase.verifyEqual(dates, datetime(2020, [1;2], [16;15], 12, 0, 0));
            testCase.verifyEqual(halo, [mean(9.75+(0:30)/40); mean(9.75+(31:59)/40)], AbsTol=1e-12);
            [~, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', days(10), Interpolation='mean', BeQuiet=true);
            epochs = temporalmeanepochs(testCase.fixtureData.dates, days(10), []);
            reference = averageepochs(testCase.fixtureData.dates, reshape(9.75+(0:59)/40, 1, 1, []), epochs);
            testCase.verifyEqual(dates, testCase.fixtureData.dates(1:10:end));
            testCase.verifyEqual(halo, reference(:), AbsTol=1e-12);
        end
        function incompleteDerivedMeansRebuiltTogether(testCase)
            sixth('EN4c14', days(10));
            path = stericfilepath('EN4c14', days(10), [], [], 'linear');
            data = load(path); data.halostericGmsl(:) = -999;
            data = rmfield(data, 'deepNonHalostericGmsl');
            save(path, '-struct', 'data', '-v7.3');
            actual = sixth('EN4c14', days(10));
            testCase.verifyEqual(actual, 9.75+(0:10:59)'/40, AbsTol=1e-12);
            testCase.verifyTrue(ismember('deepNonHalostericGmsl', who('-file', path)));
        end
        function noShallowAliasOrLayerReconstruction(testCase)
            data = rmfield(testCase.fixtureData, {'deepHalostericSl'});
            save(testCase.nativePath, '-struct', 'data');
            testCase.verifyError(@() sixth('EN4c14'), 'ULMO:steric2lonlatt:InvalidSchema');
            testCase.verifyError(@() stericproductinfo('EN4c14'), 'ULMO:stericproductinfo:InvalidSchema');
            names = string(fieldnames(testCase.fixtureData));
            data = rmfield(testCase.fixtureData, cellstr(names(startsWith(names, ["shallow", "deep"]))));
            path = stericfilepath('SIO'); mkdir(fileparts(path)); save(path, '-struct', 'data', '-v7.3');
            testCase.verifyError(@() steric2lonlatt('SIO', Depth='shallow'), 'ULMO:steric2lonlatt:UnsupportedField');
            testCase.verifyEqual(sixth('SIO'), 9.75+(0:59)'/40, AbsTol=1e-12);
        end
    end
    methods
        function writeNative(testCase)
            data = testCase.fixtureData;
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
        end
    end
end
function halo = sixth(varargin)
    [~, ~, ~, ~, ~, halo] = steric2lonlatt(varargin{:}, BeQuiet=true);
end
function value = bytes(path)
    file = fopen(path, 'r'); cleanup = onCleanup(@() fclose(file));
    value = fread(file, Inf, '*uint8');
end
