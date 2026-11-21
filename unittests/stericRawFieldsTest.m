%% STERICRAWFIELDSTEST - Raw loader, temporal GMSL and cache contract tests.
% Uses temporary native grids and isolated IFILES. No production data is read.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/20, En-Chi Lee (williameclee@gmail.com)

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
            base = [1, 2, 3; 4, 5, 6; 7, 8, 9; 10, 11, 12] + reshape((0:59) / 60, 1, 1, []);
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
                Unit = 'mm', TimeFormat = 'datenum', OutputFormat = 'ndgrid', BeQuiet = true);
            testCase.verifyEqual(raw, 1000 * permute(testCase.fixtureData.stericSl, [2, 1, 3]));
            testCase.verifyEqual(sigma, zeros(size(raw)));
            testCase.verifyEqual(dates, datenum(testCase.fixtureData.dates)); %#ok<DATNM>
            testCase.verifyEqual(lon, testCase.fixtureData.lon');
            testCase.verifyEqual(lat, testCase.fixtureData.lat');
            testCase.verifyEqual(halo, 1000 * (9.75 + (0:59)' / 40), AbsTol = 1e-11);
            testCase.verifyEqual(steric2lonlatt('EN4c14', Component = 'total', ...
                Depth = 'full-depth', BeQuiet = true), testCase.fixtureData.stericSl);
        end

        function removedComponents(testCase)

            for component = ["default", "full", "full-depth", "no halosteric drift"]
                testCase.verifyError(@() steric2lonlatt('EN4c14', Component = component), 'MATLAB:validators:mustBeMember');
            end

        end

        function componentLabelsAreNotDepthLabels(testCase)

            for depth = ["full", "total"]
                testCase.verifyError(@() steric2lonlatt('EN4c14', Depth = depth), ...
                'MATLAB:validators:mustBeMember');
                testCase.verifyError(@() steric2gmsl('EN4c14', 'total', depth), ...
                'MATLAB:validators:mustBeMember');
            end

        end

        function upperOceanProductsRejectFullDepthDespiteCachedFields(testCase)
            % Field names and cached values must not override product coverage.
            for product = ["SIO", "NOAA", "JAMSTEC", "ISAS", "IPRC"]
                data = testCase.fixtureData;
                fields = string(fieldnames(data));
                data = rmfield(data, cellstr(fields(startsWith(fields, ["shallow", "deep"]))));
                data.stericGmsl = -999 * ones(size(data.dates));
                path = stericfilepath(product);
                mkdir(fileparts(path));
                save(path, '-struct', 'data', '-v7.3');
                data.steric = data.stericSl;
                cache = stericfilepath(product, 'midmonth', [], [], 'linear');
                save(cache, '-struct', 'data', '-v7.3');
                info = stericproductinfo(product);
                testCase.verifyFalse(info.hasFullDepth);
                testCase.verifyEqual(info.nativeDepth, "shallow");
                testCase.verifyTrue(info.hasHalosteric);
                testCase.verifyEqual(steric2lonlatt(product, BeQuiet = true), data.stericSl);
                testCase.verifyError(@() steric2lonlatt(product, Depth = 'full-depth', BeQuiet = true), ...
                'ULMO:steric2lonlatt:UnsupportedField');
                testCase.verifyError(@() steric2lonlatt(product, 'midmonth', Depth = 'full-depth', BeQuiet = true), ...
                'ULMO:steric2lonlatt:UnsupportedField');
                testCase.verifyError(@() sixth(product, [], [], [], 'Depth', 'full-depth'), ...
                'ULMO:steric2lonlatt:UnsupportedField');
                testCase.verifyError(@() steric2gmsl(product, 'total', 'full-depth', BeQuiet = true), ...
                'ULMO:steric2gmsl:UnsupportedField');
                testCase.verifyEqual(steric2lonlatt(product, Depth = 'shallow', BeQuiet = true), ...
                    data.stericSl);
                testCase.verifyError(@() steric2lonlatt(product, Depth = 'deep', BeQuiet = true), ...
                'ULMO:steric2lonlatt:UnsupportedField');
                testCase.verifyError(@() steric2gmsl(product, 'halosteric', 'deep'), ...
                'ULMO:steric2gmsl:UnsupportedField');
                shallow = steric2gmsl(product, 'halosteric', 'shallow', BeQuiet = true);
                testCase.verifyEqual(shallow, 9.75 + (0:59)' / 40, AbsTol = 1e-12);
                testCase.verifyEqual(steric2gmsl(product, 'halosteric', BeQuiet = true), shallow);
                testCase.verifyEqual(sixth(product), shallow);
                testCase.verifyEqual(sixth(product, [], [], [], 'Depth', 'shallow'), ...
                    shallow, AbsTol = 1e-12);
            end

        end

        function preparedFullDepthProductsAcceptExplicitDepth(testCase)

            for product = ["EN4c14", "CMEMS", "ORAS5", "IAPv4"]
                data = testCase.fixtureData;
                path = stericfilepath(product);
                if ~isfolder(fileparts(path)), mkdir(fileparts(path)); end
                save(path, '-struct', 'data', '-v7.3');
                info = stericproductinfo(product);
                testCase.verifyTrue(info.hasFullDepth);
                testCase.verifyEqual(info.nativeDepth, "full-depth");
                testCase.verifyTrue(info.hasHalosteric);
                testCase.verifyEqual(steric2lonlatt(product, Component = 'halosteric', ...
                    Depth = "full-depth", BeQuiet = true), data.halostericSl);
                actual = steric2gmsl(product, 'halosteric', "full-depth", SaveData = false, BeQuiet = true);
                testCase.verifyEqual(actual, 9.75 + (0:59)' / 40, AbsTol = 1e-12);
            end

        end

        function nativeAndShallowShareEpochCache(testCase)
            data = testCase.fixtureData;
            fields = string(fieldnames(data));
            data = rmfield(data, cellstr(fields(startsWith(fields, ["shallow", "deep"]))));
            path = stericfilepath('SIO'); mkdir(fileparts(path));
            save(path, '-struct', 'data', '-v7.3');
            [raw, ~, ~, ~, ~, nativeHalo] = steric2lonlatt('SIO', 'midmonth', BeQuiet = true);
            cache = stericfilepath('SIO', 'midmonth', [], [], 'linear');
            stored = load(cache);
            testCase.verifyEqual(raw, stored.shallowStericSl);
            testCase.verifyEqual(nativeHalo, stored.shallowHalostericGmsl);
            testCase.verifyFalse(isfield(stored, 'halostericGmsl'));
            shallowStericSl = raw + 123;
            save(cache, 'shallowStericSl', '-append');
            delete(path);
            testCase.verifyEqual(steric2lonlatt('SIO', 'midmonth', Depth = 'shallow', BeQuiet = true), shallowStericSl);
            testCase.verifyEqual(steric2lonlatt('SIO', 'midmonth', Depth = 'native', BeQuiet = true), shallowStericSl);
        end

        function rawOnlyWithoutHalo(testCase)
            data = rmfield(testCase.fixtureData, {'halostericSl', 'shallowHalostericSl', 'deepHalostericSl'});
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
            raw = steric2lonlatt('EN4c14', BeQuiet = true);
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
            [raw, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', 'midmonth', BeQuiet = true);
            nativeHalo = steric2gmsl('EN4c14', 'halosteric', BeQuiet = true);
            testCase.verifyEqual(dates, datetime(2020, [1; 2], [16; 15], 12, 0, 0));
            testCase.verifyEqual(raw(:, :, 1), mean(testCase.fixtureData.stericSl(:, :, 1:31), 3, 'omitmissing'));
            testCase.verifyEqual(halo(1), mean(nativeHalo(1:31)), AbsTol = 1e-12);
            before = mean(testCase.fixtureData.stericSl(1, 1, 1:31) - reshape(nativeHalo(1:31), 1, 1, []), 3, 'omitmissing');
            testCase.verifyGreaterThan(abs(raw(1, 1, 1) - halo(1) - before), 0.01);
        end

        function allMeansAndFieldsInEpochCache(testCase)
            sixth('EN4c14', 'midmonth');
            path = stericfilepath('EN4c14', 'midmonth', [], [], 'linear');
            stored = load(path);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Gmsl')), 9);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Sl')), 9);
            testCase.verifyFalse(isfield(stored, 'steric'));
            testCase.verifyFalse(any(contains(lower(string(fieldnames(stored))), ...
                ["nonthermosteric", "nonhalosteric"])));
            testCase.verifySize(stored.halostericGmsl, [2, 1]);
            testCase.verifyEqual(stored.halostericGmsl, stored.shallowHalostericGmsl + stored.deepHalostericGmsl);
        end

        function selectionsShareCompleteRecordCache(testCase)
            limits = datetime(2020, 1, [1, 31]);
            [~, ~, dates] = steric2lonlatt('EN4c14', 'midmonth', [], limits, ...
                Component = 'thermosteric', Depth = 'shallow', BeQuiet = true);
            testCase.verifyEqual(dates, datetime(2020, 1, 16, 12, 0, 0));
            path = stericfilepath('EN4c14', 'midmonth', [], [], 'linear');
            stored = load(path);
            testCase.verifySize(stored.dates, [2, 1]);
            before = bytes(path);
            delete(testCase.nativePath);
            schema = stericfieldschema();

            for iDepth = 1:3
                grids = {testCase.fixtureData.(schema.gridNames(1, iDepth)), ...
                             testCase.fixtureData.(schema.gridNames(2, iDepth)), ...
                             testCase.fixtureData.(schema.gridNames(3, iDepth))};
                grids{4} = grids{1} - grids{2};
                grids{5} = grids{1} - grids{3};

                for iComponent = 1:5
                    [actual, ~, returnedDates] = steric2lonlatt('EN4c14', 'midmonth', ...
                        Component = schema.components(iComponent), Depth = schema.depths(iDepth), BeQuiet = true);
                    expected = cat(3, mean(grids{iComponent}(:, :, 1:31), 3), ...
                        mean(grids{iComponent}(:, :, 32:60), 3));
                    testCase.verifyEqual(actual, expected, AbsTol = 1e-12);
                    testCase.verifyEqual(returnedDates, stored.dates);
                end

            end

            testCase.verifyEqual(bytes(path), before);
            testCase.verifyNumElements(dir(fullfile(fileparts(path), '*.mat')), 1);
        end

        function differenceUsesProcessedEpochMeans(testCase)

            for field = ["halostericSl", "shallowHalostericSl", "deepHalostericSl"]
                testCase.fixtureData.(field)(1, 1, 1:10) = NaN;
            end

            testCase.writeNative;
            totalMean = steric2lonlatt('EN4c14', 'midmonth', Depth = 'deep', BeQuiet = true);
            actual = steric2lonlatt('EN4c14', 'midmonth', ...
                Component = 'non-halosteric', Depth = 'deep', BeQuiet = true);
            total = testCase.fixtureData.deepStericSl(:, :, 1:31);
            halo = testCase.fixtureData.deepHalostericSl(:, :, 1:31);
            expected = mean(total - halo, 3, 'omitmissing');
            separate = totalMean(:, :, 1) - mean(halo, 3, 'omitmissing');
            testCase.verifyEqual(actual(:, :, 1), separate, AbsTol = 1e-12);
            testCase.verifyGreaterThan(abs(expected(1, 1) - separate(1, 1)), 0.01);
        end

        function missingGridRebuildsCompleteCache(testCase)
            original = steric2lonlatt('EN4c14', 'midmonth', BeQuiet = true);
            path = stericfilepath('EN4c14', 'midmonth', [], [], 'linear');
            stored = load(path);
            stored.stericSl(:) = -999;
            stored = rmfield(stored, 'deepHalostericSl');
            save(path, '-struct', 'stored', '-v7.3');
            before = bytes(path);
            actual = steric2lonlatt('EN4c14', 'midmonth', SaveData = false, BeQuiet = true);
            testCase.verifyEqual(actual, original);
            testCase.verifyEqual(bytes(path), before);
            actual = steric2lonlatt('EN4c14', 'midmonth', BeQuiet = true);
            testCase.verifyEqual(actual, original);
            stored = load(path);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Sl')), 9);
        end

        function cachedMeansRemainAvailableWithoutNativeFile(testCase)
            sixth('EN4c14', 'midmonth');
            path = stericfilepath('EN4c14', 'midmonth', [], [], 'linear');
            stored = load(path);
            delete(testCase.nativePath);
            limits = datetime(2020, 2, [1, 29]);
            [raw, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', 'midmonth', [], limits, ...
                Component = 'non-thermosteric', Depth = 'deep', BeQuiet = true);
            testCase.verifyEqual(raw, stored.deepStericSl(:, :, 2) - stored.deepThermostericSl(:, :, 2));
            testCase.verifyEqual(dates, stored.dates(2));
            testCase.verifyEqual(halo, stored.deepHalostericGmsl(2));
        end

        function differencesUseProcessedInterpolationAndIgnoreOldFields(testCase)
            % Nonlinear interpolation must finish before component subtraction.
            samples = [1, 16, 31, 46, 60];
            data = testCase.fixtureData;
            data.dates = data.dates(samples);
            fields = string(fieldnames(data));
            data = rmfield(data, cellstr(fields(startsWith(fields, ["shallow", "deep"]))));
            data.stericSl = ones(4, 3) .* reshape([10, 20, 30, 40, 50], 1, 1, []);
            data.thermostericSl = ones(4, 3) .* reshape([1, 8, 2, 9, 3], 1, 1, []);
            data.halostericSl = ones(4, 3) .* reshape([1, 4, 2, 8, 3], 1, 1, []);
            testCase.fixtureData = data;
            testCase.writeNative;
            args = {'EN4c14', days(5), [], [], 'Interpolation', 'pchip', 'BeQuiet', true};
            total = steric2lonlatt(args{:});
            thermo = steric2lonlatt(args{:}, Component = 'thermosteric');
            actual = steric2lonlatt(args{:}, Component = 'non-thermosteric');
            testCase.verifyEqual(actual, total - thermo);
            before = interp1(data.dates, squeeze(data.stericSl(1, 1, :) - ...
                data.thermostericSl(1, 1, :)), data.dates(1):days(5):data.dates(end), 'pchip');
            testCase.verifyGreaterThan(max(abs(squeeze(actual(1, 1, :)) - before(:))), 0.01);
            path = stericfilepath('EN4c14', days(5), [], [], 'pchip');
            testCase.verifyFalse(any(contains(lower(string(who('-file', path))), ...
                ["nonthermosteric", "nonhalosteric"])));
            nonThermostericSl = -999 * ones(size(actual));
            save(path, 'nonThermostericSl', '-append');
            delete(testCase.nativePath);
            testCase.verifyEqual(steric2lonlatt(args{:}, Component = 'non-thermosteric'), total - thermo);
        end

        function nonlinearDepthAdditivity(testCase)
            dates = testCase.fixtureData.dates([1, 16, 31, 46, 60]);
            testCase.fixtureData.dates = dates;
            names = fieldnames(testCase.fixtureData);

            for iName = 1:numel(names)
                if endsWith(names{iName}, 'Sl'), testCase.fixtureData.(names{iName}) = testCase.fixtureData.(names{iName})(:, :, [1, 16, 31, 46, 60]); end
            end

            testCase.fixtureData.shallowHalostericSl = ones(4, 3) .* reshape([1, 4, 2, 8, 3], 1, 1, []);
            testCase.fixtureData.deepHalostericSl = ones(4, 3) .* reshape([5, 2, 7, 1, 9], 1, 1, []);
            testCase.fixtureData.halostericSl = testCase.fixtureData.shallowHalostericSl + testCase.fixtureData.deepHalostericSl;
            testCase.writeNative;
            args = {'EN4c14', days(5), [], [], 'Interpolation', 'pchip', 'SaveData', false};
            fullDepth = sixth(args{:});
            shallow = sixth(args{:}, 'Depth', 'shallow');
            deep = sixth(args{:}, 'Depth', 'deep');
            testCase.verifyEqual(fullDepth, shallow + deep);
            expected = interp1(dates, [1; 4; 2; 8; 3], dates(1):days(5):dates(end), 'pchip') + ...
                interp1(dates, [5; 2; 7; 1; 9], dates(1):days(5):dates(end), 'pchip');
            testCase.verifyEqual(fullDepth, expected(:), AbsTol = 1e-12);
        end

        function spatialGridDoesNotChangeMean(testCase)
            original = sixth('EN4c14');
            spatial = sixth('EN4c14', [], 60);
            testCase.verifyEqual(spatial, original);
        end

        function filteredDatesUseCompleteCache(testCase)
            limits = datetime(2020, 1, [10, 30]);
            [~, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', days(7), [], limits, ...
                Unit = 'mm', BeQuiet = true);
            testCase.verifyEqual(dates, datetime(2020, 1, [15; 22; 29]));
            testCase.verifyEqual(halo, 1000 * (9.75 + [14; 21; 28] / 40), AbsTol = 1e-11);
        end

        function augmentRawCacheWithoutGridRebuild(testCase)
            raw = steric2lonlatt('EN4c14', days(10), BeQuiet = true); %#ok<NASGU>
            path = stericfilepath('EN4c14', days(10), [], [], 'linear');
            cached = load(path);
            % A sentinel proves cached raw grids remain in use during augmentation.
            stericSl = cached.stericSl + 123;
            save(path, 'stericSl', '-append');
            [raw, ~, ~, ~, ~, halo] = steric2lonlatt('EN4c14', days(10), BeQuiet = true);
            testCase.verifyEqual(raw, stericSl);
            testCase.verifyEqual(halo, 9.75 + (0:10:59)' / 40, AbsTol = 1e-12);
            stored = load(path);
            testCase.verifyEqual(stored.stericSl, stericSl);
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Gmsl')), 9);
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
            [raw, ~, ~, ~, ~] = steric2lonlatt('EN4c14', 'midmonth', ForceNew = true, BeQuiet = true);
            stored = load(stericfilepath('EN4c14', 'midmonth', [], [], 'linear'));
            native = load(testCase.nativePath);
            testCase.verifyEqual(native.halostericGmsl, 2 * old.halostericGmsl, AbsTol = 1e-12);
            testCase.verifyEqual(native.deepThermostericGmsl, 2 * old.deepThermostericGmsl, AbsTol = 1e-12);
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
            cache = stericfilepath('EN4c14', 'midmonth', [], [], 'linear');
            beforeNative = bytes(testCase.nativePath); beforeCache = bytes(cache);
            actual = sixth('EN4c14', 'midmonth', [], [], 'ForceNew', true, 'SaveData', false);
            expected = 2 * [mean(9.75 + (0:30) / 40); mean(9.75 + (31:59) / 40)];
            testCase.verifyEqual(actual, expected, AbsTol = 1e-12);
            testCase.verifyEqual(bytes(testCase.nativePath), beforeNative);
            testCase.verifyEqual(bytes(cache), beforeCache);
        end

        function selectedFieldCacheIsRebuilt(testCase)
            data = testCase.fixtureData;
            data.steric = -999 * ones(4, 3, 2);
            data = rmfield(data, 'stericSl');
            data.dates = datetime(2020, [1; 2], [16; 15], 12, 0, 0);
            legacy = stericfilepath('EN4c14', 'midmonth', [], [], 'linear');
            save(legacy, '-struct', 'data');
            raw = steric2lonlatt('EN4c14', 'midmonth', SaveData = false, BeQuiet = true);
            testCase.verifyEqual(raw(:, :, 1), mean(testCase.fixtureData.stericSl(:, :, 1:31), 3));
        end

        function rawComponentsAndRequestedDepth(testCase)
            components = ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"];
            expected = {testCase.fixtureData.deepStericSl, testCase.fixtureData.deepThermostericSl, ...
                            testCase.fixtureData.deepHalostericSl, ...
                            testCase.fixtureData.deepStericSl - testCase.fixtureData.deepThermostericSl, ...
                            testCase.fixtureData.deepStericSl - testCase.fixtureData.deepHalostericSl};

            for iComponent = 1:numel(components)
                [raw, ~, ~, ~, ~, halo] = steric2lonlatt('EN4c14', ...
                    Component = components(iComponent), Depth = 'deep', BeQuiet = true);
                testCase.verifyEqual(raw, expected{iComponent});
                testCase.verifyEqual(halo, 3.25 + (0:59)' / 120, AbsTol = 1e-12);
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
            testCase.verifyEqual(native, 9.75 + (0:59)' / 40, AbsTol = 1e-12);
        end

        function forceWithoutHaloRebuildsAvailableMeans(testCase)
            fields = string(fieldnames(testCase.fixtureData));
            data = rmfield(testCase.fixtureData, cellstr(fields(contains(lower(fields), 'halosteric'))));
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
            raw = steric2lonlatt('EN4c14', days(10), ForceNew = true, BeQuiet = true);
            stored = load(testCase.nativePath);
            testCase.verifyEqual(raw, data.stericSl(:, :, 1:10:end));
            testCase.verifyEqual(nnz(endsWith(fieldnames(stored), 'Gmsl')), 6);
            testCase.verifyFalse(isfield(stored, 'halostericGmsl'));
        end

        function graceAndFixedStepMeanDiagnostics(testCase)
            folder = fileparts(fileparts(fileparts(testCase.nativePath)));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', folder));
            reference = fullfile(folder, 'RL06', 'CSR'); mkdir(reference);
            writelines(["# End of YAML header"; ...
                        "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], ...
                fullfile(reference, 'GSM-2_test_BA01_0600'));
            [~, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', 'GRACE', BeQuiet = true);
            testCase.verifyEqual(dates, datetime(2020, [1; 2], [16; 15], 12, 0, 0));
            testCase.verifyEqual(halo, [mean(9.75 + (0:30) / 40); mean(9.75 + (31:59) / 40)], AbsTol = 1e-12);
            [~, ~, dates, ~, ~, halo] = steric2lonlatt('EN4c14', days(10), Interpolation = 'mean', BeQuiet = true);
            epochs = temporalmeanepochs(testCase.fixtureData.dates, days(10), []);
            reference = averageepochs(testCase.fixtureData.dates, reshape(9.75 + (0:59) / 40, 1, 1, []), epochs);
            testCase.verifyEqual(dates, testCase.fixtureData.dates(1:10:end));
            testCase.verifyEqual(halo, reference(:), AbsTol = 1e-12);
        end

        function incompleteDerivedMeansRebuiltTogether(testCase)
            sixth('EN4c14', days(10));
            path = stericfilepath('EN4c14', days(10), [], [], 'linear');
            data = load(path); data.halostericGmsl(:) = -999;
            data = rmfield(data, 'deepHalostericGmsl');
            save(path, '-struct', 'data', '-v7.3');
            actual = sixth('EN4c14', days(10));
            testCase.verifyEqual(actual, 9.75 + (0:10:59)' / 40, AbsTol = 1e-12);
            testCase.verifyTrue(ismember('deepHalostericGmsl', who('-file', path)));
        end

        function nativeDepthResolutionDoesNotReconstructLayers(testCase)
            data = rmfield(testCase.fixtureData, {'deepHalostericSl'});
            save(testCase.nativePath, '-struct', 'data');
            testCase.verifyError(@() sixth('EN4c14'), 'ULMO:steric2lonlatt:InvalidSchema');
            testCase.verifyError(@() stericproductinfo('EN4c14'), 'ULMO:stericproductinfo:InvalidSchema');
            names = string(fieldnames(testCase.fixtureData));
            data = rmfield(testCase.fixtureData, cellstr(names(startsWith(names, ["shallow", "deep"]))));
            path = stericfilepath('SIO'); mkdir(fileparts(path)); save(path, '-struct', 'data', '-v7.3');
            testCase.verifyEqual(steric2lonlatt('SIO', Depth = 'shallow', BeQuiet = true), data.stericSl);
            testCase.verifyEqual(sixth('SIO'), 9.75 + (0:59)' / 40, AbsTol = 1e-12);
            testCase.verifyError(@() steric2lonlatt('SIO', Depth = 'full-depth'), 'ULMO:steric2lonlatt:UnsupportedField');
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
    [~, ~, ~, ~, ~, halo] = steric2lonlatt(varargin{:}, BeQuiet = true);
end

function value = bytes(path)
    file = fopen(path, 'r'); cleanup = onCleanup(@() fclose(file));
    value = fread(file, Inf, '*uint8');
end
