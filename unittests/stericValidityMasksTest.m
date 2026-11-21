%% STERICVALIDITYMASKSTEST - Tests shared native steric validity masks.
% Uses temporary grids covering observed shallow water, intermittent deep
% observations, finite zeros and cells outside the product coverage.
%
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

classdef stericValidityMasksTest < matlab.unittest.TestCase

    properties
        data
        nativePath
        deepPath
    end

    methods (TestMethodSetup)

        function prepare(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.nativePath = stericfilepath('EN4c14');
            [parent, name] = fileparts(testCase.nativePath);
            testCase.deepPath = fullfile(parent, [name, '-deep.mat']);
            mkdir(parent);

            data.lon = [0, 90, 180, 270];
            data.lat = [-45; 45];
            data.dates = datetime(2020, 1, (1:2:9)');
            data.shallowStericSl = repmat([2, 4, 8, NaN; 3, 4, 5, 6], 1, 1, 5);
            data.deepStericSl = repmat([NaN, 6, 0, NaN; 2, 2, 2, 2], 1, 1, 5);
            data.stericSl = data.shallowStericSl + data.deepStericSl;
            data.stericSl(1, 1, :) = data.shallowStericSl(1, 1, :);
            data.deepStericSl(1, 2, 2) = NaN;
            data.shallowHalostericSl = data.shallowStericSl / 2;
            data.deepHalostericSl = data.deepStericSl / 2;
            data.halostericSl = data.stericSl / 2;

            testCase.data = data;
            testCase.writeNative();
        end

    end

    methods (Test)

        function nativeAndResidualDeepGrids(testCase)
            % A finite zero in any component establishes deep coverage.
            data = testCase.data;
            data.deepStericSl(1, 3, 1:4) = NaN;
            data.deepHalostericSl(1, 3, :) = NaN;
            data.deepStericSl(1, 4, 1) = 9; % Deep value outside the product ocean mask.
            data.deepStericSl(2, 1, 1) = Inf;
            data.deepHalostericSl(2, 2, 1) = NaN;
            data.shallowStericSl(1, 1, 3) = NaN;
            data.shallowStericSl(1, 2, 3) = NaN;
            testCase.data = data;
            testCase.writeNative();

            deep = steric2lonlatt('EN4c14', Depth = 'deep', BeQuiet = true);
            testCase.verifyEqual(reshape(deep(1, 1, :), [], 1), [0; 0; NaN; 0; 0]);
            testCase.verifyEqual(reshape(deep(1, 2, :), [], 1), [6; NaN; 6; 6; 6]);
            testCase.verifyEqual(reshape(deep(1, 3, :), [], 1), [NaN; NaN; NaN; NaN; 0]);
            testCase.verifyTrue(all(isnan(deep(1, 4, :)), 'all'));
            testCase.verifyEqual(deep(2, 1, 1), Inf);
            residual = steric2lonlatt('EN4c14', Depth = 'deep', Component = 'non-halosteric', BeQuiet = true);
            testCase.verifyEqual(residual(1, 1, 1), 0);
            testCase.verifyEqual(residual(1, 2, 1), 3);
            testCase.verifyTrue(isnan(residual(2, 2, 1)));

            cache = load(testCase.deepPath);
            testCase.verifyEqual(cache.oceanMask, logical([1, 1, 1, 0; 1, 1, 1, 1]));
            testCase.verifyEqual(cache.deepOceanMask, logical([0, 1, 1, 1; 1, 1, 1, 1]));
            testCase.verifyEqual(cache.nativeLon, data.lon);
            testCase.verifyEqual(cache.nativeLat, data.lat);
            testCase.verifyEqual(cache.deepStericSl, deep);
            testCase.verifyEqual(load(testCase.nativePath), data);
        end

        function commonAreaForAllDepthMeans(testCase)
            % All positive-weight cells have the same area in this fixture.
            % A shallow-water zero retains its weight while the intermittent
            % deep value removes the whole cell from the second epoch.
            request = ["full-depth", "shallow", "deep"];
            means = steric2gmsl('EN4c14', repmat("total", size(request)), request, BeQuiet = true);
            weights = [7; 6; 7; 7; 7];
            sums = [32, 28, 32, 32, 32; 14, 8, 14, 14, 14]';
            testCase.verifyEqual(means(:, 2:3), sums ./ weights, AbsTol = 1e-12);
            testCase.verifyEqual(means(:, 1), sum(means(:, 2:3), 2));
            cached = steric2gmsl('EN4c14', repmat("total", size(request)), request, BeQuiet = true);
            testCase.verifyEqual(cached, means);

            stored = load(testCase.nativePath);
            original = fieldnames(testCase.data);
            for iField = 1:numel(original)
                name = original{iField};
                testCase.verifyEqual(stored.(name), testCase.data.(name));
            end
        end

        function fillingPrecedesEpochAveraging(testCase)
            % At this cell each source field has monthly observations, but
            % they occur on different days. The deep grid has no valid zero.
            data = testCase.data;
            data.dates(end) = datetime(2020, 1, 20);
            data.stericSl(1, 1, 2:5) = NaN;
            data.shallowStericSl(1, 1, 1) = NaN;
            testCase.data = data;
            testCase.writeNative();

            deep = steric2lonlatt('EN4c14', 'midmonth', Depth = 'deep', BeQuiet = true);
            full = steric2lonlatt('EN4c14', 'midmonth', Depth = 'full-depth', BeQuiet = true);
            shallow = steric2lonlatt('EN4c14', 'midmonth', Depth = 'shallow', BeQuiet = true);
            testCase.verifyEqual(full(1, 1), 2);
            testCase.verifyEqual(shallow(1, 1), 2);
            testCase.verifyTrue(isnan(deep(1, 1)));
            testCase.verifyEqual(deep(1, 2), 6);

            delete(testCase.nativePath);
            cached = steric2lonlatt('EN4c14', 'midmonth', Depth = 'deep', BeQuiet = true);
            testCase.verifyEqual(cached, deep);
        end

        function fillingPrecedesInterpolation(testCase)
            % Interpolate between a shallow-water zero and an observed deep
            % value, and interpolate in time over the constant zero series.
            [spatial, ~, ~, lon, lat] = steric2lonlatt('EN4c14', [], 45, Depth = 'deep', BeQuiet = true);
            testCase.verifyEqual(spatial(lat == -45, lon == 45, 1), 3, AbsTol = 1e-12);

            data = testCase.data;
            data.deepStericSl(2, 1, 2:5) = NaN;
            testCase.data = data;
            testCase.writeNative();

            for method = ["linear", "pchip", "spline", "makima"]
                [temporal, ~, dates] = steric2lonlatt('EN4c14', days(1), ...
                    Depth = 'deep', Interpolation = method, BeQuiet = true);
                testCase.verifyEqual(reshape(temporal(1, 1, :), [], 1), zeros(numel(dates), 1));
                testCase.verifyEqual(dates(2), testCase.data.dates(1) + days(1));
                testCase.verifyTrue(all(isnan(temporal(1, 4, :)), 'all'));
                testCase.verifyEqual(temporal(2, 1, 1), 2);
                testCase.verifyTrue(all(isnan(temporal(2, 1, 2:end)), 'all'));
                processedPath = stericfilepath('EN4c14', days(1), [], [], method);
                cache = load(processedPath, 'oceanMask', 'deepOceanMask', 'creationDate');
                testCase.verifyTrue(isdatetime(cache.creationDate));
                testCase.verifyFalse(cache.deepOceanMask(1, 1));
            end
        end

        function nativeCacheReadsAndForcedRefresh(testCase)
            first = steric2lonlatt('EN4c14', Depth = 'deep', SaveData = false, BeQuiet = true);
            testCase.verifyFalse(isfile(testCase.deepPath));
            saved = steric2lonlatt('EN4c14', Depth = 'deep', BeQuiet = true);
            testCase.verifyEqual(saved, first);

            creationDate = datetime(2000, 1, 1);
            save(testCase.deepPath, 'creationDate', '-append');
            oldCache = load(testCase.deepPath);
            data = testCase.data;
            data.deepStericSl(1, 1, :) = 1;
            data.stericSl(1, 1, :) = 3;
            testCase.data = data;
            testCase.writeNative();
            testCase.verifyEqual(steric2lonlatt('EN4c14', Depth = 'deep', BeQuiet = true), saved);

            fresh = steric2lonlatt('EN4c14', Depth = 'deep', ...
                ForceNew = true, SaveData = false, BeQuiet = true);
            testCase.verifyEqual(fresh(1, 1, 1), 1);
            testCase.verifyEqual(load(testCase.deepPath), oldCache);
            testCase.verifyEqual(load(testCase.nativePath), data);

            refreshed = steric2lonlatt('EN4c14', Depth = 'deep', ForceNew = true, BeQuiet = true);
            cache = load(testCase.deepPath);
            testCase.verifyEqual(refreshed, fresh);
            testCase.verifyGreaterThan(cache.creationDate, oldCache.creationDate);
            testCase.verifyTrue(cache.deepOceanMask(1, 1));
            testCase.verifyTrue(all(isnan(cache.deepHalostericSl(1, 1, :)), 'all'));
        end

    end

    methods (Access = private)

        function writeNative(testCase)
            data = testCase.data;
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
        end

    end
end
