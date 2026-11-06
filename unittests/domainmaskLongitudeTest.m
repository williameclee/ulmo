%% DOMAINMASKLONGITUDETEST - Verify masks for equivalent geographic coordinates.
% Tests use temporary mask caches and do not require steric observations.
%
% Created by
%   2026/11/05, OpenAI Codex

classdef domainmaskLongitudeTest < matlab.unittest.TestCase
    properties (TestParameter)
        buffer = {0, 1.5}
        shift = {-720, -360, 0, 360, 720}
        basin = {'oceans', 'npacific', 'spacific', 'natlantic', 'satlantic', 'indian'}
    end
    methods (TestMethodSetup)
        function temporaryRoot(t)
            folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
        end
    end

    methods (Test)
        function shiftedSouthAtlanticGrid(t, buffer)
            % Issue #101: the same native grid must select the same ocean cells.
            domain = GeoDomain('satlantic', Buffer=0.5, ...
                DefaultParams=true, Latlim=[-60, 60]);
            lon = 20.5:379.5;
            lat = (-64.5:79.5).';
            shifted = domainmask(domain, lon, lat, ...
                exteriorBuffer=buffer, Forcenew=true, SaveData=false, BeQuiet=true);
            [wrappedLon, order] = sort(mod(lon + 180, 360) - 180);
            wrapped = domainmask(domain, wrappedLon, lat, ...
                exteriorBuffer=buffer, Forcenew=true, SaveData=false, BeQuiet=true);

            t.verifyEqual(nnz(xor(shifted(:, order), wrapped)), 0, ...
                'Equivalent longitude grids must select identical physical cells.');
            t.verifyTrue(shifted(lat == -30.5, lon == 370.5));
        end

        function basinGrids(t, basin, buffer)
            % All six domains select the same physical points in both windows.
            domain = GeoDomain(basin, Buffer=0.5, DefaultParams=true);
            lon = 20.5:4:376.5;
            lat = -59.5:4:76.5;
            actual = domainmask(domain, lon, lat, exteriorBuffer=buffer, ...
                SaveData=false, BeQuiet=true);
            [wrapped, order] = sort(mod(lon + 180, 360) - 180);
            expected = domainmask(domain, wrapped, lat, exteriorBuffer=buffer, ...
                SaveData=false, BeQuiet=true);
            t.verifyEqual(actual(:, order), expected);
        end

        function quarterDegreeGrid(t)
            % The shifted 0.25-degree grid must include the eastern Atlantic.
            domain = GeoDomain('satlantic', Buffer=0.5, DefaultParams=true);
            lon = 350:0.25:380;
            lat = -50:0.25:-30;
            [actual, x, y] = domainmask(domain, lon, lat, SaveData=false, BeQuiet=true);
            expected = domainmask(domain, lon-360, lat, SaveData=false, BeQuiet=true);
            t.verifyEqual(actual, expected);
            t.verifyTrue(actual(lat == -40, lon == 370));
            [expectedX, expectedY] = meshgrid(lon, lat);
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);
        end

        function wholeCycleShifts(t, shift)
            % A rectangle provides an independent reference including boundaries.
            polygon = [350 -2;360 -2;360 2;350 2];
            lon = [-10, -5, 0, 5, 180, 350, 355, 360] + shift;
            lat = [-3, -2, 0, 2, 3];
            actual = domainmask(polygon, lon, lat, SaveData=false, BeQuiet=true);
            wrapped = mod(lon, 360);
            expected = abs(lat(:)) <= 2 & (wrapped >= 350 | wrapped == 0);
            t.verifyEqual(actual, expected);
        end

        function holesAndDisconnectedRegions(t, buffer)
            % Seam splitting must preserve a hole and a separate region.
            outer = polyshape([-4, 4, 4, -4], [-4, -4, 4, 4]);
            hole = polyshape([-1, 1, 1, -1], [-1, -1, 1, 1]);
            island = polyshape([20, 24, 24, 20], [-2, -2, 2, 2]);
            polygon = union(subtract(outer, hole), island);
            lon = [-5, -3, 0, 3, 5, 22, 26, 180];
            lat = [-3, 0, 3];
            actual = domainmask(polygon, lon, lat, exteriorBuffer=buffer, ...
                SaveData=false, BeQuiet=true);
            shifted = domainmask(translate(polygon, [360, 0]), lon+720, lat, ...
                exteriorBuffer=buffer, SaveData=false, BeQuiet=true);
            t.verifyEqual(actual, shifted);
            if buffer == 0
                expected = abs(lon) <= 4 & abs(lat(:)) <= 4 & ...
                    ~(abs(lon) < 1 & abs(lat(:)) < 1);
                expected = expected | (lon >= 20 & lon <= 24 & abs(lat(:)) <= 2);
                t.verifyEqual(actual, expected);
            else
                t.verifyTrue(all(actual(:, abs(lon) == 5), 'all'));
                t.verifyTrue(actual(lat == 0, lon == 0));
                t.verifyFalse(any(actual(:, lon == 180)));
            end
        end

        function fullLongitudeBand(t)
            % A full-width polygon must not collapse when its seam is changed.
            polygon = [-180 -2;180 -2;180 2;-180 2];
            lon = -720:90:720;
            lat = -3:3;
            actual = domainmask(polygon, lon, lat, SaveData=false, BeQuiet=true);
            t.verifyEqual(actual, repmat(abs(lat(:)) <= 2, 1, numel(lon)));
        end

        function coordinateOrderAndLayouts(t)
            % Permuted coordinates stay aligned with the observations and cache.
            polygon = [0 -1;2 -1;2 1;0 1];
            lon = [3, 0, 1, -1, 2];
            lat = [2, 1, 0, -1, -2];
            expected = lon >= 0 & lon <= 2 & abs(lat(:)) <= 1;
            [actual, x, y] = domainmask(polygon, lon, lat, BeQuiet=true);
            [expectedX, expectedY] = meshgrid(lon, lat);
            t.verifyEqual(actual, expected);
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);
            [actual, x, y] = domainmask(polygon, lon, lat, Format='ndgrid', BeQuiet=true);
            t.verifyEqual(actual, expected.');
            t.verifyEqual(x, expectedX.');
            t.verifyEqual(y, expectedY.');
        end

        function coordinateMatrices(t)
            % Paired coordinate matrices need no regular-grid reconstruction.
            polygon = [-2 -1;2 -1;2 1;-2 1];
            x = [359, 1, 3; 0, 720, -1];
            y = [0, 0, 0; 2, 0, -2];
            [actual, returnedX, returnedY] = domainmask(polygon, x, y, ...
                SaveData=false, BeQuiet=true);
            t.verifyEqual(actual, logical([1, 1, 0; 0, 1, 0]));
            t.verifyEqual(returnedX, x);
            t.verifyEqual(returnedY, y);
        end

        function maskCacheReadAndCreationDate(t)
            % Trust existing masks and timestamp only saved computations.
            polygon = [-1 -1;1 -1;1 1;-1 1];
            lon = -3:3;
            lat = (-3:3).';
            beforeSave = datetime('now', 'TimeZone', 'UTC');
            expected = domainmask(polygon, lon, lat, exteriorBuffer=0.5, BeQuiet=true);
            afterSave = datetime('now', 'TimeZone', 'UTC');
            files = dir(fullfile(getenv('IFILES'), 'DOMAINMASK', '*.mat'));
            t.assertNumElements(files, 1);
            path = fullfile(files.folder, files.name);
            data = load(path);
            t.verifyGreaterThanOrEqual(data.creationDate, beforeSave);
            t.verifyLessThanOrEqual(data.creationDate, afterSave);
            t.verifyFalse(isfield(data, 'maskVersion'));
            % An old timestamp and sentinel distinguish reads from rebuilds.
            savedDate = datetime(2000, 1, 1, 'TimeZone', 'UTC');
            data.creationDate = savedDate;
            data.mask = ~expected;
            save(path, '-struct', 'data');
            t.verifyEqual(domainmask(polygon, lon, lat, exteriorBuffer=0.5, ...
                BeQuiet=true), ~expected);
            cached = load(path);
            t.verifyEqual(cached.creationDate, savedDate);
            % Reordered vectors bypass the bounds-based cache entirely.
            t.verifyEqual(domainmask(polygon, fliplr(lon), lat, exteriorBuffer=0.5, ...
                SaveData=false, BeQuiet=true), fliplr(expected));
            t.verifyEqual(domainmask(polygon, lon, lat, exteriorBuffer=0.5, ...
                Forcenew=true, SaveData=false, BeQuiet=true), expected);
            cached = load(path);
            t.verifyEqual(cached.creationDate, savedDate);
            t.verifyEqual(cached.mask, ~expected);
            beforeSave = datetime('now', 'TimeZone', 'UTC');
            domainmask(polygon, lon, lat, exteriorBuffer=0.5, ...
                Forcenew=true, BeQuiet=true);
            afterSave = datetime('now', 'TimeZone', 'UTC');
            saved = load(path);
            t.verifyEqual(saved.mask, expected);
            t.verifyGreaterThanOrEqual(saved.creationDate, beforeSave);
            t.verifyLessThanOrEqual(saved.creationDate, afterSave);
            % Older files without a creation timestamp remain readable.
            saved = rmfield(saved, 'creationDate');
            saved.mask = ~expected;
            save(path, '-struct', 'saved');
            t.verifyEqual(domainmask(polygon, lon, lat, exteriorBuffer=0.5, ...
                SaveData=false, BeQuiet=true), ~expected);
            t.verifyFalse(ismember('creationDate', who('-file', path)));
        end
    end
end
