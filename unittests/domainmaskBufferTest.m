%% DOMAINMASKBUFFERTEST - Checks exterior mask topology and cache widths.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

classdef domainmaskBufferTest < matlab.unittest.TestCase

    properties
        Polygon = [0 -11; 32 -11; 32 11; 0 11; 0 -11]
    end

    properties (TestParameter)
        hemisphere = {1, -1}
        lonOrigin = {0, 180}
    end

    methods (TestMethodSetup)

        function temporaryRoot(t)
            % Keep generated masks separate from production data.
            folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
        end

    end

    methods (Test)

        function highLatitudeBufferAccuracy(t)
            % Spherical endpoint distances include the high-latitude case in #91.
            xy = [-4 73; 4 73; 4 77; -4 77; -4 73];
            lon = -180:180;
            lat = -90:90;
            extBuf = 10;
            mask = domainmask(xy, lon, lat, exteriorBuffer = extBuf, ...
                Forcenew = true, SaveData = false, BeQuiet = true);

            % This point is within the buffer because a vertex is only 8.895 deg away.
            vtxDist = distance(85, 32, 77, 4);
            t.assertLessThan(vtxDist, extBuf);
            t.verifyTrue(mask(lat == 85, lon == 32));
        end

        function fullLongitudeBand(t)
            % Longitude cuts must not disrupt the global band's exterior buffer.
            p = [0 -66; 360 -66; 360 66; 0 66; 0 -66];
            lon = 0.25:0.5:359.75;
            lat = (-70:0.5:70) + 0.25;
            actual = domainmask(p, lon, lat, exteriorBuffer = 1.5, ...
                SaveData = false, BeQuiet = true);
            expected = repmat(abs(lat(:)) < 67.5, 1, numel(lon));
            t.verifyEqual(actual, expected);
        end

        function segmentInteriorAndEndpoints(t)
            % Include an interior perpendicular foot and exclude the line extension.
            p = [-1 -1; 1 -1; 1 1; -1 1; -1 -1];
            lon = [-2, -1, -0.5, 0, 0.5, 1, 2];
            lat = [-2, -1.2, -1, 0, 1, 1.2, 2];
            mask = domainmask(p, lon, lat, exteriorBuffer = 0.25, ...
                SaveData = false, BeQuiet = true);
            % At half a degree between densified vertices, endpoints are too far.
            t.verifyGreaterThan(distance(1.2, 0.5, 1, 0), 0.25);
            t.verifyTrue(mask(lat == 1.2, lon == 0.5));
            t.verifyFalse(mask(lat == 2, lon == 0.5));
            t.verifyFalse(mask(lat == 1.2, lon == 2));
            t.verifyTrue(mask(lat == 0, lon == 0));
        end

        function hemisphereWidth(t)
            % A 90-degree buffer includes poles but not antipodal line extensions.
            p = [-1 -1; 1 -1; 1 1; -1 1; -1 -1];
            lon = -180:180;
            lat = -90:90;
            mask = domainmask(p, lon, lat, exteriorBuffer = 90, ...
                SaveData = false, BeQuiet = true);
            t.verifyTrue(all(mask(abs(lat) == 90, :), 'all'));
            t.verifyTrue(mask(lat == 0, lon == 90));
            t.verifyFalse(mask(lat == 0, lon == 92));
            t.verifyFalse(mask(lat == 0, lon == 180));
        end

        function matchesMeridianDistances(t, lonOrigin)
            % Analytic cross-track distances give an independent segment reference.
            p = [0 -20; 20 -20; 20 20; 0 20; 0 -20];
            p(:, 1) = p(:, 1) + lonOrigin;
            lon = (-15:0.5:-1) + lonOrigin;
            lat = -10:0.5:10;
            [x, y] = meshgrid(lon - lonOrigin, lat);
            % Include the exact cutoff despite inverse-trigonometric roundoff.
            expected = asind(abs(cosd(y) .* sind(x))) <= 7 +1e-12;
            mask = domainmask(p, lon, lat, exteriorBuffer = 7, ...
                SaveData = false, BeQuiet = true);
            t.verifyEqual(mask, expected);

            % Beyond the northern end, use vertex distance, not the infinite meridian.
            lat = 30:0.5:40;
            [x, y] = meshgrid(lon - lonOrigin, lat);
            expected = distance(y, x, 20 * ones(size(y)), zeros(size(x))) <= 15 +1e-12;
            mask = domainmask(p, lon, lat, exteriorBuffer = 15, ...
                SaveData = false, BeQuiet = true);
            t.verifyEqual(mask, expected);
        end

        function repeatedPolarVertices(t)
            % Distinct longitudes at a pole represent zero-length spherical edges.
            p = [-30 89; 30 89; 30 90; -30 90; -30 89];
            mask = domainmask(p, -180:180, [88, 89, 89.5, 90], ...
                exteriorBuffer = 0.6, SaveData = false, BeQuiet = true);
            t.verifyTrue(all(mask(4, :)));
            t.verifyTrue(all(mask(3, :)));
            t.verifyFalse(any(mask(1, :)));
        end

        function bufferBeyondGridLongitudeBounds(t)
            % The target may lie outside the grid while its exterior crosses into it.
            p = [0 -1; 1 -1; 1 1; 0 1; 0 -1];
            mask = domainmask(p, 358:359, -1:1, exteriorBuffer = 1.5, ...
                SaveData = false, BeQuiet = true);
            t.verifyFalse(mask(2, 1));
            t.verifyTrue(mask(2, 2));
        end

        function bufferTouchesPole(t, hemisphere, lonOrigin)
            % A buffer touching either pole must survive in both longitude frames.
            xy = [-4 87.5; 4 87.5; 4 88.5; -4 88.5; -4 87.5];
            xy(:, 1) = xy(:, 1) + lonOrigin;
            xy(:, 2) = hemisphere * xy(:, 2);
            lon = (-180:180) + lonOrigin;
            lat = sort(hemisphere * (70:90));
            actual = domainmask(xy, lon, lat, exteriorBuffer = 1.5, BeQuiet = true);

            % Check the original interior, the added buffer and distant points.
            t.verifyTrue(actual(lat == hemisphere * 88, lon == lonOrigin));
            t.verifyTrue(actual(lat == hemisphere * 89, lon == lonOrigin));
            t.verifyFalse(actual(lat == hemisphere * 85, lon == lonOrigin));
            t.verifyFalse(actual(lat == hemisphere * 88, lon == lonOrigin + 120));

            [cached, x, y] = domainmask(xy, lon, lat, exteriorBuffer = 1.5, ...
                Format = 'ndgrid', BeQuiet = true);
            [expectedX, expectedY] = ndgrid(lon, lat);
            t.verifyEqual(cached, actual.');
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);

            % A sentinel distinguishes a cache read from a fresh computation.
            file = dir(fullfile(getenv('IFILES'), 'DOMAINMASK', '*.mat'));
            t.assertNumElements(file, 1);
            mask = ~actual;
            save(fullfile(file.folder, file.name), 'mask', '-append');
            t.verifyEqual(domainmask(xy, lon, lat, exteriorBuffer = 1.5, ...
                BeQuiet = true), mask);
        end

        function exteriorCacheAndTopology(t)
            % A ring crossing the grid seam and a separate region with a hole.
            xy = [0 -4; 4 -4; 4 4; 0 4; 0 -4; NaN NaN; ...
                      50 -7; 65 -7; 65 7; 50 7; 50 -7; NaN NaN; ...
                      53 -4; 53 4; 62 4; 62 -4; 53 -4];
            lon = 0:359; lat = -9:9;
            m = domainmask(xy, lon, lat, exteriorBuffer = 1.5, BeQuiet = true);
            at = @(x) m(lat == 0, lon == x);
            t.verifyTrue(at(359) && at(5) && at(54) && at(61));
            t.verifyFalse(at(358) || at(6) || at(57) || at(30));
            [cached, x, y] = domainmask(xy, lon, lat, exteriorBuffer = 1.5, Format = 'ndgrid', BeQuiet = true);
            t.verifyEqual(cached, m.');
            [expectedX, expectedY] = ndgrid(lon, lat);
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);
            narrow = domainmask(xy, lon, lat, exteriorBuffer = 0.4, BeQuiet = true);
            t.verifyFalse(narrow(lat == 0, lon == 5));
            file = dir(fullfile(getenv('IFILES'), 'DOMAINMASK', '*-E0.4-H*.mat'));
            t.assertNumElements(file, 1);
        end

        function roundedMaskNameChecksExactWidth(t)
            % Both widths round to one filename, but each needs its own mask.
            for width = [1/3, 1/3 +1e-8]
                domainmask(t.Polygon, 0:33, -10:10, exteriorBuffer = width, BeQuiet = true);
                mask = dir(fullfile(getenv('IFILES'), 'DOMAINMASK', ['*-E', num2str(width), '-*.mat']));
                t.assertNumElements(mask, 1);
                cached = load(fullfile(mask.folder, mask.name), 'extBuf');
                t.verifyEqual(cached.extBuf, width);
            end

        end

    end

end
