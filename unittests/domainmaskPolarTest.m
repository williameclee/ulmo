%% DOMAINMASKPOLARTEST - Checks geographic polar masks and cache reuse.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

classdef domainmaskPolarTest < matlab.unittest.TestCase

    properties
        Root
        Polygon = [10 60; 20 60; 20 65; 10 65; 10 60]
    end

    properties (TestParameter)
        buf = {0, 1.5}
    end

    methods (TestMethodSetup)

        function arcticCoastlines(t)
            % Exercise the real Arctic reader with isolated geographic and rotated caches.
            folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Root = folder.Folder;
            coasts = fullfile(t.Root, 'COASTS');
            mkdir(coasts);
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', t.Root));
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('COASTS', coasts));
            fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fixtures));
            writeArcticCoastlineFixture(coasts, t.Polygon);
        end

    end

    methods (Test)

        function geographicMaskAndCache(t, buf)
            % A GeoDomain and its geographic polygon must select the same cells.
            lon = 0.25:0.5:359.75;
            lat = -34.75:0.5:74.75;
            expected = domainmask(t.Polygon, lon, lat, exteriorBuffer = buf, ...
                SaveData = false, BeQuiet = true);
            actual = domainmask(GeoDomain('arctic'), lon, lat, ...
                exteriorBuffer = buf, BeQuiet = true);
            t.verifyEqual(actual, expected);
            [cached, x, y] = domainmask(GeoDomain('arctic'), lon, lat, ...
                exteriorBuffer = buf, Format = 'ndgrid', BeQuiet = true);
            [expectedX, expectedY] = ndgrid(lon, lat);
            t.verifyEqual(cached, expected.');
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);
            t.verifyEqual(domainmask('arctic', lon, lat, exteriorBuffer = buf, ...
                BeQuiet = true), expected);
        end

        function ignoresRotatedMaskCache(t, buf)
            % Legacy mask files contain the old equatorial selection.
            lon = 0:2:358;
            lat = -40:2:80;
            mask = false(numel(lat), numel(lon));
            mask(lat == 0, lon == 180) = true;
            extBuf = buf;
            folder = fullfile(t.Root, 'DOMAINMASK');
            mkdir(folder);
            name = 'Arctic';

            if buf > 0
                name = sprintf('%s-E%s', name, num2str(buf));
            end

            path = fullfile(folder, sprintf('DOMAINMASK-%s-H2-Box0_358_-40_80.mat', name));
            save(path, 'mask', 'lon', 'lat', 'extBuf');
            expected = domainmask(t.Polygon, lon, lat, exteriorBuffer = buf, ...
                SaveData = false, BeQuiet = true);
            actual = domainmask(GeoDomain('arctic'), lon, lat, ...
                exteriorBuffer = buf, BeQuiet = true);
            t.verifyEqual(actual, expected);
        end

    end

end
