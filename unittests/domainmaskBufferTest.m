%% DOMAINMASKBUFFERTEST - Check exterior mask topology and cache widths.
%
% Created by
%   2026/11/04, En-Chi Lee

classdef domainmaskBufferTest < matlab.unittest.TestCase
    properties
        Polygon = [0 -11; 32 -11; 32 11; 0 11; 0 -11]
    end
    properties (TestParameter)
        hemisphere = {1, -1}
        longitudeOrigin = {0, 180}
    end
    methods (TestMethodSetup)
        function temporaryRoot(t)
            % Keep generated masks separate from production data.
            folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
        end
    end
    methods (Test)
        function bufferTouchesPole(t, hemisphere, longitudeOrigin)
            % A buffer touching either pole must survive in both longitude frames.
            xy = [-4 87.5; 4 87.5; 4 88.5; -4 88.5; -4 87.5];
            xy(:, 1) = xy(:, 1) + longitudeOrigin;
            xy(:, 2) = hemisphere * xy(:, 2);
            lon = (-180:180) + longitudeOrigin;
            lat = sort(hemisphere * (70:90));
            actual = domainmask(xy, lon, lat, exteriorBuffer=1.5, BeQuiet=true);

            % Check the original interior, the added buffer and distant points.
            t.verifyTrue(actual(lat == hemisphere*88, lon == longitudeOrigin));
            t.verifyTrue(actual(lat == hemisphere*89, lon == longitudeOrigin));
            t.verifyFalse(actual(lat == hemisphere*85, lon == longitudeOrigin));
            t.verifyFalse(actual(lat == hemisphere*88, lon == longitudeOrigin+120));

            [cached, x, y] = domainmask(xy, lon, lat, exteriorBuffer=1.5, ...
                Format='ndgrid', BeQuiet=true);
            [expectedX, expectedY] = ndgrid(lon, lat);
            t.verifyEqual(cached, actual.');
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);

            % A sentinel distinguishes a cache read from a fresh computation.
            file = dir(fullfile(getenv('IFILES'), 'DOMAINMASK', '*.mat'));
            t.assertNumElements(file, 1);
            mask = ~actual;
            save(fullfile(file.folder, file.name), 'mask', '-append');
            t.verifyEqual(domainmask(xy, lon, lat, exteriorBuffer=1.5, ...
                BeQuiet=true), mask);
        end

        function exteriorCacheAndTopology(t)
            % A ring crossing the grid seam and a separate region with a hole.
            xy = [0 -4;4 -4;4 4;0 4;0 -4;NaN NaN; ...
                50 -7;65 -7;65 7;50 7;50 -7;NaN NaN; ...
                53 -4;53 4;62 4;62 -4;53 -4];
            lon = 0:359; lat = -9:9;
            m = domainmask(xy,lon,lat,exteriorBuffer=1.5,BeQuiet=true);
            at = @(x) m(lat==0,lon==x);
            t.verifyTrue(at(359) && at(5) && at(54) && at(61));
            t.verifyFalse(at(358) || at(6) || at(57) || at(30));
            [cached,x,y] = domainmask(xy,lon,lat,exteriorBuffer=1.5,Format='ndgrid',BeQuiet=true);
            t.verifyEqual(cached,m.');
            [expectedX,expectedY] = ndgrid(lon,lat);
            t.verifyEqual(x,expectedX);
            t.verifyEqual(y,expectedY);
            narrow = domainmask(xy,lon,lat,exteriorBuffer=0.4,BeQuiet=true);
            t.verifyFalse(narrow(lat==0,lon==5));
            file = dir(fullfile(getenv('IFILES'),'DOMAINMASK','*-E0.4-H*.mat'));
            t.assertNumElements(file,1);
        end

        function roundedMaskNameChecksExactWidth(t)
            % Both widths round to one filename, but each needs its own mask.
            for width = [1/3, 1/3+1e-8]
                domainmask(t.Polygon,0:33,-10:10,exteriorBuffer=width,BeQuiet=true);
                mask = dir(fullfile(getenv('IFILES'),'DOMAINMASK',['*-E',num2str(width),'-*.mat']));
                t.assertNumElements(mask,1);
                cached = load(fullfile(mask.folder,mask.name),'exteriorBuffer');
                t.verifyEqual(cached.exteriorBuffer,width);
            end
        end
    end
end
