%% POLARDOMAINMASKTEST - Check geographic polar masks and SSH observations.
%
% Created by
%   2026/11/04, En-Chi Lee

classdef polarDomainMaskTest < matlab.unittest.TestCase
    properties
        Root
        Polygon = [10 60; 20 60; 20 65; 10 65; 10 60]
    end
    properties (TestParameter)
        buffer = {0, 1.5}
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
            XY = t.Polygon;
            p = polyshape(XY);
            XY = poly2xy(p);
            save(fullfile(coasts, 'Arctic-rotb.mat'), 'XY', 'p');
            [colatitude, longitude] = rottp(deg2rad(90-XY(:,2)), ...
                deg2rad(XY(:,1)), 0, -pi/2, 0);
            XY = [mod(rad2deg(longitude),360), 90-rad2deg(colatitude)];
            p = polyshape(XY);
            XY = poly2xy(p);
            save(fullfile(coasts, 'Arctic.mat'), 'XY', 'p');
        end
    end
    methods (Test)
        function geographicMaskAndCache(t, buffer)
            % A GeoDomain and its geographic polygon must select the same cells.
            lon = 0.25:0.5:359.75;
            lat = -34.75:0.5:74.75;
            expected = domainmask(t.Polygon, lon, lat, exteriorBuffer=buffer, ...
                SaveData=false, BeQuiet=true);
            actual = domainmask(GeoDomain('arctic'), lon, lat, ...
                exteriorBuffer=buffer, BeQuiet=true);
            t.verifyEqual(actual, expected);
            [cached, x, y] = domainmask(GeoDomain('arctic'), lon, lat, ...
                exteriorBuffer=buffer, Format='ndgrid', BeQuiet=true);
            [expectedX, expectedY] = ndgrid(lon, lat);
            t.verifyEqual(cached, expected.');
            t.verifyEqual(x, expectedX);
            t.verifyEqual(y, expectedY);
            t.verifyEqual(domainmask('arctic', lon, lat, exteriorBuffer=buffer, ...
                BeQuiet=true), expected);
        end

        function ignoresRotatedMaskCache(t, buffer)
            % Legacy mask files contain the old equatorial selection.
            lon = 0:2:358;
            lat = -40:2:80;
            mask = false(numel(lat), numel(lon));
            mask(lat == 0, lon == 180) = true;
            exteriorBuffer = buffer;
            folder = fullfile(t.Root, 'DOMAINMASK');
            mkdir(folder);
            name = 'Arctic';
            if buffer > 0
                name = sprintf('%s-E%s', name, num2str(buffer));
            end
            path = fullfile(folder, sprintf('DOMAINMASK-%s-H2-Box0_358_-40_80.mat', name));
            save(path, 'mask', 'lon', 'lat', 'exteriorBuffer');
            expected = domainmask(t.Polygon, lon, lat, exteriorBuffer=buffer, ...
                SaveData=false, BeQuiet=true);
            actual = domainmask(GeoDomain('arctic'), lon, lat, ...
                exteriorBuffer=buffer, BeQuiet=true);
            t.verifyEqual(actual, expected);
        end

        function sshUsesPolarObservations(t, buffer)
            % Polar SSH is 2 m. Equatorial SSH is 99 m, including an empty polar month.
            t.writeSSH();
            degree = 0;
            if buffer > 0
                degree = 60;
            end
            G = zeros((degree+1)^2,1);
            G(1) = 1;
            [coefficients, ~, fit] = ssh2slept('NASASSHMonthly', 'arctic', degree, ...
                GVN={G,1,1}, unit='m');
            t.verifyEqual(coefficients, [2;NaN], 'AbsTol', 1e-10);
            t.verifyEqual(fit.areaFraction, [1;0], 'AbsTol', 1e-12);
            t.verifyEqual(fit.nObs(2), 0);
        end

        function ignoresRotatedProjectionCache(t)
            % An existing Arctic projection must be rebuilt using geographic observations.
            t.writeSSH();
            folder = fullfile(t.Root, 'GLMALPHA');
            mkdir(folder);
            G = 1; V = 0.01; N = V; EL = 0; EM = 0;
            save(fullfile(folder,'glmalpha-Arctic-0-1.mat'), 'G','V','N','EL','EM');
            folder = fullfile(t.Root, 'SSH', 'SlepianExpansions');
            mkdir(folder);
            path = fullfile(folder, 'ssh2slept-NASASSHMonthly-Arctic-0-1-m.nc');
            write2nc(path, 'coefficients', [99 99], Dimensions={'mode','time'});
            ncwriteatt(path, '/', 'defaultTruncation', 1);
            dates = datetime(2020,[1;2],15);
            nccreate(path, 'time', 'Dimensions', {'time',2});
            ncwrite(path, 'time', seconds(dates-datetime(2000,1,1)));
            ncwriteatt(path, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
            ncwriteatt(path, 'time', 'calendar', 'proleptic_gregorian');
            actual = ssh2slept('NASASSHMonthly', 'arctic', 0, unit='m');
            t.verifyEqual(actual, [2;NaN], 'AbsTol', 1e-12);
            delete(fullfile(t.Root, 'SSH', 'NASASSHMonthly.mat'));
            t.verifyEqual(ssh2slept('NASASSHMonthly', 'arctic', 0, ...
                truncation=1, unit='m'), actual);
        end
    end
    methods (Access=private)
        function writeSSH(t)
            % Supply a native grid spanning the rotated and geographic polygons.
            folder = fullfile(t.Root, 'SSH');
            mkdir(folder);
            lon = (0.25:0.5:359.75)';
            lat = (-34.75:0.5:74.75)';
            dates = datetime(2020,[1;2],15);
            partialMonths = false(2,1);
            sshs = 99*ones(numel(lon),numel(lat),2);
            sshs(:,lat>50,1) = 2;
            sshs(:,lat>50,2) = NaN;
            save(fullfile(folder,'NASASSHMonthly.mat'), ...
                'sshs','lon','lat','dates','partialMonths','-v7.3');
            folder = fullfile(t.Root,'EARTHMODELS','CONSTANTS');
            mkdir(folder);
            Earth.Radius = 6371000;
            save(fullfile(folder,'Earth.mat'), 'Earth');
        end
    end
end
