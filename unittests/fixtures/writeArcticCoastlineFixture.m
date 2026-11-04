%% WRITEARCTICCOASTLINEFIXTURE - Save a synthetic Arctic boundary in both frames.
%
% Input arguments
%   coasts - Temporary coastline directory.
%   polygon - Geographic longitude and latitude vertices, in degrees.
%       Size: (n, 2).
%
% Created by
%   2026/11/04, En-Chi Lee

function writeArcticCoastlineFixture(coasts, polygon)
    XY = polygon;
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
