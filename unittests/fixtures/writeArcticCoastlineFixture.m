%% WRITEARCTICCOASTLINEFIXTURE - Saves a synthetic Arctic boundary in both frames.
%
% Syntax
%   WRITEARCTICCOASTLINEFIXTURE(coasts, polygon)
%
% Input arguments
%   coasts - Temporary coastline directory.
%   polygon - Geographic longitude and latitude vertices.
%       Size: (n, 2).
%       Unit: degrees.
%
% Output file content
%   XY - Closed boundary vertices in Arctic-rotb.mat and Arctic.mat.
%       The files use geographic and rotated coordinates, respectively.
%       Size: (n, 2).
%       Unit: degrees.
%   p - Polyshape boundary in the same coordinate frame as XY.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

function writeArcticCoastlineFixture(coasts, polygon)
    XY = polygon;
    p = polyshape(XY);
    XY = poly2xy(p);
    save(fullfile(coasts, 'Arctic-rotb.mat'), 'XY', 'p');
    [col, lon] = rottp(deg2rad(90 - XY(:, 2)), ...
        deg2rad(XY(:, 1)), 0, -pi / 2, 0);
    XY = [mod(rad2deg(lon), 360), 90 - rad2deg(col)];
    p = polyshape(XY);
    XY = poly2xy(p);
    save(fullfile(coasts, 'Arctic.mat'), 'XY', 'p');
end
