%% PREPAREGRIDCOORDINATES - Validate and order a rectangular geographic grid.
%
% Syntax
%   [lon, lat, iLon, iLat, weights, weighting] = PREPAREGRIDCOORDINATES(lon, lat)
%   [__] = PREPAREGRIDCOORDINATES(__, RequireGlobalLongitude=true)
%
% Input arguments
%   lon, lat - Finite real coordinate vectors with at least two elements.
%       Longitude must be regular; latitude must be strictly monotone and
%       within [-90, 90]. Either axis may descend. Spacing tolerance is
%       3e-5 degrees, including single-precision coordinate rounding.
%       Unit: degrees east and north, respectively.
%   RequireGlobalLongitude (name-value) - Whether to require one complete
%       periodic longitude grid after removing a duplicated endpoint.
%       Default flag: false.
%
% Output arguments
%   lon, lat - Ascending coordinate columns, with a cyclic endpoint removed.
%       Unit: degrees east and north, respectively.
%   iLon, iLat - Indices into the original coordinate vectors. Apply them to
%       fields already arranged as (longitude, latitude, time). No field
%       orientation is inferred, and longitude values are not wrapped.
%   weights - Relative cell areas in (longitude, latitude) order.
%       Size: (numel(lon), numel(lat)); dimensionless.
%   weighting - 'cosine-latitude' or 'spherical-midpoint-bands'.
%       Uniform latitude uses cosd(lat). Nonuniform latitude uses spherical
%       bands between arithmetic midpoints, with half-step outer boundaries
%       clipped at the poles. Regular longitude spacing cancels in ratios.
%
% See also
%   fitgriddedslepian, domainmask
%
% Created by
%   2026/11/06, En-Chi Lee

function [lon, lat, iLon, iLat, weights, weighting] = preparegridcoordinates(lon, lat, options)
    arguments (Input)
        lon
        lat
        options.RequireGlobalLongitude (1, 1) logical = false
    end
    if ~isnumeric(lon) || ~isnumeric(lat) || ~isvector(lon) || ~isvector(lat) || ...
            ~isreal(lon) || ~isreal(lat) || any(~isfinite(lon(:))) || ...
            any(~isfinite(lat(:))) || numel(lon) < 2 || numel(lat) < 2 || any(abs(lat) > 90)
        error('ULMO:preparegridcoordinates:InvalidGrid', 'Expected finite longitude/latitude coordinate vectors.');
    end
    lon = double(lon(:)); lat = double(lat(:));
    tolerance = 3e-5;
    longitudeStep = lon(2) - lon(1);
    latitudeSteps = diff(lat);
    if longitudeStep == 0 || any(abs(diff(lon) - longitudeStep) > tolerance) || ...
            ~(all(diff(lon) > 0) || all(diff(lon) < 0)) || ...
            ~(all(latitudeSteps > 0) || all(latitudeSteps < 0))
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
            'Longitude must be regularly spaced and latitude must be strictly monotone.');
    end
    iLon = (1:numel(lon)).';
    if abs(abs(lon(end) - lon(1)) - 360) <= tolerance, iLon(end) = []; end
    if abs(lon(end) - lon(1)) > 360 + tolerance
        error('ULMO:preparegridcoordinates:InvalidGrid', 'Longitude cannot span more than one cycle.');
    end
    if options.RequireGlobalLongitude && ...
            abs(numel(iLon)*abs(longitudeStep) - 360) > tolerance*numel(lon)
        error('ULMO:preparegridcoordinates:InvalidGrid', 'Longitude must cover one complete periodic grid.');
    end
    [lon, order] = sort(lon(iLon)); iLon = iLon(order);
    [lat, iLat] = sort(lat);
    latitudeSteps = diff(lat);
    if all(abs(latitudeSteps - latitudeSteps(1)) <= tolerance)
        latitudeWeights = cosd(lat);
        weighting = 'cosine-latitude';
    else
        edges = [lat(1)-latitudeSteps(1)/2; (lat(1:end-1)+lat(2:end))/2; ...
            lat(end)+latitudeSteps(end)/2];
        latitudeWeights = diff(sind(max(-90, min(90, edges))));
        weighting = 'spherical-midpoint-bands';
    end
    weights = repmat(latitudeWeights.', numel(lon), 1);
end
