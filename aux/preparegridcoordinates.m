%% PREPAREGRIDCOORDINATES - Validates and order a rectangular geographic grid.
%
% Syntax
%   [lon, lat, iLon, iLat, weights, wghtMthd] = PREPAREGRIDCOORDINATES(lon, lat)
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
%   wghtMthd - 'cosine-latitude' or 'spherical-midpoint-bands'.
%       Uniform latitude uses cosd(lat). Nonuniform latitude uses spherical
%       bands between arithmetic midpoints, with half-step outer boundaries
%       clipped at the poles. Regular longitude spacing cancels in ratios.
%
% See also
%   fitgriddedslepian, domainmask
%
% Created by
%   2026/11/06, En-Chi Lee

function [lon, lat, iLon, iLat, wghts, wghtMthd] = preparegridcoordinates(lon, lat, options)

    arguments (Input)
        lon {mustBeReal, mustBeFinite, mustBeVector}
        lat {mustBeReal, mustBeFinite, mustBeVector, mustBeBetween(lat, -90, 90)}
        options.RequireGlobalLongitude (1, 1) logical = false
    end

    if numel(lon) < 2 || numel(lat) < 2
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Expected longitude/latitude coordinate vectors of at least length 2.');
    end

    lon = double(lon(:));
    lat = double(lat(:));
    tol = 3e-5;
    lonStep = lon(2) - lon(1);
    latSteps = diff(lat);

    if lonStep == 0 || any(abs(diff(lon) - lonStep) > tol) || ...
            ~(all(diff(lon) > 0) || all(diff(lon) < 0)) || ...
            ~(all(latSteps > 0) || all(latSteps < 0))
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Longitude must be regularly spaced and latitude must be strictly monotone.');
    end

    iLon = (1:numel(lon)).';

    if abs(abs(lon(end) - lon(1)) - 360) <= tol
        % The first and last longitude ar the same (e.g. 0° and 360°)
        iLon(end) = [];
    elseif abs(lon(end) - lon(1)) > 360 + tol
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Longitude cannot span more than one cycle.');
    elseif options.RequireGlobalLongitude && ...
            abs(numel(iLon) * abs(lonStep) - 360) > tol * numel(lon)
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Longitude must cover one complete periodic grid.');
    end

    [lon, order] = sort(lon(iLon));
    iLon = iLon(order);
    [lat, iLat] = sort(lat);
    latSteps = diff(lat);

    if all(abs(latSteps - latSteps(1)) <= tol)
        latWghts = cosd(lat);
        wghtMthd = 'cosine-latitude';
    else
        edges = [lat(1) - latSteps(1) / 2; ...
                     (lat(1:end - 1) + lat(2:end)) / 2; ...
                     lat(end) + latSteps(end) / 2];
        edges = max(-90, min(90, edges));
        latWghts = diff(sind(edges));
        wghtMthd = 'spherical-midpoint-bands';
    end

    wghts = repmat(latWghts.', numel(lon), 1);
end
