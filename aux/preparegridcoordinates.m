%% PREPAREGRIDCOORDINATES - Validates and orders a rectangular geographic grid.
%
% Syntax
%   [lon, lat, iLon, iLat, wghts, wghtMthd] = PREPAREGRIDCOORDINATES(lon, lat)
%   [__] = PREPAREGRIDCOORDINATES(__, RequireGlobalLongitude=true)
%
% Input arguments
%   lon, lat - Finite, strictly monotone coordinate vectors with at least two entries each.
%       Longitude spacings must agree with the mean spacing to 3e-5 degrees
%       (or one input precision step for coarser single-precision coordinates).
%       Latitude must lie in [-90, 90] and may be unevenly spaced.
%       Unit: degrees east and north, respectively.
%   RequireGlobalLongitude (name-value) - Whether to require one complete
%       periodic longitude grid after removing a duplicated endpoint. In this
%       mode each offset from the first meridian must match the regular
%       360-degree grid within the coordinate tolerance, capped at 90/nLon.
%       Default flag: false.
%
% Output arguments
%   lon, lat - Ascending coordinate columns. If the original longitude
%       endpoints differ by one cycle, the last input endpoint is dropped.
%       Unit: degrees east and north, respectively.
%   iLon, iLat - Indices into the original longitude and latitude vectors.
%       Apply the same indices to the spatial dimensions of each field.
%   wghts - Relative cell areas in longitude-by-latitude order.
%       Size: (numel(lon), numel(lat)).
%   wghtMthd - 'cosine-latitude' or 'spherical-midpoint-bands'.
%       Uniform latitude uses cosd(lat). For nonuniform latitude, the
%       band edges are the latitude midpoints, with half-step outer edges
%       clipped at the poles. Longitude spacing cancels in area ratios.
%
% See also
%   fitgriddedslepian, domainmask
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function [lon, lat, iLon, iLat, wghts, wghtMthd] = preparegridcoordinates(lon, lat, options)

    arguments (Input)
        lon {mustBeCoordinateVector}
        lat {mustBeCoordinateVector}
        options.RequireGlobalLongitude (1, 1) logical = false
    end

    if any(abs(lat) > 90)
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Latitude must lie between -90 and 90 degrees.');
    end

    lonTol = 3e-5;

    if isa(lon, 'single')
        lonTol = max(lonTol, double(eps(max(abs(lon)))));
    end

    lon = double(lon(:));
    lat = double(lat(:));
    lonSteps = diff(lon);
    lonSpan = abs(lon(end) - lon(1));
    lonStep = (lon(end) - lon(1)) / (numel(lon) - 1);
    latSteps = diff(lat);

    if (~options.RequireGlobalLongitude && any(abs(lonSteps - lonStep) > lonTol)) || ...
            ~(all(lonSteps > 0) || all(lonSteps < 0)) || ...
            ~(all(latSteps > 0) || all(latSteps < 0))
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Longitude must be regularly spaced and latitude must be strictly monotone.');
    end

    iLon = (1:numel(lon)).';
    % The seam tolerance must distinguish a missing meridian from round-off.
    seamTol = min(lonTol, abs(lonStep) / 2);

    if options.RequireGlobalLongitude
        lonTol = min(lonTol, 90 / numel(lon));
        if abs(lonSpan - 360) <= lonTol
            iLon(end) = [];
        end
        step = sign(lonStep) * 360 / numel(iLon);
        offsets = (0:numel(lon) - 1)' * step;
        if any(abs((lon - lon(1)) - offsets) > lonTol)
            error('ULMO:preparegridcoordinates:InvalidGrid', ...
                'Longitude must cover one complete regular periodic grid.');
        end
    elseif abs(lonSpan - 360) <= seamTol
        iLon(end) = [];
    elseif lonSpan > 360
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
            'Longitude cannot span more than one cycle.');
    end

    [lon, order] = sort(lon(iLon));
    iLon = iLon(order);
    [lat, iLat] = sort(lat);
    latSteps = diff(lat);

    if all(abs(latSteps - latSteps(1)) <= 3e-5)
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

%% Subfunctions
function mustBeCoordinateVector(coords)

    if ~isnumeric(coords) || ~isreal(coords) || ~isvector(coords) || ...
            numel(coords) < 2 || any(~isfinite(coords))
        error('ULMO:preparegridcoordinates:InvalidGrid', ...
        'Supply a finite real numeric coordinate vector with at least two elements.');
    end

end
