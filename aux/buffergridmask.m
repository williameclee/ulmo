%% BUFFERGRIDMASK - Expand a polygon's grid mask using spherical boundary distances.
%
% Syntax
%   mask = BUFFERGRIDMASK(polygon, lon, lat, mask, width)
%
% Input arguments
%   polygon - Scalar polyshape describing the original target region.
%       All boundary rings participate, including holes and separate regions.
%       Vertices contain longitude and latitude in degrees.
%   lon, lat - Paired coordinates of the grid points to classify.
%       Regular spacing is not required. Latitudes must lie in [-90, 90].
%       Size: same two-dimensional size as mask.
%       Unit: degrees east and degrees north, respectively.
%   mask - Logical membership mask for the original target region.
%       The caller supplies this mask in the same coordinate layout as lon and
%       lat. True entries are preserved; only false entries are examined.
%       Size: same as lon and lat.
%   width - Positive exterior-buffer radius, at most 90 degrees.
%       A point is added when its shortest distance to any boundary segment,
%       including the segment endpoints, is within this radius.
%       Unit: degrees of surface arc on a sphere.
%
% Output arguments
%   mask - Original mask combined with the spherical boundary buffer.
%       Expansion adds points outside outer boundaries and into holes.
%       Size: same as the input mask.
%
% Notes
%   Why this helper exists:
%       BUFFERM failed on a full-longitude ocean polygon in MATLAB R2026a
%       Update 2 during longitude-seam reconstruction, with the identifier
%       map:gluingPolygon:expectedIncreasingY. The failure occurred both for
%       separate boundary arcs and for intact rings with 'outPlusInterior'.
%       This helper computes grid membership directly, avoiding construction
%       and seam stitching of a buffered polygon.
%
%   Boundary representation:
%       INTERPM first densifies edges with a maximum spacing of 1 degree by
%       linear interpolation in longitude and latitude. Consecutive vertices then
%       define minor great-circle segments. Distances are spherical distances
%       to these segments; the densification approximates the input polygon's
%       edges in longitude-latitude coordinates.
%
%   Distance calculation:
%       Unit vectors make the distance calculation independent of longitude
%       wrapping. A perpendicular great-circle distance is accepted only when
%       its foot lies between the segment endpoints. Endpoint distances cover
%       corners and zero-length spherical segments, including repeated poles.
%       Comparisons include a 64*eps tolerance in unit-vector coordinates.
%
%   Candidate search:
%       Exterior points are sorted into 1-degree latitude bins. A segment's
%       spherical midpoint cap, with radius width plus half the segment length,
%       bounds all points that can qualify. Candidates are examined one segment
%       at a time, without allocating a point-by-segment distance matrix.
%
% See also
%   DOMAINMASK
%
% Created by
%   2026/11/04

function mask = buffergridmask(polygon, lon, lat, mask, width)
    arguments (Input)
        polygon (1, 1) polyshape
        lon (:, :) double {mustBeReal, mustBeFinite}
        lat (:, :) double {mustBeReal, mustBeFinite}
        mask (:, :) logical
        width (1, 1) double {mustBeFinite, mustBePositive, mustBeLessThanOrEqual(width, 90)}
    end
    assert(isequal(size(lon), size(lat), size(mask)), ...
        'ULMO:buffergridmask:GridSize', 'Longitude, latitude and mask must have the same size.');
    assert(all(abs(lat(:)) <= 90), ...
        'ULMO:buffergridmask:Latitude', 'Latitudes must lie between -90 and 90 degrees.');

    % Query only exterior points; latitude bins bound the candidate search.
    outside = find(~mask);
    if isempty(outside) || isempty(polygon.Vertices)
        return
    end
    [pointLat, order] = sort(lat(outside));
    outside = outside(order);
    points = unitvectors(lon(outside), pointLat);
    counts = histcounts(pointLat, -90:91);
    offsets = [1, cumsum(counts)+1];
    selected = false(numel(outside), 1);
    [x, y] = boundary(polygon);
    [y, x] = interpm(y, x, 1);
    edges = find(isfinite(x(1:end-1)) & isfinite(x(2:end)));
    vertices = unitvectors(x, y);
    tolerance = 64*eps;
    chordWidth = 2*sind(width/2);

    for iEdge = 1:numel(edges)
        first = vertices(edges(iEdge), :);
        last = vertices(edges(iEdge)+1, :);
        normal = cross(first, last);
        normalLength = norm(normal);
        halfLength = atan2d(normalLength, dot(first, last))/2;
        midpoint = first + last;
        midpoint = midpoint/norm(midpoint);
        midpointLat = atan2d(midpoint(3), hypot(midpoint(1), midpoint(2)));
        % Every buffered point is within width+halfLength of the midpoint.
        radius = width + halfLength;
        firstBin = max(1, floor(midpointLat-radius-1e-10)+91);
        lastBin = min(181, floor(midpointLat+radius+1e-10)+91);
        rows = (offsets(firstBin):offsets(lastBin+1)-1).';
        rows = rows(~selected(rows));
        candidates = points(rows, :);
        nearby = vecnorm(candidates-midpoint, 2, 2) <= 2*sind(radius/2)+tolerance;
        rows = rows(nearby);
        candidates = candidates(nearby, :);
        if isempty(rows)
            continue
        end
        within = vecnorm(candidates-first, 2, 2) <= chordWidth+tolerance | ...
            vecnorm(candidates-last, 2, 2) <= chordWidth+tolerance;
        if normalLength > tolerance
            normal = normal/normalLength;
            % The perpendicular foot must lie between the segment endpoints.
            between = candidates*cross(normal, first).' >= -tolerance & ...
                candidates*cross(last, normal).' >= -tolerance;
            within = within | (between & abs(candidates*normal.') <= sind(width)+tolerance);
        end
        selected(rows(within)) = true;
    end
    mask(outside(selected)) = true;
end

function xyz = unitvectors(lon, lat)
    xyz = [cosd(lat(:)).*cosd(lon(:)), cosd(lat(:)).*sind(lon(:)), sind(lat(:))];
end

