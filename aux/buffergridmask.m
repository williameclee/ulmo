%% BUFFERGRIDMASK - Expands a polygon's grid mask using spherical boundary distances.
%
% Syntax
%   mask = BUFFERGRIDMASK(p, lon, lat, mask, width)
%
% Input arguments
%   p - Polyshape describing the original target region.
%       Vertices contain longitude and latitude in degrees.
%   lon, lat - Paired coordinates of the grid points to classify.
%       Size: same two-dimensional size as mask.
%       Unit: degrees east and degrees north, respectively.
%   mask - Logical mask for the original target region.
%       True entries are preserved; only false entries are examined.
%       Size: same as lon and lat.
%   width - Positive exterior-buffer radius, at most 90 degrees.
%       A point is added to the mask when its shortest distance to any boundary
%       segment is within this radius.
%       Unit: degrees of surface arc on a sphere.
%
% Output arguments
%   mask - Original mask combined with the spherical boundary buffer.
%       Size: same as the input mask.
%
% Notes
%   BUFFERM failed on a full-longitude ocean polygon in MATLAB R2026a
%   Update 2 during longitude-seam reconstruction, with the identifier
%   map:gluingPolygon:expectedIncreasingY. The failure occurred both for
%   separate boundary arcs and for intact rings with 'outPlusInterior'.
%   This helper computes grid membership directly, avoiding construction
%   and seam stitching of a buffered polygon.
%
%   Exterior points are sorted into 1-degree latitude bins. A segment's
%   spherical midpoint cap, with radius width plus half the segment length,
%   bounds all points that can qualify. Candidates are examined one segment
%   at a time, without allocating a point-by-segment distance matrix.
%
% See also
%   DOMAINMASK
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

function mask = buffergridmask(p, lon, lat, mask, width)

    arguments (Input)
        p (1, 1) polyshape
        lon (:, :) double {mustBeReal, mustBeFinite}
        lat (:, :) double {mustBeReal, mustBeFinite}
        mask (:, :) logical
        width (1, 1) double {mustBeFinite, mustBePositive, mustBeLessThanOrEqual(width, 90)}
    end

    assert(isequal(size(lon), size(lat), size(mask)), ...
        'ULMO:buffergridmask:GridSize', ...
    'Longitude, latitude and mask must have the same size.');
    assert(all(abs(lat(:)) <= 90), ...
        'ULMO:buffergridmask:Latitude', ...
    'Latitudes must lie between -90 and 90 degrees.');

    % Query only exterior points; latitude bins bound the candidate search
    notMask = find(~mask);

    if isempty(notMask) || isempty(p.Vertices)
        return
    end

    [ptLat, order] = sort(lat(notMask));
    notMask = notMask(order);
    pts = unitvectors(lon(notMask), ptLat); % Points outside the mask sorted by latitudes
    cnts = histcounts(ptLat, -90:91);
    offsets = [1, cumsum(cnts) + 1]; % Where each latitude bin starts
    selected = false(numel(notMask), 1);
    [x, y] = boundary(p);
    [y, x] = interpm(y, x, 1); % Make the polygon vertices denser
    edges = find(isfinite(x(1:end - 1)) & isfinite(x(2:end)));
    vtxs = unitvectors(x, y);
    tol = 64 * eps;
    chordWidth = 2 * sind(width / 2);

    for iEdge = 1:numel(edges)
        vtx1 = vtxs(edges(iEdge), :);
        vtx2 = vtxs(edges(iEdge) + 1, :);
        edgeNormal = cross(vtx1, vtx2);
        edgeNormalLen = norm(edgeNormal);
        edgeHalfLen = atan2d(edgeNormalLen, dot(vtx1, vtx2)) / 2;
        edgeMidpt = vtx1 + vtx2;
        edgeMidpt = edgeMidpt / norm(edgeMidpt);
        midptLat = atan2d(edgeMidpt(3), hypot(edgeMidpt(1), edgeMidpt(2)));
        % Every buffered point is within (buffer width + edge's half width) of the midpoint
        radius = width + edgeHalfLen;
        firstBin = max(1, floor(midptLat - radius -1e-10) + 91);
        lastBin = min(181, floor(midptLat + radius +1e-10) + 91);
        rows = (offsets(firstBin):offsets(lastBin + 1) - 1).';
        rows = rows(~selected(rows));
        candPts = pts(rows, :);
        nearby = vecnorm(candPts - edgeMidpt, 2, 2) <= 2 * sind(radius / 2) + tol;
        rows = rows(nearby);
        candPts = candPts(nearby, :);

        if isempty(rows)
            continue
        end

        within = vecnorm(candPts - vtx1, 2, 2) <= chordWidth + tol | ...
            vecnorm(candPts - vtx2, 2, 2) <= chordWidth + tol;

        if edgeNormalLen > tol
            edgeNormal = edgeNormal / edgeNormalLen;
            % The perpendicular foot must lie between the segment endpoints
            between = candPts * cross(edgeNormal, vtx1).' >= -tol & ...
                candPts * cross(vtx2, edgeNormal).' >= -tol;
            within = within | (between & abs(candPts * edgeNormal.') <= sind(width) + tol);
        end

        selected(rows(within)) = true;
    end

    mask(notMask(selected)) = true;
end

%% Subfunctions
function xyz = unitvectors(lon, lat)
    xyz = [cosd(lat(:)) .* cosd(lon(:)), cosd(lat(:)) .* sind(lon(:)), sind(lat(:))];
end
