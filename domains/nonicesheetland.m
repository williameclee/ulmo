%% NONICESHEETLAND
% Returns land outside the Greenland and Antarctic ice-sheet domains.
%
% Syntax
%   XY = NONICESHEETLAND(Upscale, Buffer)
%   [XY, polygon] = NONICESHEETLAND(__, 'Name', Value)
%   flag = NONICESHEETLAND('rotated')
%
% Input arguments
%   Upscale, Buffer, Latlim, MoreBuffers, LonOrigin - Common domain options.
%       Defaults: 0, 0, [-90,90], {}, and 180 degrees respectively.
%       Positive Buffer expands the land mask and the Greenland exclusion.
%       Latlim clips the resulting land; it does not turn polar oceans into land.
%   ForceNew, SaveData, BeQuiet - Forwarded to ALLOCEANS for its coastline cache.
%       Defaults: false, true, false. This function has no separate cache.
%
% Output arguments
%   XY - Closed longitude/latitude boundaries, separated by NaN; DOUBLE, N-by-2.
%   polygon - Planar POLYSHAPE in the requested longitude window.
%   flag - False: coordinates are already in their geographic orientation.
%
% Notes
%   Subtracts the ocean mask and buffered Greenland domain from the globe.
%   The Antarctic exclusion uses ICESHEETPOLY('AIS'): land south of 60 S,
%   including nearby islands. This avoids rotated Antarctic coordinates. The name
%   excludes the two ice-sheet domains, not all ice: mountain glaciers and
%   other land ice remain included. Use as a LandDomain in GRACE2FINGERPRINT.
%   Independent buffered domains need not form an exact additive partition.
%
% See also
%   GEODOMAIN, ALLOCEANS, GREENLAND, ICESHEETPOLY, GRACE2FINGERPRINT
%
% Created by
%   2026/10/01, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/01, En-Chi Lee (williameclee@arizona.edu)

function varargout = nonicesheetland(varargin)
    if nargin == 1 && strcmpi(varargin{1}, 'rotated')
        varargout = {false};
        return
    end
    % Adapt the legacy positional domain interface to validated named options.
    names = {'Upscale', 'Buffer', 'Latlim', 'MoreBuffers', 'LonOrigin', 'RotateBack'};
    positional = {};
    while ~isempty(varargin) && ~ischar(varargin{1}) && ~isstring(varargin{1})
        if numel(positional) / 2 == numel(names)
            error('ULMO:nonicesheetland:TooManyInputs', 'Too many positional inputs.');
        end
        positional = [positional, names(numel(positional) / 2 + 1), varargin(1)];
        varargin(1) = [];
    end
    [xy, polygon] = landpolygon(positional{:}, varargin{:});
    varargout = returncoastoutputs(nargout, xy, polygon);
end

function [xy, polygon] = landpolygon(options)
    arguments (Input)
        options.Upscale {mustBeNumeric, mustBeScalarOrEmpty, mustBeFinite, mustBeNonnegative} = 0
        options.Buffer {mustBeNumeric, mustBeScalarOrEmpty, mustBeFinite} = 0
        options.Latlim {mustBeNumeric, mustBeReal, mustBeValidLatlim} = [-90,90]
        options.MoreBuffers cell = {}
        options.LonOrigin {mustBeNumeric, mustBeScalarOrEmpty, mustBeFinite} = 180
        options.ForceNew (1,1) logical = false
        options.SaveData (1,1) logical = true
        options.BeQuiet (1,1) logical = false
        % GeoDomain passes these common options; this domain is unrotated.
        options.RotateBack (1,1) logical = false
        options.NearBy = []
    end
    upscale = options.Upscale;
    if isempty(upscale) || upscale == 1, upscale = 0; end
    buf = options.Buffer;
    if isempty(buf), buf = 0; end
    latlim = options.Latlim;
    if isempty(latlim) || any(isnan(latlim)), latlim = [-90,90]; end
    lonOrigin = options.LonOrigin;
    if isempty(lonOrigin), lonOrigin = 180; end
    world = polyshape(lonOrigin + [-180,180,180,-180], [-90,-90,90,90]);
    [~, ocean] = alloceans('Upscale', upscale, 'Buffer', buf, ...
        'Latlim', [-90,90], 'MoreBuffers', options.MoreBuffers, 'LonOrigin', lonOrigin, ...
        'ForceNew', options.ForceNew, 'SaveData', options.SaveData, 'BeQuiet', options.BeQuiet);
    polygon = subtract(world, ocean);
    greenlandDomain = GeoDomain('greenland', 'Buffer', buf, 'Upscale', upscale);
    xy = greenlandDomain.Lonlat([]);
    polygon = subtract(polygon, geographicpolygon(xy, world));
    antarctic = translate(icesheetPoly('AIS'), [lonOrigin-180,0]);
    polygon = subtract(polygon, antarctic);
    if isscalar(latlim), latlim = [-abs(latlim),abs(latlim)]; end
    clip = polyshape(lonOrigin + [-180,180,180,-180], latlim([1,1,2,2]));
    polygon = intersect(polygon, clip);
    [lon,lat] = boundary(polygon);
    xy = [lon,lat];
end

function polygon = geographicpolygon(xy, world)
% Unwrap Greenland before planar clipping at the requested longitude seam.
    lon = rad2deg(unwrap(deg2rad(xy(:,1))));
    lat = xy(:,2);
    % Centre the unwrapped ring on the requested window before replication.
    centre = mean(world.Vertices(:,1));
    lon = lon + 360*round((centre-mean(lon,'omitmissing'))/360);
    base = polyshape(lon,lat);
    polygon = union(base,translate(base,[-360,0]));
    polygon = union(polygon,translate(base,[360,0]));
    polygon = intersect(polygon,world);
end

function mustBeValidLatlim(value)
    if isempty(value) || any(isnan(value)), return; end
    mustBeVector(value);
    mustBeFinite(value);
    mustBeInRange(value, -90, 90);
    if numel(value) > 2 || (numel(value) == 2 && value(1) > value(2))
        error('ULMO:nonicesheetland:InvalidLatlim', ...
            'Latlim must be a scalar or an increasing pair of latitudes.');
    end
end
