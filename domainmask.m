%% DOMAINMASK - Compute a geographic domain mask on a longitude-latitude grid.
%
% Syntax
%   [mask, lonn, latt] = DOMAINMASK(domain, lon, lat)
%   [mask, lonn, latt] = DOMAINMASK(domain, [], [], lonlim, latlim, meshsize)
%   [mask, lonn, latt] = DOMAINMASK(__, "Name", value)
%
% Input arguments
%   domain - Geographic domain to compute the mask for.
%       Accepts a GeoDomain object, an (n, 2) numeric array of longitude and
%       latitude vertices, a character vector naming a boundary function,
%       or a polyshape object. Arctic and Antarctic GeoDomain objects and
%       function names use geographic coordinates; their mask cache names
%       include -rotb to distinguish them from rotated masks.
%   lon - Longitude coordinates, or [] when specifying limits and spacing.
%       Any longitude convention is accepted. Input order is preserved.
%       Unit: degrees east.
%   lat - Latitude coordinates, or [] when specifying limits and spacing.
%       Unit: degrees north.
%   lonlim (optional) - Longitude limits [minimum, maximum].
%       Used with meshsize to construct the grid.
%       Size: (1, 2).
%       Unit: degrees east.
%       Default limits: [NaN, NaN]; inferred from regular coordinate vectors.
%   latlim (optional) - Latitude limits [minimum, maximum].
%       Used with meshsize to construct the grid.
%       Size: (1, 2).
%       Unit: degrees north.
%       Default limits: [NaN, NaN]; inferred from regular coordinate vectors.
%   meshsize (optional) - Grid spacing in both longitude and latitude.
%       Supply both lonlim and latlim when setting this argument.
%       Unit: degrees.
%       Default spacing: NaN; inferred from regular coordinate vectors.
%   Format (name-value) - Output array layout.
%       'meshgrid' uses (latitude, longitude) ordering; 'ndgrid' uses
%       (longitude, latitude) ordering.
%       Default format: 'meshgrid'.
%   exteriorBuffer (name-value) - Spherical buffer around the finished polygon.
%       Includes grid points within this angular distance of the boundary.
%       Boundary edges are densified linearly at 1-degree spacing, then treated
%       as minor great-circle segments. Distances include segment endpoints.
%       Must be between 0 and 90. Positive widths use separate mask caches.
%       Unit: degrees of spherical arc.
%       Default width: 0.
%   Forcenew (name-value) - Whether to recompute an existing mask.
%       Default value: false.
%   SaveData (name-value) - Whether to save a cache for a regular grid.
%       Named domains are cacheable; numeric polygons are cacheable when
%       exteriorBuffer is positive. Polyshape inputs are not cached.
%       Default value: true.
%   BeQuiet (name-value) - Whether to suppress progress messages.
%       Default value: false.
%   CallChain (name-value) - Calling function names for progress messages.
%       Default value: {}.
%
% Output arguments
%   mask - Logical array identifying points in the buffered domain.
%       Size: (nLat, nLon) for 'meshgrid', or (nLon, nLat) for 'ndgrid'.
%   lonn - Longitude coordinates corresponding to the mask.
%       Size: same as mask.
%       Unit: degrees east.
%   latt - Latitude coordinates corresponding to the mask.
%       Size: same as mask.
%       Unit: degrees north.
%
% Notes
%   Polygon membership uses a fixed 0-to-360-degree longitude window and
%   temporary modulo-360 query coordinates. Returned coordinates retain the
%   input convention and order. Both longitude seam endpoints are included
%   when either is on the polygon boundary. Averaging callers must count a
%   repeated cyclic endpoint only once.
%   Supplied polygon rings retain their planar topology and are split at the
%   longitude seam, rather than wrapping their vertices independently.
%   Mask caches written before this longitude convention are ignored. Cached
%   projections or regional means derived from old masks need separate rebuilds.
%
% See also
%   GEODOMAIN, BUFFERGRIDMASK, SSH2LONLATT, STERIC2LONLATT
%
% Created by
%   2025/07/25, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/05, OpenAI Codex

function [mask, lonn, latt] = domainmask(domain, lon, lat, lonlim, latlim, meshsize, options)

    arguments (Input)
        domain {mustBeGeographicDomain}
        lon (:, :) double {mustBeFinite}
        lat (:, :) double {mustBeFinite}
        lonlim (1, 2) double {mustBeNumeric} = nan
        latlim (1, 2) double {mustBeNumeric} = nan
        meshsize (1, 1) double {mustBeNumeric} = nan
        options.Format (1, :) char {mustBeMember(options.Format, {'meshgrid', 'ndgrid'})} = 'meshgrid'
        options.exteriorBuffer (1, 1) double {mustBeFinite, mustBeNonnegative, mustBeLessThanOrEqual(options.exteriorBuffer, 90)} = 0
        options.Forcenew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.CallChain cell = {};
    end

    arguments (Output)
        mask (:, :) logical
        lonn {mustBeFinite, mustBeMatrix}
        latt {mustBeFinite, mustBeMatrix}
    end

    fmt = lower(options.Format);
    forceNew = options.Forcenew;
    beQuiet = options.BeQuiet;
    saveData = options.SaveData;
    callChain = [options.CallChain, {mfilename}];

    if (ischar(domain) || isstring(domain)) && ...
            ismember(lower(string(domain)), ["arctic", "antarctica"])
        domain = GeoDomain(domain);
    end
    rotateBack = isa(domain, 'GeoDomain') && ...
        ismember(domain.Domain, ["arctic", "antarctica"]);

    if ~isnan(meshsize)

        if isnan(lonlim)
            lonlim = [-180, 180] + lonlim;
        end

        if isnan(latlim)
            latlim = [-1, 1] * abs(latlim);
        end

        lon = lonlim(1):meshsize:lonlim(2);
        lat = latlim(1):meshsize:latlim(2);
    elseif isempty(lon) || isempty(lat)
        error('Both LON and LAT must be provided');
    elseif isvector(lon) && isvector(lat)
        if isscalar(unique([diff(unique(lon(:))); diff(unique(lat(:)))]))
            meshsize = abs(unique([diff(unique(lon(:))); diff(unique(lat(:)))]));
            lonlim = [min(lon(:)), max(lon(:))];
            latlim = [min(lat(:)), max(lat(:))];
        end

    end

    if isvector(lon) && isvector(lat)
        lon = lon(:)';
        lat = lat(:);
        [lonn, latt] = meshgrid(lon, lat);
    elseif ~isequal(size(lon), size(lat))
        error('LON and LAT must be vectors or matrices of the same size');
    else
        lonn = lon;
        latt = lat;
    end

    %% Check for existing data
    dataPath = datapath(domain, meshsize, lonlim, latlim, options.exteriorBuffer);
    exteriorBuffer = options.exteriorBuffer;
    maskVersion = 1;
    cacheVars = {'mask', 'maskVersion', 'lon', 'lat'};
    if exteriorBuffer > 0
        cacheVars{end+1} = 'exteriorBuffer';
    end
    vars = cacheVars;

    if ~forceNew && ~(isscalar(dataPath) && isnan(dataPath)) && ...
            exist(dataPath, 'file') && all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, cacheVars{:});
    else
        data = [];
    end

    % Different exact widths can share the same num2str filename.
    if ~isempty(data) && isequal(data.maskVersion, maskVersion) && ...
            isequal(data.lon, lon) && isequal(data.lat, lat) && ...
            isequal(size(data.mask), size(lonn)) && ...
            (exteriorBuffer == 0 || data.exteriorBuffer == exteriorBuffer)

        if ~beQuiet
            fprintf('[ULMO>%s] Loaded %s.\n', ...
                callchaintext(callChain), filehref(dataPath, 'domain mask'));
        end

        mask = data.mask;

        if nargout > 0

            if strcmpi(fmt, 'ndgrid')
                mask = mask';
                lonn = lonn';
                latt = latt';
            end

        else
            plotdomainmask(mask, lonn, latt)
        end

        return
    end

    %% Main
    if isa(domain, 'GeoDomain')
        % Native grids require geographic coordinates for rotated polar domains.
        domainPoly = ...
            domain.Lonlat(180, ...
                "RotateBack", rotateBack, "OutputFormat", 'polyshape');
    elseif ismatrix(domain) && isnumeric(domain) && size(domain, 2) == 2
        domainPoly = polyshape(domain);
    elseif ischar(domain) || isstring(domain)
        domainPoly = polyshape(feval(domain));
    elseif iscell(domain)
        domainPoly = polyshape(feval(domain{:}));
    elseif isa(domain, 'polyshape')
        domainPoly = domain;
    else
        error('Unsupported domain type: %s', class(domain));
    end

    domainPoly = wrapPolygon(domainPoly);
    queryLon = mod(lonn, 360);

    if numel(lonn) < 1e4
        % Run all at once for small grids
        mask = isinterior(domainPoly, queryLon(:), latt(:));
        mask = reshape(mask, size(lonn));
    else
        % For large grids, use arrayfyn (can be parallelised)
        if ~beQuiet
            t = tic;
            msg = 'this may take a while ... ';
            fprintf('[ULMO>%s] Computing domain mask in parts, %s\n', ...
                callchaintext(callChain), msg);
        end

        mask = false(size(lonn));
        mask = arrayfun(@(i) isinterior(domainPoly, queryLon(:, i), latt(:, i)), 1:size(mask, 2), ...
            "UniformOutput", false);
        mask = reshape(cell2mat(mask), size(lonn));

        if ~beQuiet
            fprintf(repmat('\b', 1, length(msg) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

    end

    % A boundary at 360 degrees is physically the same as one at zero.
    seam = queryLon == 0;
    if any(seam, 'all')
        mask(seam) = mask(seam) | isinterior(domainPoly, 360*ones(nnz(seam), 1), latt(seam));
    end

    if exteriorBuffer > 0
        mask = buffergridmask(domainPoly, queryLon, latt, mask, exteriorBuffer);
    end

    if ischar(dataPath) && saveData
        folder = fileparts(dataPath);
        if ~isfolder(folder)
            mkdir(folder);
        end
        tmp = [tempname(folder), '.mat'];
        cleanup = onCleanup(@() deletePartial(tmp));
        creationDate = datetime('now', 'TimeZone', 'UTC');
        save(tmp, vars{:}, 'creationDate', '-v7.3');
        movefile(tmp, dataPath, 'f');

        if ~beQuiet
            fprintf('[ULMO>%s] Saved %s.\n', ...
                callchaintext(callChain), filehref(dataPath, 'domainMask'));
        end

    end

    if nargout > 0

        if strcmpi(fmt, 'ndgrid')
            mask = mask';
            lonn = lonn';
            latt = latt';
        end

    else
        plotdomainmask(mask, lonn, latt)
    end

end

%% Subfunctions
function wrapped = wrapPolygon(polygon)
    % Clip translated copies as polygons so holes and full-width bands survive.
    vertices = polygon.Vertices;
    if isempty(vertices)
        wrapped = polygon;
        return
    end
    west = min(vertices(:, 1));
    east = max(vertices(:, 1));
    if west >= 0 && east <= 360
        wrapped = polygon;
        return
    end
    south = min(vertices(:, 2));
    north = max(vertices(:, 2));
    window = polyshape([0, 360, 360, 0], [south, south, north, north]);
    wrapped = polyshape();
    for iCopy = (floor(-east/360)+1):(ceil((360-west)/360)-1)
        piece = intersect(translate(polygon, [360*iCopy, 0]), window);
        wrapped = union(wrapped, piece);
    end
end

function dataPath = datapath(domain, meshSize, lonLim, latLim, exteriorBuffer)

    if isnan(meshSize)
        % Data cannot be saved if irregular
        dataPath = nan;
        return
    end

    dataFolder = fullfile(getenv('IFILES'), 'DOMAINMASK');

    if isa(domain, 'GeoDomain')
        domainName = domain.Id;
        if ismember(domain.Domain, ["arctic", "antarctica"])
            domainName = [domainName, '-rotb'];
        end
    elseif iscell(domain)
        domainName = strjoin(string(domain), '_');
    elseif ischar(domain) || isstring(domain)
        domainName = char(domain);
    elseif isnumeric(domain) && exteriorBuffer > 0
        domainName = hash(domain, 'sha1');
    else
        dataPath = nan;
        return
    end
    if exteriorBuffer > 0
        domainName = sprintf('%s-E%s', domainName, num2str(exteriorBuffer));
    end

    dataFile = sprintf('DOMAINMASK-%s-H%s-Box%s_%s_%s_%s.mat', ...
        domainName, num2str(meshSize), ...
        num2str(lonLim(1)), num2str(lonLim(2)), num2str(latLim(1)), num2str(latLim(2)));
    dataPath = fullfile(dataFolder, dataFile);
end

function deletePartial(path)
    if isfile(path)
        delete(path);
    end
end

function plotdomainmask(mask, lonn, latt)
    figure()

    set(gcf, "Name", 'Domain mask', 'NumberTitle', 'off');

    imagesc(lonn(1, :), latt(:, 1), double(mask));

    set(gca, "YDir", 'normal')
    axis equal tight

end
