%% DOMAINMASK
% Computes the mask for a given domain defined by a polygon over a grid.
%
% Syntax
%   [mask, lonn, latt] = DOMAINMASK(domain, lon, lat)
%   [mask, lonn, latt] = DOMAINMASK(domain, lonlim, latlim, meshsize)
%   [mask, lonn, latt] = DOMAINMASK(__, "Name", Value)
%
%Input arguments
%   domain - Geographic domain to compute the mask for
%       - GeoDomain object
%       - Nx2 numeric array of (lon, lat) vertices
%       - String function name that returns a Nx2 numeric array of
%           (lon, lat) vertices
%       - Cell array with first element as string function name and
%           subsequent elements as additional arguments to the function
%       - PolyShape object
%       See mustBeGeographicDomain and KERNELCP for more details
%       Arctic and Antarctic masks use geographic coordinates. Their cache
%       names include -rotb to distinguish masks made in rotated coordinates.
%   lon - Vector or 2D array of longitudes (degrees East)
%   lat - Vector or 2D array of latitudes (degrees North)
%   lonlim - 1x2 array specifying the longitude limits [minLon, maxLon]
%       (degrees East)
%   latlim - 1x2 array specifying the latitude limits [minLat, maxLat]
%       (degrees North)
%   meshsize - Scalar specifying the grid spacing (degrees)
%   Format - Format of the output mesh
%       - 'meshgrid' or 'ndgrid' for (lat, lon, t) or (lon, lat, t)
%           ordering.
%       The default format is 'meshgrid'.
%   exteriorBuffer - Expand the finished polygon by this spherical arc in degrees
%       using BUFFERM. Default 0. Positive widths have separate mask caches.
%	ForceNew - Logical flag to force reprocess of the data
%		The default option is false.
%	SaveData - Logical flag to save the data
%		The default option is true.
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is 'soft true', only printing important
%       messages.
%
% Output arguments
%   mask - Logical array where true indicates points inside the domain
%   lonn - 2D array of longitudes corresponding to the mask
%   latt - 2D array of latitudes corresponding to the mask
%
% See also
%   SSH2LONLATT, STERIC2LONLATT
%
% Authored by
%	2025/07/25, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%	2026/11/04 (@williameclee)

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
    else
        if isscalar(unique([diff(unique(lon(:))); diff(unique(lat(:)))]))
            meshsize = abs(unique([diff(unique(lon(:))); diff(unique(lat(:)))]));
            lonlim = [min(lon), max(lon)];
            latlim = [min(lat), max(lat)];
            lon = lonlim(1):meshsize:lonlim(2);
            lat = latlim(1):meshsize:latlim(2);
        end

    end

    if isvector(lon) && isvector(lat)
        lon = lon(:)';
        lat = lat(:);
        [lonn, latt] = meshgrid(lon, lat);
    elseif ~(ismatrix(lon) && ismatrix(lat)) || size(lon, 1) ~= size(lat, 1)
        error('LON and LAT must be vectors or matrices of the same size');
    end

    %% Check for existing data
    dataPath = datapath(domain, meshsize, lonlim, latlim, options.exteriorBuffer);
    exteriorBuffer = options.exteriorBuffer;
    cacheVars = {'mask'};
    if exteriorBuffer > 0
        cacheVars{end+1} = 'exteriorBuffer';
    end
    vars = [cacheVars, {'lon', 'lat'}];

    if ~forceNew && ~(isscalar(dataPath) && isnan(dataPath)) && ...
            exist(dataPath, 'file') && all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, cacheVars{:});
    else
        data = [];
    end

    % Different exact widths can share the same num2str filename.
    if ~isempty(data) && (exteriorBuffer == 0 || data.exteriorBuffer == exteriorBuffer)

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
            domain.Lonlat(round(mean(lonlim) / 180) * 180, ...
                "RotateBack", rotateBack, "OutputFormat", 'polyshape');
    elseif ismatrix(domain) && isnumeric(domain) && size(domain, 2) == 2
        domainPoly = polyshape(domain);
    elseif ischar(domain)
        domainPoly = polyshape(feval(domain));
    elseif iscell(domain)
        domainPoly = polyshape(feval(domain{:}));
    elseif isa(domain, 'polyshape')
        domainPoly = domain;
    else
        error('Unsupported domain type: %s', class(domain));
    end

    if options.exteriorBuffer > 0
        % Buffer the boundary as open arcs. This treats holes and separate
        % components alike and preserves the original polygon interior.
        origin = (min(lonn(:)) + max(lonn(:))) / 2;
        [x, y] = boundary(domainPoly);
        [y, x] = interpm(y, x-origin, 1);
        edges = find(isfinite(x(1:end-1)) & isfinite(x(2:end)));
        arcs = [edges, edges+1].';
        arcLon = [x(arcs); nan(1,numel(edges))];
        arcLat = [y(arcs); nan(1,numel(edges))];
        [y, x] = bufferm(arcLat(:), arcLon(:), options.exteriorBuffer);
        % Cut the spherical result back into the grid's longitude frame.
        [y, x] = maptrimp(y, x, [-90 90], [-180 180]);
        domainPoly = union(domainPoly, polyshape(x+origin, y));
    end

    if min(domainPoly.Vertices(:, 1)) > max(lon)
        error('ULMO:domainmask:InvalidInput', ...
            'Domain longitude range (min: %f) is outside the provided LON range (max: %f).', ...
            min(domainPoly.Vertices(:, 1)), max(lon));
    elseif max(domainPoly.Vertices(:, 1)) < min(lon)
        error('ULMO:domainmask:InvalidInput', ...
            'Domain longitude range (max: %f) is outside the provided LON range (min: %f).', ...
            max(domainPoly.Vertices(:, 1)), min(lon));
    end

    if numel(lonn) < 1e4
        % Run all at once for small grids
        mask = isinterior(domainPoly, lonn(:), latt(:));
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
        mask = arrayfun(@(i) isinterior(domainPoly, lonn(:, i), latt(:, i)), 1:size(mask, 2), ...
            "UniformOutput", false);
        mask = reshape(cell2mat(mask), size(lonn));

        if ~beQuiet
            fprintf(repmat('\b', 1, length(msg) + 1));
            fprintf('took %.1f seconds.\n', toc(t));
        end

    end

    if ischar(dataPath) && saveData
        folder = fileparts(dataPath);
        if ~isfolder(folder)
            mkdir(folder);
        end
        tmp = [tempname(folder), '.mat'];
        cleanup = onCleanup(@() deletePartial(tmp));
        save(tmp, vars{:}, '-v7.3');
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
    elseif ischar(domain)
        domainName = domain;
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
