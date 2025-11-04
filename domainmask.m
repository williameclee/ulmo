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
%	2025/11/03, williameclee@arizona.edu (@williameclee)

function [mask, lonn, latt] = domainmask(domain, lon, lat, lonlim, latlim, meshsize, options)

    arguments (Input)
        domain {mustBeGeographicDomain}
        lon (:, :) double {mustBeFinite}
        lat (:, :) double {mustBeFinite}
        lonlim (1, 2) double {mustBeNumeric} = nan
        latlim (1, 2) double {mustBeNumeric} = nan
        meshsize (1, 1) double {mustBeNumeric} = nan
        options.Format (1, :) char {mustBeMember(options.Format, {'meshgrid', 'ndgrid'})} = 'meshgrid'
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

    if isa(domain, 'GeoDomain')
        domainPoly = domain.Lonlat(mean(lonlim), "OutputFormat", 'polyshape');
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

    if isvector(lon) && isvector(lat)
        lon = lon(:)';
        lat = lat(:);
        [lonn, latt] = meshgrid(lon, lat);
    elseif ~(ismatrix(lon) && ismatrix(lat)) || size(lon, 1) ~= size(lat, 1)
        error('LON and LAT must be vectors or matrices of the same size');
    end

    %% Main
    dataPath = datapath(domain, meshsize, lonlim, latlim);
    vars = {'mask' 'lon' 'lat'};

    if ~forceNew && ~(isscalar(dataPath) && isnan(dataPath)) && exist(dataPath, 'file') && ...
            all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, vars{:});

        if ~beQuiet
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">domain mask</a>.\n', ...
                callchaintext(callChain), dataPath, dataPath);
        end

        mask = data.mask;
        lonn = data.lon;
        latt = data.lat;

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

    if ~beQuiet
        t = tic;
        msg = 'computing domain mask, this may take a while ... ';
        fprintf('[ULMO>%s] %s\n', callchaintext(callChain), msg);
    end

    try
        mask = isinterior(domainPoly, lonn(:), latt(:));
        mask = reshape(mask, size(lonn));
    catch
        mask = false(size(lonn));
        nDiv = 64;
        divStart = floor((0:nDiv - 1) * numel(lonn) / nDiv) + 1;
        divEnd = [divStart(2:end), numel(lonn)];

        for i = 1:nDiv
            iRange = divStart(i):divEnd(i);
            mask(iRange) = ...
                isinterior(domainPoly, lonn(iRange), latt(iRange));
        end

    end

    if ischar(dataPath) && saveData
        save(dataPath, vars{:}, '-v7.3');

        if ~beQuiet
            fprintf(repmat('\b', 1, length(msg) + 1));
            fprintf('saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">domain mask</a>, took %.1f seconds.\n', ...
                dataPath, dataPath, toc(t));
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
function dataPath = datapath(domain, meshSize, lonLim, latLim)

    if ~(isa(domain, 'GeoDomain') || isa(domain, 'char') || isa(domain, 'cell')) ...
            || isnan(meshSize)
        % Data cannot be saved if irregular
        dataPath = nan;
        return
    end

    dataFolder = fullfile(getenv('IFILES'), 'DOMAINMASK');

    if ~exist(dataFolder, 'dir')
        mkdir(dataFolder);
    end

    if isa(domain, 'GeoDomain')
        domainName = domain.Id;
    elseif iscell(domain)
        domainName = strjoin(string(domain), '_');
    end

    dataFile = sprintf('DOMAINMASK-%s-H%s-Box%s_%s_%s_%s.mat', ...
        domainName, num2str(meshSize), ...
        num2str(lonLim(1)), num2str(lonLim(2)), num2str(latLim(1)), num2str(latLim(2)));
    dataPath = fullfile(dataFolder, dataFile);
end

function plotdomainmask(mask, lonn, latt)
    figure()

    set(gcf, "Name", 'Domain mask', 'NumberTitle', 'off');

    imagesc(lonn(1, :), latt(:, 1), double(mask));

    set(gca, "YDir", 'normal')
    axis equal tight

end
