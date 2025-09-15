%% DOMAINMASK
% Computes the mask for a given domain defined by a polygon over a grid.
%
% See also
%   SSH2LONLATT, STERIC2LONLATT
%
% Authored by
%	2025/07/25, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%	2025/09/15, williameclee@arizona.edu (@williameclee)

function [mask, lonn, latt] = domainmask(varargin)

    arguments (Output)
        mask (:, :) logical
        lonn {mustBeFinite, mustBeMatrix}
        latt {mustBeFinite, mustBeMatrix}
    end

    [domain, domainPoly, lonn, latt, meshSize, lonLim, latLim, fmt, ...
         forceNew, beQuiet, saveData] = ...
        parseinputs(varargin);
    dataPath = datapath(domain, meshSize, lonLim, latLim);

    vars = {'mask' 'lon' 'lat'};

    if ~forceNew && ~(isscalar(dataPath) && isnan(dataPath)) && exist(dataPath, 'file') && ...
            all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, vars{:});

        if ~beQuiet
            fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">domain mask</a>.\n', ...
                mfilename("fullpath"), mfilename, dataPath, dataPath);
        end

        mask = data.mask;
        lonn = data.lon;
        latt = data.lat;

        if strcmpi(fmt, 'ndgrid')
            mask = mask';
            lonn = lonn';
            latt = latt';
        end

        return
    end

    if ~beQuiet
        t = tic;
        msg = 'computing domain mask, this may take a while ... ';
        fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] %s\n', ...
            mfilename("fullpath"), mfilename, msg);
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
            fprintf('saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">domain mask</a>, took %.1f seconds.', ...
                dataPath, dataPath, toc(t));
        end

    end

    if strcmpi(fmt, 'ndgrid')
        mask = mask';
        lonn = lonn';
        latt = latt';
    end

end

% Parse input arguments
function varargout = parseinputs(inputs)
    ip = inputParser;
    addRequired(ip, 'Domain', ...
        @(x) (ismatrix(x) && isnumeric(x) && size(x, 2) == 2) || ...
        isa(x, 'GeoDomain') || isa(x, 'polyshape') || ...
        ischar(x) || iscell(x));
    addOptional(ip, 'Lon', [], @(x) isnumeric(x) && ndims(x) <= 2);
    addOptional(ip, 'Lat', [], @(x) isnumeric(x) && ndims(x) <= 2);
    addParameter(ip, 'MeshSize', nan, @(x) isnumeric(x) && isscalar(x));
    addParameter(ip, 'LonLim', nan, @(x) isnumeric(x) && length(x) <= 2);
    addParameter(ip, 'LatLim', nan, @(x) isnumeric(x) && length(x) <= 2);
    addParameter(ip, 'Format', 'meshgrid', ...
        @(x) ischar(x) && ismember(lower(x), {'meshgrid', 'ndgrid'}));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    parse(ip, inputs{:});

    domain = ip.Results.Domain;
    lon = ip.Results.Lon;
    lat = ip.Results.Lat;
    meshSize = ip.Results.MeshSize;
    lonLim = ip.Results.LonLim;
    latLim = ip.Results.LatLim;
    fmt = lower(ip.Results.Format);
    forceNew = logical(ip.Results.ForceNew);
    beQuiet = logical(ip.Results.BeQuiet);
    saveData = logical(ip.Results.SaveData);

    if ~isnan(meshSize)

        if isnan(lonLim)
            lonLim = [-180, 180] + lonLim;
        end

        if isnan(latLim)
            latLim = [-1, 1] * abs(latLim);
        end

        lon = lonLim(1):meshSize:lonLim(2);
        lat = latLim(1):meshSize:latLim(2);

    elseif isempty(lon) || isempty(lat)
        error('Both LON and LAT must be provided');
    else

        if isscalar(unique([diff(unique(lon(:))); diff(unique(lat(:)))]))
            meshSize = abs(unique([diff(unique(lon(:))); diff(unique(lat(:)))]));
            lonLim = [min(lon), max(lon)];
            latLim = [min(lat), max(lat)];
        end

    end

    if isa(domain, 'GeoDomain')
        domainPoly = domain.Lonlat("OutputFormat", 'polyshape');
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
        [lon, lat] = meshgrid(lon, lat);
    elseif ~(ismatrix(lon) && ismatrix(lat)) || size(lon, 1) ~= size(lat, 1)
        error('LON and LAT must be vectors or matrices of the same size');
    end

    varargout = ...
        {domain, domainPoly, lon, lat, meshSize, lonLim, latLim, fmt, ...
         forceNew, beQuiet, saveData};

end

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
