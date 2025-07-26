%% DOMAINMASK
% Computes the mask for a given domain defined by a polygon over a grid.
% Authored by
%	2025/07/25, williameclee@arizona.edu (@williameclee)

function [mask, lon, lat] = domainmask(varargin)
    [domain, domainPoly, lon, lat, meshSize, lonLim, latLim, ...
         forceNew, beQuiet, saveData] = ...
        parseinputs(varargin);
    dataPath = datapath(domain, meshSize, lonLim, latLim);

    if exist(dataPath, 'file') && ~forceNew
        load(dataPath, 'mask', 'lon', 'lat');

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), dataPath);
        end

        return
    end

    try
        mask = isinterior(domainPoly, lon(:), lat(:));
        mask = reshape(mask, size(lon));
    catch
        mask = false(size(lon));
        nDiv = 64;
        divStart = floor((0:nDiv - 1) * numel(lon) / nDiv) + 1;
        divEnd = [divStart(2:end), numel(lon)];

        for i = 1:nDiv
            mask(divStart(i):divEnd(i)) = ...
                isinterior(domainPoly, lon(divStart(i):divEnd(i)), lat(divStart(i):divEnd(i)));
            % fprintf('Processing division %d / %d (%.0f%%)\n', i, nDiv, i / nDiv * 100);
        end

    end

    if ischar(dataPath) && saveData
        save(dataPath, 'mask', 'lon', 'lat', '-v7.3');

        if ~beQuiet
            fprintf('%s saved %s\n', upper(mfilename), dataPath);
        end

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
        {domain, domainPoly, lon, lat, meshSize, lonLim, latLim, ...
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
