%% readStericCoordinates - Extracts lat/lon/z coordinates from .nc files by matching potential variable names
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function [lon, lat, z] = readStericCoordinates(path)

    arguments (Input)
        path {mustBeFile}
    end

    arguments (Output)
        lon {mustBeNumeric}
        lat {mustBeNumeric}
        z {mustBeNumeric}
    end

    info = ncinfo(path);
    names = {info.Variables.Name};
    coordNames = {'longitude', 'latitude', 'vertical'};
    coordAliases = ...
        {{'lon', 'LONGITUDE', 'longitude'}, ...
          {'lat', 'LATITUDE', 'latitude'}, ...
          {'depth_std', 'LEVEL', 'PRES', 'DEPH', 'depth', 'DEPTH', 'lev'}};
    coords = cell(1, 3);

    % Find matching variables for each coordinate axis
    for k = 1:3
        selected = find(ismember(names, coordAliases{k}), 1);
        assert(~isempty(selected), 'ULMO:readStericCoordinates:VariableNotFound', ...
            'Could not find the coordinate %s in file %s.', coordNames{k}, path);
        coords{k} = double(ncread(path, names{selected}));
        assert(isvector(coords{k}), 'ULMO:readStericCoordinates:InvalidInputSize', ...
            ['Expected gridded source coordinate for coordinates. ', ...
         'But the coordinate %s is not a vector.'], ...
            coordNames{k});
    end

    lon = coords{1}(:)';
    lat = coords{2}(:);
    z = coords{3}(:);
end
