%% READSTERICSOURCEFIELD - Reads a NetCDF field in latitude/longitude/level order.
% Selects one time index, identifies spatial axes by dimension names,
% and removes other singleton dimensions.
%
% Syntax
%   val = readStericSourceField(path, name)
%   val = readStericSourceField(path, name, tstart)
%
% Input arguments
%   path - Path to an existing NetCDF file.
%   name - NetCDF variable name.
%   tstart (optional) - 1-based time index. 
%       Ignored if there is no recognised time dimension.
%       Default value: 1.
%
% Output arguments
%   val - Numeric field with 
%       Size: (nLat, nLon, nLevels).
%
% Notes
%   Spatial dimension names are case-sensitive:
%   - latitude: lat, LATITUDE, latitude.
%   - longitude: lon, LONGITUDE, longitude.
%   - vertical: depth_std, LEVEL, PRES, DEPH, depth, DEPTH, lev.
%   Exactly one dimension must match each spatial axis. 
%   Time dimensions match time or t case-insensitively.
%   Each recognised time dimension is sliced at tstart. All remaining 
%   dimensions must have length one.
%
% See also
%   readStericCoordinates, saveSourceStericMonth
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function val = readStericSourceField(path, name, tstart)

    arguments (Input)
        path {mustBeFile}
        name {mustBeTextScalar}
        tstart (1, 1) {mustBeInteger, mustBePositive} = 1
    end

    arguments (Output)
        val {mustBeNumeric}
    end

    fieldInfo = ncinfo(path, name);
    names = {fieldInfo.Dimensions.Name};
    coordNames = {'longitude', 'latitude', 'vertical'};
    coordAliases = ...
        {{'lat', 'LATITUDE', 'latitude'}, ...
          {'lon', 'LONGITUDE', 'longitude'}, ...
          {'depth_std', 'LEVEL', 'PRES', 'DEPH', 'depth', 'DEPTH', 'lev'}};

    start = ones(1, numel(names));
    cnt = [fieldInfo.Dimensions.Length];
    time = find(ismember(lower(string(names)), ["time", "t"]));
    cnt(time) = 1;
    start(time) = tstart;
    order = zeros(1, 3);

    for k = 1:3
        axisLoc = find(ismember(names, coordAliases{k}));
        assert(isscalar(axisLoc), 'ULMO:readStericSourceField:VariableNotFound', ...
            'Could not find the coordinate %s, or found multiple potentially matching fields in file %s.', ...
            coordNames{k}, path);
        order(k) = axisLoc; % Order fields as [lat, lon, level].
    end

    other = setdiff(1:numel(names), order, 'stable');
    assert(all(cnt(other) == 1), 'Unexpected nonsingleton source dimension.');
    val = ncread(path, name, start, cnt);
    val = reshape(val, cnt);
    val = permute(val, [order, other]);
    val = reshape(val, cnt(order));
end
