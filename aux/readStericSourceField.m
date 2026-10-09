%% readStericSourceField - Reads native NetCDF dimensions into (lat, lon, depth) order.
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
        order(k) = axisLoc; % Order the fields in the [lon, lat, z] dimension
    end

    other = setdiff(1:numel(names), order, 'stable');
    assert(all(cnt(other) == 1), 'Unexpected nonsingleton source dimension.');
    val = ncread(path, name, start, cnt);
    val = reshape(val, cnt);
    val = permute(val, [order, other]);
    val = reshape(val, cnt(order));
end
