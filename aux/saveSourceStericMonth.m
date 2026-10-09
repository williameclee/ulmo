%% saveSourceStericMonth - Saves native monthly source fields.
% Later stages append scientific variables.
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function saveSourceStericMonth(path, T, S, z, lon, lat, date, ttype, stype, ztype, options)

    arguments (Input)
        path {mustBeFile}
        T (:, :, :) {mustBeNumeric}
        S (:, :, :) {mustBeNumeric}
        z {mustBeNumeric}
        lon {mustBeNumeric, mustBeVector}
        lat {mustBeNumeric, mustBeVector}
        date (1, 1) duration
        ttype {mustBeTextScalar, mustBeMember(ttype, {'T', 'CT', 'PT'})}
        stype {mustBeTextScalar, mustBeMember(stype, {'SP', 'SA'})} = 'SP'
        ztype {mustBeTextScalar, mustBeMember(ztype, {'z', 'p'})} = 'z'
        options.ForceNew (1, 1) logical = false
    end

    if ~options.ForceNew && matFileHasVariables(path, ...
            {'date', 'lat', 'lon', 'z', 'T', 'S', 'ttype', 'ztype', 'stype'})
        data = load(path, 'date', 'lon', 'lat', 'z', 'stype');

        if isequal(data.date, date) && isequal(data.lon, lon(:)') && ...
                isequal(data.lat, lat(:)) && isequal(data.z, z(:)) && ...
                strcmp(data.stype, stype)
            return
        end

    end

    lon = lon(:)';
    lat = lat(:);
    z = z(:);
    assert(isequal(size(T), size(S)), ...
        'ULMO:saveStericSourceMonth:InvalidInputSize', ...
        ['Temperature and salinity arrays must have the same sizes. ', ...
     'Got (%d, %d, %d) for temperature and (%d, %d, %d) for salinity instead.'], ...
        size(T), size(S));
    assert(size(T, 1) == numel(lat) && size(T, 2) == numel(lon) && ...
        size(T, 3) == numel(z), 'ULMO:saveStericSourceMonth:InvalidInputSize', ...
        ['Source field and coordinate dimensions differ. ', ...
     'Expected dimension is (lat: %d, lon: %d, z: %d), but got (%d, %d, %d).'], ...
        numel(lat), numel(lon), numel(z), size(T));
    T(T < -5 | T > 50) = NaN;
    S(S < 0 | S > 50) = NaN;
    save(path, 'T', 'S', 'lon', 'lat', 'z', 'date', 'ttype', 'stype', 'ztype', '-v7.3');
end
