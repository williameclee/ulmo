%% SAVESOURCESTERICMONTH - Writes native monthly T/S fields to a MAT file.
% Stores source values and interpretation tags before conversion, using
% latitude/longitude/level field order and a shared vertical coordinate.
%
% Syntax
%   saveSourceStericMonth(path, T, S, z, lon, lat, date, ttype)
%   saveSourceStericMonth(path, T, S, z, lon, lat, date, ttype, stype, ztype)
%   saveSourceStericMonth(__, ForceNew = true)
%
% Input arguments
%   path - Output MAT-file path as a text scalar.
%       The file may be new, but its parent folder must already exist.
%   T - Temperature array.
%       Size: (nLat, nLon, nLevels).
%       Unit: °C or K.
%   S - Salinity array.
%       Size: same as T.
%       Unit: dimensionless Practical Salinity or Absolute Salinity in g/kg.
%   z - Vertical levels.
%       Size: nLevels elements, row or column.
%       Unit: m for depth, or dbar for sea pressure.
%   lon - Longitude vector
%       Size: nLon elements, row or column.
%       Unit: °E.
%   lat - Latitude vector.
%       Size: nLat elements, row or column.
%       Unit: °N.
%   date - Scalar datetime identifying the source month.
%   ttype - Text scalar describing T.
%       'T' for in-situ temperature, 'PT' for potential temperature,
%       or 'CT' for conservative temperature.
%   stype (optional) - Text scalar describing S.
%       'SP' for practical salinity, or 'SA' for absolute salinity.
%       Default value: 'SP'.
%   ztype (optional) - Text scalar describing z.
%       'z' for positive depth, or 'p' for sea pressure.
%       Default value: 'z'.
%   ForceNew (name-value) - Whether to write even when matching metadata exists.
%       Default value: false.
%
% Input file content
%   Optional existing file at path, inspected only when ForceNew=false.
%   No input file is required for an initial write.
%   T, S - Previously stored native temperature and salinity fields.
%       Both names must exist for reuse; their values and sizes are not checked.
%   lon - Longitude vector.
%       Size: (1, nLon). 
%       Unit: °E.
%   lat - Latitude vector.
%       Size: (nLat, 1). 
%       Unit: °N.
%   z - Vertical levels.
%       Size: (nLevels, 1). 
%       Unit: m or dbar according to ztype.
%   date - Stored scalar datetime, compared with the supplied date.
%   ttype, stype, ztype - Stored text-scalar interpretation tags.
%       Each must match the corresponding input argument for reuse.
%
% Output arguments
%   None.
%
% Output file content
%   T - Native temperature field, saved without conversion or type casting.
%       Size: (nLat, nLon, nLevels).
%       Unit: degrees Celsius or kelvin, as supplied.
%   S - Native salinity field, saved without conversion or type casting.
%       Size: same as T.
%       Unit: dimensionless Practical Salinity or Absolute Salinity in g/kg.
%   lon - Longitude row vector.
%       Size: (1, nLon). Unit: degrees east.
%   lat - Latitude column vector.
%       Size: (nLat, 1). Unit: degrees north.
%   z - Shared vertical levels, flattened to a column.
%       Size: (nLevels, 1). Unit: m or dbar according to ztype.
%   date - Scalar datetime saved as supplied, without midpoint adjustment.
%   ttype - Temperature tag saved as supplied: 'T', 'PT', or 'CT'.
%   stype - Salinity tag saved as supplied or defaulted: 'SP' or 'SA'.
%   ztype - Vertical tag saved as supplied or defaulted: 'z' or 'p'.
%
% Notes
%   When metadata matches, the existing file is left unchanged. T/S values
%   are not compared on reuse; use ForceNew=true when source values change.
%   A write creates a version 7.3 MAT file and replaces any existing file,
%   including derived variables appended by later processing stages.
%
% See also
%   matFileHasVariables, convertStericSourceTS, extractStericMonthDate
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function saveSourceStericMonth(path, T, S, z, lon, lat, date, ttype, stype, ztype, options)

    arguments (Input)
        path {mustBeTextScalar}
        T (:, :, :) {mustBeNumeric}
        S (:, :, :) {mustBeNumeric}
        z {mustBeNumeric}
        lon {mustBeNumeric, mustBeVector}
        lat {mustBeNumeric, mustBeVector}
        date (1, 1) datetime
        ttype {mustBeTextScalar, mustBeMember(ttype, {'T', 'CT', 'PT'})}
        stype {mustBeTextScalar, mustBeMember(stype, {'SP', 'SA'})} = 'SP'
        ztype {mustBeTextScalar, mustBeMember(ztype, {'z', 'p'})} = 'z'
        options.ForceNew (1, 1) logical = false
    end

    if ~options.ForceNew && matFileHasVariables(path, ...
            {'date', 'lat', 'lon', 'z', 'T', 'S', 'ttype', 'stype', 'ztype'})
        data = load(path, 'date', 'lat', 'lon', 'z', 'ttype', 'stype', 'ztype');

        if isequal(data.date, date) && isequal(data.lon, lon(:)') && ...
                isequal(data.lat, lat(:)) && isequal(data.z, z(:)) && ...
                strcmp(data.ttype, ttype) && strcmp(data.stype, stype) && strcmp(data.ztype, ztype)
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
    save(path, 'T', 'S', 'lon', 'lat', 'z', 'date', 'ttype', 'stype', 'ztype', '-v7.3');
end
