%% CONVERTSTERICSOURCETS - Converts stored native T/S to CT, SA, pressure, and depth.
% Updates a monthly MAT file in place using the GSW toolbox. Coordinates
% are normalised to latitude-by-level matrices and field type tags are
% updated to describe the converted values.
%
% Syntax
%   convertStericSourceTS(path)
%   convertStericSourceTS(path, ForceNew = true)
%
% Input arguments
%   path - Existing monthly MAT-file path as a character row vector.
%   ForceNew (name-value) - Whether to recompute when converted variables exist.
%       Default value: false.
%
% Input file content
%   T - Temperature field interpreted according to ttype.
%       Size: (nLat, nLon, nLevels), matching S.
%       Unit: °C or K.
%       If the maximum non-NaN value exceeds 200, the whole array is treated as kelvin.
%   S - Salinity field interpreted according to stype.
%       Size: same as T.
%       Unit: psu or g/kg.
%   lon - Longitude vector.
%       Size: (nLon).
%       Unit: °E.
%   lat - Latitude vector.
%       Size: (nLat).
%       Unit: °N.
%   z - Shared level vector or latitude-dependent vertical coordinates.
%       Size: (nLevels) or (nLat, nLevels).
%       Unit: m or dbar.
%   ttype - Temperature interpretation as a text scalar.
%       'T': in-situ; 'PT': potential; 'CT': Conservative Temperature.
%   stype (optional) - Salinity interpretation as a text scalar.
%       'SP': Practical Salinity; 'SA': Absolute Salinity.
%       Default value: 'SP'.
%   ztype - Vertical coordinate interpretation as a text scalar.
%       'p': sea pressure; 'z': positive depth.
%
% Output file content
%   T - Conservative Temperature, stored as single.
%       Size: (nLat, nLon, nLevels).
%       Unit: °C.
%   S - Absolute Salinity, stored as single.
%       Size: same as T.
%       Unit: g/kg.
%   z - Positive depth at each latitude and level, stored as double.
%       Size: (nLat, nLevels).
%       Unit: m.
%   p - Sea pressure at each latitude and level, stored as double.
%       Size: (nLat, nLevels).
%       Unit: dbar.
%   bottom - Deepest supplied level at each latitude, equal to z(:, end).
%       Size: (nLat, 1).
%       Unit: m.
%   ttype - Character row vector set to 'CT'.
%   stype - Character row vector set to 'SA'.
%   ztype - Character row vector set to 'z'.
%
% See also
%   saveSourceStericMonth, matFileHasVariables
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function convertStericSourceTS(path, options)

    arguments (Input)
        path (1, :) char
        options.ForceNew (1, 1) logical = false
    end

    if ~options.ForceNew && ...
            matFileHasVariables(path, {'S', 'T', 'p', 'z', 'bottom'})
        return
    end

    data = load(path);

    if ~isfield(data, 'stype')
        data.stype = 'SP'; % Assumes practical salinity
    end

    [T, S, p, z] = computeStandardTSVars(data.T, data.S, ...
        data.z, data.lon, data.lat, data.ttype, data.stype, data.ztype);
    T = single(T);
    S = single(S);

    bottom = z(:, end);

    % Update field types
    ttype = 'CT';
    stype = 'SA';
    ztype = 'z';

    save(path, 'T', 'S', 'p', 'z', 'bottom', 'ttype', 'stype', 'ztype', '-append');
end

%% Subfunctions
function [CT, SA, p, z] = computeStandardTSVars(T, S, z, lon, lat, ttype, stype, ztype)

    arguments (Input)
        T (:, :, :) {mustBeNumeric}
        S (:, :, :) {mustBeNumeric}
        z (:, :) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonnegative}
        lon {mustBeVector, mustBeNumeric}
        lat {mustBeVector, mustBeNumeric}
        ttype {mustBeTextScalar, mustBeMember(ttype, {'T', 'CT', 'PT'})}
        stype {mustBeTextScalar, mustBeMember(stype, {'SP', 'SA'})}
        ztype {mustBeTextScalar, mustBeMember(ztype, {'z', 'p'})} = 'z'
    end

    arguments (Output)
        CT (:, :, :) {mustBeNumeric}
        SA (:, :, :) {mustBeNumeric}
        p {mustBeNumeric}
        z (:, :) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    %% Validation and array pre-processing
    assert(isequal(size(T), size(S)), 'ULMO:convertStericSourceTS:InvalidInputSize', ...
        ['Temperature and salinity arrays must have the same sizes. ', ...
     'Got (%d, %d, %d) for T and (%d, %d, %d) for salinity instead.'], ...
        size(T), size(S));

    lon = lon(:)';
    lat = lat(:);

    nLvls = size(T, 3);

    if isequal(size(z), [numel(lat), nLvls])
        % Already latitude by level, including a single-level column.
    elseif isvector(z) && numel(z) == nLvls
        z = repmat(z(:)', numel(lat), 1);
    else
        error('ULMO:convertStericSourceTS:InvalidInputSize', ...
            ['Expected z to contain %d shared levels or have size (%d, %d). ', ...
         'Got (%d, %d) instead.'], ...
            nLvls, numel(lat), nLvls, size(z, 1), size(z, 2));
    end

    assert(all(diff(z, 1, 2) > 0, 'all'), ...
        'ULMO:convertStericSourceTS:InvalidVerticalCoordinate', ...
    'Vertical levels must increase strictly at each latitude.');
    z = double(z);

    if size(T, 1) == numel(lat) && size(T, 2) == numel(lon)
        % Correct dimensions, do nothing
    elseif size(T, 2) == numel(lat) && size(T, 1) == numel(lon)
        T = permute(T, [2, 1, 3]);
        S = permute(S, [2, 1, 3]);
    else
        error('ULMO:convertStericSourceTS:InvalidInputSize', ...
            ['T/S first two dimensions must match length of lat/lon. ', ...
         'Got (%d, %d) for the fields but (%d) for latitude and (%d) for longitude.'], ...
            size(T, 1:2), numel(lat), numel(lon));
    end

    %% Variable conversion
    if max(T, [], "all", "omitmissing") > 200 % assume K
        T = T - 273.15; % K -> °C
    end

    switch ztype
        case 'z'
            p = gsw_p_from_z(-z, lat);
        case 'p'
            p = z;
            z = -gsw_z_from_p(p, lat);
    end

    switch stype
        case 'SA'
            SA = S;
        case 'SP'
            SA = nan(size(S), "like", S);

            for k = 1:nLvls
                SA(:, :, k) = gsw_SA_from_SP(squeeze(S(:, :, k)), p(:, k), mod(lon, 360), lat);
            end

    end

    switch ttype
        case 'CT'
            CT = T;
        case 'PT'
            CT = gsw_CT_from_pt(SA, T);
        case 'T'
            CT = nan(size(T), "like", T);

            for k = 1:nLvls
                CT(:, :, k) = gsw_CT_from_t( ...
                    squeeze(SA(:, :, k)), squeeze(T(:, :, k)), p(:, k));
            end

    end

    % Simple validation
    if any(CT < -10, "all")
        warning('ULMO:convertStericSourceTS:NonphysicalData', ...
            ['Non-physical temperature below -10°C is detected and will be masked. ', ...
         'Check the input data.'])
    elseif any(CT > 50, "all")
        warning('ULMO:convertStericSourceTS:NonphysicalData', ...
            ['Non-physical temperature above 50°C is detected and will be masked. ', ...
         'Check the input data.'])
    elseif any(SA < 0, "all")
        warning('ULMO:convertStericSourceTS:NonphysicalData', ...
            ['Non-physical negative salinity is detected and will be masked. ', ...
         'Check the input data.'])
    elseif any(SA > 50, "all")
        warning('ULMO:convertStericSourceTS:NonphysicalData', ...
            ['Non-physical salinity above 50 is detected and will be masked. ', ...
         'Check the input data.'])
    end

    CT(CT < -10 | CT > 50) = NaN;
    SA(SA < 0 | SA > 50) = NaN;
end
