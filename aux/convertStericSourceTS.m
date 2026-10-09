%% convertStericSourceTS - Converts T/S units to standard conservative T and absolute salinity.
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function convertStericSourceTS(path, options)

    arguments (Input)
        path (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    if ~options.ForceNew && ...
            matFileHasVariables(path, {'salinity', 'consTemp', 'pres', 'depth', 'bottom'})
        return
    end

    data = load(path);

    if ~isfield(data, 'stype')
        data.stype = 'SP'; % Assumes practical salinity
    end

    [T, S, p] = computeStandardTSVars(data.T, data.S, ...
        data.z, data.lon, data.lat, data.stype, data.stype, data.ztype);
    T = single(T);
    S = single(S);

    if strcmp(data.ztype, 'p')
        z = -gsw_z_from_p(repmat(data.z(:)', numel(data.lat), 1), data.lat(:));
        bottom = z(:, end);
    else
        z = data.z(:);
        bottom = z(end);
    end

    save(path, 'T', 'S', 'p', 'z', 'bottom', '-append');
end

%% Subfunctions
function [CT, SA, p] = computeStandardTSVars(T, S, z, lon, lat, ttype, stype, ztype)

    arguments (Input)
        T (:, :, :) {mustBeNumeric}
        S (:, :, :) {mustBeNumeric}
        z {mustBeVector}
        lon {mustBeVector}
        lat {mustBeVector}
        ttype {mustBeTextScalar, mustBeMember(ttype, {'T', 'CT', 'PT'})}
        stype {mustBeTextScalar, mustBeMember(stype, {'SP', 'SA'})}
        ztype {mustBeTextScalar, mustBeMember(ztype, {'z', 'p'})} = 'z'
    end

    arguments (Output)
        CT (:, :, :) {mustBeNumeric}
        SA (:, :, :) {mustBeNumeric}
        p {mustBeNumeric}
    end

    %% Validation
    assert(isequal(size(T), size(S)), 'ULMO:convertStericSourceTS:InvalidInputSize', ...
        ['Temperature and salinity arrays must have the same sizes. ', ...
     'Got (%d, %d, %d) for T and (%d, %d, %d) for salinity instead.'], ...
        size(T), size(S));

    lon = lon(:)';
    lat = lat(:);

    if size(T, 3) ~= numel(z)
        error('ULMO:convertStericSourceTS:InvalidInputSize', ...
            ['Third dimension of T and S arrays must match the layers of depth. ', ...
         'Got %d and %d.'], ...
            size(T, 3), numel(z));
    elseif size(T, 1) == numel(lat) && size(T, 2) == numel(lon)
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
            p = gsw_p_from_z(repmat(-z(:)', [numel(lat), 1]), lat);
        case 'p'
            p = repmat(z(:)', numel(lat), 1);
    end

    switch stype
        case 'SA'
            SA = S;
        case 'SP'
            SA = nan(size(S), "like", S);

            for k = 1:numel(z)
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

            for k = 1:numel(z)
                CT(:, :, k) = gsw_CT_from_t( ...
                    squeeze(SA(:, :, k)), squeeze(T(:, :, k)), p(:, k));
            end

    end

end
