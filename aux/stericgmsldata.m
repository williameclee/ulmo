%% STERICGMSLDATA - Read or calculate all native steric global means together.
% Internal implementation shared by the scalar reader and gridded loaders.
%
% Syntax
%   [means, dates, schema] = STERICGMSLDATA(product, SaveData=false)
%
% Input arguments
%   product - Prepared StericProduct identifier.
%   Type, Depth (name-value) - Optional requested field to validate before
%       reading coordinates, preserving STERIC2GMSL validation order.
%       Default selection: empty Type and full Depth (all supported fields).
%   ForceNew (name-value) - Whether to recalculate all native means.
%       Default flag: false.
%   SaveData (name-value) - Whether to append calculated means to the native file.
%       Default flag: true.
%   BeQuiet (name-value) - Whether to suppress messages.
%       Default level: 0.5.
%
% Output arguments
%   means - Structure of supported named GMSL vectors in metres, (nTime, 1).
%       Contains only mean variables, excluding the cache creation timestamp.
%   dates - Native datetime column, (nTime, 1).
%   schema - Field availability from STERICFIELDSCHEMA.
%
% Notes
%   Preserves STERIC2GMSL's native support, area weights, layer closure and
%   validation. A single calculation supplies all means even without writes.
%
% See also
%   steric2gmsl, stericfieldschema
%
% Created by
%   2026/11/05, En-Chi Lee

function [means, dates, schema] = stericgmsldata(product, options)
    arguments (Input)
        product (1, 1) StericProduct
        options.Type (1, 1) string = ""
        options.Depth (1, 1) string = "full"
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) double = 0.5
    end
    inputPath = stericfilepath(product);
    if ~isfile(inputPath)
        error('ULMO:steric2gmsl:FileNotFound', 'Native product file not found: %s', inputPath);
    end
    metadata = whos('-file', inputPath);
    [schema, precision] = fieldschema(metadata);
    gridNames = schema.gridNames;
    available = schema.available;
    cacheNames = schema.meanNames;
    supported = schema.supported;
    types = schema.types;
    depths = schema.depths;
    iType = find(types == options.Type);
    iDepth = find(depths == options.Depth);
    if strlength(options.Type) > 0 && (isempty(iType) || isempty(iDepth) || ~supported(iType, iDepth))
        error('ULMO:steric2gmsl:UnsupportedField', ...
            'Product %s does not supply Type=%s, Depth=%s.', char(product), options.Type, options.Depth);
    end
    if ~all(ismember({'lon', 'lat', 'dates'}, {metadata.name}))
        error('ULMO:steric2gmsl:InvalidSchema', 'Native lon, lat and dates are required.');
    end
    coord = load(inputPath, 'lon', 'lat', 'dates');
    [dates, weights, keepLon, transposeGrid] = coordinates(coord, metadata, gridNames(available));
    names = cellstr(cacheNames(supported));
    cacheHit = ~options.ForceNew && validcache(metadata, names, numel(dates));
    if cacheHit
        means = load(inputPath, names{:});
        cacheHit = all(structfun(@(x) isreal(x) && ~any(isinf(x)), means));
    end
    if ~cacheHit
        if options.BeQuiet < 1
            fprintf('[ULMO>steric2gmsl] Calculating native means for %s.\n', char(product));
        end
        means = calculatemeans(inputPath, gridNames, available, cacheNames, supported, ...
            numel(dates), weights, keepLon, transposeGrid, precision);
        if options.SaveData
            means.gmslCreationDate = datetime('now');
            save(inputPath, '-struct', 'means', '-append');
            if options.BeQuiet < 1
                fprintf('[ULMO>steric2gmsl] Saved native means for %s.\n', char(product));
            end
        end
    end
    if isfield(means, 'gmslCreationDate')
        means = rmfield(means, 'gmslCreationDate');
    end
end

function [schema, precision] = fieldschema(metadata)
    schema = stericfieldschema(metadata);
    % Numerical requirements and source precision belong to the estimator.
    precision = eps('double');
    for iField = find(schema.available(:))'
        name = schema.gridNames(iField);
        entry = metadata(strcmp({metadata.name}, name));
        if ~ismember(entry.class, {'double', 'single'}) || entry.complex || entry.sparse
            error('ULMO:steric2gmsl:InvalidSchema', '%s must be a real floating-point grid.', name);
        end
        precision = max(precision, double(eps(entry.class)));
    end
    if schema.incompleteSplit
        error('ULMO:steric2gmsl:InvalidSchema', 'A depth split requires explicit full, shallow and deep grids.');
    end
    if schema.incompatibleSplits
        error('ULMO:steric2gmsl:InvalidSchema', 'Supplied total and halosteric depth splits must agree.');
    end
end

function [dates, weights, keepLon, transposeGrid] = coordinates(coord, metadata, names)
    lon = double(coord.lon(:));
    lat = double(coord.lat(:));
    if ~isnumeric(coord.lon) || ~isnumeric(coord.lat) || ...
            ~isvector(coord.lon) || ~isvector(coord.lat) || ...
            ~isreal(lon) || ~isreal(lat) || any(~isfinite([lon; lat])) || ...
            numel(lon) < 2 || numel(lat) < 2 || any(abs(lat) > 90)
        error('ULMO:steric2gmsl:InvalidGrid', 'Expected finite longitude/latitude coordinate vectors.');
    end
    steps = [lon(2) - lon(1), lat(2) - lat(1)];
    tolerance = 3e-5; % Allows single-precision coordinate rounding, in degrees.
    latitudeSteps = diff(lat);
    if steps(1) == 0 || any(abs(diff(lon) - steps(1)) > tolerance) || ...
            ~(all(latitudeSteps > 0) || all(latitudeSteps < 0))
        error('ULMO:steric2gmsl:InvalidGrid', ...
            'Longitude must be regularly spaced and latitude must be strictly monotone.');
    end
    keepLon = 1:numel(lon);
    if abs(abs(lon(end) - lon(1)) - 360) <= tolerance
        keepLon(end) = [];
    end
    if abs(numel(keepLon) * abs(steps(1)) - 360) > tolerance * numel(lon)
        error('ULMO:steric2gmsl:InvalidGrid', 'Longitude must cover one complete periodic grid.');
    end
    dates = coord.dates;
    if isnumeric(dates) && isreal(dates) && all(isfinite(dates(:)))
        dates = datetime(dates, 'ConvertFrom', 'datenum');
    end
    if ~isdatetime(dates) || ~isvector(dates) || isempty(dates) || ...
            any(isnat(dates(:))) || any(diff(dates(:)) <= seconds(0))
        error('ULMO:steric2gmsl:InvalidDates', 'Native dates must be finite and strictly increasing.');
    end
    dates = dates(:);
    transposeGrid = false;
    for iField = 1:numel(names)
        entry = metadata(strcmp({metadata.name}, names(iField)));
        shape = entry.size;
        shape(end + 1:3) = 1;
        if numel(shape) ~= 3 || shape(3) ~= numel(dates) || ...
                ~(isequal(shape(1:2), [numel(lat), numel(lon)]) || ...
                  isequal(shape(1:2), [numel(lon), numel(lat)]))
            error('ULMO:steric2gmsl:InvalidSchema', '%s dimensions do not match the coordinates/dates.', names(iField));
        end
        thisTranspose = ~isequal(shape(1:2), [numel(lat), numel(lon)]);
        if iField > 1 && thisTranspose ~= transposeGrid
            error('ULMO:steric2gmsl:InvalidSchema', 'All source grids must use the same orientation.');
        end
        transposeGrid = thisTranspose;
    end
    if all(abs(latitudeSteps - steps(2)) <= tolerance)
        latitudeWeights = cosd(lat);
    else
        boundaries = [lat(1) - latitudeSteps(1) / 2; ...
                      (lat(1:end - 1) + lat(2:end)) / 2; ...
                      lat(end) + latitudeSteps(end) / 2];
        boundaries = max(-90, min(90, boundaries));
        latitudeWeights = abs(diff(sind(boundaries)));
    end
    weights = repmat(latitudeWeights, 1, numel(keepLon));
end

function valid = validcache(metadata, names, nDates)
    valid = true;
    for iName = 1:numel(names)
        entry = metadata(strcmp({metadata.name}, names{iName}));
        if isempty(entry) || ~ismember(entry.class, {'double', 'single'}) || ...
                entry.complex || entry.sparse || ~isequal(entry.size, [nDates, 1])
            valid = false;
            return
        end
    end
end

function means = calculatemeans(inputPath, gridNames, available, cacheNames, supported, ...
        nDates, weights, keepLon, transposeGrid, precision)
    source = matfile(inputPath);
    % First pass uses only the reference field to determine full-record support.
    reference = find(available([3, 1, 2], 1), 1);
    priority = [3, 1, 2];
    reference = gridNames(priority(reference), 1);
    finiteCount = zeros(size(weights));
    nonzeroCount = zeros(size(weights));
    for iDate = 1:nDates
        field = readslice(source, reference, iDate, keepLon, transposeGrid);
        finiteCount = finiteCount + isfinite(field);
        nonzeroCount = nonzeroCount + (isfinite(field) & field ~= 0);
    end
    support = finiteCount >= 0.8 * nDates & nonzeroCount >= 0.8 * nDates & weights > 0;
    values = nan(5, 3, nDates);
    for iDate = 1:nDates
        fields = cell(3);
        valid = support;
        for iField = find(available(:))'
            fields{iField} = readslice(source, gridNames(iField), iDate, keepLon, transposeGrid);
            valid = valid & isfinite(fields{iField});
        end
        denominator = sum(weights(valid));
        if denominator == 0
            continue
        end
        w = weights(valid) / denominator;
        for iType = 1:3
            for iDepth = find(available(iType, :))
                field = fields{iType, iDepth};
                values(iType, iDepth, iDate) = sum(w .* field(valid));
            end
            if available(iType, 2)
                direct = values(iType, 1, iDate);
                combined = sum(values(iType, 2:3, iDate));
                scale = 1;
                for iDepth = 1:3
                    field = fields{iType, iDepth};
                    scale = max(scale, sum(w .* abs(field(valid))));
                end
                if abs(direct - combined) > 32 * precision * scale
                    error('ULMO:steric2gmsl:DepthMismatch', ...
                        'Full and shallow+deep means disagree for %s at native sample %d.', gridNames(iType, 1), iDate);
                end
                values(iType, 1, iDate) = combined;
            end
        end
        values(4, :, iDate) = values(1, :, iDate) - values(2, :, iDate);
        values(5, :, iDate) = values(1, :, iDate) - values(3, :, iDate);
    end
    means = struct;
    for iType = 1:5
        for iDepth = find(supported(iType, :))
            means.(cacheNames(iType, iDepth)) = reshape(values(iType, iDepth, :), [], 1);
        end
    end
end

function field = readslice(source, name, iDate, keepLon, transposeGrid)
    field = double(source.(char(name))(:, :, iDate));
    if transposeGrid
        field = field';
    end
    field = field(:, keepLon);
end
