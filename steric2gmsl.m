%% STERIC2GMSL - Read or cache one native global-mean steric time series.
%
% Syntax
%   [gmsl, dates] = STERIC2GMSL(product, Type="halosteric")
%   [gmsl, dates] = STERIC2GMSL(product, Type="thermosteric", Depth="shallow")
%   [gmsl, dates] = STERIC2GMSL(__, SaveData=false, ForceNew=true)
%
% Input arguments
%   product - Prepared StericProduct identifier; required.
%   Type (name-value) - Required field: total, thermosteric, halosteric,
%       non-thermosteric (total minus thermosteric), or non-halosteric
%       (total minus halosteric). No correction is applied.
%   Depth (name-value) - full, shallow, or deep. Default depth: full.
%   Unit (name-value) - m or mm. Default unit: m.
%   TimeFormat (name-value) - datetime or datenum. Default format: datetime.
%   SaveData (name-value) - Whether to append all calculated means to the
%       native MAT file. Default flag: true.
%   ForceNew (name-value) - Whether to recalculate all applicable means.
%       Default flag: false. With SaveData=false, no file is changed.
%   BeQuiet (name-value) - Whether to suppress messages. Default flag: 0.5,
%       which prints calculation/save messages; true suppresses all messages.
%
% Input file content
%   lon - Regular longitude coordinate vector spanning one complete cycle.
%   lat - Strictly monotone latitude vector; spacing may be nonuniform.
%       Both coordinate vectors are in degrees.
%   dates - Native datetime or datenum vector, strictly increasing.
%   stericSl, thermostericSl, halostericSl - Available native grids in metres,
%       (lat, lon, time) or (lon, lat, time). Square spatial grids are
%       interpreted as (lat, lon, time), as in the native aggregate format.
%       Supported depth splits supply both shallow<Component>Sl and
%       deep<Component>Sl alongside the full grid; missing layers are errors.
%       Missing or empty components are unavailable, not zero observations.
%
% Output arguments
%   gmsl - Selected uncorrected global mean. Size: (N, 1). Unit: Unit.
%   dates - Native timestamps. Size: (N, 1). Format: TimeFormat.
%
% Output file content
%   <component>Gmsl - Individual numeric (N, 1) vectors in metres. Components
%       are steric, thermosteric, halosteric, nonThermosteric and
%       nonHalosteric. Supported shallow/deep names prefix and capitalise
%       the component, e.g. shallowHalostericGmsl. Original variables remain.
%   gmslCreationDate - Scalar datetime recording when the means were processed.
%
% Notes
%   All applicable means are calculated together, then only one is returned.
%   The fixed support reference is full halosteric, otherwise full total,
%   otherwise full thermosteric. Cells need >=80% finite and >=80% finite
%   nonzero reference samples. Each date uses common finite support across
%   all supplied components/layers. Uniform latitude grids use cosine-latitude
%   weights. Nonuniform grids use spherical band areas between arithmetic
%   latitude midpoints: abs(diff(sind(boundaries))). Outer boundaries extend
%   half the adjacent spacing, clipped to [-90, 90]. Spacing differences
%   within 2e-5 degrees are treated as coordinate rounding. A duplicated
%   longitude endpoint is excluded. No supported area gives NaN. Valid zero
%   layers remain in the denominator. Residuals are differences on this
%   common support; total need not equal thermosteric plus halosteric.
%   Supplied layers are averaged directly; full means are their sum and
%   checked against the full grid with tolerance 32*eps(source precision)
%   times max(1, the weighted mean absolute full/shallow/deep values).
%   Complete caches bypass grid reads. Source/convention changes require
%   explicit ForceNew=true and manual rebuilding of affected downstream
%   caches. This function does not interpolate, remove drift or rebuild grids.
%
% See also
%   steric2lonlatt, stericfilepath
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

function [gmsl, dates] = steric2gmsl(product, options)
    arguments (Input)
        product (1, 1) StericProduct
        options.Type (1, 1) string {mustBeMember(options.Type, ...
            ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"])}
        options.Depth (1, 1) string {mustBeMember(options.Depth, ["full", "shallow", "deep"])} = "full"
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, ["m", "mm"])} = "m"
        options.TimeFormat {mustBeTextScalar, mustBeMember(options.TimeFormat, ["datetime", "datenum"])} = "datetime"
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) double {mustBeMember(options.BeQuiet, [0, 0.5, 1])} = 0.5
    end

    if ~isfield(options, 'Type')
        error('ULMO:steric2gmsl:TypeRequired', 'Specify Type explicitly.');
    end

    inputPath = stericfilepath(product);
    if ~isfile(inputPath)
        error('ULMO:steric2gmsl:FileNotFound', 'Native product file not found: %s', inputPath);
    end
    metadata = whos('-file', inputPath);
    [gridNames, available, precision] = fieldschema(metadata);
    [cacheNames, supported] = cachenames(available);
    types = ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"];
    depths = ["full", "shallow", "deep"];
    iType = find(types == options.Type);
    iDepth = find(depths == options.Depth);
    if ~supported(iType, iDepth)
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
    gmsl = double(means.(cacheNames(iType, iDepth)));
    if strcmp(options.Unit, 'mm')
        gmsl = 1e3 * gmsl;
    end
    if strcmp(options.TimeFormat, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end
end

function [names, available, precision] = fieldschema(metadata)
    names = ["stericSl", "shallowStericSl", "deepStericSl"; ...
             "thermostericSl", "shallowThermostericSl", "deepThermostericSl"; ...
             "halostericSl", "shallowHalostericSl", "deepHalostericSl"];
    available = false(3);
    precision = eps('double');
    for iField = 1:numel(names)
        entry = metadata(strcmp({metadata.name}, names(iField)));
        if isempty(entry) || prod(entry.size) == 0
            continue
        end
        if ~ismember(entry.class, {'double', 'single'}) || entry.complex || entry.sparse
            error('ULMO:steric2gmsl:InvalidSchema', '%s must be a real floating-point grid.', names(iField));
        end
        available(iField) = true;
        precision = max(precision, double(eps(entry.class)));
    end
    split = any(available(:, 2:3), 2);
    if any(split & ~all(available, 2))
        error('ULMO:steric2gmsl:InvalidSchema', 'A depth split requires explicit full, shallow and deep grids.');
    end
    if available(1, 1) && available(3, 1) && split(1) ~= split(3)
        error('ULMO:steric2gmsl:InvalidSchema', 'Supplied total and halosteric depth splits must agree.');
    end
end

function [names, supported] = cachenames(available)
    names = ["stericGmsl", "shallowStericGmsl", "deepStericGmsl"; ...
             "thermostericGmsl", "shallowThermostericGmsl", "deepThermostericGmsl"; ...
             "halostericGmsl", "shallowHalostericGmsl", "deepHalostericGmsl"; ...
             "nonThermostericGmsl", "shallowNonThermostericGmsl", "deepNonThermostericGmsl"; ...
             "nonHalostericGmsl", "shallowNonHalostericGmsl", "deepNonHalostericGmsl"];
    supported = [available; available(1, :) & available(2, :); available(1, :) & available(3, :)];
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
    tolerance = 2e-5; % Allows single-precision coordinate rounding, in degrees.
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
