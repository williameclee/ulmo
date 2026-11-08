%% STERIC2GMSL - Returns a native global-mean steric sea level series.
%
% Syntax
%   [gmsl, dates] = STERIC2GMSL(product, type)
%   [gmsl, dates] = STERIC2GMSL(product, type, depth)
%   [gmsl, dates] = STERIC2GMSL(__, "Name", value)
%
% Input arguments
%   product - Prepared product identifier accepted by StericProduct.
%       The native MAT file must already exist at STERICFILEPATH(product).
%   type - Selected component.
%       total, thermosteric, halosteric, non-thermosteric (total minus
%       thermosteric), or non-halosteric (total minus halosteric).
%       The returned mean uses uncorrected native fields.
%   depth (optional) - full, shallow, or deep, as supplied by the product.
%       Unsupported component/depth combinations raise an error.
%       Default depth: full.
%   Unit (name-value) - m or mm, for the returned series only.
%       Default unit: m.
%   TimeFormat (name-value) - datetime or datenum, for returned dates only.
%       Default format: datetime.
%   SaveData (name-value) - Whether to append all calculated means to the
%       native MAT file. Existing complete caches can still be read if false.
%       Default flag: true.
%   ForceNew (name-value) - Whether to recalculate all applicable means.
%       Set true after changing source grids or the averaging convention.
%       Combine with SaveData=false to recalculate in memory.
%       Default flag: false.
%   BeQuiet (name-value) - Whether to suppress messages.
%       Levels 0 and 0.5 report calculation and saving. Level 1 is silent.
%       Default level: 0.5.
%
% Input file content
%   lon - Regular longitude coordinate vector spanning one complete cycle.
%       Coordinates may ascend or descend and may repeat the cyclic endpoint.
%       Offsets from the first meridian must match a regular 360-degree grid.
%       The tolerance is 3e-5 degrees or one single-precision spacing at the
%       largest absolute longitude, whichever is larger for single inputs.
%       It is capped at 90 degrees divided by the stored longitude count.
%       Size: vector with nLon elements.
%       Unit: degrees east.
%   lat - Strictly monotone latitude coordinate vector in [-90, 90].
%       Spacing may be uniform or nonuniform.
%       Size: vector with nLat elements.
%       Unit: degrees north.
%   dates - Native datetime or datenum vector, strictly increasing.
%       Size: vector with N elements.
%       Unit: days for datenum values.
%   stericSl, thermostericSl, halostericSl - Available native component grids.
%       All grids must share an orientation. Square spatial grids are
%       interpreted as (lat, lon, time), as in the native aggregate format.
%       Missing or empty components are unavailable, not zero observations.
%       Size: (nLat, nLon, N) or (nLon, nLat, N).
%       Unit: m.
%   shallow<Component>Sl, deep<Component>Sl - Supplied depth-layer grids.
%       Component is Steric, Thermosteric, or Halosteric. Each split requires
%       full, shallow and deep grids. Total and halosteric splits must agree
%       when both components are present; missing layers are errors.
%       Size: same as the corresponding full grid.
%       Unit: m.
%
% Output arguments
%   gmsl - Selected uncorrected global mean on common observed support.
%       Size: (N, 1).
%       Unit: m or mm, as selected by Unit.
%   dates - Native timestamps in the selected TimeFormat.
%       Size: (N, 1).
%       Unit: days for datenum output.
%
% Output file content
%   <component>Gmsl - Individual cached mean vectors.
%       Components are steric, thermosteric, halosteric, nonThermosteric and
%       nonHalosteric. Supported shallow/deep names prefix and capitalise
%       the component, e.g. shallowHalostericGmsl. Original variables remain.
%       Size: (N, 1).
%       Unit: m.
%   gmslCreationDate - Datetime when newly calculated means were saved.
%       Size: (1, 1).
%       Reusing a cache or setting SaveData=false leaves this value unchanged.
%
% Notes
%   The fixed mask uses the first available full-depth field in priority
%   order: halosteric, total, then thermosteric. A cell must have finite,
%   nonzero reference values at 80% or more of the native timestamps.
%   Each timestamp uses cells that also have finite values in every supplied
%   component and layer. Zero layer values contribute on this common support.
%   A timestamp with no positive-weight support returns NaN.
%
% See also
%   steric2lonlatt, stericfilepath, stericfieldschema
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/08, OpenAI Codex

function [gmsl, dates] = steric2gmsl(product, type, depth, options)

    arguments (Input)
        product (1, 1) StericProduct
        type {mustBeTextScalar, mustBeMember(type, ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"])}
        depth {mustBeTextScalar, mustBeMember(depth, ["full", "shallow", "deep"])} = "full"
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, ["m", "mm"])} = "m"
        options.TimeFormat {mustBeTextScalar, mustBeMember(options.TimeFormat, ["datetime", "datenum"])} = "datetime"
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) double {mustBeMember(options.BeQuiet, [0, 0.5, 1])} = 0.5
    end

    inPath = stericfilepath(product);

    if ~isfile(inPath)
        error('ULMO:steric2gmsl:FileNotFound', 'Native product file not found: %s', inPath);
    end

    metadata = whos('-file', inPath);
    [schema, precision] = fieldschema(metadata);
    gridNames = schema.gridNames;
    available = schema.available;
    cacheNames = schema.meanNames;
    supported = schema.supported;
    types = schema.types;
    depths = schema.depths;
    iType = find(strcmpi(types, type));
    iDepth = find(strcmpi(depths, depth));

    if ~supported(iType, iDepth)
        error('ULMO:steric2gmsl:UnsupportedField', ...
            'Product %s does not supply Type=%s, Depth=%s.', char(product), type, depth);
    end

    if ~all(ismember({'lon', 'lat', 'dates'}, {metadata.name}))
        error('ULMO:steric2gmsl:InvalidSchema', 'Native lon, lat and dates are required.');
    end

    coord = load(inPath, 'lon', 'lat', 'dates');
    [dates, wghts, keepLon, transposeGrid] = coordinates(coord, metadata, gridNames(available));
    names = cellstr(cacheNames(supported));
    cacheHit = ~options.ForceNew && validcache(metadata, names, numel(dates));

    if cacheHit
        means = load(inPath, names{:});
        cacheHit = all(structfun(@(x) isreal(x) && ~any(isinf(x)), means));
    end

    if ~cacheHit

        if options.BeQuiet < 1
            fprintf('[ULMO>%s] Calculating native means for %s.\n', mfilename, char(product));
        end

        means = steric2gmslCore(inPath, gridNames, available, cacheNames, supported, ...
            numel(dates), wghts, keepLon, transposeGrid, precision);

        if options.SaveData
            means.gmslCreationDate = datetime('now');
            save(inPath, '-struct', 'means', '-append');

            if options.BeQuiet < 1
                fprintf('[ULMO>%s] Saved native means for %s.\n', mfilename, filehref(inPath, char(product)));
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

function means = steric2gmslCore(inputPath, gridNames, available, cacheNames, supported, ...
        nDates, wghts, keepLon, transposeGrid, precision)
    source = matfile(inputPath);
    % First pass uses only the reference field to determine full-record support.
    reference = find(available([3, 1, 2], 1), 1);
    priority = [3, 1, 2];
    reference = gridNames(priority(reference), 1);
    finiteCount = zeros(size(wghts));
    nonzeroCount = zeros(size(wghts));

    for iDate = 1:nDates
        field = readslice(source, reference, iDate, keepLon, transposeGrid);
        finiteCount = finiteCount + isfinite(field);
        nonzeroCount = nonzeroCount + (isfinite(field) & field ~= 0);
    end

    support = finiteCount >= 0.8 * nDates & nonzeroCount >= 0.8 * nDates & wghts > 0;
    values = nan(5, 3, nDates);

    for iDate = 1:nDates
        fields = cell(3);
        valid = support;

        for iField = find(available(:))'
            fields{iField} = readslice(source, gridNames(iField), iDate, keepLon, transposeGrid);
            valid = valid & isfinite(fields{iField});
        end

        totalWghts = sum(wghts(valid));

        if totalWghts == 0
            continue
        end

        w = wghts(valid) / totalWghts;

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

%% Subfunctions
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

function [dates, wghts, keepLon, transposeGrid] = coordinates(coord, metadata, names)
    lon = double(coord.lon(:));
    lat = double(coord.lat(:));

    if ~isnumeric(coord.lon) || ~isnumeric(coord.lat) || ...
            ~isvector(coord.lon) || ~isvector(coord.lat) || ...
            ~isreal(lon) || ~isreal(lat) || any(~isfinite([lon; lat])) || ...
            numel(lon) < 2 || numel(lat) < 2 || any(abs(lat) > 90)
        error('ULMO:steric2gmsl:InvalidGrid', 'Expected finite longitude/latitude coordinate vectors.');
    end

    tol = 3e-5; % Coordinate tolerance in degrees.
    latSteps = diff(lat);
    if ~(all(latSteps > 0) || all(latSteps < 0))
        error('ULMO:steric2gmsl:InvalidGrid', 'Latitude must be strictly monotone.');
    end

    % Both endpoints of an offset can carry single-precision rounding.
    lonTolerance = tol;
    if isa(coord.lon, 'single')
        lonTolerance = max(lonTolerance, double(eps(max(abs(coord.lon(:))))));
    end
    % Keep adjacent meridians distinguishable even on very fine grids.
    lonTolerance = min(lonTolerance, 90 / numel(lon));
    lonSpan = lon(end) - lon(1);
    keepLon = 1:numel(lon);
    if abs(abs(lonSpan) - 360) <= lonTolerance
        keepLon(end) = [];
    end

    % A complete cycle fixes the spacing independently of the first interval.
    lonStep = sign(lonSpan) * 360 / numel(keepLon);
    expectedOffsets = (0:numel(lon) - 1)' * lonStep;
    if lonSpan == 0 || any(abs((lon - lon(1)) - expectedOffsets) > lonTolerance)
        error('ULMO:steric2gmsl:InvalidGrid', 'Longitude must cover one complete regular periodic grid.');
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

    if all(abs(latSteps - latSteps(1)) <= tol)
        latWghts = cosd(lat);
    else
        bdrys = [lat(1) - latSteps(1) / 2; ...
                     (lat(1:end - 1) + lat(2:end)) / 2; ...
                     lat(end) + latSteps(end) / 2];
        bdrys = max(-90, min(90, bdrys));
        latWghts = abs(diff(sind(bdrys)));
    end

    wghts = repmat(latWghts, 1, numel(keepLon));
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

function field = readslice(source, name, iDate, keepLon, transposeGrid)
    field = double(source.(char(name))(:, :, iDate));

    if transposeGrid
        field = field';
    end

    field = field(:, keepLon);
end
