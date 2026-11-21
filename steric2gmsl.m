%% STERIC2GMSL - Returns a native global-mean steric sea level series.
%
% Syntax
%   [gmsl, dates] = STERIC2GMSL(product, component)
%   [gmsl, dates] = STERIC2GMSL(product, component, depth)
%   [gmsl, dates] = STERIC2GMSL(__, "Name", value)
%
% Input arguments
%   product - Prepared product identifier accepted by StericProduct.
%       The native MAT file must already exist at STERICFILEPATH(product).
%   component - Component or array of components that make up total steric sea level.
%       Explicit component and depth inputs have equal element counts.
%       Pairs follow MATLAB linear indexing order; duplicates are retained.
%       total, thermosteric, halosteric, non-thermosteric (total minus
%       thermosteric), or non-halosteric (total minus halosteric).
%   depth (optional) - native, full-depth, shallow, or deep.
%       - 'native' resolves to the prepared product's depth range: full-depth for
%       EN4c14, CMEMS, ORAS5 and IAPv4; shallow for SIO, NOAA, JAMSTEC, ISAS and IPRC.
%       - 'full-depth' means integration from the sea surface to the seafloor.
%       Default depth: native for each component.
%   Unit (name-value) - Unit of the returned time series, m or mm.
%       Default unit: m.
%   TimeFormat (name-value) - datetime or datenum.
%       Default format: datetime.
%   SaveData (name-value) - Whether to append all calculated means to the native MAT file.
%       Default flag: true.
%   ForceNew (name-value) - Whether to recalculate all applicable means.
%       Default flag: false.
%   BeQuiet (name-value) - Whether to suppress messages.
%       Levels 0 and 0.5 report calculation and saving. Level 1 is silent.
%       Default level: 0.5.
%
% Input file content
%   lon - Regular longitude coordinate vector spanning one complete cycle.
%       Size: vector with nLon elements.
%       Unit: degrees east.
%   lat - Strictly monotone latitude coordinate vector in [-90, 90].
%       Spacing may be uniform or nonuniform.
%       Size: vector with nLat elements.
%       Unit: degrees north.
%   dates - Datetime or datenum vector, strictly increasing.
%       Size: vector with nDate elements.
%       Unit: days for datenum values.
%   stericSl, thermostericSl, halostericSl - Component rasters.
%       All grids must share an orientation. Square spatial grids are
%       interpreted as (lat, lon, time), as in the native aggregate format.
%       Missing or empty components are unavailable, not zero observations.
%       Size: (nLat, nLon, nDate) or (nLon, nLat, nDate).
%       Unit: m.
%   shallow<Component>Sl, deep<Component>Sl - Depth-layer rasters.
%       Component is Steric, Thermosteric, or Halosteric. Each split requires
%       full-depth, shallow and deep grids.
%       Size: same as the corresponding full-depth grid.
%       Unit: m.
%
% Output arguments
%   gmsl - Selected uncorrected global means on common observed support.
%       Size: (N, nSeries), one column per requested series.
%       Unit: m or mm, as selected by Unit.
%   dates - Native timestamps in the selected TimeFormat.
%       Size: (N, 1).
%       Unit: days for datenum output.
%
% Output file content
%   <component>Gmsl - Individual cached mean vectors.
%       Components are steric, thermosteric and halosteric, with shallow/deep
%       prefixes where supported. Difference components are calculated from
%       the finished means on request and are never saved. Original variables remain.
%       Size: (N, 1).
%       Unit: m.
%   gmslCreationDate - Datetime when newly calculated means were saved.
%       Size: (1, 1).
%
% Notes
%   The fixed mask uses the first available unprefixed native field in priority
%   order: halosteric, total, then thermosteric. A cell must have finite,
%   nonzero reference values at 80% or more of the native timestamps.
%   Deep-ocean NaNs at shallow-water cells are filled using STERICCOVERAGE
%   and FILLABSENTDEEPSTERIC before selecting the common finite cells across
%   supplied components and layers. Shallow-water zeros retain their area
%   weight. A timestamp with no positive-weight valid cells returns NaN.
%   Use ForceNew to refresh native means saved under earlier conventions.
%
% See also
%   steric2lonlatt
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function [gmsl, dates] = steric2gmsl(product, component, depth, options)

    arguments (Input)
        product (1, 1) StericProduct
        component string {mustBeNonempty, mustBeMember(component, ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"])}
        depth string {mustBeMember(depth, ["native", "full-depth", "shallow", "deep"])} = repmat("native", size(component))
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["m", "mm"])} = "m"
        options.TimeFormat (1, 1) DateFormat = "datetime"
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) double {mustBeMember(options.BeQuiet, [0, 0.5, 1])} = 0.5
    end

    if numel(component) ~= numel(depth)
        error('ULMO:steric2gmsl:PairCountMismatch', 'Component and depth must have the same number of elements.');
    end

    inPath = stericfilepath(product);

    if ~isfile(inPath)
        error('ULMO:steric2gmsl:FileNotFound', 'Native product file not found: %s', inPath);
    end

    nativDepth = product.nativeDepth();

    depth(depth == "native") = nativDepth;
    invalidDepth = find(ismember(depth, ["full-depth", "deep"]) & ~product.supportsFullDepth(), 1);

    if ~isempty(invalidDepth)
        error('ULMO:steric2gmsl:UnsupportedField', ...
            'Product %s does not supply Depth=%s.', char(product), depth(invalidDepth));
    end

    metadata = whos('-file', inPath);
    [schema, precision] = fieldschema(metadata, nativDepth);
    rasterNames = schema.gridNames;
    available = schema.available;
    meanNames = schema.meanNames;
    supported = schema.supported;
    components = schema.components;
    depths = schema.depths;
    [~, iComp] = ismember(component(:), components);
    [~, iDepth] = ismember(depth(:), depths);
    series = sub2ind(size(supported), iComp, iDepth);
    unsupported = find(~supported(series), 1);

    if ~isempty(unsupported)
        error('ULMO:steric2gmsl:UnsupportedField', ...
            'Product %s does not supply Component=%s, Depth=%s.', ...
            char(product), component(unsupported), depth(unsupported));
    end

    if ~all(ismember({'lon', 'lat', 'dates'}, {metadata.name}))
        error('ULMO:steric2gmsl:InvalidSchema', 'Native lon, lat and dates are required.');
    end

    coord = load(inPath, 'lon', 'lat', 'dates');
    [dates, wghts, keepLon, transGrid] = coordinates(coord, metadata, rasterNames(available));
    cachedComponents = supported;
    cachedComponents(4:5, :) = false;
    names = cellstr(meanNames(cachedComponents));
    cacheHit = ~options.ForceNew && validcache(metadata, names, numel(dates));

    if cacheHit
        means = load(inPath, names{:});
        cacheHit = all(structfun(@(x) isreal(x) && ~any(isinf(x)), means));
    end

    if ~cacheHit

        if options.BeQuiet < 1
            fprintf('[ULMO>%s] Calculating GMSLs for %s.\n', mfilename, char(product));
        end

        means = steric2gmslCore(inPath, schema, ...
            numel(dates), wghts, keepLon, transGrid, precision, find(depths == nativDepth));

        if options.SaveData
            means.gmslCreationDate = datetime('now');
            save(inPath, '-struct', 'means', '-append');

            if options.BeQuiet < 1
                fprintf('[ULMO>%s] Saved GMSLs for %s.\n', mfilename, filehref(inPath, char(product)));
            end

        end

    end

    % Derive each requested difference from the finished native means.
    gmsl = zeros(numel(dates), numel(series));

    for iSeries = 1:numel(series)

        if iComp(iSeries) > 3
            gmsl(:, iSeries) = double(means.(meanNames(1, iDepth(iSeries)))) - ...
                double(means.(meanNames(iComp(iSeries) - 2, iDepth(iSeries))));
        else
            gmsl(:, iSeries) = double(means.(meanNames(iComp(iSeries), iDepth(iSeries))));
        end

    end

    if options.Unit == "mm"
        gmsl = 1e3 * gmsl;
    end

    if options.TimeFormat == DateFormat.datenum
        dates = datenum(dates); %#ok<DATNM>
    end

end

function gmslObj = steric2gmslCore(inPath, schema, ...
        nDates, wghts, keepLon, transGrid, precision, nativDepthIdx)
    inFileObj = matfile(inPath);
    rasterNames = schema.gridNames;
    fieldAvail = schema.available;
    meanNames = schema.meanNames;
    supported = schema.supported;
    splitComponents = find(schema.split)';
    fieldIndices = find(fieldAvail(:))';

    if ~isempty(splitComponents)
        masks = stericcoverage(inFileObj, rasterNames, fieldAvail, nDates);
        oceanMask = masks.oceanMask;
        deepOceanMask = masks.deepOceanMask;

        if transGrid
            oceanMask = oceanMask';
            deepOceanMask = deepOceanMask';
        end

        oceanMask = oceanMask(:, keepLon);
        deepOceanMask = deepOceanMask(:, keepLon);
    end

    % Determine full-record validity from the existing reference field.
    priority = [3, 1, 2];
    refField = find(fieldAvail(priority, nativDepthIdx), 1);
    refField = rasterNames(priority(refField), nativDepthIdx);
    finiteCnt = zeros(size(wghts));
    nonzeroCnt = zeros(size(wghts));

    % Counts which grid points have enough valid points to be used
    for iDate = 1:nDates
        field = readtimeslice(inFileObj, refField, iDate, keepLon, transGrid);
        finiteCnt = finiteCnt + isfinite(field);
        nonzeroCnt = nonzeroCnt + (isfinite(field) & field ~= 0);
    end

    support = (finiteCnt >= 0.8 * nDates) & (nonzeroCnt >= 0.8 * nDates) & wghts > 0;
    gmsls = nan(3, 3, nDates);

    for iDate = 1:nDates
        fields = cell(3);
        valid = support;

        for iField = fieldIndices
            fields{iField} = readtimeslice(inFileObj, rasterNames(iField), iDate, keepLon, transGrid);
        end

        for iComp = splitComponents
            fields{iComp, 3} = fillabsentdeepsteric( ...
                fields{iComp, 1}, fields{iComp, 2}, fields{iComp, 3}, oceanMask, deepOceanMask);
        end

        for iField = fieldIndices
            valid = valid & isfinite(fields{iField});
        end

        totalWghts = sum(wghts(valid));

        if totalWghts == 0
            continue
        end

        w = wghts(valid) / totalWghts;

        for iComp = 1:3

            for iDepth = find(fieldAvail(iComp, :))
                field = fields{iComp, iDepth};
                gmsls(iComp, iDepth, iDate) = sum(w .* field(valid));
            end

            % If shallow and deep fields are available, make sure they sum to the full-depth value
            if all(fieldAvail(iComp, :))
                gmsls(iComp, 1, iDate) = checkDepthPartsSumToFullDepth( ...
                    gmsls, iComp, iDate, fields, w, valid, precision, rasterNames(iComp, 1));
            end

        end

    end

    gmslObj = struct;

    for iComp = 1:3

        for iDepth = find(supported(iComp, :))
            gmslObj.(meanNames(iComp, iDepth)) = reshape(gmsls(iComp, iDepth, :), [], 1);
        end

    end

end

%% Subfunctions
function [schema, precision] = fieldschema(metadata, nativeDepth)
    schema = stericfieldschema(metadata, nativeDepth);
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
        error('ULMO:steric2gmsl:InvalidSchema', 'A depth split requires explicit full-depth, shallow and deep grids.');
    end

    if schema.incompatibleSplits
        error('ULMO:steric2gmsl:InvalidSchema', 'Supplied total and halosteric depth splits must agree.');
    end

end

%% COORDINATES - Validate time and spatial coordinates of the field grids
function [dates, wghts, keepLon, trans] = coordinates(coord, metadata, names)
    lon = double(coord.lon(:));
    lat = double(coord.lat(:));

    mustBeCoordinateVector = @(x) isreal(x) && all(isfinite(x(:))) && isvector(x) && numel(x) >= 2;

    if ~mustBeCoordinateVector(coord.lon) || ~mustBeCoordinateVector(coord.lat) ...
            || any(abs(lat) > 90)
        error('ULMO:steric2gmsl:InvalidGrid', 'Expected finite longitude/latitude coordinate vectors.');
    end

    tol = 3e-5; % Coordinate tolerance in degrees.
    latSteps = diff(lat);

    if ~(all(latSteps > 0) || all(latSteps < 0))
        error('ULMO:steric2gmsl:InvalidGrid', 'Latitude must be strictly increasing or decreasing.');
    end

    % Both endpoints of an offset can carry single-precision rounding.
    lonTol = tol;

    if isa(coord.lon, 'single')
        lonTol = max(lonTol, double(eps(max(abs(coord.lon(:))))));
    end

    lonTol = min(lonTol, 90 / numel(lon));
    lonSpan = lon(end) - lon(1);

    % Drop the last longitude slice if it is the same as the first longitude (e.g. 0° and 360°)
    % to avoid double counting its contribution
    keepLon = 1:(numel(lon) - (abs(abs(lonSpan) - 360) <= lonTol));

    % Check the longitude vector is equally spaced
    lonStep = sign(lonSpan) * 360 / numel(keepLon);
    expdLons = (0:numel(lon) - 1)' * lonStep;

    if lonSpan == 0 || any(abs((lon - lon(1)) - expdLons) > lonTol)
        error('ULMO:steric2gmsl:InvalidGrid', 'Longitude must cover one complete regular periodic grid.');
    end

    dates = coord.dates;

    if isnumeric(dates) && isreal(dates) && all(isfinite(dates(:)))
        dates = datetime(dates, 'ConvertFrom', 'datenum');
    end

    if ~isdatetime(dates) || ~isvector(dates) || isempty(dates) || ...
            any(isnat(dates(:))) || any(diff(dates(:)) <= seconds(0))
        error('ULMO:steric2gmsl:InvalidDates', ...
        'Native dates must be finite and strictly increasing.');
    end

    dates = dates(:);

    % Checks the size of the fields, and whether it needs to be transposed
    % The expected shape is (nLat, nLon, nTime)
    trans = false;

    for iField = 1:numel(names)
        field = metadata(strcmp({metadata.name}, names(iField)));
        shape = field.size;
        shape(end + 1:3) = 1;

        if numel(shape) ~= 3 || shape(3) ~= numel(dates) || ...
                ~(isequal(shape(1:2), [numel(lat), numel(lon)]) || ...
                isequal(shape(1:2), [numel(lon), numel(lat)]))
            error('ULMO:steric2gmsl:InvalidSchema', ...
                '%s dimensions do not match the coordinates/dates.', names(iField));
        end

        thisTranspose = ~isequal(shape(1:2), [numel(lat), numel(lon)]);

        if iField > 1 && thisTranspose ~= trans
            error('ULMO:steric2gmsl:InvalidSchema', ...
            'All source grids must use the same orientation.');
        end

        trans = thisTranspose;
    end

    % Construct the grid weights for GMSL
    if all(abs(latSteps - latSteps(1)) <= tol)
        % Regularly spaced grid
        latWghts = cosd(lat);
    else
        % Irregularly spaced grid: Compute the weighting be differencing cell boundaries
        bdrys = [lat(1) - latSteps(1) / 2; ...
                     (lat(1:end - 1) + lat(2:end)) / 2; ...
                     lat(end) + latSteps(end) / 2];
        bdrys = max(-90, min(90, bdrys));
        latWghts = abs(diff(sind(bdrys)));
    end

    wghts = repmat(latWghts, 1, numel(keepLon));
    % Don't normalise yet since we haven't check which grids are actually valid yet
end

%% VALIDCACHE - Validates that the stored cache exists and and has the correct properties
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

%% READTIMESLICE - Reads a time slice of a field saved in a .mat file
function field = readtimeslice(matObj, field, iDate, keepLon, trans)
    field = double(matObj.(char(field))(:, :, iDate));

    if trans
        field = field';
    end

    field = field(:, keepLon);
end

%% CHECKDEPTHPARTSSUMTOFULLDEPTH - Checks that the full-depth value agrees with the sum of shallow and deep parts
function combined = checkDepthPartsSumToFullDepth(values, iComponent, iDate, fields, w, valid, precision, fieldName)
    direct = values(iComponent, 1, iDate);
    combined = sum(values(iComponent, 2:3, iDate));
    scale = 1;

    for iDepth = 1:3
        field = fields{iComponent, iDepth};
        scale = max(scale, sum(w .* abs(field(valid))));
    end

    if abs(direct - combined) > 32 * precision * scale
        error('ULMO:steric2gmsl:DepthMismatch', ...
            'Full-depth and shallow + deep means disagree for field %s at time step number %d.', ...
            fieldName, iDate);
    end

end
