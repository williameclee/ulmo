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
%   component - Component that make up the total steric sea level.
%       total, thermosteric, halosteric, non-thermosteric (total minus
%       thermosteric), or non-halosteric (total minus halosteric).
%   depth (optional) - native, full-depth, shallow, or deep.
%       - 'native' resolves to the prepared product's depth range: full-depth for
%       EN4c14, CMEMS, ORAS5 and IAPv4; shallow for SIO, NOAA, JAMSTEC, ISAS and IPRC.
%       - 'full-depth' means integration from the sea surface to the seafloor.
%       Default depth: native.
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
%   gmsl - Selected uncorrected global mean on common observed support.
%       Size: (N, 1).
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
%   Each timestamp uses cells that also have finite values in every supplied
%   component and layer. Zero layer values contribute on this common support.
%   A timestamp with no positive-weight support returns NaN.
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
        component (1, 1) string {mustBeStericComponent}
        depth (1, 1) string {mustBeStericDepth} = "native"
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["m", "mm"])} = "m"
        options.TimeFormat (1, 1) DateFormat = "datetime"
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) double {mustBeMember(options.BeQuiet, [0, 0.5, 1])} = 0.5
    end

    inPath = stericfilepath(product);

    if ~isfile(inPath)
        error('ULMO:steric2gmsl:FileNotFound', 'Native product file not found: %s', inPath);
    end

    nativDepth = product.nativeDepth();

    if depth == "native"
        depth = nativDepth;
    elseif ismember(depth, ["full-depth", "deep"]) && ~product.supportsFullDepth()
        error('ULMO:steric2gmsl:UnsupportedField', ...
            'Product %s does not supply steric sea level at Depth=%s.', char(product), depth);
    end

    metadata = whos('-file', inPath);
    [schema, precision] = fieldschema(metadata, nativDepth);
    rasterNames = schema.gridNames;
    available = schema.available;
    meanNames = schema.meanNames;
    supported = schema.supported;
    components = schema.components;
    depths = schema.depths;
    iComp = find(strcmpi(components, component));
    iDepth = find(strcmpi(depths, depth));

    if ~supported(iComp, iDepth)
        error('ULMO:steric2gmsl:UnsupportedField', ...
            'Product %s does not supply Component=%s, Depth=%s.', char(product), component, depth);
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

        means = steric2gmslCore(inPath, rasterNames, available, meanNames, supported, ...
            numel(dates), wghts, keepLon, transGrid, precision, find(depths == nativDepth));

        if options.SaveData
            means.gmslCreationDate = datetime('now');
            save(inPath, '-struct', 'means', '-append');

            if options.BeQuiet < 1
                fprintf('[ULMO>%s] Saved GMSLs for %s.\n', mfilename, filehref(inPath, char(product)));
            end

        end

    end

    % Derive difference components only after the native means are complete.
    if iComp > 3
        gmsl = double(means.(meanNames(1, iDepth))) - ...
            double(means.(meanNames(iComp - 2, iDepth)));
    else
        gmsl = double(means.(meanNames(iComp, iDepth)));
    end

    if options.Unit == "mm"
        gmsl = 1e3 * gmsl;
    end

    if options.TimeFormat == DateFormat.datenum
        dates = datenum(dates); %#ok<DATNM>
    end

end

function gmslObj = steric2gmslCore(inPath, rasterNames, fieldAvail, meanNames, supported, ...
        nDates, wghts, keepLon, transGrid, precision, nativDepthIdx)
    inFileObj = matfile(inPath);
    % First pass uses only the reference field to determine full-record support
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

        for iField = find(fieldAvail(:))'
            fields{iField} = readtimeslice(inFileObj, rasterNames(iField), iDate, keepLon, transGrid);
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
function [dates, wghts, keepLon, transposeGrid] = coordinates(coord, metadata, names)
    try
        [~, ~, iLon, iLat, orderedWeights] = preparegridcoordinates( ...
            coord.lon, coord.lat, RequireGlobalLongitude=true);
    catch cause
        if strcmp(cause.identifier, 'ULMO:preparegridcoordinates:InvalidGrid')
            throwAsCaller(addCause(MException('ULMO:steric2gmsl:InvalidGrid', '%s', cause.message), cause));
        end
        rethrow(cause)
    end
    lon = coord.lon; lat = coord.lat;
    % Keep the native reader's original traversal and summation order.
    [keepLon, orderLon] = sort(iLon);
    keepLon = keepLon.';
    wghts = zeros(numel(lat), numel(keepLon));
    wghts(iLat, :) = orderedWeights(orderLon, :).';
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
