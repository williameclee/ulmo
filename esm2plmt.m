%% ESM2PLMT - Averages ESAESM2 fields over specified GRACE intervals.
% Reads one component of the six-hour, degree-180 ESAESM2 archive. Returns
% the same coefficient layout, units and date conventions as GRACE2PLMT_NEW.
%
% Syntax
%   [esmPlmt, esmStdPlmt, dates, metadata] = esm2plmt(folder, intervals)
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt(folder, intervals, ...
%       "Component", "H", "Unit", "GRAV", "Loutput", 60)
%
% Input arguments
%   folder - Folder containing extracted mtmshc_COMPONENT_YYYYMMDD_HH.gfc
%       files. Subdirectories are searched recursively.
%   intervals - Start and exclusive end of each averaging interval.
%       Use the GRCOF2 interval returned by PARSEGRACESOURCEFILE.
%       Size: (N, 2), datetime, unzoned UTC.
%       Boundaries must coincide with the six-hour sampling grid.
%   Component (name-value) - Individual component identifier.
%       Options: 'A', 'O', 'H', 'I', 'S', 'AOHIS'. Default component: 'O'.
%   Unit (name-value) - Output coefficient unit.
%       'GRAV': dimensionless Stokes; 'POT': radius times Stokes, following
%       CONVERTGRAVITY; 'SD': kg/m^2, equivalent to mm of fresh water.
%       Default unit: 'SD'.
%   Loutput (name-value) - Maximum output degree. Empty keeps source degree.
%       Larger values pad coefficient columns with NaN. Default degree: [].
%   LoveNumSource (name-value) - Load Love numbers used by CONVERTGRAVITY.
%       Default source: 'ISSM', in the CF frame.
%   OutputFormat (name-value) - 'timefirst' or 'traditional'.
%       Default format: 'timefirst'.
%   TimeFormat (name-value) - 'datenum' or 'datetime'.
%       Default format: 'datenum'.
%
% Output arguments
%   esmPlmt - Degree, order, cosine, sine coefficient columns.
%       Size: (N, M, 4), or (M, 4, N) for 'traditional'.
%       M = (Loutput + 1) * (Loutput + 2) / 2.
%   esmStdPlmt - Same layout, with NaN uncertainties and valid degree/order.
%       ESAESM2 supplies no formal uncertainties; NaN does not mean zero error.
%   dates - Interval midpoints, size (N, 1), in TimeFormat.
%   metadata - Source files, timestamps, intervals, reference parameters,
%       sampling, normalization and conversion settings for reproducibility.
%
% Notes
%   Uses an arithmetic mean of all six-hour samples in [start, end).
%   Missing or duplicate samples, incomplete harmonics and mixed reference
%   parameters fail. This is an interval match, not GRACE orbit weighting.
%   ESM fields are anomalies. No WGS84 subtraction, demeaning, low-degree
%   replacement, de-aliasing, GIA correction or filtering is applied.
%   Component C00 and degree 1 are preserved. Source GM/radius are retained;
%   matching the layout does not imply the same physical signal as GSM.
%   SD uses the same density (5517 kg/m^3) and conversion as GRACE2PLMT_NEW.
%   Selecting ISSM ensures pipeline consistency, not an exact inversion of
%   the Love-number spectrum originally used to construct ESAESM2.
%   For S, SD is equivalent surface mass, not a physical surface load.
%
% Data sources
%   ESAESM2. https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/
%
% See also
%   GRACE2PLMT_NEW, PARSEGRACESOURCEFILE, CONVERTGRAVITY
%
% Created by
%   2026/10/29

function [esmPlmt, esmStdPlmt, dates, metadata] = esm2plmt(folder, intervals, options)
    arguments (Input)
        folder (1, 1) string {mustBeFolder}
        intervals (:, 2) datetime
        options.Component (1, 1) string {mustBeMember(options.Component, ...
            ["A", "O", "H", "I", "S", "AOHIS"])} = "O"
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["GRAV", "POT", "SD"])} = "SD"
        options.Loutput double = []
        options.LoveNumSource (1, 1) string = "ISSM"
        options.OutputFormat (1, 1) string {mustBeMember(options.OutputFormat, ...
            ["timefirst", "traditional"])} = "timefirst"
        options.TimeFormat (1, 1) string {mustBeMember(options.TimeFormat, ...
            ["datenum", "datetime"])} = "datenum"
    end

    assert(~isempty(intervals) && ~any(isnat(intervals), 'all') && ...
        isempty(intervals.TimeZone) && all(intervals(:, 2) > intervals(:, 1)), ...
        'ULMO:esm2plmt:InvalidIntervals', 'Provide increasing, unzoned UTC intervals.');
    assert(all(mod(hour(intervals), 6) == 0 & minute(intervals) == 0 & ...
        second(intervals) == 0, 'all'), 'ULMO:esm2plmt:InvalidIntervals', ...
        'Interval boundaries must be on the six-hour sampling grid.');
    if ~isempty(options.Loutput)
        validateattributes(options.Loutput, {'double'}, ...
            {'scalar', 'integer', 'nonnegative', 'finite'});
    end

    component = options.Component;
    files = dir(fullfile(folder, '**', "mtmshc_" + component + "_*.gfc"));
    fileDates = NaT(numel(files), 1);
    for k = 1:numel(files)
        token = regexp(files(k).name, '^mtmshc_' + component + ...
            '_(\d{8}_\d{2})\.gfc$', 'tokens', 'once');
        assert(~isempty(token), 'ULMO:esm2plmt:InvalidName', ...
            'Invalid ESM filename: %s', files(k).name);
        fileDates(k) = datetime(token{1}, 'InputFormat', 'yyyyMMdd_HH');
    end

    nEpochs = size(intervals, 1);
    metadata = struct('dataset', "ESAESM2", 'component', component, ...
        'intervals', intervals, 'endExclusive', true, 'samplingHours', 6, ...
        'normalization', "fully_normalized", 'frame', "CF", ...
        'unit', options.Unit, 'loveNumSource', options.LoveNumSource, ...
        'earthDensity', 5517, 'outputFormat', options.OutputFormat, ...
        'timeFormat', options.TimeFormat, 'uncertaintyAvailable', false);
    metadata.files = cell(nEpochs, 1);
    metadata.sampleDates = cell(nEpochs, 1);
    metadata.sampleCounts = zeros(nEpochs, 1);
    reference = [];
    for iEpoch = 1:nEpochs
        expected = (intervals(iEpoch, 1):hours(6):intervals(iEpoch, 2) - hours(6))';
        ids = find(fileDates >= intervals(iEpoch, 1) & fileDates < intervals(iEpoch, 2));
        [actual, order] = sort(fileDates(ids));
        ids = ids(order);
        assert(isequal(actual, expected), 'ULMO:esm2plmt:IncompleteSampling', ...
            'Missing, duplicate or off-grid %s samples in interval %d.', component, iEpoch);
        paths = string(fullfile({files(ids).folder}, {files(ids).name}))';
        total = [];
        for k = 1:numel(ids)
            [plm, header] = readesmfile(paths(k));
            if isempty(reference)
                reference = header;
                nativeDegree = header.maxDegree;
                outDegree = options.Loutput;
                if isempty(outDegree)
                    outDegree = nativeDegree;
                end
                labels = harmoniclabels(outDegree);
                nRows = size(labels, 1);
                esmPlmt = nan(nEpochs, nRows, 4);
                esmPlmt(:, :, 1:2) = repmat(reshape(labels, 1, nRows, 2), nEpochs, 1, 1);
                esmStdPlmt = esmPlmt;
                metadata.sourceDegree = nativeDegree;
                metadata.outputDegree = outDegree;
                metadata.gravityParam = header.gravityParam;
                metadata.equatorRadius = header.equatorRadius;
            else
                assert(isequal(header, reference), 'ULMO:esm2plmt:MixedMetadata', ...
                    'ESM degree, normalization or reference parameters differ: %s', paths(k));
            end
            if isempty(total)
                total = zeros(size(plm, 1), 2);
            end
            total = total + plm(:, 3:4);
        end
        plm(:, 3:4) = total / numel(ids);
        if options.Unit ~= "GRAV"
            plm = convertgravity(plm, 'GRAV', options.Unit, ...
                reference.gravityParam, reference.equatorRadius, ...
                'LoveNumSource', options.LoveNumSource);
        end
        nCopy = min(size(plm, 1), size(esmPlmt, 2));
        esmPlmt(iEpoch, 1:nCopy, 3:4) = reshape(plm(1:nCopy, 3:4), 1, nCopy, 2);
        metadata.files{iEpoch} = paths;
        metadata.sampleDates{iEpoch} = actual;
        metadata.sampleCounts(iEpoch) = numel(ids);
    end
    dates = intervals(:, 1) + (intervals(:, 2) - intervals(:, 1)) / 2;
    if options.TimeFormat == "datenum"
        dates = datenum(dates);
    end
    if options.OutputFormat == "traditional"
        esmPlmt = permute(esmPlmt, [2, 3, 1]);
        esmStdPlmt = permute(esmStdPlmt, [2, 3, 1]);
    end
end

function [plm, header] = readesmfile(path)
    contents = fileread(path);
    endHeader = regexp(contents, '(?m)^end_of_head[^\r\n]*\r?\n', 'end', 'once');
    assert(~isempty(endHeader), 'ULMO:esm2plmt:InvalidHeader', ...
        'Missing end_of_head in %s', path);
    head = contents(1:endHeader);
    normalization = headertext(head, 'norm');
    assert(strcmp(normalization, 'fully_normalized') && ...
        strcmp(headertext(head, 'product_type'), 'gravity_field') && ...
        strcmp(headertext(head, 'error'), 'no'), 'ULMO:esm2plmt:InvalidHeader', ...
        'Expected fully normalized ESAESM2 gravity coefficients without formal errors.');
    [~, stem] = fileparts(path);
    assert(strcmp(headertext(head, 'modelname'), stem), ...
        'ULMO:esm2plmt:InvalidHeader', 'Model name does not match filename: %s', path);
    header.gravityParam = headernumber(head, 'earth_gravity_constant');
    header.equatorRadius = headernumber(head, 'radius');
    header.maxDegree = headernumber(head, 'max_degree');
    header.normalization = normalization;
    assert(header.gravityParam > 0 && header.equatorRadius > 0 && ...
        header.maxDegree >= 0 && header.maxDegree == fix(header.maxDegree), ...
        'ULMO:esm2plmt:InvalidHeader', 'Invalid reference constants or degree.');
    data = regexprep(contents(endHeader + 1:end), '[dD]', 'E');
    columns = textscan(data, '%s%f%f%f%f', 'ReturnOnError', false);
    plm = [columns{2:5}];
    plm = sortrows(plm, [1, 2]);
    assert(all(strcmp(columns{1}, 'gfc')) && all(isfinite(plm), 'all') && ...
        isequal(plm(:, 1:2), harmoniclabels(header.maxDegree)), ...
        'ULMO:esm2plmt:InvalidCoefficients', 'Incomplete or invalid harmonics in %s', path);
end

function value = headertext(head, key)
    match = regexp(head, ['(?m)^', key, '[ \t]+(\S+)[ \t]*\r?$'], 'tokens', 'once');
    assert(~isempty(match), 'ULMO:esm2plmt:InvalidHeader', 'Missing header key %s.', key);
    value = match{1};
end

function value = headernumber(head, key)
    value = str2double(regexprep(headertext(head, key), '[dD]', 'E'));
    assert(isfinite(value), 'ULMO:esm2plmt:InvalidHeader', 'Invalid header value %s.', key);
end

function labels = harmoniclabels(degree)
    labels = zeros((degree + 1) * (degree + 2) / 2, 2);
    for l = 0:degree
        rows = l * (l + 1) / 2 + (1:l + 1);
        labels(rows, :) = [repmat(l, l + 1, 1), (0:l)'];
    end
end
