%% ESM2PLMT - Loads cached six-hour or monthly ESAESM2 fields.
% Reads one component of the six-hour, degree-180 ESAESM2 archive. Returns
% the same coefficient layout, units and date conventions as GRACE2PLMT_NEW.
%
% Syntax
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt()
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt(timestep, timelim)
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt("midmonth", timelim, ...
%       "Component", "H", "Unit", "GRAV", "Loutput", 60)
%
% Input arguments
%   timestep (optional) - Temporal averaging step.
%       [] returns raw six-hour fields; 'midmonth' returns calendar-month
%       arithmetic means at the exact midpoint of each month.
%       Default step: [].
%   timelim (optional) - Inclusive time range applied to returned timestamps.
%       Two-element datetime or datenum vector. Empty returns all dates.
%       Default range: [].
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
%   ForceNew (name-value) - Whether to rebuild from the source GFC files.
%       Use after adding or changing source data. Default option: false.
%   SaveData (name-value) - Whether to save raw and monthly MAT caches.
%       Default option: true.
%   BeQuiet (name-value) - Whether to suppress cache progress messages.
%       Default option: 0.5 (report cache reads and writes).
%
% Output arguments
%   esmPlmt - Degree, order, cosine, sine coefficient columns.
%       Size: (N, M, 4), or (M, 4, N) for 'traditional'.
%       M = (Loutput + 1) * (Loutput + 2) / 2.
%   esmStdPlmt - Same layout, with NaN uncertainties and valid degree/order.
%       ESAESM2 supplies no formal uncertainties; NaN does not mean zero error.
%   dates - Raw timestamps or calendar-month midpoints, size (N, 1), in TimeFormat.
%
% Notes
%   Reads mtmshc_COMPONENT_YYYYMMDD_HH.gfc files from
%   $IFILES/ESM/ESAESM2/raw/COMPONENT, including its subdirectories.
%   IFILES must name the configured data root; there is no folder argument.
%   Caches contain all available dates at native degree, before time selection,
%   degree truncation/padding and output formatting. Cache names follow SSH2LONLATT:
%   ESAESM2_COMPONENT_UNIT[_LOVESOURCE].mat and
%   ESAESM2_COMPONENT_UNIT[_LOVESOURCE]-Tmidmonth-mean.mat, in
%   $IFILES/ESM/ESAESM2. Love-number source is included for SD only.
%   Cache variables are esmPlmt, esmStdPlmt and dates (datetime, timefirst).
%   Rebuilding a saved raw cache invalidates its monthly cache. Cache hits
%   do not scan GFC files; ForceNew refreshes both levels when needed.
%   Monthly means use [month start, next month start). A month missing any
%   six-hour sample has NaN coefficients and emits an IncompleteMonth warning.
%   Raw mode returns available samples, allowing gaps. Duplicate/off-grid
%   timestamps and incomplete harmonic rows fail when reading source files.
%   Calendar-month means do not reproduce irregular GRACE reporting windows
%   or orbit-dependent observation weights.
%   ESM fields are anomalies. No WGS84 subtraction, demeaning, low-degree
%   replacement, de-aliasing, GIA correction or filtering is applied.
%   Component C00 and degree 1 are preserved. The source radius is read for
%   unit conversion; header metadata are not checked or returned. Inputs are
%   assumed to follow the fully normalized ESAESM2 coefficient convention.
%   Matching the layout does not imply the same physical signal as GSM.
%   SD uses the same density (5517 kg/m^3) and conversion as GRACE2PLMT_NEW.
%   Selecting ISSM ensures pipeline consistency, not an exact inversion of
%   the Love-number spectrum originally used to construct ESAESM2.
%   For S, SD is equivalent surface mass, not a physical surface load.
%
% Data sources
%   ESAESM2. https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/
%
% See also
%   GRACE2PLMT_NEW, SSH2LONLATT, CONVERTGRAVITY
%
% Created by
%   2026/10/29

function [esmPlmt, esmStdPlmt, dates] = esm2plmt(timestep, timelim, options)
    arguments (Input)
        timestep = []
        timelim {mustBeTimeRange} = []
        options.Component (1, 1) string {mustBeMember(options.Component, ...
            ["A", "O", "H", "I", "S", "AOHIS"])} = "O"
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["GRAV", "POT", "SD"])} = "SD"
        options.Loutput double = []
        options.LoveNumSource (1, 1) string = "ISSM"
        options.OutputFormat (1, 1) string {mustBeMember(options.OutputFormat, ...
            ["timefirst", "traditional"])} = "timefirst"
        options.TimeFormat (1, 1) string {mustBeMember(options.TimeFormat, ...
            ["datenum", "datetime"])} = "datenum"
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
    end

    assert(isempty(timestep) || ((ischar(timestep) || isstring(timestep)) && ...
        isscalar(string(timestep)) && strcmp(timestep, 'midmonth')), ...
        'ULMO:esm2plmt:InvalidTimeStep', 'Use [] or ''midmonth'' for timestep.');
    if ~isempty(timelim)
        if isnumeric(timelim)
            timelim = datetime(timelim, 'ConvertFrom', 'datenum');
        end
        assert(~any(isnat(timelim)) && isempty(timelim.TimeZone), ...
            'ULMO:esm2plmt:InvalidTimeRange', 'Use finite, unzoned UTC time limits.');
    end
    if ~isempty(options.Loutput)
        validateattributes(options.Loutput, {'double'}, ...
            {'scalar', 'integer', 'nonnegative', 'finite'});
    end

    dataRoot = getenv('IFILES');
    assert(~isempty(dataRoot), 'ULMO:esm2plmt:MissingDataRoot', ...
        'Set IFILES to the configured data root before reading ESM fields.');
    [rawPath, monthlyPath] = outputpaths(dataRoot, options.Component, ...
        options.Unit, options.LoveNumSource);
    outputPath = rawPath;
    if ~isempty(timestep)
        outputPath = monthlyPath;
    end
    vars = {'esmPlmt', 'esmStdPlmt', 'dates'};

    if ~options.ForceNew && validcache(outputPath, vars)
        data = load(outputPath, vars{:});
        reportcache('Loaded', outputPath, options.BeQuiet);
    else
        if ~options.ForceNew && validcache(rawPath, vars)
            data = load(rawPath, vars{:});
            reportcache('Loaded', rawPath, options.BeQuiet);
        else
            folder = fullfile(dataRoot, 'ESM', 'ESAESM2', 'raw', options.Component);
            [data.esmPlmt, data.esmStdPlmt, data.dates] = ...
                readraw(folder, options.Component, options.Unit, options.LoveNumSource);
            if options.SaveData
                save(rawPath, '-struct', 'data', vars{:}, '-v7.3');
                % Derived averages must not survive a raw-cache rebuild.
                if isfile(monthlyPath)
                    delete(monthlyPath);
                end
                reportcache('Saved', rawPath, options.BeQuiet);
            end
        end
        if ~isempty(timestep)
            [data.esmPlmt, data.esmStdPlmt, data.dates] = ...
                monthlymean(data.esmPlmt, data.dates);
            if options.SaveData
                save(monthlyPath, '-struct', 'data', vars{:}, '-v7.3');
                reportcache('Saved', monthlyPath, options.BeQuiet);
            end
        end
    end
    [esmPlmt, esmStdPlmt, dates] = formatoutput(data.esmPlmt, ...
        data.esmStdPlmt, data.dates, timelim, options.Loutput, ...
        options.OutputFormat, options.TimeFormat);
end

function [rawPath, monthlyPath] = outputpaths(dataRoot, component, unit, loveSource)
    product = "ESAESM2_" + component + "_" + unit;
    if unit == "SD"
        product = product + "_" + replace(loveSource, " ", "");
    end
    folder = fullfile(dataRoot, 'ESM', 'ESAESM2');
    rawPath = fullfile(folder, product + ".mat");
    monthlyPath = fullfile(folder, product + "-Tmidmonth-mean.mat");
end

function valid = validcache(path, vars)
    valid = isfile(path) && all(ismember(vars, who('-file', path)));
end

function reportcache(action, path, beQuiet)
    if beQuiet < 1
        fprintf('[ULMO>esm2plmt] %s %s.\n', action, path);
    end
end

function [esmPlmt, esmStdPlmt, dates] = readraw(folder, component, unit, loveSource)
    assert(isfolder(folder), 'ULMO:esm2plmt:MissingDataFolder', ...
        'ESM component folder does not exist: %s', folder);
    files = dir(fullfile(folder, '**', "mtmshc_" + component + "_*.gfc"));
    assert(~isempty(files), 'ULMO:esm2plmt:MissingData', ...
        'No ESM files found in %s.', folder);
    dates = NaT(numel(files), 1);
    for k = 1:numel(files)
        token = regexp(files(k).name, '^mtmshc_' + component + ...
            '_(\d{8}_\d{2})\.gfc$', 'tokens', 'once');
        assert(~isempty(token), 'ULMO:esm2plmt:InvalidName', ...
            'Invalid ESM filename: %s', files(k).name);
        dates(k) = datetime(token{1}, 'InputFormat', 'yyyyMMdd_HH');
    end
    [dates, order] = sort(dates);
    files = files(order);
    assert(numel(unique(dates)) == numel(dates) && all(mod(hour(dates), 6) == 0), ...
        'ULMO:esm2plmt:InvalidSampling', 'Duplicate or off-grid %s samples.', component);
    for k = 1:numel(files)
        path = fullfile(files(k).folder, files(k).name);
        [plm, equatorRadius] = readesmfile(path);
        if k == 1
            labels = plm(:, 1:2);
            esmPlmt = emptyplmt(numel(files), labels);
            esmStdPlmt = esmPlmt;
        else
            assert(isequal(plm(:, 1:2), labels), ...
                'ULMO:esm2plmt:InvalidCoefficients', 'Harmonic rows differ: %s', path);
        end
        if unit ~= "GRAV"
            plm = convertgravity(plm, 'GRAV', unit, ...
                'EquatorRadius', equatorRadius, 'LoveNumSource', loveSource);
        end
        esmPlmt(k, :, 3:4) = reshape(plm(:, 3:4), 1, size(plm, 1), 2);
    end
end

function [monthly, uncertainty, dates] = monthlymean(raw, rawDates)
    first = dateshift(rawDates(1), 'start', 'month');
    last = dateshift(rawDates(end), 'start', 'month');
    starts = (first:calmonths(1):last)';
    ends = starts + calmonths(1);
    dates = starts + (ends - starts) / 2;
    labels = reshape(raw(1, :, 1:2), [], 2);
    monthly = emptyplmt(numel(dates), labels);
    uncertainty = monthly;
    incomplete = false(numel(dates), 1);
    for k = 1:numel(dates)
        ids = rawDates >= starts(k) & rawDates < ends(k);
        expected = (starts(k):hours(6):ends(k) - hours(6))';
        if ~isequal(rawDates(ids), expected)
            incomplete(k) = true;
            continue
        end
        monthly(k, :, 3:4) = mean(raw(ids, :, 3:4), 1);
    end
    if any(incomplete)
        warning('ULMO:esm2plmt:IncompleteMonth', ...
            '%d incomplete calendar months have NaN coefficients.', nnz(incomplete));
    end
end

function plmt = emptyplmt(nDates, labels)
    plmt = nan(nDates, size(labels, 1), 4);
    plmt(:, :, 1:2) = repmat(reshape(labels, 1, size(labels, 1), 2), nDates, 1, 1);
end

function [plmt, uncertainty, dates] = ...
        formatoutput(plmt, uncertainty, dates, timelim, degree, outputFormat, timeFormat)
    if ~isempty(timelim)
        keep = dates >= timelim(1) & dates <= timelim(2);
        plmt = plmt(keep, :, :);
        uncertainty = uncertainty(keep, :, :);
        dates = dates(keep);
    end
    if ~isempty(degree)
        labels = harmoniclabels(degree);
        nRows = size(labels, 1);
        if nRows <= size(plmt, 2)
            plmt = plmt(:, 1:nRows, :);
            uncertainty = uncertainty(:, 1:nRows, :);
        else
            padded = emptyplmt(numel(dates), labels);
            paddedStd = padded;
            padded(:, 1:size(plmt, 2), :) = plmt;
            paddedStd(:, 1:size(uncertainty, 2), :) = uncertainty;
            plmt = padded;
            uncertainty = paddedStd;
        end
    end
    if timeFormat == "datenum"
        dates = datenum(dates); %#ok<DATNM>
    end
    dates = dates(:);
    if outputFormat == "traditional"
        plmt = permute(plmt, [2, 3, 1]);
        uncertainty = permute(uncertainty, [2, 3, 1]);
    end
end

function [plm, equatorRadius] = readesmfile(path)
    contents = fileread(path);
    endHeader = regexp(contents, '(?m)^end_of_head[^\r\n]*\r?\n', 'end', 'once');
    assert(~isempty(endHeader), 'ULMO:esm2plmt:InvalidHeader', ...
        'Missing end_of_head in %s', path);
    radiusText = regexp(contents(1:endHeader), '(?m)^radius[ \t]+(\S+)', 'tokens', 'once');
    equatorRadius = [];
    if ~isempty(radiusText)
        equatorRadius = str2double(regexprep(radiusText{1}, '[dD]', 'E'));
    end
    data = regexprep(contents(endHeader + 1:end), '[dD]', 'E');
    columns = textscan(data, '%s%f%f%f%f', 'ReturnOnError', false);
    plm = [columns{2:5}];
    plm = sortrows(plm, [1, 2]);
    assert(~isempty(plm) && all(strcmp(columns{1}, 'gfc')) && all(isfinite(plm), 'all') && ...
        all(plm(:, 1:2) >= 0 & plm(:, 1:2) == fix(plm(:, 1:2)), 'all') && ...
        isequal(plm(:, 1:2), harmoniclabels(plm(end, 1))), ...
        'ULMO:esm2plmt:InvalidCoefficients', 'Incomplete or invalid harmonics in %s', path);
end

function labels = harmoniclabels(degree)
    labels = zeros((degree + 1) * (degree + 2) / 2, 2);
    for l = 0:degree
        rows = l * (l + 1) / 2 + (1:l + 1);
        labels(rows, :) = [repmat(l, l + 1, 1), (0:l)'];
    end
end
