%% ESM2PLMT - Loads cached six-hour or monthly ESAESM2 fields.
% Loads a prepared native lmcosi MAT file; use PREPAREESM to import GFC files. Returns
% the same coefficient layout, units and date conventions as GRACE2PLMT_NEW.
%
% Syntax
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt()
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt(component, timelim, timestep)
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt("H", timelim, "midmonth", ...
%       "Unit", "GRAV", "Loutput", 60)
%
% Input arguments
%   component (optional) - Individual component identifier.
%       Options: 'A', 'O', 'H', 'I', 'S', 'AOHIS'. Default component: 'O'.
%   timelim (optional) - Inclusive time range applied to returned timestamps.
%       Two-element datetime or datenum vector. Empty returns all dates.
%       Default range: [].
%   timestep (optional) - Temporal averaging step as char or string scalar.
%       '' returns six-hour fields; 'midmonth' returns calendar-month means.
%       Default step: ''.
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
%   ForceNew (name-value) - Whether to rebuild averaged caches from the prepared MAT file.
%       Never reads GFC files. Default option: false.
%   SaveData (name-value) - Whether to save averaged MAT caches.
%       Default option: true.
%   BeQuiet (name-value) - Whether to suppress cache progress messages.
%       Default option: 0.5 (report cache reads and writes).
%
% Input file content
%   $IFILES/ESM/ESAESM2/ESAESM2_COMPONENT.mat, written by PREPAREESM:
%   esmPlmt - Native-degree dimensionless Stokes lmcosi array.
%       Size: (N, M, 4), timefirst, columns degree, order, C and S.
%   dates - Sorted unique six-hour UTC timestamps. Size: (N, 1), datetime.
%   equatorRadius - Radius used for unit conversion. Size: (N, 1). Unit: m.
%   esmStdPlmt (optional) - Formal one-sigma Stokes errors with the same layout
%       as esmPlmt. Absent for ESAESM2; missing errors are represented by NaN.
%
% Output arguments
%   esmPlmt - Degree, order, cosine, sine coefficient columns.
%       Size: (N, M, 4), or (M, 4, N) for 'traditional'.
%       M = (Loutput + 1) * (Loutput + 2) / 2.
%   esmStdPlmt - Uncertainties in the same units/layout, with valid labels.
%       Raw errors are NaN unless supplied. Means propagate independent formal
%       errors when all are known; otherwise use std(samples)/sqrt(n), with
%       sample normalization n-1. Fewer than two samples give NaN for this
%       empirical estimate. No temporal correlation correction is applied.
%   dates - Sample timestamps or epoch midpoints, size (N, 1), in TimeFormat.
%
% Notes
%   IFILES names the fixed data root; ESM2PLMT never reads a GFC folder.
%   Raw requests read only the selected dates/degrees from the prepared v7.3
%   MAT file. Conversion is in memory; no unit-specific raw copies are saved.
%   Averaged cache hits return immediately. Otherwise every supplied epoch is
%   averaged at native degree, one window at a time, then saved before output
%   selection/truncation. PREPAREESM invalidates these caches on reimport.
%   Cache names follow SSH2LONLATT: ESAESM2_COMPONENT_UNIT[_LOVESOURCE]-TSTEP.mat
%   in $IFILES/ESM/ESAESM2, with variables esmPlmt, esmStdPlmt and dates.
%   Love-number source is included for SD only. ForceNew bypasses this cache;
%   SaveData=false suppresses writes without modifying the prepared input.
%   CALENDARMONTHEPOCHS supplies [start, end) windows to AVERAGEEPOCHS.
%   A missing six-hour sample invalidates the month with NaNs and a warning.
%   Calendar months do not reproduce irregular GRACE windows/observation weights.
%   No WGS84 subtraction, demeaning, low-degree replacement, de-aliasing,
%   GIA correction or filtering is applied. C00 and degree 1 are preserved.
%   SD uses CONVERTGRAVITY with density 5517 kg/m^3 and the CF frame.
%   ISSM matches the pipeline convention, not necessarily ESM's Earth model.
%   For S, SD is equivalent mass rather than a physical surface load.
%
% Data sources
%   ESAESM2. https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/
%
% See also
%   PREPAREESM, GRACE2PLMT_NEW, SSH2LONLATT, CALENDARMONTHEPOCHS, AVERAGEEPOCHS, CONVERTGRAVITY
%
% Created by
%   2026/10/29
%
% Last modified
%   2026/10/30

function [esmPlmt, esmStdPlmt, dates] = esm2plmt(component, varargin)
    % Consume the component before parsing optional/name-value arguments:
    % MATLAB otherwise interprets 'O' as an abbreviation of OutputFormat.
    if nargin == 0
        component = 'O';
    end
    [esmPlmt, esmStdPlmt, dates] = loadesm(component, varargin{:});
end

function [esmPlmt, esmStdPlmt, dates] = loadesm(component, timelim, timestep, options)
    arguments (Input)
        component {mustBeTextScalar, mustBeMember(component, ...
            ["A", "O", "H", "I", "S", "AOHIS"])}
        timelim {mustBeTimeRange} = []
        timestep {mustBeTextScalar, mustBeMember(timestep, ["", "midmonth"])} = ''
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, ["GRAV", "POT", "SD"])} = 'SD'
        options.Loutput double = []
        options.LoveNumSource {mustBeTextScalar} = 'ISSM'
        options.OutputFormat {mustBeTextScalar, mustBeMember(options.OutputFormat, ...
            ["timefirst", "traditional"])} = 'timefirst'
        options.TimeFormat {mustBeTextScalar, mustBeMember(options.TimeFormat, ...
            ["datenum", "datetime"])} = 'datenum'
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
    end
    if ~isempty(timelim)
        if isnumeric(timelim)
            timelim = datetime(timelim, 'ConvertFrom', 'datenum');
        end
        assert(~any(isnat(timelim) | isinf(timelim)) && isempty(timelim.TimeZone), ...
            'ULMO:esm2plmt:InvalidTimeRange', 'Use finite, unzoned UTC time limits.');
    end
    if ~isempty(options.Loutput)
        validateattributes(options.Loutput, {'double'}, {'scalar', 'integer', 'nonnegative', 'finite'});
    end
    dataRoot = getenv('IFILES');
    assert(~isempty(dataRoot), 'ULMO:esm2plmt:MissingDataRoot', 'Set IFILES before reading ESM fields.');
    folder = fullfile(dataRoot, 'ESM', 'ESAESM2');
    inputPath = fullfile(folder, "ESAESM2_" + string(component) + ".mat");
    vars = {'esmPlmt', 'esmStdPlmt', 'dates'};
    if strlength(timestep) > 0
        outputPath = cachepath(folder, component, options.Unit, options.LoveNumSource, timestep);
        if ~options.ForceNew && isfile(outputPath) && all(ismember(vars, who('-file', outputPath)))
            data = load(outputPath, vars{:});
            reportcache('Loaded', outputPath, options.BeQuiet);
            [esmPlmt, esmStdPlmt, dates] = formatoutput(data.esmPlmt, data.esmStdPlmt, ...
                data.dates, timelim, options.Loutput, options.OutputFormat, options.TimeFormat);
            return
        end
    end
    assert(isfile(inputPath), 'ULMO:esm2plmt:MissingPreparedData', ...
        'Prepared ESM file is missing: %s. Import the source folder with prepareesm first.', inputPath);
    source = matfile(inputPath);
    variables = who(source);
    assert(all(ismember({'esmPlmt', 'dates', 'equatorRadius'}, variables)), ...
        'ULMO:esm2plmt:InvalidPreparedData', 'Prepared ESM file lacks required variables.');
    nativeDates = source.dates;
    shape = size(source, 'esmPlmt');
    assert(isdatetime(nativeDates) && iscolumn(nativeDates) && ~isempty(nativeDates) && ...
        isempty(nativeDates.TimeZone) && all(isfinite(nativeDates)) && ...
        all(diff(nativeDates) > seconds(0)) && ...
        all(mod(hours(nativeDates - datetime(2000, 1, 1)), 6) == 0) && ...
        numel(shape) == 3 && shape(1) == numel(nativeDates) && shape(3) == 4 && ...
        isequal(size(source, 'equatorRadius'), [numel(nativeDates), 1]), ...
        'ULMO:esm2plmt:InvalidPreparedData', 'Invalid prepared dimensions or six-hour timestamps.');
    labels = reshape(source.esmPlmt(1, :, 1:2), [], 2);
    assert(isequal(labels, harmoniclabels(labels(end, 1))), ...
        'ULMO:esm2plmt:InvalidPreparedData', 'Prepared harmonics must use complete ascending degree/order rows.');
    hasStd = ismember('esmStdPlmt', variables);
    assert(~hasStd || isequal(size(source, 'esmStdPlmt'), shape), ...
        'ULMO:esm2plmt:InvalidPreparedData', 'Prepared uncertainties must match coefficient dimensions.');
    if strlength(timestep) == 0
        ids = true(size(nativeDates));
        if ~isempty(timelim)
            ids = nativeDates >= timelim(1) & nativeDates <= timelim(2);
        end
        nRows = size(labels, 1);
        if ~isempty(options.Loutput)
            nRows = min(nRows, (options.Loutput + 1) * (options.Loutput + 2) / 2);
        end
        [data.esmPlmt, data.esmStdPlmt] = readnative(source, find(ids), labels(1:nRows, :), ...
            hasStd, options.Unit, options.LoveNumSource);
        data.dates = nativeDates(ids);
    else
        epochs = averagingepochs(timestep, nativeDates);
        data.dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;
        data.esmPlmt = emptyplmt(size(epochs, 1), labels);
        data.esmStdPlmt = data.esmPlmt;
        for k = 1:size(epochs, 1)
            ids = find(nativeDates >= epochs(k, 1) & nativeDates < epochs(k, 2));
            if isempty(ids)
                warning('ULMO:averageepochs:IncompleteEpoch', 'No samples in epoch starting %s.', string(epochs(k, 1)));
                continue
            end
            [raw, rawStd] = readnative(source, ids, labels, hasStd, options.Unit, options.LoveNumSource);
            [coefficients, ~, coefficientStd] = averageepochs(nativeDates(ids), ...
                permute(raw(:, :, 3:4), [2, 3, 1]), epochs(k, :), hours(6), ...
                'CheckCompleteness', true, 'SamplingOrigin', datetime(2000, 1, 1), ...
                'SampleStd', permute(rawStd(:, :, 3:4), [2, 3, 1]));
            data.esmPlmt(k, :, 3:4) = permute(coefficients, [3, 1, 2]);
            data.esmStdPlmt(k, :, 3:4) = permute(coefficientStd, [3, 1, 2]);
        end
        if options.SaveData
            save(outputPath, '-struct', 'data', vars{:}, '-v7.3');
            reportcache('Saved', outputPath, options.BeQuiet);
        end
    end
    [esmPlmt, esmStdPlmt, dates] = formatoutput(data.esmPlmt, data.esmStdPlmt, ...
        data.dates, timelim, options.Loutput, options.OutputFormat, options.TimeFormat);
end

function [plmt, uncertainty] = readnative(source, ids, labels, hasStd, unit, loveSource)
    plmt = emptyplmt(numel(ids), labels);
    uncertainty = plmt;
    if isempty(ids)
        return
    end
    plmt = source.esmPlmt(ids, 1:size(labels, 1), 1:4);
    if hasStd
        uncertainty(:, :, 3:4) = source.esmStdPlmt(ids, 1:size(labels, 1), 3:4);
    end
    if ~strcmp(unit, 'GRAV')
        % Conversion is linear in radius. Evaluate degree factors once and
        % apply each sample's source radius to coefficients and errors alike.
        factors = convertgravity([labels, ones(size(labels))], 'GRAV', unit, ...
            'EquatorRadius', 1, 'LoveNumSource', loveSource);
        radii = source.equatorRadius(ids, 1);
        assert(all(isfinite(radii) & radii > 0), 'ULMO:esm2plmt:InvalidPreparedData', ...
            'Source radii must be positive and finite for conversion.');
        scale = radii .* reshape(factors(:, 3:4), 1, size(labels, 1), 2);
        plmt(:, :, 3:4) = plmt(:, :, 3:4) .* scale;
        uncertainty(:, :, 3:4) = uncertainty(:, :, 3:4) .* abs(scale);
    end
end

function epochs = averagingepochs(timestep, nativeDates)
    switch string(timestep)
        case "midmonth"
            epochs = calendarmonthepochs(nativeDates);
        otherwise
            error('ULMO:esm2plmt:InvalidTimeStep', 'Unsupported averaging step: %s.', timestep);
    end
end

function path = cachepath(folder, component, unit, loveSource, timestep)
    product = "ESAESM2_" + string(component) + "_" + string(unit);
    if strcmp(unit, 'SD')
        product = product + "_" + replace(string(loveSource), " ", "");
    end
    path = fullfile(folder, product + "-T" + string(timestep) + ".mat");
end

function reportcache(action, path, beQuiet)
    if beQuiet < 1
        fprintf('[ULMO>esm2plmt] %s %s.\n', action, path);
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
    if strcmp(timeFormat, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end
    dates = dates(:);
    if strcmp(outputFormat, 'traditional')
        plmt = permute(plmt, [2, 3, 1]);
        uncertainty = permute(uncertainty, [2, 3, 1]);
    end
end

function labels = harmoniclabels(degree)
    labels = zeros((degree + 1) * (degree + 2) / 2, 2);
    for l = 0:degree
        rows = l * (l + 1) / 2 + (1:l + 1);
        labels(rows, :) = [repmat(l, l + 1, 1), (0:l)'];
    end
end
