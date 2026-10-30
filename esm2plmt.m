%% ESM2PLMT - Loads cached six-hour or monthly ESAESM2 fields.
% Loads a prepared compact harmonic NetCDF file; use PREPAREESM to import GFC files. Returns
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
%   ForceNew (name-value) - Whether to rebuild averaged caches from the prepared NetCDF file.
%       Never reads GFC files. Default option: false.
%   SaveData (name-value) - Whether to save averaged NetCDF caches.
%       Default option: true.
%   BeQuiet (name-value) - Whether to suppress cache progress messages.
%       Default option: 0.5 (report cache reads and writes).
%
% Input file content
%   $IFILES/ESM/ESAESM2/ESAESM2_COMPONENT.nc, written by PREPAREESM:
%   degree, order - uint16 vectors, stored once on the lm dimension.
%   coefficients - Native-degree dimensionless Stokes values. Size: (M, 2, N).
%       Dimensions: lm, cs (cosine then sine), time. Precision: double.
%   time - Double seconds since 2000-01-01 00:00:00 UTC. Size: (N, 1).
%   equatorRadius (global attribute) - Constant reference radius, metres.
%   uncertainty (optional) - Formal one-sigma Stokes errors, same shape as
%       coefficients. Absent for ESAESM2; missing errors are returned as NaN.
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
%   Native and cached requests read time first, then only selected date/degree
%   slices from NetCDF. Conversion is in memory; no unit-specific raw copies are saved.
%   Averaged cache hits return immediately. Otherwise every supplied epoch is
%   averaged at native degree, one window at a time, then saved before output
%   selection/truncation. PREPAREESM invalidates these caches on reimport.
%   Cache names follow SSH2LONLATT: ESAESM2_COMPONENT_UNIT[_LOVESOURCE]-TSTEP.nc
%   in $IFILES/ESM/ESAESM2, using the same compact coefficient layout.
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
    inputPath = fullfile(folder, "ESAESM2_" + string(component) + ".nc");
    if strlength(timestep) > 0
        outputPath = cachepath(folder, component, options.Unit, options.LoveNumSource, timestep);
        if ~options.ForceNew && isfile(outputPath)
            [coefficients, uncertainty, dates, labels] = readshnc(outputPath, timelim, options.Loutput, ...
                'OutputFormat', 'compact');
            reportcache('Loaded', outputPath, options.BeQuiet);
            [esmPlmt, esmStdPlmt, dates] = formatoutput(coefficients, uncertainty, labels, ...
                dates, options.Loutput, options.OutputFormat, options.TimeFormat);
            return
        end
    end
    assert(isfile(inputPath), 'ULMO:esm2plmt:MissingPreparedData', ...
        'Prepared ESM file is missing: %s. Import the source folder with prepareesm first.', inputPath);
    if strlength(timestep) == 0
        [coefficients, uncertainty, dates, labels] = readnative(inputPath, timelim, ...
            options.Loutput, options.Unit, options.LoveNumSource, false);
    else
        nativeDates = readshnctime(inputPath);
        epochs = averagingepochs(timestep, nativeDates);
        [coefficients, uncertainty, dates, labels] = averagefile(inputPath, epochs, ...
            options.Unit, options.LoveNumSource);
        if options.SaveData
            writeshnc(outputPath, coefficients, uncertainty, dates, options.Unit, 'Labels', labels);
            reportcache('Saved', outputPath, options.BeQuiet);
        end
        % Build the complete epoch product before selecting the returned range.
        if ~isempty(timelim)
            keep = dates >= timelim(1) & dates <= timelim(2);
            coefficients = coefficients(:, :, keep);
            uncertainty = uncertainty(:, :, keep);
            dates = dates(keep);
        end
    end
    [esmPlmt, esmStdPlmt, dates] = formatoutput(coefficients, uncertainty, labels, dates, ...
        options.Loutput, options.OutputFormat, options.TimeFormat);
end

function [coefficients, uncertainty, dates, labels] = averagefile(path, epochs, unit, loveSource)
    labels = [double(ncread(path, 'degree')), double(ncread(path, 'order'))];
    dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;
    coefficients = nan(size(labels, 1), 2, size(epochs, 1));
    uncertainty = coefficients;
    for k = 1:size(epochs, 1)
        [raw, rawStd, rawDates] = readnative(path, epochs(k, :), [], unit, loveSource, true);
        if isempty(rawDates)
            warning('ULMO:averageepochs:IncompleteEpoch', 'No samples in epoch starting %s.', string(epochs(k, 1)));
            continue
        end
        [coefficients(:, :, k), ~, uncertainty(:, :, k)] = averageepochs(rawDates, ...
            raw, epochs(k, :), hours(6), 'CheckCompleteness', true, ...
            'SamplingOrigin', datetime(2000, 1, 1), 'SampleStd', rawStd);
    end
end

function [coefficients, uncertainty, dates, labels] = readnative(path, timelim, degree, unit, loveSource, exclusive)
    [coefficients, uncertainty, dates, labels] = readshnc(path, timelim, degree, ...
        'OutputFormat', 'compact', 'EndExclusive', exclusive);
    assert(all(mod(hours(dates - datetime(2000, 1, 1)), 6) == 0), ...
        'ULMO:esm2plmt:InvalidPreparedData', 'Prepared samples must lie on the six-hour grid.');
    if ~strcmp(unit, 'GRAV') && ~isempty(dates)
        radius = ncreadatt(path, '/', 'equatorRadius');
        validateattributes(radius, {'double'}, {'scalar', 'positive', 'finite'});
        factors = convertgravity([labels, ones(size(labels))], 'GRAV', unit, ...
            'EquatorRadius', radius, 'LoveNumSource', loveSource);
        coefficients = coefficients .* factors(:, 3:4);
        uncertainty = uncertainty .* abs(factors(:, 3:4));
    end
end

function epochs = averagingepochs(timestep, nativeDates)
    switch string(timestep)
        case "midmonth"
            epochs = calendarmonthepochs(nativeDates([1, end]));
        otherwise
            error('ULMO:esm2plmt:InvalidTimeStep', 'Unsupported averaging step: %s.', timestep);
    end
end

function path = cachepath(folder, component, unit, loveSource, timestep)
    product = "ESAESM2_" + string(component) + "_" + string(unit);
    if strcmp(unit, 'SD')
        product = product + "_" + replace(string(loveSource), " ", "");
    end
    path = fullfile(folder, product + "-T" + string(timestep) + ".nc");
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
        formatoutput(coefficients, coefficientStd, labels, dates, degree, outputFormat, timeFormat)
    % Apply degree selection while compact, then expand labels for lmcosi.
    nRows = size(labels, 1);
    if ~isempty(degree)
        labels = harmoniclabels(degree);
        nRows = min(nRows, size(labels, 1));
    end
    plmt = emptyplmt(numel(dates), labels);
    uncertainty = plmt;
    plmt(:, 1:nRows, 3:4) = permute(coefficients(1:nRows, :, :), [3, 1, 2]);
    uncertainty(:, 1:nRows, 3:4) = permute(coefficientStd(1:nRows, :, :), [3, 1, 2]);
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
