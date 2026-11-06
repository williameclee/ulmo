%% ESM2PLMT - Loads cached six-hour or monthly ESAESM2 fields.
% Loads a prepared compact harmonic NetCDF file; use PREPAREESM to import GFC files.
%
% Syntax
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt()
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt(comp, timelim, timestep)
%   [esmPlmt, esmStdPlmt, dates] = esm2plmt("H", timelim, "midmonth", ...
%       "Unit", "GRAV", "Loutput", 60)
%
% Input arguments
%   comp (optional) - Individual component identifier.
%       Options: 'A', 'O', 'H', 'I', 'S', 'AOHIS'.
%       Default component: 'O'.
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
%   OutputFormat (name-value) - 'timefirst' or 'traditional'.
%       Default format: 'timefirst'.
%   TimeFormat (name-value) - 'datenum' or 'datetime'.
%       Default format: 'datetime'.
%   ForceNew (name-value) - Whether to rebuild averaged caches from the prepared NetCDF file.
%       Default option: false.
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
%       sample normalisation n-1. Fewer than two samples give NaN for this
%       empirical estimate. No temporal correlation correction is applied.
%   dates - Sample timestamps or epoch midpoints, size (N, 1), in TimeFormat.
%
% Notes
%   Conversion is in memory; no unit-specific raw copies are saved.
%   A missing six-hour sample invalidates the month with NaNs and a warning.
%   No WGS84 subtraction, demeaning, low-degree replacement, de-aliasing,
%   GIA correction or filtering is applied. C00 and degree 1 are preserved.
%
% Data sources
%   ESAESM2. https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/
%
% See also
%   PREPAREESM, CALENDARMONTHEPOCHS, AVERAGEEPOCHS, CONVERTGRAVITY
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/03, En-Chi Lee (williameclee@gmail.com)

function [esmPlmt, esmStdPlmt, dates] = esm2plmt(comp, varargin)
    % Consume the component before parsing optional/name-value arguments:
    % MATLAB otherwise interprets 'O' as an abbreviation of OutputFormat.
    if nargin == 0
        comp = 'O';
    end

    [esmPlmt, esmStdPlmt, dates] = loadesm(comp, varargin{:});
end

% Parse options and choose a cached, native, or newly averaged product.
function [esmPlmt, esmStdPlmt, dates] = loadesm(comp, timelim, timestep, options)

    arguments (Input)
        comp {mustBeTextScalar, mustBeMember(comp, ["A", "O", "H", "I", "S", "AOHIS"])}
        timelim {mustBeTimeRange} = []
        timestep {mustBeTextScalar, mustBeMember(timestep, ["", "midmonth"])} = ''
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, ["GRAV", "POT", "SD"])} = 'SD'
        options.Loutput double = []
        options.OutputFormat {mustBeTextScalar, mustBeMember(options.OutputFormat, ...
                                  ["timefirst", "traditional"])} = 'timefirst'
        options.TimeFormat {mustBeTextScalar, mustBeMember(options.TimeFormat, ...
                                ["datenum", "datetime"])} = 'datetime'
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
    end

    % Normalise selection limits before cache lookup or native-data reads.
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

    % Resolve prepared inputs and derived products under the configured data root.
    dataRoot = getenv('IFILES');
    assert(~isempty(dataRoot), 'ULMO:esm2plmt:MissingDataRoot', 'Set IFILES before reading ESM fields.');
    dataFolder = fullfile(dataRoot, 'ESM', 'ESAESM2');
    inPath = fullfile(dataFolder, "ESAESM2_" + string(comp) + ".nc");

    % An averaged cache already has the requested units; read only the selected
    % time/degree slab and skip native conversion and averaging on a cache hit.
    if strlength(timestep) > 0
        outPath = getoutputpath(dataFolder, comp, options.Unit, timestep);

        if ~options.ForceNew && isfile(outPath)
            [coeffs, uncs, dates, labels] = readshnc(outPath, timelim, options.Loutput, ...
                'OutputFormat', 'compact');
            reportcache('Loaded', outPath, options.BeQuiet);
            [esmPlmt, esmStdPlmt, dates] = formatoutput(coeffs, uncs, labels, ...
                dates, options.Loutput, options.OutputFormat, options.TimeFormat);
            return
        end

    end

    assert(isfile(inPath), 'ULMO:esm2plmt:MissingPreparedData', ...
        'Prepared ESM file is missing: %s. Import the source folder with prepareesm first.', inPath);

    % Native requests convert only the requested six-hour samples and degree rows.
    if strlength(timestep) == 0
        [coeffs, uncs, dates, labels] = readnative(inPath, timelim, ...
            options.Loutput, options.Unit, false);
        [esmPlmt, esmStdPlmt, dates] = formatoutput(coeffs, uncs, labels, dates, ...
            options.Loutput, options.OutputFormat, options.TimeFormat);
        return
    end

    % On a cache miss, average the entire record at native degree so later
    % requests can reuse the cache regardless of their time or degree limits.
    nativeDates = readshnctime(inPath);
    epochs = averagingepochs(timestep, nativeDates);
    [coeffs, uncs, dates, labels] = averagefile(inPath, epochs, ...
        options.Unit);

    if options.SaveData
        writeshnc(outPath, coeffs, uncs, dates, options.Unit, 'Labels', labels);
        reportcache('Saved', outPath, options.BeQuiet);
    end

    % Build the complete epoch product before selecting the returned range.
    if ~isempty(timelim)
        keep = dates >= timelim(1) & dates <= timelim(2);
        coeffs = coeffs(:, :, keep);
        uncs = uncs(:, :, keep);
        dates = dates(keep);
    end

    [esmPlmt, esmStdPlmt, dates] = formatoutput(coeffs, uncs, labels, dates, ...
        options.Loutput, options.OutputFormat, options.TimeFormat);
end

% Average one window at a time to bound the memory used by native samples.
function [coeffs, uncs, dates, labels] = averagefile(path, epochs, unit)
    labels = [double(ncread(path, 'degree')), double(ncread(path, 'order'))];
    dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;
    coeffs = nan(size(labels, 1), 2, size(epochs, 1));
    uncs = coeffs;

    % Read all harmonic rows but only this half-open time window.
    for k = 1:size(epochs, 1)
        [raw, rawStd, rawDates] = readnative(path, epochs(k, :), [], unit, true);

        if isempty(rawDates)
            warning('ULMO:averageepochs:IncompleteEpoch', 'No samples in epoch starting %s.', string(epochs(k, 1)));
            continue
        end

        % Require the complete six-hour grid; missing samples invalidate the mean.
        % Formal errors are propagated when known, otherwise sample variability is used.
        [coeffs(:, :, k), ~, uncs(:, :, k)] = averageepochs(rawDates, ...
            raw, epochs(k, :), hours(6), 'CheckCompleteness', true, ...
            'SamplingOrigin', datetime(2000, 1, 1), 'SampleStd', rawStd);
    end

end

% Read compact native Stokes coefficients and validate the ESAESM2 cadence.
function [coeffs, uncs, dates, labels] = readnative(path, timelim, degree, unit, exclusive)
    [coeffs, uncs, dates, labels] = readshnc(path, timelim, degree, ...
        'OutputFormat', 'compact', 'EndExclusive', exclusive);
    assert(all(mod(hours(dates - datetime(2000, 1, 1)), 6) == 0), ...
        'ULMO:esm2plmt:InvalidPreparedData', 'Prepared samples must lie on the six-hour grid.');

    % Convert a unit coefficient to obtain each harmonic's scale factor.
    % Standard deviations use its magnitude so uncertainties remain nonnegative.
    if ~strcmp(unit, 'GRAV') && ~isempty(dates)
        radius = ncreadatt(path, '/', 'equatorRadius');
        validateattributes(radius, {'double'}, {'scalar', 'positive', 'finite'});
        factors = convertgravity([labels, ones(size(labels))], 'GRAV', unit, ...
            'EquatorRadius', radius, 'LoveNumSource', 'ISSM');
        coeffs = coeffs .* factors(:, 3:4);
        uncs = uncs .* abs(factors(:, 3:4));
    end

end

% Keep calendar-window generation separate from the averaging calculation.
function epochs = averagingepochs(timestep, nativeDates)

    switch string(timestep)
        case "midmonth"
            epochs = calendarmonthepochs(nativeDates([1, end]));
        otherwise
            error('ULMO:esm2plmt:InvalidTimeStep', 'Unsupported averaging step: %s.', timestep);
    end

end

% Identify derived caches by component, unit and averaging step.
function path = getoutputpath(folder, component, unit, timestep)
    product = "ESAESM2_" + string(component) + "-" + string(unit);

    path = fullfile(folder, product + "-T" + string(timestep) + ".nc");
end

% Report cache activity for both verbose and default intermediate verbosity.
function reportcache(action, path, beQuiet)

    if beQuiet < 1
        fprintf('[ULMO>esm2plmt] %s %s.\n', action, path);
    end

end

% Allocate time-first lmcosi output with valid labels and NaN coefficient slots.
function plmt = emptyplmt(nDates, labels)
    plmt = nan(nDates, size(labels, 1), 4);
    plmt(:, :, 1:2) = repmat(reshape(labels, 1, size(labels, 1), 2), nDates, 1, 1);
end

function [varPlmt, uncPlmt, dates] = ...
        formatoutput(varCoeffs, uncCoeffs, labels, dates, degree, outFmt, timeFmt)
    % Apply degree selection while compact, then expand labels for lmcosi.
    nRows = size(labels, 1);

    if ~isempty(degree)
        labels = harmoniclabels(degree);
        nRows = min(nRows, size(labels, 1));
    end

    % Copy only available rows; requested degrees above the source retain NaNs.
    varPlmt = emptyplmt(numel(dates), labels);
    uncPlmt = varPlmt;
    varPlmt(:, 1:nRows, 3:4) = permute(varCoeffs(1:nRows, :, :), [3, 1, 2]);
    uncPlmt(:, 1:nRows, 3:4) = permute(uncCoeffs(1:nRows, :, :), [3, 1, 2]);

    % Apply date and array-layout conventions only after constructing the output.
    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    dates = dates(:);

    if strcmp(outFmt, 'traditional')
        varPlmt = permute(varPlmt, [2, 3, 1]);
        uncPlmt = permute(uncPlmt, [2, 3, 1]);
    end

end

% Enumerate degree-major (l,m) pairs using triangular row offsets.
function labels = harmoniclabels(degree)
    labels = zeros((degree + 1) * (degree + 2) / 2, 2);

    for l = 0:degree
        rows = l * (l + 1) / 2 + (1:l + 1);
        labels(rows, :) = [repmat(l, l + 1, 1), (0:l)'];
    end

end
