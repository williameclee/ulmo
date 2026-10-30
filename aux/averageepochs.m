%% AVERAGEEPOCHS - Averages data and estimates mean uncertainties over epoch windows.
% Overlapping windows share samples. Time is the third array dimension,
%
% Syntax
%   [avgs, dates] = averageepochs(nativeDates, data, epochs)
%   [averages, dates, uncs] = averageepochs(__, expTimestep, ...
%       "CheckCompleteness", true, "SampleStd", sampleStd)
%
% Input arguments
%   nativeDates - Sorted, unique timestamps or MATLAB datenums.
%       Datetime inputs must be unzoned.
%       Size: (N, 1) or (1, N).
%   data - Real floating-point samples, with NaN denoting missing values.
%       No product-specific or harmonic-label columns.
%       Size: (A, B, N).
%   epochs - Start and exclusive end of each window, in datetime or datenum.
%       Size: (K, 2).
%   expTimestep (optional) - Positive sampling interval as duration or days.
%       Required only when CheckCompleteness is true.
%       Default interval: NaT.
%   CheckCompleteness (name-value) - Whether to require the entire expected
%       timestamp grid and finite data for every element in each epoch.
%       Default option: false (average available nonmissing values).
%   SamplingOrigin (name-value) - Timestamp defining the sampling grid phase.
%       Used only for checking.
%       Default origin: [] (first native timestamp).
%   SampleStd (name-value) - Formal one-sigma uncertainties for input samples.
%       NaN means unknown.
%       Size: same as data.
%       Default uncertainties: [].
%
% Output arguments
%   avgs - Arithmetic means within each [start, end) window.
%       Size: (A, B, K).
%   dates - Exact epoch midpoints, size (K, 1), datetime.
%   uncs - One-sigma uncertainty estimate of each mean, in data units.
%       With known formal errors for every included sample, returns
%       sqrt(sum(sigma_i^2))/n. Otherwise estimates std(samples,0)/sqrt(n)
%       from the same n samples (sample variance normalization n-1).
%       The empirical estimate requires at least two samples; otherwise NaN.
%       Rejected/empty means have NaN uncertainty. A single sample retains its
%       supplied formal sigma. Both estimates assume independent samples;
%       neither accounts for temporal correlation or systematic model error.
%
% See also
%   CALENDARMONTHEPOCHS
%
% Created by
%   2026/10/30, En-Chi Lee (williameclee@gmail.com)

function [avgs, dates, stds] = averageepochs(nativeDates, data, epochs, expTimestep, options)

    arguments (Input)
        nativeDates {mustBeA(nativeDates, {'datetime', 'numeric'}), mustBeVector}
        data (:, :, :) {mustBeFloat, mustBeReal}
        epochs (:, 2) {mustBeA(epochs, {'datetime', 'numeric'})}
        expTimestep (1, 1) ...
            {mustBeA(expTimestep, {'duration', 'numeric'}), mustBeFinite, mustBePositive} = NaT
        options.CheckCompleteness (1, 1) logical = false
        options.SamplingOrigin = []
        options.SampleStd {mustBeFloat, mustBeReal} = []
    end

    arguments (Output)
        avgs (:, :, :) {mustBeFloat}
        dates (:, 1) datetime
        stds (:, :, :) {mustBeFloat}
    end

    nativeDates = parsedates(nativeDates);
    nativeDates = nativeDates(:);
    epochs = parsedates(epochs);
    assert(~isempty(nativeDates) && size(data, 3) == numel(nativeDates), ...
        'ULMO:averageepochs:DimensionMismatch', ...
    'Data time dimension must match nonempty native dates.');
    assert(all(diff(nativeDates) > seconds(0)), ...
        'ULMO:averageepochs:InvalidDates', ...
    'Input dates must be sorted and unique.');
    assert(all(epochs(:, 2) > epochs(:, 1)), ...
        'ULMO:averageepochs:InvalidEpochs', ...
    'Epoch ends must be later than their starts.');
    assert(~any(isinf(data), 'all'), ...
        'ULMO:averageepochs:InvalidData', ...
    'Inifinite value(s) detected; data must be finite or NaN.');
    sampleStd = options.SampleStd;

    if ~isempty(sampleStd)
        assert(isequal(size(sampleStd), size(data)), ...
            'ULMO:averageepochs:DimensionMismatch', ...
        'SampleStd must have the same dimensions as data.');
        assert(all(isnan(sampleStd) | (isfinite(sampleStd) & sampleStd >= 0), 'all'), ...
            'ULMO:averageepochs:InvalidUncertainty', ...
        'SampleStd must be nonnegative and finite, or NaN.');
    end

    if options.CheckCompleteness

        if isnumeric(expTimestep)
            expTimestep = days(expTimestep);
        end

        assert(~isnat(expTimestep) && ~isnan(expTimestep), ...
            'ULMO:averageepochs:InvalidTimeStep', ...
        'Checking data completeness requires supplying a valid expected time step.');
        origin = options.SamplingOrigin;

        if isempty(origin)
            origin = nativeDates(1);
        else
            origin = parsedates(origin);
            assert(isscalar(origin), 'ULMO:averageepochs:InvalidOrigin', ...
            'SamplingOrigin must be a single timestamp.');
        end

    end

    dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;

    % Preallocate arrays
    avgs = nan(size(data, 1), size(data, 2), size(epochs, 1), 'like', data);
    stds = avgs;
    incomplete = false(size(epochs, 1), 1);

    % Populate each epoch's value
    for k = 1:size(epochs, 1)
        start = epochs(k, 1);
        finish = epochs(k, 2);
        isInEpoch = nativeDates >= start & nativeDates < finish;

        if ~any(isInEpoch)
            incomplete(k) = true;
            continue
        end

        if options.CheckCompleteness
            nExpData = floor((finish - max(start, origin)) / expTimestep);

            if nExpData ~= sum(isInEpoch)
                incomplete(k) = true;
                continue
            end

        end

        samples = data(:, :, isInEpoch);
        valid = ~isnan(samples);
        validCnt = sum(valid, 3);
        rejected = validCnt == 0;

        if options.CheckCompleteness
            rejected = validCnt ~= nnz(isInEpoch);
        end

        means = mean(samples, 3, 'omitmissing');
        means(rejected) = NaN;
        avgs(:, :, k) = means;
        incomplete(k) = any(rejected, 'all');
        uncs = std(samples, 0, 3, 'omitmissing') ./ sqrt(validCnt);
        uncs(validCnt < 2) = NaN;

        if ~isempty(sampleStd)
            sigmas = sampleStd(:, :, isInEpoch);
            sigmas(~valid) = 0;
            formal = sqrt(sum(sigmas .^ 2, 3)) ./ validCnt;
            known = all(~isnan(sigmas), 3);
            uncs(known) = formal(known);
        end

        uncs(rejected) = NaN;
        stds(:, :, k) = uncs;
    end

    if any(incomplete)
        warning('ULMO:averageepochs:IncompleteEpoch', ...
            '%d epochs contain empty or incomplete averages.', nnz(incomplete));
    end

end

%% Subfunctions
function dates = parsedates(dates)

    if isnumeric(dates)
        assert(all(isfinite(dates), 'all'), ...
            'ULMO:averageepochs:InvalidDates', 'Dates must be finite.');
        dates = datetime(dates, 'ConvertFrom', 'datenum');
    end

    assert(isdatetime(dates) && ~any(isnat(dates) | isinf(dates), 'all') && ...
        isempty(dates.TimeZone), ...
        'ULMO:averageepochs:InvalidDates', ...
    'Use finite, unzoned timestamps or datenums.');
end
