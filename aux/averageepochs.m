%% AVERAGEEPOCHS - Averages data and propagates uncertainties over epoch windows.
% Overlapping windows share samples. Time is the third array dimension,
% following the gridded-data convention used by SSH2LONLATT and STERIC2LONLATT.
%
% Syntax
%   [averages, dates] = averageepochs(nativeDates, data, epochs)
%   [averages, dates, averageStd] = averageepochs(__, expectedStep, ...
%       "CheckCompleteness", true, "SampleStd", sampleStd)
%
% Input arguments
%   nativeDates - Sorted, unique timestamps or MATLAB datenums.
%       Size: (N, 1) or (1, N). Datetime inputs must be unzoned.
%   data - Real floating-point samples, with NaN denoting missing values.
%       Size: (A, B, N). No product-specific or harmonic-label columns.
%   epochs - Start and exclusive end of each window, in datetime or datenum.
%       Size: (K, 2). Windows may overlap; their input order is preserved.
%   expectedStep (optional) - Positive sampling interval as duration or days.
%       Required only when CheckCompleteness is true. Default interval: [].
%   CheckCompleteness (name-value) - Whether to require the entire expected
%       timestamp grid and finite data for every element in each epoch.
%       Default option: false (average available nonmissing values).
%   SamplingOrigin (name-value) - Timestamp defining the sampling grid phase.
%       Default origin: [] (first native timestamp). Used only for checking.
%   SampleStd (name-value) - Formal one-sigma uncertainties for input samples.
%       Size: same as data. NaN means unknown. Default uncertainties: [].
%
% Output arguments
%   averages - Arithmetic means within each [start, end) window.
%       Size: (A, B, K), same floating-point class as data.
%   dates - Exact epoch midpoints, size (K, 1), datetime.
%   averageStd - Formal one-sigma uncertainty of each mean, in data units.
%       Assuming independent sample errors, sigma(mean) = sqrt(sum(sigma_i^2))/n,
%       using exactly the n samples included in the mean. Missing input sigma,
%       absent samples or rejected epochs give NaN. A single sample retains
%       its supplied sigma. This is neither sample variability nor an empirical
%       standard error estimated from variability. Correlated input errors
%       require their covariance and are not represented by SampleStd alone.
%
% Notes
%   Without checking, NaN data are omitted separately for each array element;
%   sigma values are selected with the same mask. Unknown sigma is never
%   omitted for an otherwise valid data sample. With checking, an incomplete
%   timestamp grid invalidates the entire epoch; missing data invalidate only
%   the affected elements. Empty or incomplete results emit IncompleteEpoch.
%   No interpolation or calendar assumptions are made. The first two outputs
%   retain the AVERAGEEPOCHS interface used by the GRACE epoch averaging branch.
%
% See also
%   CALENDARMONTHEPOCHS, SSH2LONLATT, STERIC2LONLATT
%
% Created by
%   2026/10/30

function [averages, dates, averageStd] = averageepochs(nativeDates, data, epochs, expectedStep, options)
    arguments (Input)
        nativeDates {mustBeA(nativeDates, {'datetime', 'numeric'}), mustBeVector}
        data (:, :, :) {mustBeFloat, mustBeReal}
        epochs (:, 2) {mustBeA(epochs, {'datetime', 'numeric'})}
        expectedStep = []
        options.CheckCompleteness (1, 1) logical = false
        options.SamplingOrigin = []
        options.SampleStd {mustBeFloat, mustBeReal} = []
    end
    nativeDates = parsedates(nativeDates);
    epochs = parsedates(epochs);
    nativeDates = nativeDates(:);
    assert(~isempty(nativeDates) && size(data, 3) == numel(nativeDates), ...
        'ULMO:averageepochs:DimensionMismatch', 'Data time dimension must match nonempty native dates.');
    assert(all(diff(nativeDates) > seconds(0)), 'ULMO:averageepochs:InvalidDates', ...
        'Native dates must be sorted and unique.');
    assert(all(epochs(:, 2) > epochs(:, 1)), 'ULMO:averageepochs:InvalidEpochs', ...
        'Epoch ends must be later than their starts.');
    assert(~any(isinf(data), 'all'), 'ULMO:averageepochs:InvalidData', ...
        'Data must be finite or NaN.');
    sampleStd = options.SampleStd;
    if ~isempty(sampleStd)
        assert(isequal(size(sampleStd), size(data)), 'ULMO:averageepochs:DimensionMismatch', ...
            'SampleStd must have the same dimensions as data.');
        assert(all(isnan(sampleStd) | (isfinite(sampleStd) & sampleStd >= 0), 'all'), ...
            'ULMO:averageepochs:InvalidUncertainty', 'SampleStd must be nonnegative and finite, or NaN.');
    end
    if options.CheckCompleteness
        if isnumeric(expectedStep)
            expectedStep = days(expectedStep);
        end
        assert(isduration(expectedStep) && isscalar(expectedStep) && ...
            isfinite(seconds(expectedStep)) && expectedStep > seconds(0), ...
            'ULMO:averageepochs:InvalidTimeStep', 'Checking requires a positive finite expected time step.');
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
    averages = nan(size(data, 1), size(data, 2), size(epochs, 1), 'like', data);
    averageStd = averages;
    incomplete = false(size(epochs, 1), 1);
    for k = 1:size(epochs, 1)
        start = epochs(k, 1);
        finish = epochs(k, 2);
        ids = nativeDates >= start & nativeDates < finish;
        if ~any(ids)
            incomplete(k) = true;
            continue
        end
        if options.CheckCompleteness
            firstSample = origin + ceil(seconds(start - origin) / seconds(expectedStep)) * expectedStep;
            expected = (firstSample:expectedStep:finish)';
            expected = expected(expected < finish);
            if ~isequal(nativeDates(ids), expected)
                incomplete(k) = true;
                continue
            end
        end
        samples = data(:, :, ids);
        valid = ~isnan(samples);
        count = sum(valid, 3);
        rejected = count == 0;
        if options.CheckCompleteness
            rejected = count ~= nnz(ids);
        end
        means = mean(samples, 3, 'omitmissing');
        means(rejected) = NaN;
        averages(:, :, k) = means;
        incomplete(k) = any(rejected, 'all');
        if ~isempty(sampleStd)
            sigmas = sampleStd(:, :, ids);
            sigmas(~valid) = 0;
            uncertainty = sqrt(sum(sigmas .^ 2, 3)) ./ count;
            uncertainty(rejected) = NaN;
            averageStd(:, :, k) = uncertainty;
        end
    end
    if any(incomplete)
        warning('ULMO:averageepochs:IncompleteEpoch', ...
            '%d epochs contain empty or incomplete averages.', nnz(incomplete));
    end
end

function dates = parsedates(dates)
    if isnumeric(dates)
        assert(all(isfinite(dates), 'all'), 'ULMO:averageepochs:InvalidDates', 'Dates must be finite.');
        dates = datetime(dates, 'ConvertFrom', 'datenum');
    end
    assert(isdatetime(dates) && ~any(isnat(dates) | isinf(dates), 'all') && ...
        isempty(dates.TimeZone), 'ULMO:averageepochs:InvalidDates', ...
        'Use finite, unzoned timestamps or datenums.');
end
