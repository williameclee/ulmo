%% AVERAGEEPOCHS - Averages data and estimates mean uncertainties over epoch windows.
% Overlapping windows share samples. Time is the third array dimension,
%
% Syntax
%   [avgs, dates] = averageepochs(inDates, data, epochs)
%   [avgs, dates, uncs] = averageepochs(__, expTimestep, ...
%       "CheckCompleteness", true, "SampleStd", sampleStd)
%
% Input arguments
%   inDates - Sorted, unique timestamps or MATLAB datenums.
%       Datetime inputs must be unzoned.
%       Size: (N, 1) or (1, N).
%   data - Real floating-point samples, with NaN denoting missing values.
%       No product-specific or harmonic-label columns.
%       Size: (A, B, N).
%   epochs - Start and exclusive end of each window, in datetime or datenum.
%       Size: (K, 2).
%   expTimestep (optional) - Positive sampling interval as duration or days.
%       Required only when CheckCompleteness is true.
%       Default interval: NaN.
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
%       from the same n samples (sample variance normalisation n-1).
%       Rejected/empty means have NaN uncertainty. A single sample retains its
%       supplied formal sigma. Both estimates assume independent samples.
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)

function [avgs, dates, stds] = averageepochs(inDates, data, epochs, expTimestep, options)

    arguments (Input)
        inDates {mustBeA(inDates, {'datetime', 'numeric'}), mustBeVector}
        data (:, :, :) {mustBeFloat, mustBeReal}
        epochs (:, 2) {mustBeA(epochs, {'datetime', 'numeric'})}
        expTimestep (1, 1) ...
            {mustBeA(expTimestep, {'duration', 'numeric'})} = NaN
        options.CheckCompleteness (1, 1) logical = false
        options.SamplingOrigin = []
        options.SampleStd {mustBeFloat, mustBeReal} = []
    end

    arguments (Output)
        avgs (:, :, :) {mustBeFloat}
        dates (:, 1) datetime
        stds (:, :, :) {mustBeFloat}
    end

    [inDates, epochs, dataStd, expTimestep, origin] = validatePreprocessInputs( ...
        inDates, epochs, data, options.SampleStd, ...
        options.CheckCompleteness, expTimestep, options.SamplingOrigin);

    % Use window midpoints as output timestamps, even for rejected windows.
    dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;

    % Preallocate arrays
    avgs = nan(size(data, 1), size(data, 2), size(epochs, 1), 'like', data);
    stds = avgs;
    incompl = false(size(epochs, 1), 1);

    % Populate each epoch's value
    for iEpoch = 1:size(epochs, 1)
        epochStart = epochs(iEpoch, 1);
        epochEnd = epochs(iEpoch, 2);
        % Exclude the end boundary so adjacent windows do not share that sample.
        isInEpoch = inDates >= epochStart & inDates < epochEnd;

        if ~any(isInEpoch)
            incompl(iEpoch) = true;
            continue
        end

        if options.CheckCompleteness
            % Build the exact grid inside [start, finish), including shifted origins.
            % Matching counts alone would not detect missing or off-grid timestamps.
            firstIdx = ceil((epochStart - origin) / expTimestep);
            lastIdx = ceil((epochEnd - origin) / expTimestep) - 1;
            expdDates = origin + (firstIdx:lastIdx)' * expTimestep;

            if ~isequal(inDates(isInEpoch), expdDates)
                incompl(iEpoch) = true;
                continue
            end

        end

        % Count available samples per element; missing values can differ by element
        epochData = data(:, :, isInEpoch);
        valid = ~isnan(epochData);
        validCnt = sum(valid, 3);
        rejected = validCnt == 0;

        if options.CheckCompleteness
            rejected = validCnt ~= nnz(isInEpoch);
        end

        % Strict completeness rejects any element missing a sample; otherwise omit NaNs
        means = mean(epochData, 3, 'omitmissing');
        means(rejected) = NaN;
        avgs(:, :, iEpoch) = means;
        incompl(iEpoch) = any(rejected, 'all');
        % Estimate the standard error under sample independence.
        % One sample cannot provide an empirical variance.
        epochStd = std(epochData, 0, 3, 'omitmissing') ./ sqrt(validCnt);
        epochStd(validCnt < 2) = NaN;

        if ~isempty(dataStd)
            % Use formal propagation only when every contributing sample has an error.
            % Missing data contribute neither variance nor sample count; unknown errors
            % for valid data retain the empirical estimate above.
            sigmas = dataStd(:, :, isInEpoch);
            sigmas(~valid) = 0;
            formal = sqrt(sum(sigmas .^ 2, 3)) ./ validCnt;
            known = all(~isnan(sigmas), 3);
            epochStd(known) = formal(known);
        end

        epochStd(rejected) = NaN;
        stds(:, :, iEpoch) = epochStd;
    end

    % Report rejected windows once, after all output slices have been filled
    if any(incompl)
        warning('ULMO:averageepochs:IncompleteEpoch', ...
            '%d epochs contain empty or incomplete averages.', nnz(incompl));
    end

end

%% Subfunctions
% Converts datenums and reject zoned or invalid dates before comparisons.
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

function [inDates, epochs, dataStd, expTimestep, origin] = validatePreprocessInputs( ...
        inDates, epochs, data, dataStd, checkComplet, expTimestep, origin)
    % Normalise dates before checking alignment with the third data dimension
    inDates = parsedates(inDates);
    inDates = inDates(:);
    epochs = parsedates(epochs);

    assert(~isempty(inDates) && size(data, 3) == numel(inDates), ...
        'ULMO:averageepochs:DimensionMismatch', ...
    'Data time dimension must match nonempty native dates.');
    assert(all(diff(inDates) > seconds(0)), ...
        'ULMO:averageepochs:InvalidDates', ...
    'Input dates must be sorted and unique.');
    assert(all(epochs(:, 2) > epochs(:, 1)), ...
        'ULMO:averageepochs:InvalidEpochs', ...
    'Epoch ends must be later than their starts.');
    assert(~any(isinf(data), 'all'), ...
        'ULMO:averageepochs:InvalidData', ...
    'Inifinite value(s) detected; data must be finite or NaN.');

    % Validate formal errors separately; NaN marks an unknown uncertainty
    if ~isempty(dataStd)
        assert(isequal(size(dataStd), size(data)), ...
            'ULMO:averageepochs:DimensionMismatch', ...
        'SampleStd must have the same dimensions as data.');
        assert(all(isnan(dataStd) | (isfinite(dataStd) & dataStd >= 0), 'all'), ...
            'ULMO:averageepochs:InvalidUncertainty', ...
        'SampleStd must be nonnegative and finite, or NaN.');
    end

    % Resolve the expected sampling interval and grid phase only when required
    if ~checkComplet
        expTimestep = NaN;
        origin = NaT;
        return
    end

    if isnumeric(expTimestep)
        expTimestep = days(expTimestep);
    end

    assert(isfinite(expTimestep) && expTimestep > seconds(0), ...
        'ULMO:averageepochs:InvalidTimeStep', ...
    'Checking data completeness requires supplying a valid expected time step.');

    if isempty(origin)
        origin = inDates(1);
    else
        origin = parsedates(origin);
        assert(isscalar(origin), 'ULMO:averageepochs:InvalidOrigin', ...
        'SamplingOrigin must be a single timestamp.');
    end

end
