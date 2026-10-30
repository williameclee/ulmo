%% AVERAGEESMEPOCHS - Averages six-hour ESM harmonics over supplied epochs.
% Accepts arbitrary, including overlapping, windows without calendar logic.
%
% Syntax
%   [esmPlmt, esmStdPlmt, dates] = averageesmepochs(nativeDates, raw, epochs)
%
% Input arguments
%   nativeDates - Sorted, unique six-hour UTC timestamps.
%       Size: (N, 1), unzoned datetime.
%   raw - Native coefficients, columns degree, order, cosine, sine.
%       Size: (N, M, 4), in the desired output units.
%   epochs - Start and exclusive end of each averaging window.
%       Size: (K, 2), unzoned datetime. Windows may overlap and boundaries
%       need not fall on the six-hour grid. Input order is preserved.
%
% Output arguments
%   esmPlmt - Mean coefficients for each supplied epoch.
%       Size: (K, M, 4).
%   esmStdPlmt - Matching labels with NaN coefficient uncertainties.
%   dates - Exact epoch midpoints, size (K, 1), datetime.
%
% Notes
%   Selects six-hour UTC samples in [start, end) with equal weights. Missing
%   expected samples or windows containing no samples produce NaN coefficient
%   columns and an IncompleteEpoch warning. Labels remain valid.
%   No interpolation, inferred calendar boundaries or time-range filtering
%   is performed. This is the same epoch-window contract as AVERAGEEPOCHS,
%   with ESM's complete-sampling requirement and harmonic array layout.
%
% See also
%   ESM2PLMT, ESMEPOCHS
%
% Created by
%   2026/10/30

function [esmPlmt, esmStdPlmt, dates] = averageesmepochs(nativeDates, raw, epochs)
    arguments (Input)
        nativeDates (:, 1) datetime
        raw (:, :, 4) double
        epochs (:, 2) datetime
    end
    assert(~isempty(nativeDates) && size(raw, 1) == numel(nativeDates), ...
        'ULMO:averageesmepochs:DimensionMismatch', ...
        'Raw time dimension must match nonempty native dates.');
    assert(~any(isnat(nativeDates)) && isempty(nativeDates.TimeZone) && ...
        all(diff(nativeDates) > seconds(0)) && ...
        all(mod(hour(nativeDates), 6) == 0 & minute(nativeDates) == 0 & second(nativeDates) == 0), ...
        'ULMO:averageesmepochs:InvalidDates', 'Native dates must be sorted, unique six-hour UTC timestamps.');
    assert(~any(isnat(epochs), 'all') && isempty(epochs.TimeZone) && ...
        all(epochs(:, 2) > epochs(:, 1)), 'ULMO:averageesmepochs:InvalidEpochs', ...
        'Epochs must contain valid unzoned UTC start/end times with end > start.');
    dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;
    esmPlmt = nan(size(epochs, 1), size(raw, 2), 4);
    esmPlmt(:, :, 1:2) = repmat(raw(1, :, 1:2), size(epochs, 1), 1, 1);
    esmStdPlmt = esmPlmt;
    incomplete = false(size(epochs, 1), 1);
    for k = 1:size(epochs, 1)
        start = epochs(k, 1);
        finish = epochs(k, 2);
        ids = nativeDates >= start & nativeDates < finish;
        % Align to the UTC sampling grid, independently of epoch boundaries.
        dayStart = dateshift(start, 'start', 'day');
        firstSample = dayStart + hours(6 * ceil(hours(start - dayStart) / 6));
        expected = (firstSample:hours(6):finish)';
        expected = expected(expected < finish);
        if isempty(expected) || ~isequal(nativeDates(ids), expected)
            incomplete(k) = true;
            continue
        end
        esmPlmt(k, :, 3:4) = mean(raw(ids, :, 3:4), 1);
    end
    if any(incomplete)
        warning('ULMO:esm2plmt:IncompleteEpoch', ...
            '%d epochs lack complete six-hour sampling and have NaN coefficients.', nnz(incomplete));
    end
end
