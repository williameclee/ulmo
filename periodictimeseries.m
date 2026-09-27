%% PERIODICTIMESERIES
%
% Created by
%   2025/06/03, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2026/09/27, williameclee@arizona.edu (@williameclee)

function [polyCoeffs, periodicCoeffs, xFit, polyCoeffSigmas, periodicCoeffSigmas] = ...
        periodictimeseries(t, x, sigma, p, periods, options)

    arguments (Input)
        t {mustBeA(t, {'numeric', 'datetime', 'duration'}), mustBeVector, mustBeNonempty}
        x {mustBeNumeric, mustBeReal, mustBeMatrix, mustBeNonempty}
        sigma {mustBeNumeric, mustBeReal, mustBeMatrix} = []
        p (1, 1) double {mustBeInteger, mustBeNonnegative} = 2
        periods {mustBeA(periods, {'numeric', 'duration'}), mustBeVectorOrEmpty} = []
        options.PeriodicFormat (1, 1) string ...
            {mustBeMember(options.PeriodicFormat, ["sin-cos", "cos-sin", "amp-phase"])} ...
            = "cos-sin"
        options.PolynomialFormat (1, 1) string ...
            {mustBeMember(options.PolynomialFormat, ["coefficients", "average-derivatives"])} ...
            = "coefficients"
    end

    arguments (Output)
        polyCoeffs (:, :) {mustBeNumeric}
        periodicCoeffs {mustBeNumeric}
        xFit (:, :) {mustBeNumeric}
        polyCoeffSigmas (:, :) {mustBeNumeric}
        periodicCoeffSigmas {mustBeNumeric}
    end

    %% Input sanitisation
    t = t(:);
    nTime = numel(t);
    [x, isTrans] = normalisedata(double(x), nTime, 'x');
    nSeries = size(x, 2);

    if ~isempty(sigma)
        sigma = normalisedata(double(sigma), nTime, 'sigma');

        if size(sigma, 2) == 1
            sigma = repmat(sigma, 1, nSeries);
        elseif size(sigma, 2) ~= nSeries
            error('ULMO:periodictimeseries:SizeMismatch', ...
            'sigma must have the same size as x.');
        end

        if any(isinf(sigma) | sigma <= 0, 'all')
            error('ULMO:periodictimeseries:InvalidSigma', ...
            'Nonmissing sigma values must be finite and positive.');
        end

    end

    isTime = isdatetime(t);

    if isTime
        firstValid = find(~isnat(t), 1);

        if isempty(firstValid)
            error('periodictimeseries:InvalidTime', 'At least one valid time is required.');
        end

        t = years(t - dateshift(t(firstValid), 'start', 'year'));
    elseif isduration(t)
        t = years(t);
    else

        if ~isreal(t)
            error('periodictimeseries:InvalidTime', 'Times must be real.');
        end

        t = double(t);
    end

    if isduration(periods)
        periods = years(periods);
    end

    periods = double(periods(:).');

    if ~isreal(periods) || any(~isfinite(periods) | periods <= 0)
        error('periodictimeseries:InvalidPeriods', 'Periods must be finite, real, and positive.');
    end

    polyCoeffs = nan(p + 1, nSeries);
    polyCoeffSigmas = nan(p + 1, nSeries);

    periodicCoeffs = nan(numel(periods), 2, nSeries);
    periodicCoeffSigmas = periodicCoeffs;
    xFit = nan(size(x));

    %% Fitting
    for iSeries = 1:nSeries
        isValid = isfinite(t) & isfinite(x(:, iSeries));
        sigma_tofit = [];

        if ~isempty(sigma)
            isValid = isValid & isfinite(sigma(:, iSeries));
            sigma_tofit = sigma(isValid, iSeries);
        end

        nMinValids = p + 1 + 2 * numel(periods);

        if nnz(isValid) <= nMinValids
            error('ULMO:periodictimeseries:InsufficientData', ...
                ['Series %d needs more valid observations than fitted parameters.', ...
             '%d observations is needed, but only got %d valid ones.'], ...
                iSeries, nMinValids, nnz(isValid));
        end

        [polyCoeffs(:, iSeries), periodicCoeffs(:, :, iSeries), xFit(isValid, iSeries), ...
             polyCoeffSigmas(:, iSeries), periodicCoeffSigmas(:, :, iSeries)] = ...
            fitsingleseries( ...
            t(isValid), x(isValid, iSeries), sigma_tofit, p, periods, ...
            options.PeriodicFormat, isTime);
    end

    if isTrans
        xFit = xFit.';
    end

end

%% Subfunctions
function [data, transposed] = normalisedata(data, nTimes, name)
    transposed = false;

    if size(data, 1) == nTimes
        return
    elseif size(data, 2) == nTimes
        data = data.';
        transposed = true;
    else
        error('ULMO:periodictimeseries:SizeMismatch', ...
            '%s must have a dimension matching the shared time axis (%d).', name, nTimes);
    end

end

function [polys, harmons, xFit, polySigmas, harmonSigmas] = ...
        fitsingleseries(t, x, sigma, p, periods, harmonFmt, isTime)
    N = numel(t);
    %% Fitting
    X = zeros([N, p + 1 + 2 * numel(periods)]);
    X(:, 1:p + 1) = t .^ (0:p);

    for i = 1:numel(periods)
        X(:, p + 1 + (i - 1) * 2 + [1, 2]) = ...
            [cos(2 * pi * t / periods(i)), sin(2 * pi * t / periods(i))];
    end

    if ~isempty(sigma)
        W = diag(1 ./ sigma .^ 2);
        coeffs = (X' * W * X) \ (X' * W * x);
        res = x - X * coeffs;
        wght_res_var = sum((res ./ sigma) .^ 2) / (N - size(X, 2));
        cov_matrix = wght_res_var * inv(X' * W * X);
        coeffSigmas = sqrt(diag(cov_matrix));
    else
        coeffs = X \ x;
        % Estimate uncertainties of coeffs when sigma is not provided
        res = x - X * coeffs;
        res_var = sum(res .^ 2) / (N - size(X, 2));
        cov_matrix = res_var * inv(X' * X);
        coeffSigmas = sqrt(diag(cov_matrix));
    end

    xFit = X * coeffs;

    polys = coeffs(1:p + 1);
    harmons = coeffs(p + 2:end);

    polySigmas = coeffSigmas(1:p + 1);
    harmonSigmas = coeffSigmas(p + 2:end);

    if strcmpi(harmonFmt, "sin-cos")
        cosCoeffs = harmons(1:2:end);
        sinCoeffs = harmons(2:2:end);
        harmons(1:2:end) = sinCoeffs;
        harmons(2:2:end) = cosCoeffs;

        cosCoeffs = harmonSigmas(1:2:end);
        sinCoeffs = harmonSigmas(2:2:end);
        harmonSigmas(1:2:end) = sinCoeffs;
        harmonSigmas(2:2:end) = cosCoeffs;
    elseif strcmp(harmonFmt, "amp-phase")
        harmons = reshape(harmons, [2, numel(periods)]);
        harmons = ...
            [sqrt(sum(harmons .^ 2, 1)); ...
             wrapTo2Pi(atan2(harmons(2, :), harmons(1, :)))];
        harmons = harmons';

        if isTime
            harmons(:, 2) = harmons(:, 2) / (2 * pi) * days(years(1));
        end

        harmonSigmas = reshape(harmonSigmas, [2, numel(periods)]);
        harmonSigmas = ...
            [sqrt(sum(harmonSigmas .^ 2, 1)); ...
             wrapTo2Pi(atan2(harmonSigmas(2, :), harmonSigmas(1, :)))];
        harmonSigmas = harmonSigmas';

        if isTime
            harmonSigmas(:, 2) = harmonSigmas(:, 2) / (2 * pi) * days(years(1));
        end

    end

end

function mustBeVectorOrEmpty(value)

    if ~isempty(value) && ~isvector(value)
        error('periodictimeseries:InvalidShape', 'Input must be a vector or empty.');
    end

end

function mustBeCompatibleRange(value, t)
    sameType = (isnumeric(t) && isnumeric(value)) ...
        || (isdatetime(t) && isdatetime(value)) ...
        || (isduration(t) && isduration(value));

    if ~sameType || (isnumeric(value) && ~isreal(value)) || any(ismissing(value), 'all')
        error('periodictimeseries:InvalidRange', ...
        'Ranges must contain real, nonmissing endpoints of the same time type as t.');
    end

    if value(1) > value(2)
        error('periodictimeseries:InvalidRange', ...
        'Range endpoints must be in ascending order.');
    end

end
