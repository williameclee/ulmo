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

function [polyCoeffs, periodicCoeffs, dataFit, polyCoeffSigmas, periodicCoeffSigmas] = ...
        fitsingleseries(t_tofit, x_tofit, sigma_tofit, p, periods, periodicFormat, isTime)
    N = numel(t_tofit);
    %% Fitting
    X = zeros([N, p + 1 + 2 * numel(periods)]);
    X(:, 1:p + 1) = t_tofit .^ (0:p);

    for i = 1:numel(periods)

        switch periodicFormat
            case "sin-cos"
                X(:, p + 1 + (i - 1) * 2 + [1, 2]) = ...
                    [sin(2 * pi * t_tofit / periods(i)), cos(2 * pi * t_tofit / periods(i))];
            otherwise
                X(:, p + 1 + (i - 1) * 2 + [1, 2]) = ...
                    [cos(2 * pi * t_tofit / periods(i)), sin(2 * pi * t_tofit / periods(i))];
        end

    end

    if ~isempty(sigma_tofit)
        W = diag(1 ./ sigma_tofit .^ 2);
        coeffs = (X' * W * X) \ (X' * W * x_tofit);
        residuals = x_tofit - X * coeffs;
        weighted_residual_variance = sum((residuals ./ sigma_tofit) .^ 2) / (N - size(X, 2));
        cov_matrix = weighted_residual_variance * inv(X' * W * X);
        coeffSigmas = sqrt(diag(cov_matrix));
    else
        coeffs = X \ x_tofit;
        % Estimate uncertainties of coeffs when sigma is not provided
        residuals = x_tofit - X * coeffs;
        residual_variance = sum(residuals .^ 2) / (N - size(X, 2));
        cov_matrix = residual_variance * inv(X' * X);
        coeffSigmas = sqrt(diag(cov_matrix));
    end

    dataFit = X * coeffs;

    polyCoeffs = coeffs(1:p + 1);
    periodicCoeffs = coeffs(p + 2:end);

    polyCoeffSigmas = coeffSigmas(1:p + 1);
    periodicCoeffSigmas = coeffSigmas(p + 2:end);

    if strcmp(periodicFormat, "amp-phase")
        periodicCoeffs = reshape(periodicCoeffs, [2, numel(periods)]);
        periodicCoeffs = ...
            [sqrt(sum(periodicCoeffs .^ 2, 1)); ...
             wrapTo2Pi(atan2(periodicCoeffs(2, :), periodicCoeffs(1, :)))];
        periodicCoeffs = periodicCoeffs';

        if isTime
            periodicCoeffs(:, 2) = periodicCoeffs(:, 2) / (2 * pi) * days(years(1));
        end

        periodicCoeffSigmas = reshape(periodicCoeffSigmas, [2, numel(periods)]);
        periodicCoeffSigmas = ...
            [sqrt(sum(periodicCoeffSigmas .^ 2, 1)); ...
             wrapTo2Pi(atan2(periodicCoeffSigmas(2, :), periodicCoeffSigmas(1, :)))];
        periodicCoeffSigmas = periodicCoeffSigmas';

        if isTime
            periodicCoeffSigmas(:, 2) = periodicCoeffSigmas(:, 2) / (2 * pi) * days(years(1));
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
