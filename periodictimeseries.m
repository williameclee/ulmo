%% PERIODICTIMESERIES
%
% Created by
%   2025/06/03, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2026/09/27, williameclee@arizona.edu (@williameclee)

function [polyCoeffs, periodicCoeffs, dataFit, polyCoeffSigmas, periodicCoeffSigmas] = ...
        periodictimeseries(t, x, sigma, p, periods, options)

    arguments (Input)
        t {mustBeA(t, {'numeric', 'datetime', 'duration'}), mustBeVector, mustBeNonempty}
        x {mustBeNumeric, mustBeVector, mustBeNonempty}
        sigma {mustBeNumeric, mustBeVectorOrEmpty} = []
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
        polyCoeffs (:, 1) {mustBeNumeric}
        periodicCoeffs (:, :) {mustBeNumeric}
        dataFit (:, 1) {mustBeNumeric}
        polyCoeffSigmas (:, 1) {mustBeNumeric}
        periodicCoeffSigmas (:, :) {mustBeNumeric}

    end

    N = numel(t);

    if N ~= numel(x)
        error('dates and data must have the same number of elements');
    end

    t = t(:);
    x = x(:);
    sigma = sigma(:);

    if ~isempty(sigma) && any(sigma <= 0)
        error('sigma must be positive');
    end

    isTime = false;

    if isdatetime(t)
        isTime = true;
        t = years(t - datetime(year(t(1)), 1, 1));
    end

    isValid = ~isnan(t) & ~isnan(x);
    t_tofit = t(isValid);
    x_tofit = x(isValid);

    if ~isempty(sigma)
        sigma_tofit = sigma(isValid);
    else
        sigma_tofit = [];
    end

    N = numel(t_tofit);

    if ~isempty(periods) && isduration(periods)
        periods = years(periods);
    end

    X = zeros([N, p + 1 + 2 * numel(periods)]);
    X(:, 1:p + 1) = t_tofit .^ (0:p);

    for i = 1:numel(periods)

        switch options.PeriodicFormat
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

    if strcmp(options.PeriodicFormat, "amp-phase")
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
