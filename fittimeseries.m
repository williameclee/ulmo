%% FITTIMESERIES - Fits a polynomial and periodic terms to time series sharing a time axis.
% Each series is fitted independently by ordinary or weighted least squares.
% The polynomial can be returned as coefficients or interval-average
% derivatives, with associated standard errors.
%
% Syntax
%   polys = fittimeseries(t, x)
%   [polys, harmons, xFit] = fittimeseries(t, x, sigma, p, periods)
%   [polys, harmons, xFit, polySigmas, harmonSigmas] = fittimeseries(__)
%   [polys, harmons, xFit, polySigmas, harmonSigmas, polyCovariances] = fittimeseries(__)
%   [__] = fittimeseries(__, "Name", Value)
%
% Input arguments
%   t - Shared observation times
%       Times may be supplied in either orientation. Numeric values are used
%       directly, without shifting their origin. Datetimes are converted to
%       years since January 1 of the first nonmissing input time's year.
%       Durations are converted to years without shifting their origin.
%       Data type: NUMERIC | DATETIME | DURATION
%       Dimension: [N x 1] | [1 x N]
%   x - Observations for M time series
%       The dimension matching N is treated as the observation dimension.
%       If x is square, rows are observations. Nonfinite observations are
%       excluded from fitting separately for each series.
%       Data type: NUMERIC (real; converted to DOUBLE)
%       Dimension: [N x M] | [M x N]
%   sigma (optional) - Observation standard deviations used for weighting
%       A vector is shared by all series; a matrix supplies individual
%       uncertainties. Nonmissing entries must be finite and positive.
%       NaN entries exclude the corresponding observations from fitting.
%       The weights are 1./sigma.^2; see Notes for covariance scaling.
%       The default value is [], for ordinary least squares.
%       Unit: same as x
%       Data type: NUMERIC (real; converted to DOUBLE)
%       Dimension: [] | [N x 1] | [1 x N] | [N x M] | [M x N]
%   p (optional) - Polynomial degree
%       A nonnegative integer. The polynomial includes powers 0 through p.
%       The default value is 2.
%       Data type: DOUBLE
%       Dimension: scalar
%   periods (optional) - Periods of the fitted harmonics
%       Each period adds both a cosine and a sine term. Values must be
%       finite and positive. Duration values are converted to years.
%       Numeric periods must use the same units as the converted time axis:
%       years for datetime/duration t, otherwise the numeric units of t.
%       The default value is [], for no periodic terms.
%       Data type: NUMERIC | DURATION
%       Dimension: [] | [K x 1] | [1 x K]
%   PeriodicFormat (name-value) - Representation of the harmonic outputs
%       - "cos-sin": cosine and sine coefficients (default).
%       - "sin-cos": sine and cosine coefficients.
%       - "amp-phase": amplitude and phase lag, using A*cos(omega*t-phi).
%       This option does not change the fitted model or reconstruction.
%       Data type: STRING | CHAR
%   PolynomialFormat (name-value) - Representation of the polynomial outputs
%       - "coefficients": coefficients in ascending power order (default).
%       - "average-derivatives": interval averages of derivatives of orders
%           0 through p over AverageRange. Order 0 is the mean level,
%           order 1 the mean trend, and order 2 the mean acceleration.
%       This option does not change the fitted model or reconstruction.
%       Data type: STRING | CHAR
%   Reconstruction (name-value) - Components included in xFit
%       - "includeharmonics": polynomial and periodic terms (default).
%       - "omitharmonics": polynomial component only.
%       - "onlyharmonics": periodic terms only; zero at valid times when
%           periods is empty, and NaN at nonfinite times.
%       All components are still included in the fit regardless of which
%       components are reconstructed in xFit.
%       Data type: STRING | CHAR
%   FitRange (name-value) - Inclusive time interval selecting observations for fitting
%       Endpoints must be nonmissing and in ascending order, with the same
%       time type and units as t. This does not restrict reconstruction.
%       The default interval is [min(t), max(t)], omitting missing times.
%       Data type: NUMERIC | DATETIME | DURATION (matching t)
%       Dimension: [1 x 2]
%   AverageRange (name-value) - Interval used for average-derivative outputs
%       Endpoints use the same time type and units as t. In
%       "average-derivatives" mode, endpoints must be finite and the
%       interval must have positive length. This does not select fit data.
%       The default interval is FitRange, independently of missing
%       observations in x or sigma or shorter input-time coverage.
%       The interval may extend beyond FitRange, implying extrapolation.
%       This option has no effect in "coefficients" mode.
%       Data type: NUMERIC | DATETIME | DURATION (matching t)
%       Dimension: [1 x 2]
%
% Output arguments
%   polys - Polynomial coefficients or average derivatives
%       Row k+1 corresponds to power/derivative order k, for k = 0:p.
%       In "average-derivatives" mode, each entry is the integral of the
%       kth derivative of the fitted polynomial divided by interval length.
%       For q(t) = b0+b1*t+b2*t^2 over [a,b], the average trend is
%       b1+b2*(a+b), and the average acceleration is 2*b2.
%       Unit: units of x / (time unit)^k; the time unit is years for
%           datetime/duration input and the units of t for numeric input
%       Data type: DOUBLE
%       Dimension: [p+1 x M]
%   harmons - Harmonic coefficients or amplitude/phase pairs
%       Each row corresponds to a supplied period, in the input order.
%       The two columns follow PeriodicFormat. For datetime input,
%       amplitude/phase output expresses the phase lag in days; otherwise
%       it is in radians. Phase is referenced to the converted time origin.
%       At exactly zero amplitude, phase is NaN.
%       Unit: same as x for coefficients/amplitudes; days or radians for phase
%       Data type: DOUBLE
%       Dimension: [K x 2 x M] (or [K x 2] for one series)
%           K = numel(periods); empty when periods is empty
%   xFit - Reconstructed values at the original input times
%       Reconstruction uses raw fitted coefficients regardless of
%       PolynomialFormat. Values are evaluated even at times outside
%       FitRange or where x/sigma was missing, provided the time is valid.
%       Unit: same as x
%       Data type: DOUBLE
%       Dimension: same as x, including its original orientation and order
%   polySigmas - Standard errors of polys
%       Average-derivative errors propagate the full fitted polynomial
%       coefficient covariance, including off-diagonal terms.
%       The units and dimensions are the same as polys.
%       Data type: DOUBLE
%   harmonSigmas - Standard errors of harmons
%       Amplitude/phase errors use first-order covariance propagation.
%       At exactly zero amplitude, both transformed errors are NaN.
%       The units and dimensions are the same as harmons.
%       Data type: DOUBLE
%
%   polyCovariances - Raw polynomial coefficient covariance for each series
%       Always in ascending coefficient order, independent of PolynomialFormat.
%       Dimension: [p+1 x p+1 x M]; scaled by residual variance like polySigmas.
%       Data type: DOUBLE
%
% Notes
%   N is the number of input times, M the number of series, and K the number
%   of harmonic periods. Invalid times and nonfinite observations are
%   excluded from each fit. Each series must retain more than p+1+2*K
%   observations; reliable coefficient estimates also require full rank.
%   Coefficient covariance is scaled by residual variance, using the
%   residual degrees of freedom. Weighted fits use weighted residual
%   variance, so supplied sigma values act as relative uncertainty weights.
%   Reported uncertainties are standard errors, not confidence intervals;
%   the fit does not model temporal correlation in residuals.
%   FitRange and AverageRange serve separate purposes: a full-record fit
%   can report an average derivative over a shorter interval, while a
%   restricted fit can report over its own or a different averaging interval.
%
% Author
%   2025/06/03, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified by
%   2026/09/28, En-Chi Lee (williameclee@arizona.edu)

function [polys, harmons, xFit, polySigmas, harmonSigmas, polyCovariances] = ...
        fittimeseries(t, x, sigma, p, periods, options)

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
        options.Reconstruction (1, 1) string ...
            {mustBeMember(options.Reconstruction, ["includeharmonics", "omitharmonics", "onlyharmonics"])} ...
            = "includeharmonics"
        options.FitRange (1, 2) {mustBeCompatibleRange(options.FitRange, t)} = ...
            [min(t(:), [], 'omitmissing'), max(t(:), [], 'omitmissing')]
        options.AverageRange (1, 2) {mustBeCompatibleRange(options.AverageRange, t)}
    end

    arguments (Output)
        polys (:, :) {mustBeNumeric}
        harmons {mustBeNumeric}
        xFit (:, :) {mustBeNumeric}
        polySigmas (:, :) {mustBeNumeric}
        harmonSigmas {mustBeNumeric}
        polyCovariances {mustBeNumeric}
    end

    %% Input sanitisation
    if ~isfield(options, 'AverageRange')
        options.AverageRange = options.FitRange;
    end

    t = t(:);
    nTime = numel(t);
    [x, isTrans] = normalisedata(double(x), nTime, 'x');
    nSeries = size(x, 2);

    if ~isempty(sigma)
        sigma = normalisedata(double(sigma), nTime, 'sigma');

        if size(sigma, 2) == 1
            sigma = repmat(sigma, 1, nSeries);
        elseif size(sigma, 2) ~= nSeries
            error('ULMO:fittimeseries:SizeMismatch', ...
            'sigma must have the same size as x.');
        end

        if any(isinf(sigma) | sigma <= 0, 'all')
            error('ULMO:fittimeseries:InvalidSigma', ...
            'Nonmissing sigma values must be finite and positive.');
        end

    end

    fitRange = options.FitRange(:);
    averageRange = options.AverageRange(:);

    isTime = isdatetime(t);

    if isTime
        firstValid = find(~isnat(t), 1);

        if isempty(firstValid)
            error('fittimeseries:InvalidTime', 'At least one valid time is required.');
        end

        tRef = dateshift(t(firstValid), 'start', 'year');
        t = years(t - tRef);
        fitRange = years(fitRange - tRef);
        averageRange = years(averageRange - tRef);
    elseif isduration(t)
        t = years(t);
        fitRange = years(fitRange);
        averageRange = years(averageRange);
    else

        if ~isreal(t)
            error('fittimeseries:InvalidTime', 'Times must be real.');
        end

        t = double(t);
        fitRange = double(fitRange);
        averageRange = double(averageRange);
    end

    if isduration(periods)
        periods = years(periods);
    end

    periods = double(periods(:).');

    if ~isreal(periods) || any(~isfinite(periods) | periods <= 0)
        error('fittimeseries:InvalidPeriods', 'Periods must be finite, real, and positive.');
    end

    polyTransform = [];

    if options.PolynomialFormat == "average-derivatives"

        if any(~isfinite(averageRange)) || averageRange(2) <= averageRange(1)
            error('ULMO:fittimeseries:InvalidAverageRange', ...
            'AverageRange must have finite endpoints and positive length.');
        end

        polyTransform = averagederivativematrix(p, averageRange);
    end

    polys = nan(p + 1, nSeries);
    polySigmas = nan(p + 1, nSeries);
    polyCovariances = nan(p + 1, p + 1, nSeries);

    harmons = nan(numel(periods), 2, nSeries);
    harmonSigmas = harmons;
    xFit = nan(size(x));

    %% Fitting
    for iSeries = 1:nSeries
        isValid = isfinite(t) & isfinite(x(:, iSeries)) & (t >= fitRange(1) & t <= fitRange(2));
        sigma_tofit = [];

        if ~isempty(sigma)
            isValid = isValid & isfinite(sigma(:, iSeries));
            sigma_tofit = sigma(isValid, iSeries);
        end

        nMinValids = p + 1 + 2 * numel(periods);

        if nnz(isValid) <= nMinValids
            error('ULMO:fittimeseries:InsufficientData', ...
                ['Series %d needs more valid observations than fitted parameters.', ...
             '%d observations is needed, but only got %d valid ones.'], ...
                iSeries, nMinValids, nnz(isValid));
        end

        [polys(:, iSeries), harmons(:, :, iSeries), xFit(:, iSeries), ...
             polySigmas(:, iSeries), harmonSigmas(:, :, iSeries), polyCovariances(:, :, iSeries)] = ...
            fitsingleseries( ...
            t(isValid), x(isValid, iSeries), sigma_tofit, p, periods, ...
            options.PeriodicFormat, isTime, t, options.Reconstruction, polyTransform);
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
        error('ULMO:fittimeseries:SizeMismatch', ...
            '%s must have a dimension matching the shared time axis (%d).', name, nTimes);
    end

end

function [polys, harmons, xFit, polySigmas, harmonSigmas, polyCovariance] = ...
        fitsingleseries(t, x, sigma, p, periods, harmonFmt, isTime, tFit, fitMethod, polyTransform)
    %% Fitting
    N = numel(t);
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

    polyCovariance = cov_matrix(1:p + 1, 1:p + 1);
    polys = coeffs(1:p + 1);
    polySigmas = coeffSigmas(1:p + 1);

    if ~isempty(polyTransform)
        polys = polyTransform * polys;
        polyCov = polyTransform * polyCovariance * polyTransform.';
        polyCov = (polyCov + polyCov.') / 2;
        polySigmas = sqrt(max(diag(polyCov), 0));
    end

    harmons = coeffs(p + 2:end);
    harmonSigmas = coeffSigmas(p + 2:end);

    harmons = reshape(harmons, [2, numel(periods)])';
    harmonSigmas = reshape(harmonSigmas, [2, numel(periods)])';

    %% Post-processing coefficients
    if strcmpi(harmonFmt, "sin-cos")
        harmons = harmons(:, [2, 1]);
        harmonSigmas = harmonSigmas(:, [2, 1]);
    elseif strcmp(harmonFmt, "amp-phase")
        nPeriods = numel(periods);
        harmons = nan(nPeriods, 2);
        harmonSigmas = nan(nPeriods, 2);

        for k = 1:nPeriods
            % Indices in the original coefficient vector and full covariance.
            idx = p + 1 + 2 * (k - 1) + [1, 2];
            c = coeffs(idx(1));
            s = coeffs(idx(2));
            Ccs = cov_matrix(idx, idx);

            A = hypot(c, s);
            harmons(k, 1) = A;

            % At zero amplitude, phase and this linearisation are undefined.
            if A == 0
                continue
            end

            phi = mod(atan2(s, c), 2 * pi);

            J = [c / A, s / A;
                 -s / A ^ 2, c / A ^ 2];

            Cap = J * Ccs * J.';
            Cap = (Cap + Cap.') / 2; % Remove roundoff asymmetry

            % Phase in days for datetime input, otherwise in radians.
            if isTime
                phaseScale = days(years(periods(k))) / (2 * pi);
                phi = phi * phaseScale;

                D = diag([1, phaseScale]);
                Cap = D * Cap * D.';
            end

            harmons(k, :) = [A, phi];

            % Clamp tiny negative variances caused by roundoff.
            harmonSigmas(k, :) = sqrt(max(diag(Cap), 0)).';
        end

    end

    %% Reconstruction
    N = numel(tFit);

    if strcmpi(fitMethod, "onlyharmonics")
        X = zeros(N, 2 * numel(periods));

        for i = 1:numel(periods)
            X(:, (i - 1) * 2 + [1, 2]) = ...
                [cos(2 * pi * tFit / periods(i)), sin(2 * pi * tFit / periods(i))];
        end

        xFit = X * coeffs(p + 2:end);
        xFit(~isfinite(tFit)) = NaN;
        return
    elseif strcmpi(fitMethod, "omitharmonics")
        X = zeros([N, p + 1]);
        X(:, 1:p + 1) = tFit .^ (0:p);
    else
        X = zeros([N, p + 1 + 2 * numel(periods)]);
        X(:, 1:p + 1) = tFit .^ (0:p);

        for i = 1:numel(periods)
            X(:, p + 1 + (i - 1) * 2 + [1, 2]) = ...
                [cos(2 * pi * tFit / periods(i)), sin(2 * pi * tFit / periods(i))];
        end

    end

    xFit = X * coeffs(1:size(X, 2));
end

function mustBeVectorOrEmpty(value)

    if ~isempty(value) && ~isvector(value)
        error('fittimeseries:InvalidShape', 'Input must be a vector or empty.');
    end

end

% Make sure the input time range is valid
function mustBeCompatibleRange(value, t)
    sameType = (isnumeric(t) && isnumeric(value)) ...
        || (isdatetime(t) && isdatetime(value)) ...
        || (isduration(t) && isduration(value));

    if ~sameType || (isnumeric(value) && ~isreal(value)) || any(ismissing(value), 'all')
        error('ULMO:fittimeseries:InvalidRange', ...
        'Ranges must contain real, nonmissing endpoints of the same time type as t.');
    end

    if value(1) > value(2)
        error('ULMO:fittimeseries:InvalidRange', ...
        'Range endpoints must be in ascending order.');
    end

end
