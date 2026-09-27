%% PERIODICTIMESERIES
%
% Created by
%   2025/06/03, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified by
%   2026/09/27, En-Chi Lee (williameclee@arizona.edu)

function [polys, harmons, xFit, polySigmas, harmonSigmas] = ...
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
        options.Reconstruction (1, 1) string ...
            {mustBeMember(options.Reconstruction, ["includeharmonics", "omitharmonics"])} ...
            = "includeharmonics"
        options.FitRange (1, 2) {mustBeCompatibleRange(options.FitRange, t)} = ...
            [min(t(:), [], 'omitmissing'), max(t(:), [], 'omitmissing')]
        options.AverageRange (1, 2) {mustBeCompatibleRange(options.AverageRange, t)} = ...
            [min(t(:), [], 'omitmissing'), max(t(:), [], 'omitmissing')]
    end

    arguments (Output)
        polys (:, :) {mustBeNumeric}
        harmons {mustBeNumeric}
        xFit (:, :) {mustBeNumeric}
        polySigmas (:, :) {mustBeNumeric}
        harmonSigmas {mustBeNumeric}
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

    fitRange = options.FitRange(:);
    averageRange = options.AverageRange(:);

    isTime = isdatetime(t);

    if isTime
        firstValid = find(~isnat(t), 1);

        if isempty(firstValid)
            error('periodictimeseries:InvalidTime', 'At least one valid time is required.');
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
            error('periodictimeseries:InvalidTime', 'Times must be real.');
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
        error('periodictimeseries:InvalidPeriods', 'Periods must be finite, real, and positive.');
    end

    polyTransform = [];

    if options.PolynomialFormat == "average-derivatives"

        if any(~isfinite(averageRange)) || averageRange(2) <= averageRange(1)
            error('ULMO:periodictimeseries:InvalidAverageRange', ...
            'AverageRange must have finite endpoints and positive length.');
        end

        polyTransform = averagederivativematrix(p, averageRange);
    end

    polys = nan(p + 1, nSeries);
    polySigmas = nan(p + 1, nSeries);

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
            error('ULMO:periodictimeseries:InsufficientData', ...
                ['Series %d needs more valid observations than fitted parameters.', ...
             '%d observations is needed, but only got %d valid ones.'], ...
                iSeries, nMinValids, nnz(isValid));
        end

        [polys(:, iSeries), harmons(:, :, iSeries), xFit(:, iSeries), ...
             polySigmas(:, iSeries), harmonSigmas(:, :, iSeries)] = ...
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
        error('ULMO:periodictimeseries:SizeMismatch', ...
            '%s must have a dimension matching the shared time axis (%d).', name, nTimes);
    end

end

function [polys, harmons, xFit, polySigmas, harmonSigmas] = ...
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

    polys = coeffs(1:p + 1);
    polySigmas = coeffSigmas(1:p + 1);

    if ~isempty(polyTransform)
        polys = polyTransform * polys;
        polyCov = polyTransform * cov_matrix(1:p + 1, 1:p + 1) * polyTransform.';
        polyCov = (polyCov + polyCov.') / 2;
        polySigmas = sqrt(max(diag(polyCov), 0));
    end

    harmons = coeffs(p + 2:end);
    harmonSigmas = coeffSigmas(p + 2:end);

    harmons = reshape(harmons, [2, numel(periods)])';
    harmonSigmas = reshape(harmonSigmas, [2, numel(periods)])';

    %% Post-processing coefficients
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

    if strcmpi(fitMethod, "omitharmonics")
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
        error('periodictimeseries:InvalidShape', 'Input must be a vector or empty.');
    end

end

% Make sure the input time range is valid
function mustBeCompatibleRange(value, t)
    sameType = (isnumeric(t) && isnumeric(value)) ...
        || (isdatetime(t) && isdatetime(value)) ...
        || (isduration(t) && isduration(value));

    if ~sameType || (isnumeric(value) && ~isreal(value)) || any(ismissing(value), 'all')
        error('ULMO:periodictimeseries:InvalidRange', ...
        'Ranges must contain real, nonmissing endpoints of the same time type as t.');
    end

    if value(1) > value(2)
        error('ULMO:periodictimeseries:InvalidRange', ...
        'Range endpoints must be in ascending order.');
    end

end

function A = averagederivativematrix(p, interval)
    % A(k+1,j+1) is the interval average of d^k(t^j)/dt^k.
    a = interval(1);
    b = interval(2);
    A = zeros(p + 1);

    for k = 0:p

        for j = k:p
            n = j - k;
            % Divided difference of t^(n+1), evaluated without subtracting
            % nearby endpoint powers or dividing by a small interval length.
            meanPower = sum(a .^ (0:n) .* b .^ (n:-1:0)) / (n + 1);
            A(k + 1, j + 1) = prod((j - k + 1):j) * meanPower;
        end

    end

end
