%% fittimeseriesTest
%
% Last modified by
%   2026/09/28, En-Chi Lee (williameclee@arizona.edu)

classdef fittimeseriesTest < matlab.unittest.TestCase
    % Run with runtests('fittimeseriesTest').
    properties (TestParameter)
        format = {'cos-sin', 'sin-cos', 'amp-phase'}
        transposeTime = {false, true}
        transposeData = {false, true}
        transposeSigma = {false, true}
        uncertaintyMode = {'none', 'column', 'row'}
        timeType = {'numeric', 'duration', 'datetime'}
        harmonicPeriods = struct('two', [1, 2], 'three', [1, 2, 0.7])
        invalidInput = {'dataSize', 'sigmaSize', 'zeroSigma', 'infiniteSigma', ...
                            'zeroPeriod', 'missingData', 'differentTimeAxes', 'repeatedTimeAxis'}
    end

    properties
        T
        X
        Sigma
    end

    methods (TestClassSetup)

        function addUlmoPath(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(fileparts(mfilename('fullpath')))));
        end

    end

    methods (TestMethodSetup)

        function createSeries(testCase)
            t = linspace(0, 6, 73)';
            testCase.T = t;
            testCase.X = [2 + 3 * t + 0.2 * t .^ 2 + cos(2 * pi * t), ...
                              -1 + 0.5 * t - 0.1 * t .^ 2 + 2 * sin(2 * pi * t)] + [0.02 * sin(5 * t), 0.03 * cos(3 * t)];
            testCase.Sigma = [0.1 + 0.01 * t, 0.3 - 0.01 * t];
        end

    end

    methods (Test)

        function averageRangeDefaultsToFitRange(testCase, timeType)
            u = testCase.T;
            x = 2 + 3 * u + 0.2 * u .^ 2;
            convert = @(v) v;

            if strcmp(timeType, 'duration')
                convert = @(v) years(v);
            elseif strcmp(timeType, 'datetime')
                convert = @(v) datetime(2005, 1, 1) + years(v);
            end

            t = convert(u);
            % Asymmetric intervals distinguish their mean derivatives from
            % the input-time mean, including when FitRange exceeds coverage.
            for limits = {[1, 4], [-1, 9]}
                range = limits{1};
                actual = fittimeseries(t, x, [], 2, [], ...
                    FitRange = convert(range), PolynomialFormat = 'average-derivatives');
                testCase.verifyEqual(actual(2), 3 + 0.2 * sum(range), AbsTol = 1e-9);
                explicit = fittimeseries(t, x, [], 2, [], ...
                    FitRange = convert(range), AverageRange = convert([2, 3]), ...
                    PolynomialFormat = 'average-derivatives');
                testCase.verifyEqual(explicit(2), 4, AbsTol = 1e-9);
            end

            default = fittimeseries(t, x, [], 2, [], PolynomialFormat = 'average-derivatives');
            testCase.verifyEqual(default(2), 4.2, AbsTol = 1e-9);
        end

        function onlyHarmonics(testCase, format, transposeData, uncertaintyMode)
            t = testCase.T;
            periodic = [cos(2 * pi * t) + 2 * sin(pi * t), ...
                            3 * sin(2 * pi * t) - cos(pi * t)];
            x = [2 + 3 * t + 0.2 * t .^ 2, -1 + 0.5 * t] + periodic;
            sigma = [];
            if strcmp(uncertaintyMode, 'column'), sigma = testCase.Sigma(:, 1); end
            if strcmp(uncertaintyMode, 'row'), sigma = testCase.Sigma(:, 1).'; end
            x(20, 1) = NaN;
            if transposeData, x = x.'; periodic = periodic.'; end
            args = {sigma, 2, [1, 2], 'PeriodicFormat', format, ...
                        'PolynomialFormat', 'average-derivatives', 'FitRange', [1, 5]};
            [a, b, full, d, e] = fittimeseries(t, x, args{:});
            [aa, bb, fitted, dd, ee] = fittimeseries(t, x, args{:}, ...
                Reconstruction = 'onlyharmonics');
            [~, ~, polynomial] = fittimeseries(t, x, args{:}, ...
                Reconstruction = 'omitharmonics');
            testCase.verifyEqual(fitted, periodic, AbsTol = 1e-9);
            testCase.verifyEqual(full, polynomial + fitted, AbsTol = 1e-9);
            testCase.verifyEqual(aa, a);
            testCase.verifyEqual(bb, b);
            testCase.verifyEqual(dd, d);
            testCase.verifyEqual(ee, e);
        end

        function onlyHarmonicsInvalidTimes(testCase)
            t = testCase.T;
            t([8, 12, 15]) = [NaN, Inf, -Inf];

            for periods = {[], 1}
                [~, ~, fitted] = fittimeseries(t, testCase.X, [], 2, periods{1}, ...
                    Reconstruction = "onlyharmonics");
                testCase.verifyTrue(all(isnan(fitted(~isfinite(t), :)), 'all'));
                testCase.verifyTrue(all(isfinite(fitted(isfinite(t), :)), 'all'));

                if isempty(periods{1})
                    testCase.verifyEqual(fitted(isfinite(t), :), zeros(nnz(isfinite(t)), 2));
                end

            end

        end

        function multiplePeriodSinCos(testCase, harmonicPeriods, uncertaintyMode)
            t = testCase.T;
            % Build the reference model directly in sine/cosine order.
            G = [ones(size(t)), t];

            for period = harmonicPeriods
                G = [G, sin(2 * pi * t / period), cos(2 * pi * t / period)]; %#ok<AGROW>
            end

            beta = [(1:size(G, 2))', - (size(G, 2):-1:1)'];
            x = G * beta + [0.02 * sin(5 * t), 0.03 * cos(3 * t)];
            sigma = [];
            if strcmp(uncertaintyMode, 'column'), sigma = testCase.Sigma(:, 1); end
            if strcmp(uncertaintyMode, 'row'), sigma = testCase.Sigma(:, 1).'; end
            [~, harmonics, fitted, ~, errors] = fittimeseries( ...
                t, x, sigma, 1, harmonicPeriods, PeriodicFormat = 'sin-cos');

            % Independent QR reference checks period order and error pairing.
            for j = 1:size(x, 2)
                design = G;
                data = x(:, j);

                if ~isempty(sigma)
                    design = design ./ sigma(:);
                    data = data ./ sigma(:);
                end

                [Q, R] = qr(design, 0);
                expected = R \ (Q' * data);
                inverseR = R \ eye(size(R));
                variance = sum((data - design * expected) .^ 2) / ...
                    (size(design, 1) - size(design, 2));
                se = sqrt(diag(variance * (inverseR * inverseR.')));
                testCase.verifyEqual(harmonics(:, :, j), ...
                    reshape(expected(3:end), 2, []).', AbsTol = 1e-9);
                testCase.verifyEqual(errors(:, :, j), ...
                    reshape(se(3:end), 2, []).', AbsTol = 1e-9);
                testCase.verifyEqual(fitted(:, j), G * expected, AbsTol = 1e-9);
            end

        end

        function matrixOrientations(testCase, format, transposeTime, transposeData, transposeSigma)
            t = testCase.T; x = testCase.X; sigma = testCase.Sigma;
            if transposeTime, t = t.'; end
            if transposeData, x = x.'; end
            if transposeSigma, sigma = sigma.'; end
            [poly, periodic, fitted, polySigma, periodicSigma] = ...
                fittimeseries(t, x, sigma, 2, 1, PeriodicFormat = format);
            testCase.verifySize(fitted, size(x));
            testCase.verifySize(periodic, [1, 2, 2]);
            testCase.verifySize(periodicSigma, [1, 2, 2]);
            if transposeData, fitted = fitted.'; end

            for j = 1:2
                [a, b, c, d, e] = fittimeseries(testCase.T, testCase.X(:, j), ...
                    testCase.Sigma(:, j), 2, 1, PeriodicFormat = format);
                testCase.verifyEqual(poly(:, j), a, AbsTol = 1e-9, RelTol = 1e-9);
                testCase.verifyEqual(periodic(:, :, j), b, AbsTol = 1e-9, RelTol = 1e-9);
                testCase.verifyEqual(fitted(:, j), c, AbsTol = 1e-9, RelTol = 1e-9);
                testCase.verifyEqual(polySigma(:, j), d, AbsTol = 1e-9, RelTol = 1e-9);
                testCase.verifyEqual(periodicSigma(:, :, j), e, AbsTol = 1e-9, RelTol = 1e-9);
            end

        end

        function sharedUncertainties(testCase, uncertaintyMode)
            t = testCase.T; x = testCase.X;
            sigma = [];
            if strcmp(uncertaintyMode, 'column'), sigma = testCase.Sigma(:, 1); end
            if strcmp(uncertaintyMode, 'row'), sigma = testCase.Sigma(:, 1).'; end
            [a, b, c, d, e] = fittimeseries(t, x, sigma, 2, 1);
            [aa, bb, cc, dd, ee] = fittimeseries(t.', x.', sigma, 2, 1);
            testCase.verifyEqual(a, aa, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(b, bb, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, cc.', AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(d, dd, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(e, ee, AbsTol = 1e-9, RelTol = 1e-9);
            % Shared weights must also agree with independent scalar fits.
            for j = 1:2
                expected = fittimeseries(t, x(:, j), sigma, 2, 1);
                testCase.verifyEqual(a(:, j), expected, AbsTol = 1e-9, RelTol = 1e-9);
            end

        end

        function singleRowSeries(testCase)
            [a, ~, c] = fittimeseries(testCase.T.', testCase.X(:, 1).', [], 2, 1);
            [b, ~, d] = fittimeseries(testCase.T, testCase.X(:, 1), [], 2, 1);
            testCase.verifySize(c, [1, numel(testCase.T)]);
            testCase.verifyEqual(a, b, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, d.', AbsTol = 1e-9, RelTol = 1e-9);
        end

        function missingObservationsAndReconstruction(testCase)
            t = testCase.T; x = testCase.X; sigma = testCase.Sigma;
            t(8) = NaN; t(12) = Inf;
            x(3, 1) = NaN; x(5, 2) = Inf; sigma(10, 2) = NaN;
            [a, h, fitted, d] = fittimeseries(t, x, sigma, 2, 1);
            valid = isfinite(t) & isfinite(x) & isfinite(sigma);
            % Missing observations/weights are excluded from fitting, but
            % predictions remain defined wherever the original time is finite.
            testCase.verifyEqual(isfinite(fitted), repmat(isfinite(t), 1, 2));

            for j = 1:2
                keep = valid(:, j);
                [aa, ~, cc, dd] = fittimeseries(t(keep), x(keep, j), sigma(keep, j), 2, 1);
                testCase.verifyEqual(a(:, j), aa, AbsTol = 1e-9, RelTol = 1e-9);
                testCase.verifyEqual(fitted(keep, j), cc, AbsTol = 1e-9, RelTol = 1e-9);
                testCase.verifyEqual(d(:, j), dd, AbsTol = 1e-9, RelTol = 1e-9);
                times = t(isfinite(t));
                expected = (times .^ (0:2)) * aa + cos(2 * pi * times) * h(1, 1, j) + sin(2 * pi * times) * h(1, 2, j);
                testCase.verifyEqual(fitted(isfinite(t), j), expected, AbsTol = 1e-9, RelTol = 1e-9);
            end

            [~, ~, transposed] = fittimeseries(t.', x.', sigma.', 2, 1);
            testCase.verifyEqual(fitted, transposed.', AbsTol = 1e-9, RelTol = 1e-9);
        end

        function timeTypesAndOrientation(testCase, timeType)
            t = testCase.T; periods = 1;

            if strcmp(timeType, 'duration')
                t = years(t); periods = years(1);
            elseif strcmp(timeType, 'datetime')
                t = datetime(2005, 1, 1) + years(t); periods = years(1);
            end

            [a, ~, c] = fittimeseries(t, testCase.X, [], 2, periods);
            [b, ~, d] = fittimeseries(t.', testCase.X.', [], 2, periods);
            [expected, ~, expectedFit] = fittimeseries(testCase.T, testCase.X, [], 2, 1);
            testCase.verifyEqual(a, b, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, d.', AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(a, expected, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, expectedFit, AbsTol = 1e-9, RelTol = 1e-9);
        end

        function missingFirstDatetime(testCase)
            dates = datetime(2005, 1, 1) + years(testCase.T); dates(1) = NaT;
            [~, ~, fitted] = fittimeseries(dates, testCase.X);
            testCase.verifyTrue(all(isnan(fitted(1, :))));
            testCase.verifyTrue(all(isfinite(fitted(2:end, :)), 'all'));
        end

        function emptyHarmonicsAndExactPolynomial(testCase, format)
            t = testCase.T; exact = [1 + 2 * t, 3 - 4 * t];
            [a, b, c, d, e] = fittimeseries(t, exact, [], 1, [], PeriodicFormat = format);
            testCase.verifyEqual(a, [1 3; 2 -4], AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, exact, AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifySize(b, [0 2 2]); testCase.verifySize(e, [0 2 2]);
            testCase.verifyLessThan(d, 1e-10);
        end

        function descendingTimes(testCase)
            t = flipud(testCase.T); exact = [1 + 2 * t, 3 - 4 * t];
            [a, ~, c] = fittimeseries(t, exact, [], 1);
            testCase.verifyEqual(a, [1 3; 2 -4], AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, exact, AbsTol = 1e-9, RelTol = 1e-9);
        end

        function squareDataUsesRowsAsObservations(testCase)
            t = (1:4)'; x = t * (1:4);
            [a, ~, c] = fittimeseries(t, x, [], 1);
            testCase.verifyEqual(a, [zeros(1, 4); 1:4], AbsTol = 1e-9, RelTol = 1e-9);
            testCase.verifyEqual(c, x, AbsTol = 1e-9, RelTol = 1e-9);
        end

        function invalidInputs(testCase, invalidInput)
            t = testCase.T; x = testCase.X;

            switch invalidInput
                case 'dataSize'
                    callback = @() fittimeseries(t, zeros(5, 2));
                    id = 'ULMO:fittimeseries:SizeMismatch';
                case 'sigmaSize'
                    callback = @() fittimeseries(t, x, ones(numel(t), 3));
                    id = 'ULMO:fittimeseries:SizeMismatch';
                case 'zeroSigma'
                    callback = @() fittimeseries(t, x, zeros(numel(t), 1));
                    id = 'ULMO:fittimeseries:InvalidSigma';
                case 'infiniteSigma'
                    callback = @() fittimeseries(t, x, inf(numel(t), 1));
                    id = 'ULMO:fittimeseries:InvalidSigma';
                case 'zeroPeriod'
                    callback = @() fittimeseries(t, x, [], 2, 0);
                    id = 'fittimeseries:InvalidPeriods';
                case 'missingData'
                    callback = @() fittimeseries(t, nan(size(x)));
                    id = 'ULMO:fittimeseries:InsufficientData';
                case 'differentTimeAxes'
                    callback = @() fittimeseries([t, t + 1], x);
                    id = 'MATLAB:validators:mustBeVector';
                case 'repeatedTimeAxis'
                    callback = @() fittimeseries(repmat(t, 1, 3), x);
                    id = 'MATLAB:validators:mustBeVector';
            end

            testCase.verifyError(callback, id);
        end

    end

end
