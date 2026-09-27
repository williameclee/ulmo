classdef periodicAverageDerivativesTest < matlab.unittest.TestCase

    properties (TestParameter)
        degree = num2cell(0:5)
        weighted = struct('unweighted', false, 'weighted', true)
        timeType = {'duration', 'datetime'}
        invalidRange = struct('zeroLength', [2, 2], 'infiniteEndpoint', [1, Inf])
    end

    methods (TestClassSetup)

        function addUlmoPath(testCase)
            ulmoPath = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(ulmoPath));
        end

    end

    methods (Test)

        function quadraticExample(testCase)
            t = linspace(0, 4, 101)';
            y = 1 + 2 * t + 3 * t .^ 2;
            [averages, ~, fitted] = fittimeseries(t, y, [], 2, [], ...
                PolynomialFormat = "average-derivatives", AverageRange = [1, 3]);

            testCase.verifyEqual(averages, [18; 14; 6], AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(fitted, y, AbsTol = 1e-7, RelTol = 1e-7);
        end

        function polynomialDegrees(testCase, degree)
            t = linspace(0, 4, 101)';
            interval = [1, 3];
            beta = (1:degree + 1)' / 10;
            y = (t .^ (0:degree)) * beta;
            [averages, ~, fitted] = fittimeseries(t, y, [], degree, [], ...
                PolynomialFormat = "average-derivatives", AverageRange = interval);

            % Integrate differentiated polynomials independently of the implementation.
            expected = testCase.integralTransform(degree, interval) * beta;
            testCase.verifyEqual(averages, expected, AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(fitted, y, AbsTol = 1e-7, RelTol = 1e-7);
        end

        function covarianceAndOrientation(testCase, weighted)
            t = linspace(0, 4, 101)';
            interval = [1, 3];
            p = 3;
            periods = [1, .6, 1.7];
            X = t .^ (0:p);

            for period = periods
                X = [X, cos(2 * pi * t / period), sin(2 * pi * t / period)]; %#ok<AGROW>
            end

            y = [X * (1:size(X, 2))' / 10, X * (size(X, 2):-1:1)' / 20];
            y = y + [0.01 * sin(11 * t), 0.02 * cos(13 * t)];
            y(30, 2) = NaN;
            A = testCase.integralTransform(p, interval);
            sigma = [];

            if weighted
                sigma = [1 + .1 * t, 2 + .2 * t];
            end

            [raw, harmonics, rawFit] = fittimeseries(t, y, sigma, p, periods, ...
                FitRange = [.4, 3.6]);
            [averages, averagedHarmonics, fitted, se] = fittimeseries(t, y, sigma, p, periods, ...
                FitRange = [.4, 3.6], PolynomialFormat = "average-derivatives", AverageRange = interval);
            testCase.verifyEqual(averages, A * raw, AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(averagedHarmonics, harmonics, AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(fitted, rawFit, AbsTol = 1e-7, RelTol = 1e-7);

            % Reference covariance from a QR solve, including different valid
            % observations and weights for each series.
            for j = 1:2
                valid = t >= .4 & t <= 3.6 & isfinite(y(:, j));
                G = X(valid, :);
                d = y(valid, j);

                if weighted
                    G = G ./ sigma(valid, j);
                    d = d ./ sigma(valid, j);
                end

                [Q, R] = qr(G, 0);
                b = R \ (Q' * d);
                inverseR = R \ eye(size(R));
                C = (sum((d - G * b) .^ 2) / (size(G, 1) - size(G, 2))) * (inverseR * inverseR.');
                expectedSE = sqrt(diag(A * C(1:p + 1, 1:p + 1) * A.'));
                testCase.verifyEqual(se(:, j), expectedSE, AbsTol = 1e-7, RelTol = 1e-7);
            end

            [transposedAverages, ~, transposedFit] = fittimeseries(t.', y.', sigma.', p, periods, ...
                FitRange = [.4, 3.6], PolynomialFormat = "average-derivatives", AverageRange = interval);
            testCase.verifyEqual(transposedAverages, averages, AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(transposedFit, fitted.', AbsTol = 1e-7, RelTol = 1e-7);
        end

        function defaultAverageRangeIsIndependentOfFitRange(testCase)
            t = linspace(0, 4, 101)';
            averages = fittimeseries(t, 1 + 2 * t + 3 * t .^ 2, [], 2, [], ...
                FitRange = [1, 2], PolynomialFormat = "average-derivatives");
            testCase.verifyEqual(averages, [21; 14; 6], AbsTol = 1e-7, RelTol = 1e-7);
        end

        function timeUnits(testCase, timeType)
            t = linspace(0, 4, 101)';
            times = years(t);
            range = years([1, 3]);
            fitRange = years([.4, 3.6]);

            if strcmp(timeType, 'datetime')
                origin = datetime(2000, 1, 1);
                times = origin + times;
                range = origin + range;
                fitRange = origin + fitRange;
            end

            y = 1 + 2 * t + 3 * t .^ 2;
            [averages, ~, fitted] = fittimeseries(times, y, [], 2, [], FitRange = fitRange, ...
                PolynomialFormat = "average-derivatives", AverageRange = range, Reconstruction = "omitharmonics");
            testCase.verifyEqual(averages, [18; 14; 6], AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(fitted, y, AbsTol = 1e-7, RelTol = 1e-7);
        end

        function averageRangeChangesSummaryOnly(testCase)
            t = linspace(0, 4, 101)';
            y = 1 + 2 * t + 3 * t .^ 2;
            [averages, ~, fitted] = fittimeseries(t, y, [], 2, [], ...
                PolynomialFormat = "average-derivatives", AverageRange = [0, 1]);
            testCase.verifyEqual(averages, [3; 5; 6], AbsTol = 1e-7, RelTol = 1e-7);
            testCase.verifyEqual(fitted, y, AbsTol = 1e-7, RelTol = 1e-7);
        end

        function rejectsInvalidAverageRange(testCase, invalidRange)
            t = linspace(0, 4, 101)';
            testCase.verifyError(@() fittimeseries(t, t, [], 1, [], ...
                PolynomialFormat = "average-derivatives", AverageRange = invalidRange), ...
            'ULMO:fittimeseries:InvalidAverageRange');
        end

    end

    methods (Static, Access = private)

        function A = integralTransform(p, range)
            A = zeros(p + 1);

            for j = 0:p
                c = [1, zeros(1, j)];

                for k = 0:p
                    A(k + 1, j + 1) = integral(@(t) polyval(c, t), range(1), range(2)) / (range(2) - range(1));
                    c = polyder(c);
                end

            end

        end

    end

end
