classdef xyzs2slepTest < matlab.unittest.TestCase

    properties
        Polygon = [340 -25; 380 -20; 375 30; 345 20; 340 -25]
        Lon
        Lat
        Basis
        Eigenvalues
        ShannonNumber
    end

    methods (TestClassSetup)

        function setup(testCase)
            fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fixtures));
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() rmdir(root, 's'));
            original = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', original));
            setenv('IFILES', root);

            for folder = ["GLMALPHA", "KERNELC", "KERNELCP", "SDWCAP", "GLMALPHAPTO", "LEGENDRE"]
                mkdir(fullfile(root, folder));
            end

            % Bravo generates unused default points even with explicit inputs.
            % Give that path its radius constant without reading production data.
            constants = fullfile(root, 'EARTHMODELS', 'CONSTANTS');
            mkdir(constants);
            Earth.Radius = 6371000;
            save(fullfile(constants, 'Earth.mat'), 'Earth');
            originalRng = rng;
            testCase.addTeardown(@() rng(originalRng));
            % Deterministic, non-gridded points including the longitude seam
            % and locations outside the concentration region.
            k = (1:80)';
            testCase.Lon = mod(round(k * 137.507764), 360);
            testCase.Lat = round(asind(-0.95 + 1.9 * (k - 0.5) / 80));
            [testCase.Basis, testCase.Eigenvalues, ~, ~, testCase.ShannonNumber] = ...
                glmalpha_new(testCase.Polygon, 3, 'BeQuiet', true);
        end

    end

    methods (Test)

        function synthesisRecovery(testCase)
            J = 6;
            c = (1:J)' / 7;
            data = testCase.synthesize(testCase.Basis(:, 1:J) * c);
            [actual, V, N, j, n, used, rankA, s, condition, fitted, residual, dof] = ...
                xyzs2slep(data, testCase.Lon, testCase.Lat, testCase.Polygon, 3, truncation = J);
            testCase.verifyEqual(actual, c, 'AbsTol', 2e-11);
            testCase.verifyEqual(V, testCase.Eigenvalues(1:J), 'AbsTol', 1e-12);
            testCase.verifyEqual(N, testCase.ShannonNumber);
            testCase.verifyEqual([j n rankA dof], [J 80 J 80 - J]);
            testCase.verifyEqual(used, (1:80)');
            testCase.verifySize(s, [J 1]);
            testCase.verifyEqual(condition, s(1) / s(end));
            testCase.verifyEqual(fitted, data, 'AbsTol', 2e-11);
            testCase.verifyEqual(residual, zeros(80, 1), 'AbsTol', 2e-11);
            % Existing coefficient projection independently checks ordering.
            lmcosi = testCase.toPlm(testCase.Basis(:, 1:J) * c);
            inverse = slep2plm_new(c, testCase.Polygon, 3, 0, 0, 0, J);
            testCase.verifyEqual(inverse, lmcosi, 'AbsTol', 2e-12);
            projected = plm2slep_new(lmcosi, testCase.Polygon, 3, 'BeQuiet', true);
            testCase.verifyEqual(actual, projected(1:J), 'AbsTol', 2e-11);
        end

        function noisyLeastSquaresAndBlocks(testCase)
            J = 5;
            A = zeros(80, J);

            for j = 1:J
                A(:, j) = testCase.synthesize(testCase.Basis(:, j));
            end

            data = A * (1:J)' + 0.02 * cos((1:80)');
            expected = A \ data;
            actual = xyzs2slep(data', testCase.Lon', testCase.Lat', testCase.Polygon, 3, ...
                truncation = J, blockSize = 7);
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-11);
            other = xyzs2slep(data, testCase.Lon + 720, testCase.Lat, testCase.Polygon, 3, ...
                truncation = J, blockSize = 200);
            testCase.verifyEqual(actual, other, 'AbsTol', 2e-11);
        end

        function analyticOffGridHarmonics(testCase)
            lon = testCase.Lon + 0.12345;
            lat = testCase.Lat + 0.23456;
            data = 2 + sqrt(3) * (0.3 * sind(lat) + 0.7 * cosd(lat) .* cosd(lon) ...
                - 0.4 * cosd(lat) .* sind(lon));
            lmcosi = testCase.toPlm(zeros(16, 1));
            lmcosi(1, 3) = 2;
            lmcosi(2, 3) = 0.3;
            lmcosi(3, 3:4) = [0.7 -0.4];
            expected = plm2slep_new(lmcosi, testCase.Polygon, 3, 'BeQuiet', true);
            actual = xyzs2slep(data, lon, lat, testCase.Polygon, 3, truncation = 16, blockSize = 1);
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-11);
        end

        function defaultTruncation(testCase)
            [~, ~, N, J] = xyzs2slep(ones(80, 1), testCase.Lon, testCase.Lat, testCase.Polygon, 3);
            testCase.verifyEqual(J, max(1, round(N)));
        end

        function missingAndInvalidInputs(testCase)
            data = ones(80, 1); data(3) = NaN;
            lon = testCase.Lon; lon(5) = NaN;
            testCase.verifyError(@() xyzs2slep(data, lon, testCase.Lat, testCase.Polygon, 3), ...
            'ULMO:xyzs2slep:MissingObservation');
            [~, ~, ~, ~, n, used] = xyzs2slep(data, lon, testCase.Lat, testCase.Polygon, 3, missingPolicy = "omit");
            testCase.verifyEqual(n, 78);
            testCase.verifyEqual(used, setdiff((1:80)', [3; 5]));
            testCase.verifyError(@() xyzs2slep([1 2], 1, [0 0], 30, 1), 'ULMO:xyzs2slep:ObservationSize');
            testCase.verifyError(@() xyzs2slep(1, 0, 91, 30, 1), 'ULMO:xyzs2slep:InvalidObservation');
            testCase.verifyError(@() xyzs2slep(Inf, 0, 0, 30, 1, missingPolicy = "omit"), 'ULMO:xyzs2slep:InvalidObservation');
            testCase.verifyError(@() xyzs2slep(NaN, 0, 0, 30, 1, missingPolicy = "omit"), 'ULMO:xyzs2slep:ObservationSize');
            testCase.verifyError(@() xyzs2slep([1 2], [0 1], [0 1], 30, 1, truncation = 3), 'ULMO:xyzs2slep:InvalidTruncation');
        end

        function rankAndExactlyDetermined(testCase)
            testCase.verifyError(@() xyzs2slep(ones(10, 1), zeros(10, 1), zeros(10, 1), ...
                testCase.Polygon, 3, truncation = 3), 'ULMO:xyzs2slep:RankDeficient');
            J = 4;
            data = testCase.synthesize(testCase.Basis(:, 1:J) * (1:J)');
            [actual, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, dof] = xyzs2slep(data(1:J), ...
                testCase.Lon(1:J), testCase.Lat(1:J), testCase.Polygon, 3, truncation = J);
            testCase.verifyEqual(actual, (1:J)', 'AbsTol', 2e-10);
            testCase.verifyEqual(dof, 0);
            testCase.verifyError(@() xyzs2slep(data, testCase.Lon, testCase.Lat, ...
                testCase.Polygon, 3, truncation = J, rankTolerance = 0.999), 'ULMO:xyzs2slep:RankDeficient');
        end

        function geographicDomainTopology(testCase)
            domain = GeoDomain('xyzs2slepregionfixture');
            xy = domain.Lonlat;
            testCase.verifyEqual(sum(isnan(xy(:, 1))), 2);
            [G, ~] = glmalpha_new(domain, 3, 'BeQuiet', true);
            c = (1:5)' / 3;
            data = testCase.synthesize(G(:, 1:5) * c);
            actual = xyzs2slep(data, testCase.Lon, testCase.Lat, domain, 3, truncation = 5);
            testCase.verifyEqual(actual, c, 'AbsTol', 2e-11);
            byName = xyzs2slep(data, testCase.Lon, testCase.Lat, 'xyzs2slepregionfixture', 3, truncation = 5);
            byCell = xyzs2slep(data, testCase.Lon, testCase.Lat, {'xyzs2slepregionfixture', 0}, 3, truncation = 5);
            testCase.verifyEqual(byName, actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(byCell, actual, 'AbsTol', 2e-11);
            % Verify that the actual kernel excludes the island and includes
            % the disconnected component, rather than only a bounding box.
            [~, ~, ~, ~, actualN] = glmalpha_new(domain, 3, 'BeQuiet', true);
            [~, ~, ~, ~, outerN] = glmalpha_new(xy(1:5, :), 3, 'BeQuiet', true);
            [~, ~, ~, ~, holeN] = glmalpha_new(xy(7:11, :), 3, 'BeQuiet', true);
            [~, ~, ~, ~, islandN] = glmalpha_new(xy(13:17, :), 3, 'BeQuiet', true);
            testCase.verifyEqual(actualN, outerN - holeN + islandN, 'AbsTol', 2e-3);
        end

        function uncertaintyAbsentPreservesOutputs(testCase)
            data = cos((1:80)');
            original = cell(1, 12);
            extended = cell(1, 15);
            [original{:}] = xyzs2slep(data, testCase.Lon, testCase.Lat, testCase.Polygon, 3, truncation = 4);
            [extended{:}] = xyzs2slep(data, testCase.Lon, testCase.Lat, testCase.Polygon, 3, truncation = 4);
            testCase.verifyEqual(extended(1:12), original);
            testCase.verifyEmpty(extended{13});
            testCase.verifyEmpty(extended{14});
            testCase.verifyEqual(extended{15}, "none");
        end

        function standardDeviationAndDiagonalCovariance(testCase)
            J = 4;
            A = testCase.designMatrix(J);
            data = A * (1:J)' + 0.1 * sin((1:80)');
            sigma = linspace(0.2, 0.8, 80)';
            [Q, R] = qr(A ./ sigma, 0);
            expected = R \ (Q.' * (data ./ sigma));
            inverseR = R \ eye(J);
            expectedCov = inverseR * inverseR.';
            [actual, covariance, se, source, fitted, residual] = ...
                testCase.fitUncertainty(data, J, dataStd = sigma');
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-11);
            testCase.verifyEqual(covariance, expectedCov, 'AbsTol', 2e-12);
            testCase.verifyEqual(se, sqrt(diag(covariance)), 'AbsTol', 1e-14);
            testCase.verifyEqual(source, "dataStd");
            [scaledFit, scaledCov] = testCase.fitUncertainty(data, J, dataStd = 3 * sigma);
            testCase.verifyEqual(scaledFit, actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(scaledCov, 9 * covariance, 'AbsTol', 2e-12);
            [scaledFit, scaledCov] = testCase.fitUncertainty(3 * data, J, dataStd = 3 * sigma);
            testCase.verifyEqual(scaledFit, 3 * actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(scaledCov, 9 * covariance, 'AbsTol', 2e-12);
            testCase.verifyEqual(fitted, A * actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(residual, data - fitted, 'AbsTol', 1e-14);
            [diagonalFit, diagonalCov, ~, diagonalSource] = ...
                testCase.fitUncertainty(data, J, dataCovariance = diag(sigma .^ 2));
            testCase.verifyEqual(diagonalFit, actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(diagonalCov, covariance, 'AbsTol', 2e-12);
            testCase.verifyEqual(diagonalSource, "dataCovariance");
            uniform = testCase.fitUncertainty(data, J, dataStd = 0.4 * ones(80, 1));
            testCase.verifyEqual(uniform, A \ data, 'AbsTol', 2e-11);
        end

        function correlatedCovarianceAndScaling(testCase)
            J = 5;
            A = testCase.designMatrix(J);
            data = A * (1:J)' + sin((1:80)');
            covarianceData = toeplitz(0.6 .^ (0:79));
            factor = chol(covarianceData, 'lower');
            B = factor \ A;
            [Q, R] = qr(B, 0);
            expected = R \ (Q.' * (factor \ data));
            inverseR = R \ eye(J);
            [actual, covariance] = testCase.fitUncertainty(data, J, dataCovariance = covarianceData);
            testCase.verifyEqual(actual, expected, 'AbsTol', 2e-11);
            testCase.verifyEqual(covariance, inverseR * inverseR.', 'AbsTol', 2e-12);
            [~, ~, ~, ~, ~, ~, rankA, s, condition] = xyzs2slep(data, testCase.Lon, ...
                testCase.Lat, testCase.Polygon, 3, truncation = J, dataCovariance = covarianceData);
            testCase.verifyEqual(rankA, J);
            testCase.verifyEqual(s, svd(B), 'AbsTol', 2e-11);
            testCase.verifyEqual(condition, cond(B), 'AbsTol', 2e-11);
            [scaledFit, scaledCov] = testCase.fitUncertainty(data, J, dataCovariance = 9 * covarianceData);
            testCase.verifyEqual(scaledFit, actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(scaledCov, 9 * covariance, 'AbsTol', 2e-12);
            [scaledFit, scaledCov] = testCase.fitUncertainty(3 * data, J, dataCovariance = 9 * covarianceData);
            testCase.verifyEqual(scaledFit, 3 * actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(scaledCov, 9 * covariance, 'AbsTol', 2e-12);
            % Roundoff-level asymmetry is explicitly averaged, not treated
            % as a request to use just one triangular half.
            nearlySymmetric = covarianceData;
            nearlySymmetric(1, 2) = nearlySymmetric(1, 2) + eps;
            [roundoffFit, roundoffCov] = testCase.fitUncertainty(data, J, dataCovariance = nearlySymmetric);
            testCase.verifyEqual(roundoffFit, actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(roundoffCov, covariance, 'AbsTol', 2e-12);
        end

        function residualVarianceIsExplicit(testCase)
            J = 3;
            A = testCase.designMatrix(J);
            data = A * (1:J)' + 0.2 * cos((1:80)');
            [actual, covariance, ~, source, ~, residual] = ...
                testCase.fitUncertainty(data, J, estimateNoiseVariance = true);
            [~, R] = qr(A, 0);
            inverseR = R \ eye(J);
            expected = sum(residual .^ 2) / (80 - J) * (inverseR * inverseR.');
            testCase.verifyEqual(actual, A \ data, 'AbsTol', 2e-11);
            testCase.verifyEqual(covariance, expected, 'AbsTol', 2e-12);
            testCase.verifyEqual(source, "estimatedIid");
            [~, knownCov] = testCase.fitUncertainty(data, J, dataStd = ones(80, 1));
            [~, differentResidualCov] = testCase.fitUncertainty(20 * data, J, dataStd = ones(80, 1));
            testCase.verifyEqual(knownCov, differentResidualCov);
            testCase.verifyEqual(knownCov, inverseR * inverseR.', 'AbsTol', 2e-12);
        end

        function correlatedMissingSubset(testCase)
            J = 3;
            A = testCase.designMatrix(J);
            data = A * (1:J)' + 0.1 * sin((1:80)');
            data([2 8]) = NaN;
            C = toeplitz(0.7 .^ (0:79));
            kept = setdiff((1:80)', [2; 8]);
            [actual, covariance] = testCase.fitUncertainty(data, J, ...
                dataCovariance = C, missingPolicy = "omit");
            factor = chol(C(kept, kept), 'lower');
            [Q, R] = qr(factor \ A(kept, :), 0);
            inverseR = R \ eye(J);
            testCase.verifyEqual(actual, R \ (Q.' * (factor \ data(kept))), 'AbsTol', 2e-11);
            testCase.verifyEqual(covariance, inverseR * inverseR.', 'AbsTol', 2e-12);
            sigma = linspace(1, 2, 80)';
            [actual, covariance] = testCase.fitUncertainty(data, J, ...
                dataStd = sigma, missingPolicy = "omit");
            [Q, R] = qr(A(kept, :) ./ sigma(kept), 0);
            inverseR = R \ eye(J);
            testCase.verifyEqual(actual, R \ (Q.' * (data(kept) ./ sigma(kept))), 'AbsTol', 2e-11);
            testCase.verifyEqual(covariance, inverseR * inverseR.', 'AbsTol', 2e-12);
        end

        function analyticConstantFieldUncertainty(testCase)
            % Degree zero is a constant, so its weighted-mean covariance can
            % be checked by hand without another matrix factorization.
            [c, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, covariance] = ...
                xyzs2slep([2 4], [0 50], [-20 10], 60, 0, dataStd = [1 2]);
            testCase.verifyEqual(c, 2.4, 'AbsTol', 1e-13);
            testCase.verifyEqual(covariance, 0.8, 'AbsTol', 1e-13);
            [c, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, covariance] = ...
                xyzs2slep([2 4], [0 50], [-20 10], 60, 0, dataCovariance = [1 .5; .5 4]);
            testCase.verifyEqual(c, 2.25, 'AbsTol', 1e-13);
            testCase.verifyEqual(covariance, 0.9375, 'AbsTol', 1e-13);
        end

        function uncertaintyValidation(testCase)
            base = {ones(80, 1), testCase.Lon, testCase.Lat, testCase.Polygon, 3, 'truncation', 3};
            testCase.verifyError(@() xyzs2slep(base{:}, dataStd = ones(80, 1), dataCovariance = eye(80)), ...
            'ULMO:xyzs2slep:ConflictingUncertainty');
            testCase.verifyError(@() xyzs2slep(base{:}, dataStd = ones(80, 1), estimateNoiseVariance = true), ...
            'ULMO:xyzs2slep:ConflictingUncertainty');
            testCase.verifyError(@() xyzs2slep(base{:}, dataCovariance = eye(80), estimateNoiseVariance = true), ...
            'ULMO:xyzs2slep:ConflictingUncertainty');

            for bad = {1, zeros(80, 1), -ones(80, 1), nan(80, 1), inf(80, 1), ones(8, 10)}
                testCase.verifyError(@() xyzs2slep(base{:}, dataStd = bad{1}), ...
                'ULMO:xyzs2slep:InvalidDataStd');
            end

            for bad = {eye(79), nan(80), inf(80), ones(80, 80, 2)}
                testCase.verifyError(@() xyzs2slep(base{:}, dataCovariance = bad{1}), ...
                'ULMO:xyzs2slep:InvalidDataCovariance');
            end

            for bad = {zeros(80), ones(80), -eye(80)}
                testCase.verifyError(@() xyzs2slep(base{:}, dataCovariance = bad{1}), ...
                'ULMO:xyzs2slep:NonPositiveCovariance');
            end

            asymmetric = eye(80); asymmetric(1, 2) = 0.01;
            testCase.verifyError(@() xyzs2slep(base{:}, dataCovariance = asymmetric), ...
            'ULMO:xyzs2slep:AsymmetricCovariance');
            % Invalid errors remain invalid even on rows omitted from data.
            base{1}(1) = NaN;
            badStd = ones(80, 1); badStd(1) = 0;
            testCase.verifyError(@() xyzs2slep(base{:}, missingPolicy = "omit", dataStd = badStd), ...
            'ULMO:xyzs2slep:InvalidDataStd');
            badCov = eye(80); badCov(1, 1) = NaN;
            testCase.verifyError(@() xyzs2slep(base{:}, missingPolicy = "omit", dataCovariance = badCov), ...
            'ULMO:xyzs2slep:InvalidDataCovariance');
            badCov(1, 1) = -1;
            testCase.verifyError(@() xyzs2slep(base{:}, missingPolicy = "omit", dataCovariance = badCov), ...
            'ULMO:xyzs2slep:NonPositiveCovariance');
        end

        function exactlyDeterminedUncertainty(testCase)
            J = 4;
            A = testCase.designMatrix(J);
            A = A(1:J, :);
            data = A * (1:J)';
            base = {data, testCase.Lon(1:J), testCase.Lat(1:J), testCase.Polygon, 3, 'truncation', J};
            [actual, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, dof, covariance] = ...
                xyzs2slep(base{:}, dataStd = ones(J, 1));
            inverseA = A \ eye(J);
            testCase.verifyEqual(actual, (1:J)', 'AbsTol', 2e-10);
            testCase.verifyEqual(covariance, inverseA * inverseA.', 'AbsTol', 2e-10);
            testCase.verifyEqual(dof, 0);
            testCase.verifyError(@() xyzs2slep(base{:}, estimateNoiseVariance = true), ...
            'ULMO:xyzs2slep:NoiseDegreesOfFreedom');
        end

        function monteCarloCorrelatedPropagation(testCase)
            rng(7123);
            J = 3;
            A = testCase.designMatrix(J);
            C = 0.04 * toeplitz(0.5 .^ (0:79));
            factor = chol(C, 'lower');
            truth = (1:J)' / 10;
            [~, predicted] = testCase.fitUncertainty(A * truth, J, dataCovariance = C);
            draws = 600;
            noise = factor * randn(80, draws);
            recovered = zeros(J, draws);

            for draw = 1:draws
                recovered(:, draw) = xyzs2slep(A * truth + noise(:, draw), ...
                    testCase.Lon, testCase.Lat, testCase.Polygon, 3, ...
                    truncation = J, dataCovariance = C);
            end

            empirical = cov(recovered.');
            % Gaussian sample covariance has this elementwise variance.
            tolerance = 6 * sqrt((predicted .^ 2 + diag(predicted) * diag(predicted).') / (draws - 1));
            testCase.verifyLessThanOrEqual(abs(empirical - predicted), tolerance);
            testCase.verifyLessThanOrEqual(abs(mean(recovered, 2) - truth), ...
                6 * sqrt(diag(predicted) / draws));
        end

        function capAndDegreeZero(testCase)
            [G, V] = glmalpha(40, 3);
            [~, order] = sort(V(:), 'descend');
            c = (1:4)';
            data = testCase.synthesize(G(:, order(1:4)) * c);
            actual = xyzs2slep(data, testCase.Lon, testCase.Lat, 40, 3, truncation = 4);
            testCase.verifyEqual(actual, c, 'AbsTol', 2e-11);
            % Bravo uses unit-normalized functions and its current galpha
            % already corrects the phase. Both calls reuse the same cached
            % cap basis, avoiding sign differences across eigensolves.
            [bravo, ~, ~, ~, ~, ~, sampled] = xyz2slep(data, deg2rad(90 - testCase.Lat), ...
                deg2rad(testCase.Lon), 40, 3, 0, 0, 0, 4);
            testCase.verifyEqual(bravo / sqrt(4 * pi), actual, 'AbsTol', 2e-11);
            testCase.verifyEqual(sampled.' * bravo, data, 'AbsTol', 2e-11);
            [c0, v0, n0, j0] = xyzs2slep([2 4], [0 50], [-20 10], 60, 0);
            testCase.verifyEqual([c0 v0 n0 j0], [3 .25 .25 1], 'AbsTol', 1e-14);
            [cg, ~, ~, jg] = xyzs2slep([2 4], [0 50], [-20 10], testCase.Polygon, 0);
            testCase.verifyEqual(abs(cg), 3, 'AbsTol', 1e-14);
            testCase.verifyEqual(jg, 1);
        end

    end

    methods (Access = private)

        function A = designMatrix(testCase, J)
            A = zeros(numel(testCase.Lon), J);

            for j = 1:J
                A(:, j) = testCase.synthesize(testCase.Basis(:, j));
            end

        end

        function [c, covariance, se, source, fitted, residual] = fitUncertainty(testCase, data, J, varargin)
            [c, ~, ~, ~, ~, ~, ~, ~, ~, fitted, residual, ~, covariance, se, source] = ...
                xyzs2slep(data, testCase.Lon, testCase.Lat, testCase.Polygon, 3, ...
                'truncation', J, varargin{:});
        end

        function data = synthesize(testCase, coefficients)
            % Independently synthesize a grid and select irregular paired
            % locations. ULMO's plm2xyz currently rejects vector coordinates.
            grid = plm2xyz(testCase.toPlm(coefficients), 1, 'BeQuiet', true);
            data = grid(sub2ind(size(grid), 91 - testCase.Lat, testCase.Lon + 1));
        end

        function lmcosi = toPlm(~, coefficients)
            [~, ~, ~, lmcosi, ~, ~, ~, ~, ~, indices] = addmon(3);
            lmcosi(2 * size(lmcosi, 1) + indices(1:16)) = coefficients;
        end

    end

end
