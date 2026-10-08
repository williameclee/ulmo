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
