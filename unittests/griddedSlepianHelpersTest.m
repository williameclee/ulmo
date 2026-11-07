%% GRIDDEDSLEPIANHELPERSTEST - Verify product-independent spatial workers.
% Synthetic grids test analytic weighted means, orientation, empty epochs,
% target/buffer diagnostics and shared designs. All caches use temporary roots.
%
% Syntax
%   results = runtests('unittests/griddedSlepianHelpersTest.m')
%
% Created by
%   2026/11/06, En-Chi Lee

classdef griddedSlepianHelpersTest < matlab.unittest.TestCase
    properties
        Polygon = [0 -60; 120 -60; 120 60; 0 60; 0 -60]
        Folder
    end
    properties (TestParameter)
        descending = {false, true}
        nonuniform = {false, true}
    end
    methods (TestMethodSetup)
        function isolate(t)
            tmp = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Folder = tmp.Folder;
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', t.Folder));
            for folder = ["GLMALPHA", "KERNELC", "KERNELCP", "LEGENDRE", "EARTHMODELS/CONSTANTS"]
                mkdir(fullfile(t.Folder, folder));
            end
            Earth.Radius = 6371000;
            save(fullfile(t.Folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
        end
    end
    methods (Test)
        function orderedCoordinatesAndAreas(t, descending, nonuniform)
            lon = (0:60:360).'; lat = [-60; -20; 20; 60];
            expected = cosd(lat);
            if nonuniform
                lat = [-70; -30; 0; 50];
                expected = diff(sind([-90; -50; -15; 25; 75]));
            end
            if descending, lon = flip(lon); lat = flip(lat); end
            [x, y, ix, iy, w, name] = preparegridcoordinates(lon, lat, RequireGlobalLongitude=true);
            t.verifyEqual(x, lon(ix)); t.verifyEqual(y, lat(iy));
            t.verifyTrue(issorted(x)); t.verifyTrue(issorted(y));
            t.verifySize(w, [6 4]);
            t.verifyEqual(w, repmat(expected.', 6, 1), AbsTol=1e-15);
            t.verifyFalse(any(ix == numel(lon)));
            if nonuniform
                t.verifyEqual(name, 'spherical-midpoint-bands');
            else
                t.verifyEqual(name, 'cosine-latitude');
            end
        end

        function regionalCoverageIsExplicit(t)
            [x, ~, ~, ~, w] = preparegridcoordinates(10:10:80, [-30 0 30]);
            t.verifyNumElements(x, 8); t.verifySize(w, [8 3]);
            t.verifyError(@() preparegridcoordinates(10:10:80, [-30 0 30], ...
                RequireGlobalLongitude=true), 'ULMO:preparegridcoordinates:InvalidGrid');
        end

        function rejectsMalformedCoordinates(t)
            bad = {{[0 1 3], [-30 0 30]}, {[0 1 2], [-30 30 0]}, ...
                {[0 1 2], [-30 0 91]}, {[0 NaN 2], [-30 0 30]}, ...
                {[0 0 0], [-30 0 30]}, {[0 180 360 540], [-30 0 30]}};
            for iCase = 1:numel(bad)
                t.verifyError(@() preparegridcoordinates(bad{iCase}{:}), ...
                    'ULMO:preparegridcoordinates:InvalidGrid');
            end
        end

        function singlePrecisionFractionalGrid(t)
            lon = single((0:4319)/12); lat = single(-89.5:89.5);
            [x, ~, ~, ~, w] = preparegridcoordinates(lon, lat, RequireGlobalLongitude=true);
            t.verifyNumElements(x, 4320);
            t.verifyEqual(w(1, :), cosd(double(lat)), AbsTol=1e-15);
        end

        function analyticMeanAndEmptyEpoch(t, descending, nonuniform)
            lon = (0:30:360).'; lat = [-50; -20; 10; 40];
            latitudeWeights = cosd(lat);
            if nonuniform
                lat = [-50; -30; 0; 40];
                latitudeWeights = diff(sind([-60; -40; -15; 20; 60]));
            end
            [x, y] = ndgrid(lon, lat);
            values = 2 + x/100 + y/70;
            values(2, 2) = NaN;
            % The duplicate endpoint must never enter the fit.
            values(end, :) = 1e6;
            selected = x >= 0 & x <= 120 & x < 360 & isfinite(values);
            w = repmat(latitudeWeights.', numel(lon), 1);
            expected = sum(w(selected).*values(selected))/sum(w(selected));
            expectedArea = sum(w(selected))/sum(w(x <= 120));
            values = cat(3, values, nan(size(values)), values+3);
            if descending
                lon = flip(lon); lat = flip(lat);
                values = values(end:-1:1, end:-1:1, :);
                % The retained cyclic copy is now the original last row.
                values(1, :, :) = values(end, :, :);
            end
            dates = datetime(2020, 1, (1:3).');
            [c, fits, buffer] = fitgriddedslepian(struct('height', values), dates, ...
                lon, lat, t.Polygon, 0, {1, 1, 1}, SaveData=false);
            t.verifyEqual(c.height([1 3]), [expected; expected+3], AbsTol=2e-12);
            t.verifyTrue(isnan(c.height(2)));
            t.verifyEqual(fits.height.nObs, [nnz(selected); 0; nnz(selected)]);
            t.verifyEqual(fits.height.areaFraction, [expectedArea; 0; expectedArea], AbsTol=1e-14);
            t.verifyTrue(isnan(fits.height.condNum(2)));
            t.verifyTrue(isnan(fits.height.residualRms(2)));
            t.verifyEqual(buffer, 0);
        end

        function laterFieldsActivateUnusedDesignRows(t)
            lon = (0:30:330).'; lat = (-45:30:45).';
            [x, y] = ndgrid(lon, lat);
            first = 2 + cosd(x).*cosd(y); first(2:2:end, :) = NaN;
            second = 4 + sind(y); second(1:2:end, :) = NaN;
            dates = datetime(2020, 1, 1);
            basis = {eye(4), ones(4,1), 4};
            [combined, fits] = fitgriddedslepian(struct('a', first, 'b', second), ...
                dates, lon, lat, t.Polygon, 1, basis, SaveData=false);
            [alone, aloneFits] = fitgriddedslepian(struct('b', second), ...
                dates, lon, lat, t.Polygon, 1, basis, SaveData=false);
            t.verifyEqual(combined.b, alone.b, AbsTol=1e-12);
            t.verifyEqual(fits.b.nObs, aloneFits.b.nObs);
            t.verifyEqual(fits.b{:,2:4}, aloneFits.b{:,2:4}, AbsTol=1e-14);
            t.verifyLessThan(fits.a.residualRms, 1e-12);
            t.verifyLessThan(fits.b.residualRms, 1e-12);
        end

        function bufferOnlyObservationsHaveZeroTargetCoverage(t)
            lon = (0:10:180).'; lat = (-50:10:50).';
            values = nan(numel(lon), numel(lat)); values(lon == 140, lat == 0) = 7;
            basis = {[1; 0; 0; 0], 1, 1};
            [c, fits, buffer] = fitgriddedslepian(struct('height', values), datetime(2020,1,1), ...
                lon, lat, t.Polygon, 1, basis, SaveData=false);
            t.verifyEqual(c.height, 7, AbsTol=1e-12);
            t.verifyEqual(fits.height.nObs, 1);
            t.verifyEqual(fits.height.areaFraction, 0);
            t.verifyEqual(buffer, 90);
        end

        function invalidAndUnderdeterminedFields(t)
            lon = (0:30:330).'; lat = [-45; -15; 15; 45];
            values = nan(12,4);
            call = @(data, basis) fitgriddedslepian(struct('height', data), datetime(2020,1,1), ...
                lon, lat, t.Polygon, 1, basis, SaveData=false);
            t.verifyError(@() call(values, {eye(4), ones(4,1), 4}), 'ULMO:fitgriddedslepian:EmptyRecord');
            values(2,2) = 1;
            t.verifyError(@() call(values, {eye(4), ones(4,1), 4}), 'ULMO:fitgriddedslepian:FitFailed');
            values(2,2) = Inf;
            t.verifyError(@() call(values, {eye(4), ones(4,1), 4}), 'ULMO:fitgriddedslepian:InvalidObservation');
            t.verifyError(@() call(values.', {eye(4), ones(4,1), 4}), 'ULMO:fitgriddedslepian:InvalidGrid');
        end

        function suppliedBasisKeepsSignsOrderAndCount(t)
            G = eye(4); G = G(:, [4 1 3 2]); G(:, 1) = -G(:, 1);
            V = [0.2; 0.9; 0.3; 0.7];
            b = prepareslepianbasis(t.Polygon, 1, [], {G, V, 2});
            t.verifyEqual(b, {G(:,1:2), V(1:2), 2});
            b = prepareslepianbasis(t.Polygon, 1, [], {G(:,1:3), V(1:3), 1});
            t.verifyEqual(b, {G(:,1:3), V(1:3), 1});
            b = prepareslepianbasis(t.Polygon, 1, 1, {G, V, 2});
            t.verifyEqual(b, {G(:,1), V(1), 2});
            t.verifyError(@() prepareslepianbasis(t.Polygon, 1, 4, {G(:,1:3),V(1:3),1}), ...
                'ULMO:prepareslepianbasis:Truncation');
        end

        function orthonormalityIsAnExplicitCallerRequirement(t)
            G = single(2*eye(4));
            b = prepareslepianbasis(t.Polygon, 1, 4, {G, ones(4,1), 4});
            t.verifyEqual(b{1}, double(G));
            t.verifyError(@() prepareslepianbasis(t.Polygon, 1, 4, {G, ones(4,1), 4}, ...
                RequireOrthonormal=true), 'ULMO:prepareslepianbasis:Basis');
            G = single(eye(4)); G(1,2) = 1e-7;
            prepareslepianbasis(t.Polygon, 1, 4, {G, ones(4,1), 4}, RequireOrthonormal=true);
        end

        function computedBasisMatchesSolverAndDoesNotSave(t)
            before = dir(fullfile(t.Folder, '**', '*'));
            b = prepareslepianbasis(t.Polygon, 1, 4, {}, SaveData=false);
            [x, y] = ndgrid(0:30:330, -45:30:45);
            [~, ~, ~, ~, design] = xyzs2slep(ones(numel(x),1), x(:), y(:), t.Polygon, 1, ...
                GVN=b, truncation=4);
            t.verifyEqual(design.GVN, b);
            after = dir(fullfile(t.Folder, '**', '*'));
            t.verifyEqual({after(~[after.isdir]).name}, {before(~[before.isdir]).name});
        end
    end
end
