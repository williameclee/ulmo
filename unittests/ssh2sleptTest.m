%% SSH2SLEPTTEST - Checks monthly SSH fits, observation masks and projection caches.
%
% Created by
%   2026/11/03, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

classdef ssh2sleptTest < matlab.unittest.TestCase

    properties
        Folder
        Dates
        Polygon = [0 -11; 32 -11; 32 11; 0 11; 0 -11]
    end

    properties (TestParameter)
        polarBuffer = {0, 1.5}
        nearby = {false, true}
        forceNew = {false, true}
    end

    methods (TestMethodSetup)

        function nativeFile(t)
            % Use four native months with known coefficients and one empty month.
            tmp = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Folder = fullfile(tmp.Folder, 'SSH');
            mkdir(t.Folder);
            original = getenv('IFILES');
            t.addTeardown(@() setenv('IFILES', original));
            setenv('IFILES', tmp.Folder);
            fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fixtures));

            for folder = ["GLMALPHA", "KERNELC", "KERNELCP", "LEGENDRE"]
                mkdir(fullfile(getenv('IFILES'), folder));
            end

            constants = fullfile(getenv('IFILES'), 'EARTHMODELS', 'CONSTANTS');
            mkdir(constants);
            Earth.Radius = 6371000;
            save(fullfile(constants, 'Earth.mat'), 'Earth');
            lon = (0.25:0.5:33.75)';
            lat = (-10.75:0.5:10.75)';
            dates = datetime(2020, 1, 1) + days([15.5; 45.5; 75.5; 105]);
            partialMonths = [true; false; false; true];
            [x, y] = ndgrid(lon, lat);
            field = sqrt(3) * (2 * cosd(y) .* cosd(x) + 3 * sind(y) + 4 * cosd(y) .* sind(x));
            sshs = field + reshape(1:4, 1, 1, []);
            sshs(:, :, 1) = NaN;
            sshs(2:5, :, 2) = NaN;
            sshs(:, 3:10, 3) = NaN;

            t.Dates = dates;
            save(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'sshs', 'dates', 'lon', 'lat', 'partialMonths', '-v7.3');
        end

    end

    methods (Test)

        function descendingNonuniformCyclicSource(t)
            % The source is lon-by-lat even when both spatial axes descend.
            lon = (360:-30:0).'; lat = [50; 15; -10; -40];
            dates = datetime(2020, 1, 15); partialMonths = false;
            [x, y] = ndgrid(lon, lat); sshs = 2 + x / 100 + y / 80;
            sshs(end, :) = 1e6; % Excluded duplicate of the first longitude.
            save(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'sshs', 'dates', 'lon', 'lat', 'partialMonths');
            polygon = [10 -60; 130 -60; 130 60; 10 60; 10 -60];
            wlat = abs(diff(sind([67.5; 32.5; 2.5; -25; -55])));
            wghts = repmat(wlat.', numel(lon), 1);
            selected = x >= 10 & x <= 130;
            expected = sum(wghts(selected) .* sshs(selected)) / sum(wghts(selected));
            [c, when, fit] = ssh2slept('NASASSHMonthly', polygon, 0, GVN = {1, 1, 1}, unit = 'm');
            t.verifyEqual(c, expected, AbsTol = 1e-12); t.verifyEqual(when, dates);
            t.verifyEqual(fit.nObs, nnz(selected)); t.verifyEqual(fit.areaFraction, 1);
            t.verifyFalse(fit.partialMonth);
        end

        function invalidCoordinatesDoNotDependOnStericValidation(t)
            % Check the public SSH error for unordered or out-of-range latitudes.
            native = fullfile(t.Folder, 'NASASSHMonthly.mat');
            data = load(native); data.lat([2 3]) = data.lat([3 2]);
            save(native, '-struct', 'data');
            t.verifyError(@() ssh2slept('NASASSHMonthly', t.Polygon, 0, GVN = {1, 1, 1}), ...
            'ULMO:ssh2slept:SourceGrid');
            data.lat = sort(data.lat);
            data.lat(end) = 91;
            save(native, '-struct', 'data');
            t.verifyError(@() ssh2slept('NASASSHMonthly', t.Polygon, 0, GVN = {1, 1, 1}), ...
            'ULMO:ssh2slept:SourceGrid');
        end

        function unsavedComputedBasisCreatesNoCaches(t)
            % Check that a computed basis honours SaveData for each cache.
            root = fileparts(t.Folder);
            before = dir(fullfile(root, '**', '*'));
            ssh2slept('NASASSHMonthly', t.Polygon, 0, SaveData = false, unit = 'm');
            after = dir(fullfile(root, '**', '*'));
            t.verifyEqual({after(~[after.isdir]).name}, {before(~[before.isdir]).name});
        end

        function outputFormats(t, forceNew)
            % Fresh fits and cache reads orient coefficients after selecting dates.
            args = {"NASASSHMonthly", t.Polygon, 1};
            [expected, expectedDates, expectedFit, expectedBasis] = ...
                ssh2slept(args{:}, truncation = 4, unit = "m");
            file = dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc'));
            t.assertNumElements(file, 1);
            path = fullfile(file.folder, file.name);
            originalBytes = readbinary(path);

            if ~forceNew
                delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            end

            [actual, dates, fit, basis] = ssh2slept(args{:}, truncation = 4, ...
                unit = "m", OutputFormat = "traditional", ForceNew = forceNew);
            t.verifyEqual(actual, expected.');
            t.verifyEqual(dates, expectedDates);
            t.verifyEqual(fit, expectedFit);
            t.verifyEqual(basis, expectedBasis);
            [selected, dates, fit] = ssh2slept(args{:}, t.Dates([2 4]), ...
                truncation = 4, OutputFormat = 'traditional', ...
                dropPartialMonths = true, ForceNew = forceNew);
            t.verifySize(selected, [4 2]);
            t.verifyEqual(selected, expected(2:3, :).' * 1000);
            t.verifyEqual(dates, expectedDates(2:3));
            t.verifyEqual(fit.partialMonth, expectedFit.partialMonth(2:3));
            t.verifyEqual(fit.residualRms, expectedFit.residualRms(2:3) * 1000);
            singleDate = ssh2slept(args{:}, t.Dates([3 3]), truncation = 4, ...
                unit = "m", OutputFormat = "traditional", ForceNew = forceNew);
            t.verifySize(singleDate, [4 1]);
            t.verifyEqual(singleDate, expected(3, :).');
            [none, dates, fit] = ssh2slept(args{:}, datetime(2010, [1 2], 1), ...
                truncation = 4, OutputFormat = "traditional", ForceNew = forceNew);
            t.verifySize(none, [4 0]);
            t.verifyEmpty(dates);
            t.verifyEqual(height(fit), 0);
            t.verifyEqual(ssh2slept(args{:}, truncation = 4, unit = "m", ...
                OutputFormat = "timefirst"), expected);

            if ~forceNew
                t.verifyEqual(readbinary(path), originalBytes);
            end

            t.verifyNumElements(dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc')), 1);
        end

        function invalidOutputFormat(t)
            % Reject an unsupported orientation before fitting or creating a cache.
            t.verifyError(@() ssh2slept("NASASSHMonthly", t.Polygon, 1, ...
                OutputFormat = "invalid"), 'MATLAB:validators:mustBeMember');
            t.verifyFalse(isfolder(fullfile(t.Folder, 'SlepianExpansions')));
        end

        function readerReturnsNativeMonthlyFlags(t)
            % Subsetting by time must keep flags aligned with the native grid.
            limits = t.Dates([3 4]);
            [ssh, sigma, dates, lon, lat, partial] = ssh2lonlatt("NASASSHMonthly", ...
                [], [], limits, OutputFormat = 'ndgrid', BeQuiet = true);
            t.verifySize(ssh, [numel(lon), numel(lat), 2]);
            t.verifyEmpty(sigma);
            t.verifyEqual(dates, limits);
            t.verifyEqual(partial, [false; true]);
            t.verifyNumElements(dir(fullfile(t.Folder, '*.mat')), 1);
        end

        function serialDateLimitsSelectProjectionRows(t, forceNew)
            % Numeric limits must select the same rows on a refit and a source-free cache hit.
            args = {'NASASSHMonthly', t.Polygon, 1};
            [full, ~, diagnostics] = ssh2slept(args{:}, truncation = 4);
            if ~forceNew
                delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            end

            limits = datenum(t.Dates([2 3])); %#ok<DATNM>
            [selected, dates, fit] = ssh2slept(args{:}, limits, truncation = 4, ForceNew = forceNew);
            t.verifyEqual(dates, t.Dates(2:3));
            t.verifyEqual(selected, full(2:3, :));
            t.verifyEqual(fit, diagnostics(2:3, :));
        end

        function enumProductSharesProjectionCache(t)
            % Fit with the enum and read the same projection with either spelling.
            product = SSHProduct.NASASSHMonthly;
            [expected, dates] = ssh2slept(product, t.Polygon, 0);
            t.verifyEqual(dates, t.Dates);
            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            t.verifyEqual(ssh2slept("NASASSHMonthly", t.Polygon, 0), expected);
            t.verifyEqual(ssh2slept(product, t.Polygon, 0), expected);
            t.verifyError(@() ssh2slept(SSHProduct.NASASSH, t.Polygon, 0), ...
                'MATLAB:validators:mustBeMember');
        end

        function readerReusesMonthlyGridCache(t)
            % A flag-free request can read a legacy cache after the source is removed.
            source = fullfile(t.Folder, 'NASASSHMonthly.mat');
            cached = load(source, 'sshs', 'dates', 'lon', 'lat');
            cached.sshs = cached.sshs + 1;
            save(fullfile(t.Folder, 'NASASSHMonthly-linear.mat'), '-struct', 'cached');
            args = {'NASASSHMonthly', 'OutputFormat', 'ndgrid', 'BeQuiet', true};
            [native, ~, ~, ~, ~, flags] = ssh2lonlatt(args{:});
            t.verifyEqual(flags, [true; false; false; true]);
            t.verifyEqual(native + 1, cached.sshs);
            % An empty method resolves to the source filename, even during a forced read.
            sourceBytes = readbinary(source);
            t.verifyEqual(ssh2lonlatt(args{:}, Interpolation = "", ForceNew = true), native);
            t.verifyEqual(readbinary(source), sourceBytes);
            t.verifyEqual(ssh2lonlatt(args{:}, ForceNew = true, SaveData = false), native);
            t.verifyEqual(ssh2lonlatt(args{:}), cached.sshs);
            t.verifyEqual(ssh2lonlatt(args{:}, ForceNew = true), native);

            delete(source);
            [ssh, sigma, dates, lon, lat] = ssh2lonlatt(args{:});
            t.verifyEqual(ssh, native);
            t.verifyEmpty(sigma);
            t.verifyEqual(dates, cached.dates);
            t.verifyEqual(lon, cached.lon(:));
            t.verifyEqual(lat, cached.lat(:).');
        end

        function readerAcceptsOpenDateRanges(t)
            % Infinite endpoints select the full record or the requested side of a date.
            args = {'NASASSHMonthly', [], []};
            [~, ~, dates] = ssh2lonlatt(args{:}, [-Inf Inf], BeQuiet = true);
            t.verifyEqual(dates, t.Dates);
            middle = datenum(t.Dates(2)); %#ok<DATNM>
            [~, ~, dates] = ssh2lonlatt(args{:}, [-Inf middle], BeQuiet = true);
            t.verifyEqual(dates, t.Dates(1:2));
            [~, ~, dates] = ssh2lonlatt(args{:}, [middle Inf], BeQuiet = true);
            t.verifyEqual(dates, t.Dates(2:end));
            t.verifyError(@() ssh2slept('NASASSHMonthly', t.Polygon, 0, [-Inf Inf]), ...
                'ULMO:notValidTimeRange');
        end

        function degreeDependentExteriorFit(t)
            % A constant basis isolates which observed cells enter the fit.
            G = zeros(61 ^ 2, 1);
            G(1) = 1;
            args = {"NASASSHMonthly", t.Polygon, 60};
            data = matfile(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'Writable', true);
            data.sshs(end - 3:end, :, 1) = 900 * ones(4, size(data, 'sshs', 2));
            [c, ~, fit] = ssh2slept(args{:}, GVN = {G, 1, 1}, unit = "m", SaveData = false);
            data = load(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            % The first month only has SSH at 32.25, 32.75, 33.25 and 33.75 E.
            t.verifyEqual(fit.nObs(1), 3 * numel(data.lat));
            t.verifyEqual(fit.areaFraction(1), 0);
            t.verifyEqual(c(1), 900, 'AbsTol', 1e-8);
            t.verifyLessThanOrEqual(fit.areaFraction, ones(4, 1) + eps);
        end

        function subsetBuildsWholeRecord(t)
            % A single-date request must cache coefficients and diagnostics for all dates.
            args = {"NASASSHMonthly", t.Polygon, 1};
            [c, d, fit, b] = ssh2slept(args{:}, t.Dates([3 3]), truncation = 4);
            t.verifyEqual(c * b{1}.', [3 2 3 4] * 1000, 'AbsTol', 1e-8);
            t.verifyEqual(d, t.Dates(3));
            t.verifySize(b{1}, [4 4]);
            t.verifyLessThan(fit.residualRms, 1e-9);
            t.verifyEqual(fit.areaFraction, 36/44, 'AbsTol', 0.005);
            file = dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc'));
            t.assertNumElements(file, 1);
            expectedName = sprintf('ssh2slept-NASASSHMonthly-%s-1-4-m.nc', hash(t.Polygon, 'sha1'));
            t.verifyEqual(file.name, expectedName);
            path = fullfile(file.folder, file.name);
            info = ncinfo(path);
            t.verifyEmpty(intersect({info.Variables.Name}, {'G', 'V', 'N'}));
            t.verifyEmpty(intersect({info.Attributes.Name}, {'G', 'V', 'N'}));
            t.verifySize(ncread(path, 'coefficients'), [4 4]);
            t.verifyEqual(ncreadatt(path, 'coefficients', 'units'), 'm');
            % Cache hits must work with no source cube, including a new unit.
            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            [allC, dates, allFit, b2] = ssh2slept(args{:}, truncation = 4, unit = "m");
            t.verifyEqual(dates, t.Dates);
            t.verifyTrue(all(isnan(allC(1, :))));
            t.verifyEqual(allC(2:4, :) * b2{1}.', [(2:4)', repmat([2 3 4], 3, 1)], 'AbsTol', 1e-10);
            t.verifyEqual(allFit.nObs, [0; 64 * 44; 68 * 36; 68 * 44]);
            t.verifyEqual(allFit.partialMonth, [true; false; false; true]);
            t.verifyEqual(allFit.areaFraction([1 4]), [0; 1], 'AbsTol', 1e-12);
            t.verifyEqual(b2, b);
            [none, when, emptyFit] = ssh2slept(args{:}, datetime(2010, [1 2], 1), truncation = 4);
            t.verifySize(none, [0 4]);
            t.verifyEmpty(when);
            t.verifyEqual(height(emptyFit), 0);
            ncwriteatt(path, 'time', 'units', 'days since 2000-01-01');
            t.verifyError(@() ssh2slept(args{:}, truncation = 4), 'ULMO:readshnc:InvalidTime');
        end

        function partialMonthSelectionPreservesCache(t, forceNew)
            % Fresh fits save every month; filtering changes only returned rows.
            args = {"NASASSHMonthly", t.Polygon, 1};
            [selected, dates, fit, basis] = ssh2slept(args{:}, truncation = 4, ...
                unit = "m", dropPartialMonths = true, ForceNew = forceNew);
            t.verifyEqual(dates, t.Dates(2:3));
            t.verifyEqual(fit.partialMonth, false(2, 1));
            t.verifyEqual(selected * basis{1}.', [(2:3)', repmat([2 3 4], 2, 1)], ...
                'AbsTol', 1e-10);
            file = dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc'));
            t.assertNumElements(file, 1);
            path = fullfile(file.folder, file.name);
            t.verifySize(ncread(path, 'coefficients'), [4 4]);
            t.verifyEqual(logical(ncread(path, 'partialMonth')), [true; false; false; true]);
            originalBytes = readbinary(path);

            % All output counts must filter correctly without the source cube.
            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            oneOutput = ssh2slept(args{:}, truncation = 4, unit = "m", dropPartialMonths = true);
            t.verifyEqual(oneOutput, selected);
            [twoOutputs, dates] = ssh2slept(args{:}, truncation = 4, dropPartialMonths = true);
            t.verifyEqual(twoOutputs, selected * 1000);
            t.verifyEqual(dates, t.Dates(2:3));
            [allC, allDates, allFit, allBasis] = ssh2slept(args{:}, truncation = 4, unit = "m");
            t.verifyEqual(allDates, t.Dates);
            t.verifyEqual(allC(2:3, :), selected);
            t.verifyEqual(allFit(2:3, :), fit);
            t.verifyEqual(allBasis, basis);
            [subset, dates, subsetFit, subsetBasis] = ssh2slept(args{:}, ...
                t.Dates([3 4]), truncation = 4, unit = "m", dropPartialMonths = true);
            t.verifyEqual(subset, selected(2, :));
            t.verifyEqual(dates, t.Dates(3));
            t.verifyEqual(subsetFit, fit(2, :));
            t.verifyEqual(subsetBasis, basis);
            [none, dates, emptyFit] = ssh2slept(args{:}, t.Dates([4 4]), ...
                truncation = 4, dropPartialMonths = true);
            t.verifySize(none, [0 4]);
            t.verifyEmpty(dates);
            t.verifyEqual(height(emptyFit), 0);
            t.verifyEqual(readbinary(path), originalBytes);
            t.verifyNumElements(dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc')), 1);
        end

        function partialMonthSelectionWithoutCache(t)
            % Custom-basis fits filter partial months and allow an empty result.
            args = {"NASASSHMonthly", t.Polygon, 1};
            basis = {eye(4), ones(4, 1), 4};
            [selected, dates, fit] = ssh2slept(args{:}, GVN = basis, unit = "m", ...
                dropPartialMonths = true);
            t.verifyEqual(selected, [(2:3)', repmat([2 3 4], 2, 1)], 'AbsTol', 1e-10);
            t.verifyEqual(dates, t.Dates(2:3));
            t.verifyEqual(fit.partialMonth, false(2, 1));
            [none, dates, fit] = ssh2slept(args{:}, t.Dates([4 4]), ...
                GVN = basis, dropPartialMonths = true);
            t.verifySize(none, [0 4]);
            t.verifyEmpty(dates);
            t.verifyEqual(height(fit), 0);
            t.verifyFalse(isfolder(fullfile(t.Folder, 'SlepianExpansions')));
        end

        function repeatedMasksKeepMonthsAligned(t)
            % Months two and four share a mask but have different coefficients.
            data = matfile(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'Writable', true);
            data.sshs(:, :, 4) = data.sshs(:, :, 2) + 2;
            [c, dates, fit] = ssh2slept("NASASSHMonthly", t.Polygon, 1, ...
                GVN = {eye(4), ones(4, 1), 4}, unit = "m");
            t.verifyEqual(dates, t.Dates);
            t.verifyTrue(all(isnan([c(1, :), fit.condNum(1), fit.residualRms(1)])));
            t.verifyEqual(c(2:4, :), [(2:4)', repmat([2 3 4], 3, 1)], 'AbsTol', 1e-10);
            t.verifyEqual(fit.nObs, [0; 64 * 44; 68 * 36; 64 * 44]);
            t.verifyEqual(fit.areaFraction(2), fit.areaFraction(4));
            t.verifyEqual(fit.condNum(2), fit.condNum(4));
            t.verifyLessThan(fit.residualRms(2:4), 1e-10);
        end

        function computedBasisAndNamedDomain(t)
            % The region name and its GeoDomain must give the same projection.
            region = 'xyzs2slepregionfixture';
            [c, dates, fit, basis] = ssh2slept("NASASSHMonthly", region, 0, ...
                t.Dates([2 4]), SaveData = false, unit = "m");
            t.verifyFalse(isfolder(fullfile(t.Folder, 'SlepianExpansions')));
            t.verifyEqual(dates, t.Dates(2:4));
            t.verifySize(c, [3 1]);
            t.verifySize(basis{1}, [1 1]);
            t.verifyGreaterThan(fit.nObs, 0);
            % ForceNew must also save the first cache when the folder is absent.
            fromObject = ssh2slept("NASASSHMonthly", GeoDomain(region), 0, ...
                t.Dates([2 4]), ForceNew = true, unit = "m");
            t.verifyEqual(fromObject, c, 'AbsTol', 1e-12);
            geo = GeoDomain(region);
            file = dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc'));
            t.verifyEqual({file.name}, {sprintf('ssh2slept-NASASSHMonthly-%s-0-1-m.nc', geo.Id)});
            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            t.verifyEqual(ssh2slept("NASASSHMonthly", region, 0, t.Dates([2 4]), ...
                truncation = 1, unit = "m"), c);
            t.verifyEqual(ssh2slept("NASASSHMonthly", region, 0, t.Dates([2 4]), unit = "m"), c);
        end

        function suppliedBasisBypassesProjectionCache(t)
            % Fit the current source with the supplied basis even when a standard cache exists.
            args = {"NASASSHMonthly", t.Polygon, 1, t.Dates([4 4])};
            custom = {eye(4), ones(4, 1), 4};
            expected = [4 2 3 4];
            folder = fullfile(t.Folder, 'SlepianExpansions');
            [c, ~, ~, basis] = ssh2slept(args{:}, GVN = custom, unit = "m");
            t.verifyEqual(basis, custom);
            t.verifyEqual(c, expected, 'AbsTol', 1e-10);
            t.verifyFalse(isfolder(folder));
            t.verifyEmpty(dir(fullfile(getenv('IFILES'), 'DOMAINMASK', '*.mat')));

            % Make the source differ from its standard cache at the same L and J.
            [cached, ~, ~, standard] = ssh2slept(args{:}, truncation = 4, unit = "m");
            data = matfile(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'Writable', true);
            data.sshs(:, :, 4) = 2 * data.sshs(:, :, 4);
            fresh = ssh2slept(args{:}, GVN = standard, truncation = 4, unit = "m", ...
                ForceNew = false, SaveData = true);
            t.verifyEqual(fresh, 2 * cached, 'AbsTol', 1e-10);
            t.verifyEqual(ssh2slept(args{:}, truncation = 4, unit = "m"), cached);
            t.verifyNumElements(dir(fullfile(folder, '*.nc')), 1);
        end

        function separateTruncationsAndRebuild(t)
            % Check default mode selection and cache replacement after successful and failed fits.
            args = {"NASASSHMonthly", t.Polygon, 1};
            original = ssh2slept(args{:}, truncation = 4);
            % An explicit four-mode file cannot satisfy the one-mode default.
            [defaultC, ~, ~, expectedBasis] = ssh2slept(args{:});
            [~, ~, ~, cachedBasis] = ssh2slept(args{:});
            t.verifyEqual(cachedBasis, expectedBasis);
            t.verifySize(cachedBasis{1}, [4 1]);
            t.verifySize(defaultC, [4 1]);
            t.verifyEqual(ssh2slept(args{:}, truncation = 1), defaultC);
            t.verifyNumElements(dir(fullfile(t.Folder, 'SlepianExpansions', '*.nc')), 2);
            data = matfile(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'Writable', true);
            data.sshs(:, :, 4) = 2 * data.sshs(:, :, 4);
            rebuilt = ssh2slept(args{:}, t.Dates([4 4]), truncation = 4, ForceNew = true);
            t.verifyEqual(rebuilt, 2 * original(4, :), 'AbsTol', 1e-8);
            % A failed refresh preserves the complete previous cache.
            bad = nan(size(data, 'sshs', 1), size(data, 'sshs', 2));
            bad(1) = 1;
            data.sshs(:, :, 4) = bad;
            t.verifyError(@() ssh2slept(args{:}, truncation = 4, ForceNew = true), 'ULMO:ssh2slept:FitFailed');
            t.verifyEqual(ssh2slept(args{:}, t.Dates([4 4]), truncation = 4), rebuilt);
        end

        function sshUsesPolarObservations(t, polarBuffer)
            % Polar SSH is 2 m. Equatorial SSH is 99 m, including an empty polar month.
            t.writePolarSSH();
            degree = 0;

            if polarBuffer > 0
                degree = 60;
            end

            G = zeros((degree + 1) ^ 2, 1);
            G(1) = 1;
            [coeffs, ~, fit] = ssh2slept('NASASSHMonthly', 'arctic', degree, ...
                GVN = {G, 1, 1}, unit = 'm');
            t.verifyEqual(coeffs, [2; NaN], 'AbsTol', 1e-10);
            t.verifyEqual(fit.areaFraction, [1; 0], 'AbsTol', 1e-12);
            t.verifyEqual(fit.nObs(2), 0);
        end

        function ignoresRotatedProjectionCache(t)
            % An existing Arctic projection must be rebuilt using geographic observations.
            t.writePolarSSH();
            folder = fullfile(fileparts(t.Folder), 'GLMALPHA');
            G = 1;
            V = 0.01;
            N = V;
            EL = 0;
            EM = 0;
            save(fullfile(folder, 'glmalpha-Arctic-0-1.mat'), 'G', 'V', 'N', 'EL', 'EM');
            folder = fullfile(fileparts(t.Folder), 'SSH', 'SlepianExpansions');
            mkdir(folder);
            path = fullfile(folder, 'ssh2slept-NASASSHMonthly-Arctic-0-1-m.nc');
            write2nc(path, 'coefficients', [99 99], Dimensions = {'mode', 'time'});
            ncwriteatt(path, '/', 'defaultTruncation', 1);
            dates = datetime(2020, [1; 2], 15);
            nccreate(path, 'time', 'Dimensions', {'time', 2});
            ncwrite(path, 'time', seconds(dates - datetime(2000, 1, 1)));
            ncwriteatt(path, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
            ncwriteatt(path, 'time', 'calendar', 'proleptic_gregorian');
            actual = ssh2slept('NASASSHMonthly', 'arctic', 0, unit = 'm');
            t.verifyEqual(actual, [2; NaN], 'AbsTol', 1e-12);
            delete(fullfile(fileparts(t.Folder), 'SSH', 'NASASSHMonthly.mat'));
            t.verifyEqual(ssh2slept('NASASSHMonthly', 'arctic', 0, ...
                truncation = 1, unit = 'm'), actual);
        end

        function nearbySelectsProjectionCache(t, nearby)
            % The same projection settings must retain separate NearBy coefficients.
            root = fileparts(t.Folder);
            folder = fullfile(t.Folder, 'SlepianExpansions');
            mkdir(folder);
            identifiers = {'Greenland0-2', 'Greenland0-2-nearby'};
            dates = datetime(2020, [1; 2], 15);

            for iDomain = 1:2
                id = identifiers{iDomain};
                path = fullfile(root, 'SSH', 'SlepianExpansions', ...
                    sprintf('ssh2slept-NASASSHMonthly-%s-0-1-m.nc', id));
                write2nc(path, 'coefficients', iDomain * [10 11], Dimensions = {'mode', 'time'});
                nccreate(path, 'time', 'Dimensions', {'time', 2});
                ncwrite(path, 'time', seconds(dates - datetime(2000, 1, 1)));
                ncwriteatt(path, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
                ncwriteatt(path, 'time', 'calendar', 'proleptic_gregorian');
                ncwriteatt(path, '/', 'defaultTruncation', 1);
            end

            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            domain = GeoDomain('greenland', 'Buffer', 2, 'NearBy', nearby);
            expectedIndex = 1 + double(nearby);
            [coeffs, actualDates] = ssh2slept('NASASSHMonthly', domain, 0, unit = 'm');
            t.verifyEqual(coeffs, expectedIndex * [10; 11]);
            t.verifyEqual(actualDates, dates);
            t.verifyEqual(ssh2slept('NASASSHMonthly', domain, 0, ...
                truncation = 1, unit = 'm'), coeffs);
        end

    end

    methods (Access = private)

        function writePolarSSH(t)
            % Supply polar coastlines and SSH spanning both coordinate frames.
            coasts = fullfile(fileparts(t.Folder), 'COASTS');
            mkdir(coasts);
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('COASTS', coasts));
            writeArcticCoastlineFixture(coasts, [10 60; 20 60; 20 65; 10 65; 10 60]);
            folder = fullfile(fileparts(t.Folder), 'SSH');
            lon = (0.25:0.5:359.75)';
            lat = (-34.75:0.5:74.75)';
            dates = datetime(2020, [1; 2], 15);
            partialMonths = false(2, 1);
            sshs = 99 * ones(numel(lon), numel(lat), 2);
            sshs(:, lat > 50, 1) = 2;
            sshs(:, lat > 50, 2) = NaN;
            save(fullfile(folder, 'NASASSHMonthly.mat'), ...
                'sshs', 'lon', 'lat', 'dates', 'partialMonths', '-v7.3');
        end

    end

end

function bytes = readbinary(path)
    file = fopen(path, 'rb');
    cleanup = onCleanup(@() fclose(file));
    bytes = fread(file, Inf, '*uint8');
end
