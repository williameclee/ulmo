%% STERIC2SLEPTTEST - Verify complete steric projections and correction metadata.
% Uses synthetic native products and isolated data roots. Tests exercise
% temporal-before-spatial fitting, independent masks, unit-C00 projection,
% all-field NetCDF replacement and matching coefficient/GMSL orientation.
%
% Created by
%   2026/11/05, En-Chi Lee

classdef steric2sleptTest < matlab.unittest.TestCase
    properties
        Folder
        NativePath
        Data
        Polygon = [0 -60; 120 -60; 120 60; 0 60; 0 -60]
    end
    properties (TestParameter)
        temporal = struct('native', {{[], 'linear'}}, ...
            'monthly', {{'midmonth', 'linear'}}, 'grace', {{'GRACE', 'linear'}}, ...
            'linear', {{days(7), 'linear'}}, 'pchip', {{days(7), 'pchip'}}, ...
            'mean', {{days(10), 'mean'}})
        refresh = {false, true}
        ordering = {'timefirst', 'traditional'}
    end
    methods (TestMethodSetup)
        function prepare(t)
            tmp = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Folder = tmp.Folder;
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', t.Folder));
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', t.Folder));
            for folder = ["GLMALPHA", "KERNELC", "KERNELCP", "LEGENDRE", "EARTHMODELS/CONSTANTS"]
                mkdir(fullfile(t.Folder, folder));
            end
            Earth.Radius = 6371000;
            save(fullfile(t.Folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
            reference = fullfile(t.Folder, 'RL06', 'CSR'); mkdir(reference);
            writelines(["# End of YAML header"; ...
                "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], ...
                fullfile(reference, 'GSM-2_test_BA01_0600'));
            data.lon = (15:30:345).';
            data.lat = (-75:30:75).';
            data.dates = datetime(2020, 1, (1:60).');
            [x, y] = meshgrid(data.lon, data.lat);
            base = 10 + cosd(y).*cosd(x) + 0.5*sind(y) + ...
                reshape(sin((1:60)/11), 1, 1, []);
            data.shallowStericSl = 2*base;
            data.deepStericSl = 0.3*base + 0.1;
            data.stericSl = data.shallowStericSl + data.deepStericSl;
            data.shallowHalostericSl = base/10;
            data.deepHalostericSl = 0.02*base;
            data.halostericSl = data.shallowHalostericSl + data.deepHalostericSl;
            data.shallowThermostericSl = 1.4*base;
            data.deepThermostericSl = 0.15*base;
            data.thermostericSl = data.shallowThermostericSl + data.deepThermostericSl;
            t.Data = data;
            t.NativePath = stericfilepath('EN4c14');
            mkdir(fileparts(t.NativePath));
            t.writeNative;
        end
    end
    methods (Test)
        function decimalDomainCacheReadsWithoutNativeInput(t)
            % Decimal IDs must remain intact during automatic cache discovery.
            domain = GeoDomain('satlantic', Buffer=.5, DefaultParams=true, Latlim=[-60, 60]);
            t.assertTrue(contains(domain.Id, '0.5'));
            expected = steric2slept('EN4c14', domain, 0, BeQuiet=1);
            path = t.cachePath;
            original = readbinary(path);
            delete(t.NativePath);
            % Both explicit and inferred truncation must read the existing file.
            explicit = steric2slept('EN4c14', domain, 0, truncation=1, SaveData=false, BeQuiet=1);
            t.verifyEqual(explicit, expected);
            [actual, dates] = steric2slept('EN4c14', domain, 0, SaveData=false, BeQuiet=1);
            t.verifyEqual(actual, expected);
            t.verifyEqual(dates, t.Data.dates);
            t.verifyEqual(readbinary(path), original);
            t.verifyNumElements(dir(fullfile(fileparts(path), '*.nc')), 1);
        end

        function agreesWithLoaderThenFit(t, temporal)
            % Epoch residuals must use paired native samples before averaging.
            t.Data.thermostericSl(2, 2, 1:9) = NaN;
            t.Data.shallowThermostericSl(2, 2, 1:9) = NaN;
            t.Data.deepThermostericSl(2, 2, 1:9) = NaN;
            t.Data.shallowStericSl(3, 3, 3:8) = NaN;
            t.Data.stericSl(3, 3, 3:8) = NaN;
            t.Data.deepStericSl(3, 3, 3:8) = NaN;
            t.writeNative;
            step = temporal{1}; method = temporal{2};
            args = {'EN4c14', t.Polygon, 1, step};
            steric2slept(args{:}, truncation=4, Interpolation=method, BeQuiet=1);
            schema = stericfieldschema(whos('-file', t.NativePath));
            for iField = find(schema.supported(:)).'
                [iType, iDepth] = ind2sub(size(schema.supported), iField);
                [c, dates, fit, basis, h, constant] = steric2slept(args{:}, ...
                    Type=schema.types(iType), Depth=schema.depths(iDepth), ...
                    truncation=4, Interpolation=method, BeQuiet=1);
                [grid, ~, expectedDates, lon, lat, expectedH] = steric2lonlatt('EN4c14', step, ...
                    Type=schema.types(iType), Depth=schema.depths(iDepth), Interpolation=method, ...
                    OutputFormat='ndgrid', SaveData=false, BeQuiet=1);
                [expected, expectedFit] = t.fitReference(grid, lon, lat, basis);
                t.verifyEqual(c, expected, AbsTol=2e-10);
                t.verifyEqual(dates, expectedDates);
                t.verifyEqual(h, expectedH, AbsTol=1e-12);
                t.verifyEqual(fit.nObs, expectedFit.nObs(:));
                t.verifyEqual(fit.residualRms, expectedFit.residualRms(:), AbsTol=2e-10);
                t.verifyEqual(constant, basis{1}(1, :).');
            end
        end

        function completeCacheAndFormats(t, refresh, ordering)
            % One selected output still creates every supported field/diagnostic.
            args = {'EN4c14', t.Polygon, 1, []};
            first = steric2slept(args{:}, truncation=4, BeQuiet=1);
            path = t.cachePath;
            info = ncinfo(path);
            names = string({info.Variables.Name});
            t.verifyEqual(nnz(endsWith(names, 'Slept')), 16); % 15 fields + constant.
            t.verifyEqual(nnz(endsWith(names, 'Gmsl')), 15);
            t.verifyEqual(nnz(endsWith(names, '_nObs')), 15);
            t.verifyEqual(string({ncinfo(path, 'stericSlept').Dimensions.Name}), ["mode", "time"]);
            t.verifyEqual(ncread(path, 'stericSlept'), first.');
            t.verifyEqual(ncreadatt(path, 'stericSlept', 'units'), 'm');
            t.verifyEqual(ncreadatt(path, 'constantSlept', 'units'), '1');
            t.verifyGreaterThan(ncread(path, 'creationDate'), 0);
            original = readbinary(path);
            if ~refresh, delete(t.NativePath); end
            [c, dates, fit, basis, h, constant] = steric2slept(args{:}, ...
                t.Data.dates([5 15]), truncation=4, ForceNew=refresh, ...
                Unit='mm', TimeFormat='datenum', OutputFormat=ordering, BeQuiet=1);
            expectedH = 1000*ncread(path, 'halostericGmsl', 5, 11);
            expected = 1000*first(5:15, :);
            if strcmp(ordering, 'traditional')
                expected = expected.'; expectedH = expectedH.';
                corrected = c - constant*h;
            else
                corrected = c - h*constant.';
            end
            t.verifySize(corrected, size(c));
            t.verifyEqual(c, expected, AbsTol=1e-9);
            t.verifyEqual(h, expectedH, AbsTol=1e-10);
            t.verifyEqual(dates, datenum(t.Data.dates(5:15))); %#ok<DATNM>
            t.verifyEqual(height(fit), 11);
            t.verifyLessThan(max(fit.residualRms), 1e-8);
            t.verifyEqual(constant, basis{1}(1, :).');
            % Empty and singleton time dimensions must survive formatting.
            for limits = {datetime(2010, [1 2], 1), t.Data.dates(5) + seconds([-1 1])}
                [c, dates, ~, ~, h] = steric2slept(args{:}, limits{1}, ...
                    truncation=4, OutputFormat=ordering, BeQuiet=1);
                if strcmp(ordering, 'traditional')
                    t.verifySize(c, [4, numel(dates)]); t.verifySize(h, [1, numel(dates)]);
                else
                    t.verifySize(c, [numel(dates), 4]); t.verifySize(h, [numel(dates), 1]);
                end
            end
            if ~refresh, t.verifyEqual(readbinary(path), original); end
            t.verifyNumElements(dir(fullfile(fileparts(path), '*.nc')), 1);
        end

        function forceRefreshUpdatesAllFieldsAndMeans(t)
            args = {'EN4c14', t.Polygon, 0, 'midmonth'};
            first = steric2slept(args{:}, BeQuiet=1);
            path = t.cachePath;
            previous = ncread(path, 'deepHalostericSlept');
            oldMean = ncread(path, 'deepHalostericGmsl');
            schema = stericfieldschema(whos('-file', t.NativePath));
            data = t.Data;
            for name = schema.gridNames(schema.available).'
                data.(name) = 2*data.(name);
            end
            % Keep saved means stale to prove ForceNew reaches the native reader.
            save(t.NativePath, '-struct', 'data', '-append');
            t.verifyEqual(steric2slept(args{:}, BeQuiet=1), first);
            actual = steric2slept(args{:}, ForceNew=true, BeQuiet=1);
            t.verifyEqual(actual, 2*first, AbsTol=1e-11);
            t.verifyEqual(ncread(path, 'deepHalostericSlept'), 2*previous, AbsTol=1e-11);
            t.verifyEqual(ncread(path, 'deepHalostericGmsl'), 2*oldMean, AbsTol=1e-12);
            native = load(t.NativePath, 'deepHalostericGmsl');
            t.verifyEqual(native.deepHalostericGmsl, 2*stericgmslExpected(t.Data), AbsTol=1e-12);
        end

        function noWritesAndSuppliedBasis(t)
            args = {'EN4c14', t.Polygon, 1, days(7)};
            original = snapshot(t.Folder);
            [c, ~, ~, basis, h, constant] = steric2slept(args{:}, ...
                GVN={eye(4), ones(4,1), 4}, SaveData=false, ForceNew=true, BeQuiet=1);
            t.verifyEqual(snapshot(t.Folder), original);
            % Reorder and flip signs; C00 must use exactly those columns.
            G = eye(4); G = G(:, [3 1 4 2]); G(:, 2) = -G(:, 2);
            [other, ~, ~, supplied, otherH, otherConstant] = steric2slept(args{:}, ...
                GVN={G, (1:4).', 4}, SaveData=false, BeQuiet=1);
            t.verifyEqual(other*supplied{1}.', c*basis{1}.', AbsTol=1e-11);
            t.verifyEqual(otherH, h);
            t.verifyEqual(otherConstant, [0;-1;0;0]);
            t.verifyEqual(constant, [1;0;0;0]);
            t.verifyEqual(snapshot(t.Folder), original);
            t.verifyError(@() steric2slept(args{:}, GVN={2*eye(4),ones(4,1),4}, ...
                SaveData=false, BeQuiet=1), 'ULMO:steric2slept:Basis');
        end

        function noWritesComputedBasisAndCachedRebuild(t)
            args = {'EN4c14', t.Polygon, 1, 'midmonth'};
            original = snapshot(t.Folder);
            c = steric2slept(args{:}, SaveData=false, truncation=4, BeQuiet=1);
            t.verifyEqual(snapshot(t.Folder), original);
            saved = steric2slept(args{:}, truncation=4, BeQuiet=1);
            t.verifyEqual(c, saved, AbsTol=1e-10);
            original = snapshot(t.Folder);
            rebuilt = steric2slept(args{:}, ForceNew=true, SaveData=false, truncation=4, BeQuiet=1);
            t.verifyEqual(rebuilt, saved, AbsTol=1e-10);
            t.verifyEqual(snapshot(t.Folder), original);
            supplied = steric2slept(args{:}, GVN={eye(4),ones(4,1),4}, SaveData=false, BeQuiet=1);
            t.verifyNotEqual(supplied, saved);
            t.verifyEqual(snapshot(t.Folder), original);
        end

        function missingDeepUsesZeroWithoutReconstruction(t)
            % Zero filling cannot repair an inconsistent full/shallow split.
            [x,y] = meshgrid(t.Data.lon, t.Data.lat);
            extra = 2*cosd(2*x).*cosd(y).^2;
            t.Data.deepStericSl = t.Data.deepStericSl + extra;
            t.Data.stericSl = t.Data.shallowStericSl + t.Data.deepStericSl;
            t.Data.deepStericSl(2:4, 1:2, 1:5) = NaN;
            t.writeNative;
            args = {'EN4c14', t.Polygon, 1, []};
            [full, ~, fullFit, basis, ~, constant] = steric2slept(args{:}, truncation=4, BeQuiet=1);
            [deep, ~, deepFit] = steric2slept(args{:}, truncation=4, Depth='deep', BeQuiet=1);
            shallow = steric2slept(args{:}, truncation=4, Depth='shallow', BeQuiet=1);
            t.verifyEqual(fullFit.nObs, deepFit.nObs);
            t.verifyGreaterThan(norm(full(1,:)-shallow(1,:)-deep(1,:)), 0.01);
            t.verifyEqual(full(6:end,:), shallow(6:end,:)+deep(6:end,:), AbsTol=1e-10);
            t.verifyEqual(constant, basis{1}(1,:).');
            [grid,~,~,lon,lat] = steric2lonlatt('EN4c14', Depth='deep', OutputFormat='ndgrid', BeQuiet=1);
            grid(1:2,2:4,1:5) = 0;
            t.verifyEqual(deep, t.fitReference(grid, lon, lat, basis), AbsTol=1e-10);
        end

        function absentDeepLayersCloseBeforeAveraging(t)
            % An absent contribution during part of a month must enter its mean
            % as zero, rather than be omitted from the temporal denominator.
            schema = stericfieldschema(whos('-file',t.NativePath));
            for iType = 1:3
                deepName = schema.gridNames(iType,3);
                shallowName = schema.gridNames(iType,2);
                fullName = schema.gridNames(iType,1);
                t.Data.(deepName)(2:4,1:2,1:10) = 0;
                t.Data.(fullName) = t.Data.(shallowName) + t.Data.(deepName);
                t.Data.(deepName)(2:4,1:2,1:10) = NaN;
            end
            t.writeNative;
            args = {'EN4c14',t.Polygon,1,'midmonth'};
            for type = schema.types
                [full,~,fullFit] = steric2slept(args{:},Type=type,truncation=4,BeQuiet=1);
                [shallow,~,shallowFit] = steric2slept(args{:},Type=type,Depth='shallow',truncation=4,BeQuiet=1);
                [deep,dates,deepFit,basis,~,constant] = steric2slept(args{:},Type=type,Depth='deep',truncation=4,BeQuiet=1);
                t.verifyEqual(full,shallow+deep,AbsTol=1e-10);
                t.verifyEqual(deepFit.nObs,fullFit.nObs);
                t.verifyEqual(deepFit.areaFraction,shallowFit.areaFraction);
                t.verifyEqual(constant,basis{1}(1,:).');
                t.verifyNumElements(dates,2);
            end
            % An independent temporal mean of the physically zero-filled deep grid.
            reference = t.Data.deepStericSl; reference(2:4,1:2,1:10) = 0;
            reference = cat(3,mean(reference(:,:,1:31),3),mean(reference(:,:,32:60),3));
            reference = permute(reference,[2 1 3]);
            actual = steric2slept(args{:},Depth='deep',truncation=4,BeQuiet=1);
            t.verifyEqual(actual,t.fitReference(reference,t.Data.lon,t.Data.lat,basis),AbsTol=1e-10);
            % Filling is a projection convention; the supplied grids remain intact.
            stored = load(t.NativePath,'deepStericSl');
            t.verifyEqual(stored.deepStericSl,t.Data.deepStericSl);
            t.verifySubstring(ncreadatt(t.cachePath,'/','deepSupport'),'before temporal processing');
        end

        function deepFitsUsePairedUpperSupport(t)
            % Preserve full/shallow observations and exclude unmatched deep cells.
            t.Data.stericSl(2,2,1:5) = NaN;
            t.Data.shallowStericSl(3,3,1:5) = NaN;
            t.writeNative;
            [actual,~,fit,basis] = steric2slept('EN4c14',t.Polygon,1,truncation=4,Depth='deep',BeQuiet=1);
            expectedGrid = t.Data.deepStericSl;
            expectedGrid(2,2,1:5) = NaN; expectedGrid(3,3,1:5) = NaN;
            [expected,expectedFit] = t.fitReference(permute(expectedGrid,[2 1 3]),t.Data.lon,t.Data.lat,basis);
            t.verifyEqual(actual,expected,AbsTol=1e-10);
            t.verifyEqual(fit.nObs,expectedFit.nObs(:));
        end

        function constantProjectionWithTruncation(t)
            % Harmonic C00 projection differs from fitting ones regionally.
            G = [1; 1; 0; 0]/sqrt(2);
            args = {'EN4c14', t.Polygon, 1};
            [~,~,~,basis,~,constant] = steric2slept(args{:}, GVN={G,1,1}, SaveData=false, BeQuiet=1);
            t.verifyEqual(constant, 1/sqrt(2), AbsTol=1e-15);
            [x,y] = ndgrid(t.Data.lon,t.Data.lat);
            fitted = xyzs2slep(ones(numel(x),1), x(:), y(:), t.Polygon, 1, ...
                GVN=basis, weight=cosd(y(:)));
            t.verifyGreaterThan(abs(fitted-constant), 0.001);
            t.Data.stericSl(2,2,1:6) = NaN;
            t.writeNative;
            [~,~,~,~,~,other] = steric2slept(args{:}, GVN=basis, SaveData=false, BeQuiet=1);
            t.verifyEqual(other, constant);
        end

        function emptyEpochAndFailurePreserveCache(t)
            schema = stericfieldschema(whos('-file', t.NativePath));
            for name = schema.gridNames(schema.available).'
                t.Data.(name)(:,:,1) = NaN;
            end
            t.writeNative;
            args = {'EN4c14', t.Polygon, 1};
            [c,~,fit] = steric2slept(args{:}, truncation=4, BeQuiet=1);
            t.verifyTrue(all(isnan(c(1,:))));
            t.verifyEqual(fit.nObs(1),0); t.verifyEqual(fit.areaFraction(1),0);
            t.verifyTrue(isnan(fit.condNum(1)) && isnan(fit.residualRms(1)));
            path = t.cachePath; original = readbinary(path);
            value = t.Data.stericSl(3,2,2);
            t.Data.stericSl(:,:,2) = NaN; t.Data.stericSl(3,2,2) = value;
            t.writeNative;
            t.verifyError(@() steric2slept(args{:}, truncation=4, ForceNew=true, BeQuiet=1), ...
                'ULMO:steric2slept:FitFailed');
            t.verifyEqual(readbinary(path), original);
            t.verifyNumElements(dir(fullfile(fileparts(path), '*.nc')),1);
        end

        function failedWritePreservesCompleteCache(t)
            % Inject a write failure after coordinate creation to exercise cleanup.
            args = {'EN4c14', t.Polygon, 0};
            steric2slept(args{:}, BeQuiet=1);
            path = t.cachePath;
            original = readbinary(path);
            folder = fullfile(t.Folder, 'writeFailure'); mkdir(folder);
            writelines(["function write2nc(varargin)"; ...
                "error('ULMO:test:WriteFailure','Injected NetCDF write failure.');"; "end"], ...
                fullfile(folder,'write2nc.m'));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
            t.verifyError(@() steric2slept(args{:}, ForceNew=true, BeQuiet=1), 'ULMO:test:WriteFailure');
            t.verifyEqual(readbinary(path), original);
            t.verifyNumElements(dir(fullfile(fileparts(path), '*.nc')),1);
        end

        function unavailableFieldsAndZeroDeep(t)
            t.Data.deepStericSl(:) = 0;
            t.Data.stericSl = t.Data.shallowStericSl;
            t.Data.deepHalostericSl(:) = 0;
            t.Data.halostericSl = t.Data.shallowHalostericSl;
            t.writeNative;
            [deep,~,fit,~,h] = steric2slept('EN4c14',t.Polygon,0, Depth='deep', BeQuiet=1);
            t.verifyEqual(deep,zeros(60,1)); t.verifyEqual(h,zeros(60,1));
            t.verifyGreaterThan(min(fit.nObs),0);
            % A distinct product with only total still projects without a halo diagnostic.
            data = struct('lon',t.Data.lon,'lat',t.Data.lat,'dates',t.Data.dates,'stericSl',t.Data.stericSl);
            path = stericfilepath('SIO'); mkdir(fileparts(path)); save(path,'-struct','data','-v7.3');
            c = steric2slept('SIO',t.Polygon,0, BeQuiet=1);
            t.verifySize(c,[60,1]);
            t.verifyError(@() requestHalo('SIO',t.Polygon), 'ULMO:steric2slept:UnsupportedField');
            t.verifyError(@() steric2slept('SIO',t.Polygon,0,Depth='deep',BeQuiet=1), ...
                'ULMO:steric2slept:UnsupportedField');
        end

        function irregularDescendingCyclicGrid(t)
            % Native means and fits must exclude a duplicated cyclic endpoint.
            schema = stericfieldschema(whos('-file',t.NativePath));
            t.Data.lat = [80;45;10;-5;-30;-80];
            t.Data.lon = [t.Data.lon; t.Data.lon(1)+360];
            for name = schema.gridNames(schema.available).'
                value = t.Data.(name);
                t.Data.(name) = cat(2,value,value(:,1,:));
            end
            t.writeNative;
            [c,~,fit] = steric2slept('EN4c14',t.Polygon,0,BeQuiet=1);
            [x,y] = ndgrid(t.Data.lon(1:end-1),flip(t.Data.lat));
            mask = domainmask(t.Polygon,x(:,1),y(1,:),Format='ndgrid',SaveData=false,BeQuiet=true);
            lat = flip(t.Data.lat); d = diff(lat);
            edges = [lat(1)-d(1)/2; (lat(1:end-1)+lat(2:end))/2; lat(end)+d(end)/2];
            w = repmat(diff(sind(max(-90,min(90,edges)))).',12,1);
            values = permute(t.Data.stericSl(end:-1:1,1:end-1,:),[2 1 3]);
            values = reshape(values,[],60);
            expected = w(mask).'*values(mask(:),:)/sum(w(mask));
            t.verifyEqual(c,expected.',AbsTol=1e-11);
            t.verifyEqual(fit.nObs,repmat(nnz(mask),60,1));
            t.verifyEqual(ncreadatt(t.cachePath,'/','weighting'),'spherical-midpoint-bands');
        end

        function textOptionsAcceptCharAndString(t)
            % Every text selector accepts either representation without conversion.
            args = {'EN4c14', t.Polygon, 0, days(7)};
            [c, dates, ~, ~, h] = steric2slept(args{:}, Type='halosteric', ...
                Depth='shallow', Interpolation='LiNeAr', OutputFormat='traditional', ...
                Unit='mm', TimeFormat='datenum', BeQuiet=true);
            [actual, actualDates, ~, ~, actualH] = steric2slept(args{:}, Type="halosteric", ...
                Depth="shallow", Interpolation="linear", OutputFormat="traditional", ...
                Unit="mm", TimeFormat="datenum", BeQuiet=1, SaveData=0);
            t.verifyEqual(actual, c);
            t.verifyEqual(actualDates, dates);
            t.verifyEqual(actualH, h);
            t.verifySize(h, [1, numel(dates)]);
            t.verifyClass(dates, 'double');
        end

        function interpolationValidationBelongsToInterpolator(t)
            % Native sampling does not use the method; fixed steps delegate it.
            args = {'EN4c14', t.Polygon, 0};
            steric2slept(args{:}, Interpolation='unused-method', BeQuiet=1);
            expected = steric2slept(args{:}, days(7), Interpolation='linear', BeQuiet=1);
            actual = steric2slept(args{:}, days(7), Interpolation='LiN', BeQuiet=1);
            t.verifyEqual(actual, expected, AbsTol=1e-12);
            expected = steric2slept(args{:}, days(7), Interpolation='mean', BeQuiet=1);
            actual = steric2slept(args{:}, days(7), Interpolation="MEAN", BeQuiet=1);
            t.verifyEqual(actual, expected);
            try
                interptemporal(t.Data.dates, t.Data.stericSl, days(7), 'invalid-method');
                t.assertFail('The interpolation function must reject this method.');
            catch exception
                expectedId = exception.identifier;
            end
            t.verifyError(@() steric2slept(args{:}, days(7), ...
                Interpolation='invalid-method', BeQuiet=1), expectedId);
        end

        function textInputsRejectNontextAndArrays(t)
            % Reject values that a typed string argument could silently convert.
            for name = {'Type', 'Depth', 'Interpolation', 'OutputFormat', 'Unit', 'TimeFormat'}
                for value = {42, ["linear", "nearest"]}
                    t.verifyError(@() steric2slept('EN4c14', t.Polygon, 0, ...
                        name{1}, value{1}), 'MATLAB:validators:mustBeTextScalar');
                end
            end
        end

        function timeRangeUsesLoaderValidator(t)
            for limits = {t.Data.dates([5 5]), t.Data.dates([6 5]), t.Data.dates(5)}
                t.verifyError(@() steric2slept('EN4c14', t.Polygon, 0, [], limits{1}), ...
                    'ULMO:notValidTimeRange');
            end
        end

        function nativeQuarterlyAndCacheIdentity(t)
            % Native dates need not be monthly; fractional steps must not collide.
            schema = stericfieldschema(whos('-file',t.NativePath));
            t.Data.dates = datetime(2020,[2;5;8;11],15);
            for name = schema.gridNames(schema.available).'
                t.Data.(name) = t.Data.(name)(:,:,1:4);
            end
            t.writeNative;
            [~,dates] = steric2slept('EN4c14',t.Polygon,0,BeQuiet=1);
            t.verifyEqual(dates,t.Data.dates);
            steric2slept('EN4c14',t.Polygon,0,days(90),BeQuiet=1);
            steric2slept('EN4c14',t.Polygon,0,90.25,BeQuiet=1);
            files = dir(fullfile(fileparts(t.NativePath),'SlepianExpansions','*.nc'));
            t.verifyNumElements(files,3);
            steric2slept('EN4c14',t.Polygon,0,days(90.25),BeQuiet=1);
            t.verifyNumElements(dir(fullfile(files(1).folder,'*.nc')),3);
        end
    end
    methods
        function writeNative(t)
            data = t.Data;
            save(t.NativePath,'-struct','data','-v7.3');
        end
        function path = cachePath(t)
            files = dir(fullfile(fileparts(t.NativePath),'SlepianExpansions','*.nc'));
            t.assertNumElements(files,1);
            path = fullfile(files.folder,files.name);
        end
        function [c,fit] = fitReference(t,grid,lon,lat,basis)
            [x,y] = ndgrid(lon,lat);
            mask = domainmask(t.Polygon,lon,lat,Format='ndgrid',exteriorBuffer=90,SaveData=false,BeQuiet=true);
            values = reshape(grid,[],size(grid,3));
            [c,~,~,fit] = xyzs2slep(values(mask(:),:),x(mask),y(mask),t.Polygon,1, ...
                GVN=basis,truncation=size(basis{1},2),weight=cosd(y(mask)),missingPolicy='omit');
            c = c.';
        end
    end
end

function bytes = readbinary(path)
    fid = fopen(path,'rb'); cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid,Inf,'*uint8');
end

function files = snapshot(folder)
    entries = dir(fullfile(folder,'**','*'));
    entries = entries(~[entries.isdir]);
    files = cell(numel(entries),2);
    for iFile = 1:numel(entries)
        path = fullfile(entries(iFile).folder,entries(iFile).name);
        files{iFile,1} = erase(path,folder);
        files{iFile,2} = readbinary(path);
    end
end

function requestHalo(product,domain)
    [~,~,~,~,h] = steric2slept(product,domain,0,BeQuiet=1); %#ok<ASGLU>
end

function expected = stericgmslExpected(data)
    weights = repmat(cosd(data.lat),1,numel(data.lon));
    expected = reshape(sum(data.deepHalostericSl.*weights,[1 2])/sum(weights,'all'),[],1);
end
