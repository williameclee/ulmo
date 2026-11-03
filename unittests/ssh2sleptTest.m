classdef ssh2sleptTest < matlab.unittest.TestCase
    % Verify native grid ordering and full-record NetCDF cache behavior.
    properties
        Folder
        Dates
        Polygon = [0 -11; 32 -11; 32 11; 0 11; 0 -11]
    end
    methods (TestMethodSetup)
        function nativeFile(t)
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
            dates = datetime(2020,1,1) + days([15.5;45.5;75.5;105]);
            partialMonths = [true;false;false;true];
            [x,y] = ndgrid(lon,lat);
            field = sqrt(3)*(2*cosd(y).*cosd(x) + 3*sind(y) + 4*cosd(y).*sind(x));
            sshs = field + reshape(1:4,1,1,[]);
            sshs(:,:,1) = NaN;
            sshs(2:5,:,2) = NaN;
            sshs(:,3:10,3) = NaN;

            t.Dates = dates;
            save(fullfile(t.Folder,'NASASSHMonthly.mat'),'sshs','dates','lon','lat','partialMonths','-v7.3');
        end
    end
    methods (Test)
        function degreeDependentExteriorFit(t)
            % A constant basis isolates which observed cells enter the fit.
            G = zeros(61^2,1);
            G(1) = 1;
            args = {"NASASSHMonthly",t.Polygon,60};
            data = matfile(fullfile(t.Folder,'NASASSHMonthly.mat'),'Writable',true);
            data.sshs(end-3:end,:,1) = 900*ones(4,size(data,'sshs',2));
            [c,~,fit] = ssh2slept(args{:}, GVN={G,1,1}, unit="m", SaveData=false);
            data = load(fullfile(t.Folder,'NASASSHMonthly.mat'));
            % The first month only has SSH at 32.25, 32.75, 33.25 and 33.75 E.
            t.verifyEqual(fit.nObs(1),3*numel(data.lat));
            t.verifyEqual(fit.areaFraction(1),0);
            t.verifyEqual(c(1),900,'AbsTol',1e-8);
            t.verifyLessThanOrEqual(fit.areaFraction,ones(4,1)+eps);
        end

        function exteriorCacheAndTopology(t)
            % A ring crossing the grid seam and a separate region with a hole.
            xy = [0 -4;4 -4;4 4;0 4;0 -4;NaN NaN; ...
                50 -7;65 -7;65 7;50 7;50 -7;NaN NaN; ...
                53 -4;53 4;62 4;62 -4;53 -4];
            lon = 0:359; lat = -9:9;
            m = domainmask(xy,lon,lat,exteriorBuffer=1.5,BeQuiet=true);
            at = @(x) m(lat==0,lon==x);
            t.verifyTrue(at(359) && at(5) && at(54) && at(61));
            t.verifyFalse(at(358) || at(6) || at(57) || at(30));
            [cached,x,y] = domainmask(xy,lon,lat,exteriorBuffer=1.5,Format='ndgrid',BeQuiet=true);
            t.verifyEqual(cached,m.');
            [expectedX,expectedY] = ndgrid(lon,lat);
            t.verifyEqual(x,expectedX);
            t.verifyEqual(y,expectedY);
            narrow = domainmask(xy,lon,lat,exteriorBuffer=0.4,BeQuiet=true);
            t.verifyFalse(narrow(lat==0,lon==5));
            file = dir(fullfile(getenv('IFILES'),'DOMAINMASK','*-E0.4-H*.mat'));
            t.assertNumElements(file,1);
        end

        function roundedMaskNameChecksExactWidth(t)
            for width = [1/3, 1/3+1e-8]
                domainmask(t.Polygon,0:33,-10:10,exteriorBuffer=width,BeQuiet=true);
                mask = dir(fullfile(getenv('IFILES'),'DOMAINMASK',['*-E',num2str(width),'-*.mat']));
                t.assertNumElements(mask,1);
                cached = load(fullfile(mask.folder,mask.name),'exteriorBuffer');
                t.verifyEqual(cached.exteriorBuffer,width);
            end
        end

        function subsetBuildsWholeRecord(t)
            args = {"NASASSHMonthly",t.Polygon,1};
            [c,d,fit,b] = ssh2slept(args{:},t.Dates([3 3]),truncation=4);
            t.verifyEqual(c*b{1}.',[3 2 3 4]*1000,'AbsTol',1e-8);
            t.verifyEqual(d,t.Dates(3));
            t.verifySize(b{1},[4 4]);
            t.verifyLessThan(fit.residualRms,1e-9);
            t.verifyEqual(fit.areaFraction,36/44,'AbsTol',0.005);
            file = dir(fullfile(t.Folder,'SlepianExpansions','*.nc'));
            t.assertNumElements(file,1);
            expectedName = sprintf('ssh2slept-NASASSHMonthly-%s-1-4-m.nc',hash(t.Polygon,'sha1'));
            t.verifyEqual(file.name,expectedName);
            path = fullfile(file.folder,file.name);
            info = ncinfo(path);
            t.verifyEmpty(intersect({info.Variables.Name}, {'G','V','N'}));
            t.verifyEmpty(intersect({info.Attributes.Name}, {'G','V','N'}));
            t.verifySize(ncread(path,'coefficients'),[4 4]);
            t.verifyEqual(ncreadatt(path,'coefficients','units'),'m');
            % Cache hits must work with no source cube, including a new unit.
            delete(fullfile(t.Folder,'NASASSHMonthly.mat'));
            [allC,dates,allFit,b2] = ssh2slept(args{:},truncation=4,unit="m");
            t.verifyEqual(dates,t.Dates);
            t.verifyTrue(all(isnan(allC(1,:))));
            t.verifyEqual(allC(2:4,:)*b2{1}.',[(2:4)',repmat([2 3 4],3,1)],'AbsTol',1e-10);
            t.verifyEqual(allFit.nObs,[0;64*44;68*36;68*44]);
            t.verifyEqual(allFit.partialMonth,[true;false;false;true]);
            t.verifyEqual(allFit.areaFraction([1 4]),[0;1],'AbsTol',1e-12);
            t.verifyEqual(b2,b);
            [none,when,emptyFit] = ssh2slept(args{:},datetime(2010,[1 2],1),truncation=4);
            t.verifySize(none,[0 4]);
            t.verifyEmpty(when);
            t.verifyEqual(height(emptyFit),0);
            ncwriteatt(path, 'time', 'units', 'days since 2000-01-01');
            t.verifyError(@() ssh2slept(args{:},truncation=4), 'ULMO:readshnc:InvalidTime');
        end
        function repeatedMasksKeepMonthsAligned(t)
            % Months two and four share a mask but have different coefficients.
            data = matfile(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'Writable', true);
            data.sshs(:,:,4) = data.sshs(:,:,2) + 2;
            [c, dates, fit] = ssh2slept("NASASSHMonthly", t.Polygon, 1, ...
                GVN={eye(4), ones(4,1), 4}, unit="m");
            t.verifyEqual(dates, t.Dates);
            t.verifyTrue(all(isnan([c(1,:), fit.condNum(1), fit.residualRms(1)])));
            t.verifyEqual(c(2:4,:), [(2:4)', repmat([2 3 4],3,1)], 'AbsTol', 1e-10);
            t.verifyEqual(fit.nObs, [0;64*44;68*36;64*44]);
            t.verifyEqual(fit.areaFraction(2), fit.areaFraction(4));
            t.verifyEqual(fit.condNum(2), fit.condNum(4));
            t.verifyLessThan(fit.residualRms(2:4), 1e-10);
        end

        function computedBasisAndNamedDomain(t)
            region = 'xyzs2slepregionfixture';
            [c, dates, fit, basis] = ssh2slept("NASASSHMonthly", region, 0, ...
                t.Dates([2 4]), SaveData=false, unit="m");
            t.verifyFalse(isfolder(fullfile(t.Folder, 'SlepianExpansions')));
            t.verifyEqual(dates, t.Dates(2:4));
            t.verifySize(c, [3 1]);
            t.verifySize(basis{1}, [1 1]);
            t.verifyGreaterThan(fit.nObs, 0);
            fromObject = ssh2slept("NASASSHMonthly", GeoDomain(region), 0, ...
                t.Dates([2 4]), unit="m");
            t.verifyEqual(fromObject, c, 'AbsTol', 1e-12);
            geo = GeoDomain(region);
            file = dir(fullfile(t.Folder,'SlepianExpansions','*.nc'));
            t.verifyEqual({file.name},{sprintf('ssh2slept-NASASSHMonthly-%s-0-1-m.nc',geo.Id)});
            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            t.verifyEqual(ssh2slept("NASASSHMonthly",region,0,t.Dates([2 4]), ...
                truncation=1,unit="m"),c);
            t.verifyEqual(ssh2slept("NASASSHMonthly", region, 0, t.Dates([2 4]), unit="m"), c);
        end

        function suppliedBasisBypassesProjectionCache(t)
            args = {"NASASSHMonthly", t.Polygon, 1, t.Dates([4 4])};
            custom = {eye(4), ones(4,1), 4};
            expected = [4 2 3 4];
            folder = fullfile(t.Folder, 'SlepianExpansions');
            [c, ~, ~, basis] = ssh2slept(args{:}, GVN=custom, unit="m");
            t.verifyEqual(basis, custom);
            t.verifyEqual(c, expected, 'AbsTol', 1e-10);
            t.verifyFalse(isfolder(folder));
            t.verifyEmpty(dir(fullfile(getenv('IFILES'), 'DOMAINMASK', '*.mat')));

            % Make the source differ from its standard cache at the same L and J.
            [cached, ~, ~, standard] = ssh2slept(args{:}, truncation=4, unit="m");
            data = matfile(fullfile(t.Folder, 'NASASSHMonthly.mat'), 'Writable', true);
            data.sshs(:,:,4) = 2 * data.sshs(:,:,4);
            fresh = ssh2slept(args{:}, GVN=standard, truncation=4, unit="m", ...
                ForceNew=false, SaveData=true);
            t.verifyEqual(fresh, 2 * cached, 'AbsTol', 1e-10);
            t.verifyEqual(ssh2slept(args{:}, truncation=4, unit="m"), cached);
            t.verifyNumElements(dir(fullfile(folder, '*.nc')), 1);
        end

        function separateTruncationsAndRebuild(t)
            args = {"NASASSHMonthly",t.Polygon,1};
            original = ssh2slept(args{:},truncation=4);
            % An explicit four-mode file cannot satisfy the one-mode default.
            [defaultC,~,~,expectedBasis] = ssh2slept(args{:});
            [~,~,~,cachedBasis] = ssh2slept(args{:});
            t.verifyEqual(cachedBasis,expectedBasis);
            t.verifySize(cachedBasis{1},[4 1]);
            t.verifySize(defaultC,[4 1]);
            t.verifyEqual(ssh2slept(args{:},truncation=1),defaultC);
            t.verifyNumElements(dir(fullfile(t.Folder,'SlepianExpansions','*.nc')),2);
            data = matfile(fullfile(t.Folder,'NASASSHMonthly.mat'),'Writable',true);
            data.sshs(:,:,4) = 2*data.sshs(:,:,4);
            rebuilt = ssh2slept(args{:},t.Dates([4 4]),truncation=4,ForceNew=true);
            t.verifyEqual(rebuilt,2*original(4,:),'AbsTol',1e-8);
            % A failed refresh preserves the complete previous cache.
            bad = nan(size(data,'sshs',1),size(data,'sshs',2));
            bad(1) = 1;
            data.sshs(:,:,4) = bad;
            t.verifyError(@() ssh2slept(args{:},truncation=4,ForceNew=true),'ULMO:ssh2slept:FitFailed');
            t.verifyEqual(ssh2slept(args{:},t.Dates([4 4]),truncation=4),rebuilt);
        end
    end
end
