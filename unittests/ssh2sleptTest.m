classdef ssh2sleptTest < matlab.unittest.TestCase
    % Verify native grid ordering and full-record NetCDF cache behavior.
    properties
        Folder
        Dates
        Polygon = [0 -11; 32 -11; 32 11; 0 11; 0 -11]
        Basis = {eye(4), ones(4,1), 4}
    end
    methods (TestMethodSetup)
        function nativeFile(t)
            tmp = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Folder = fullfile(tmp.Folder, 'SSH');
            mkdir(t.Folder);
            original = getenv('IFILES');
            t.addTeardown(@() setenv('IFILES', original));
            setenv('IFILES', tmp.Folder);
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
            sshs(x > 32) = 900; % Only outside the target polygon in month one.
            t.Dates = dates;
            save(fullfile(t.Folder,'NASASSHMonthly.mat'),'sshs','dates','lon','lat','partialMonths','-v7.3');
        end
    end
    methods (Test)
        function subsetBuildsWholeRecord(t)
            args = {"NASASSHMonthly",t.Polygon,1};
            [c,d,fit,b] = ssh2slept(args{:},t.Dates([3 3]),GVN=t.Basis);
            t.verifyEqual(c,[3 2 3 4]*1000,'AbsTol',1e-8);
            t.verifyEqual(d,t.Dates(3));
            t.verifyEqual(b,t.Basis);
            t.verifyLessThan(fit.residualRms,1e-9);
            t.verifyEqual(fit.areaFraction,36/44,'AbsTol',0.005);
            file = dir(fullfile(t.Folder,'SlepianExpansions','*.nc'));
            t.assertNumElements(file,1);
            path = fullfile(file.folder,file.name);
            t.verifySize(ncread(path,'coefficients'),[4 4]);
            t.verifyEqual(ncreadatt(path,'coefficients','units'),'m');
            % Cache hits must work with no source cube, including a new unit.
            delete(fullfile(t.Folder,'NASASSHMonthly.mat'));
            [allC,dates,allFit,b2] = ssh2slept(args{:},GVN=t.Basis,unit="m");
            t.verifyEqual(dates,t.Dates);
            t.verifyTrue(all(isnan(allC(1,:))));
            t.verifyEqual(allC(2:4,:),[(2:4)',repmat([2 3 4],3,1)],'AbsTol',1e-10);
            t.verifyEqual(allFit.nObs,[0;60*44;64*36;64*44]);
            t.verifyEqual(allFit.partialMonth,[true;false;false;true]);
            t.verifyEqual(allFit.areaFraction([1 4]),[0;1],'AbsTol',1e-12);
            t.verifyEqual(b2,b);
            [none,when,emptyFit] = ssh2slept(args{:},datetime(2010,[1 2],1),GVN=t.Basis);
            t.verifySize(none,[0 4]);
            t.verifyEmpty(when);
            t.verifyEqual(height(emptyFit),0);
            ncwriteatt(path, 'time', 'units', 'days since 2000-01-01');
            t.verifyError(@() ssh2slept(args{:},GVN=t.Basis), 'ULMO:readshnc:InvalidTime');
        end
        function computedBasisAndNamedDomain(t)
            fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fixtures));
            for folder = ["GLMALPHA", "KERNELC", "KERNELCP", "LEGENDRE"]
                mkdir(fullfile(getenv('IFILES'), folder));
            end
            constants = fullfile(getenv('IFILES'), 'EARTHMODELS', 'CONSTANTS');
            mkdir(constants);
            Earth.Radius = 6371000;
            save(fullfile(constants, 'Earth.mat'), 'Earth');
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
            delete(fullfile(t.Folder, 'NASASSHMonthly.mat'));
            t.verifyEqual(ssh2slept("NASASSHMonthly", region, 0, t.Dates([2 4]), unit="m"), c);
        end

        function basisIdentityAndRebuild(t)
            args = {"NASASSHMonthly",t.Polygon,1};
            original = ssh2slept(args{:},GVN=t.Basis);
            flipped = t.Basis;
            flipped{1}(:,1) = -flipped{1}(:,1);
            changed = ssh2slept(args{:},GVN=flipped);
            t.verifyEqual(changed(2:end,1),-original(2:end,1),'AbsTol',1e-8);
            t.verifyNumElements(dir(fullfile(t.Folder,'SlepianExpansions','*.nc')),2);
            data = matfile(fullfile(t.Folder,'NASASSHMonthly.mat'),'Writable',true);
            data.sshs(:,:,4) = 2*data.sshs(:,:,4);
            rebuilt = ssh2slept(args{:},t.Dates([4 4]),GVN=t.Basis,ForceNew=true);
            t.verifyEqual(rebuilt,2*original(4,:),'AbsTol',1e-8);
            % A failed refresh preserves the complete previous cache.
            bad = nan(size(data,'sshs',1),size(data,'sshs',2));
            bad(1) = 1;
            data.sshs(:,:,4) = bad;
            t.verifyError(@() ssh2slept(args{:},GVN=t.Basis,ForceNew=true),'ULMO:ssh2slept:FitFailed');
            t.verifyEqual(ssh2slept(args{:},t.Dates([4 4]),GVN=t.Basis),rebuilt);
        end
    end
end
