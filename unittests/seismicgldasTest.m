%% SEISMICGLDASTEST - Checks native storage integration and interval coverage.
% Synthetic Noah-format NetCDF files test analytic harmonics, unequal month
% overlaps, missing components, units, masks, coverage and save-free reads.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

classdef seismicgldasTest < matlab.unittest.TestCase
    properties
        Folder
        Files
        Epochs = [datetime(2004, 1, 1), datetime(2004, 2, 1); ...
            datetime(2004, 2, 1), datetime(2004, 3, 1)]
        Variables = {'SoilMoi0_10cm_inst', 'SoilMoi10_40cm_inst', ...
            'SoilMoi40_100cm_inst', 'SoilMoi100_200cm_inst', 'CanopInt_inst', 'SWE_inst'}
    end
    methods (TestMethodSetup)
        function setup(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root, IncludingSubfolders=true));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
            old = getenv('GLDASDATA');
            testCase.addTeardown(@() setenv('GLDASDATA', old));
            setenv('GLDASDATA', testCase.Folder);
            oldRoot = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', oldRoot));
            setenv('IFILES', testCase.Folder);
            constants = fullfile(testCase.Folder, 'EARTHMODELS', 'CONSTANTS');
            mkdir(constants);
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(constants, 'Earth.mat'), 'Earth');
            mkdir(fullfile(testCase.Folder, 'LOVENUMS'));
            loadinggravitationalpotential = zeros(3, 1);
            save(fullfile(testCase.Folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), 'loadinggravitationalpotential');
            testCase.Files = strings(2, 1);
            for iMonth = 1:2
                path = fullfile(testCase.Folder, sprintf('GLDAS_NOAH10_M.A20040%d.021.nc4', iMonth));
                testCase.Files(iMonth) = path;
                nccreate(path, 'lon', Dimensions={'lon', 360}, Format='netcdf4');
                ncwrite(path, 'lon', (-179.5:179.5).');
                nccreate(path, 'lat', Dimensions={'lat', 150});
                ncwrite(path, 'lat', (-59.5:89.5).');
                for iVariable = 1:numel(testCase.Variables)
                    name = testCase.Variables{iVariable};
                    nccreate(path, name, Dimensions={'lon', 360, 'lat', 150, 'time', 1}, ...
                        Datatype='single', FillValue=single(-9999));
                    ncwriteatt(path, name, 'units', 'kg m-2');
                    ncwrite(path, name, single(ones(360, 150) * (2 * iMonth - 1)));
                end
            end
        end
    end
    methods (Test)
        function uniformStorageMatchesAnalyticZonalCoefficient(testCase)
            p = testCase.read(testCase.Epochs);
            testCase.verifyEqual(p(:, 2, 3), [1; 3] * 6 * sqrt(3) / 16, AbsTol=1e-12);
            testCase.verifyEqual(p(:, 1, 3:4), zeros(2, 1, 2));
            testCase.verifyLessThan(max(abs(p(:, [3, 5, 6], 3:4)), [], 'all'), 1e-12);
            testCase.verifyNumElements(dir(fullfile(testCase.Folder, '*.nc4')), 2);
            testCase.verifyEmpty(dir(fullfile(testCase.Folder, '*.nc')));
        end
        function unitsUseTheRequestedGravityReference(testCase)
            sd = testCase.read(testCase.Epochs);
            radius = 6371000;
            pot = seismicgldas(testCase.Epochs, 2, 'POT', 'ISSM', radius);
            grav = seismicgldas(testCase.Epochs, 2, 'GRAV', 'ISSM', radius);
            testCase.verifyEqual(pot(:, 4, 3), sd(:, 4, 3) * 3 / (5517 * 5), RelTol=1e-12);
            testCase.verifyEqual(grav(:, :, 3:4), pot(:, :, 3:4) / radius, AbsTol=1e-20);
            testCase.verifyEqual(grav(:, :, 1:2), sd(:, :, 1:2));
        end
        function overlapUsesDaysAndLeapMonth(testCase)
            monthly = testCase.read(testCase.Epochs);
            epoch = [datetime(2004, 1, 16), datetime(2004, 2, 16)];
            actual = testCase.read(epoch);
            testCase.verifyEqual(actual(:, :, 3:4), ...
                (16 * monthly(1, :, 3:4) + 15 * monthly(2, :, 3:4)) / 31, AbsTol=1e-12);
        end
        function missingMonthFailsInsteadOfExtrapolating(testCase)
            testCase.verifyError(@() testCase.read([datetime(2004, 2, 15), datetime(2004, 3, 2)]), ...
                'ULMO:seismicgldas:Coverage');
        end
        function duplicateRequiredMonthFails(testCase)
            mkdir(fullfile(testCase.Folder, 'duplicate'));
            copyfile(testCase.Files(1), fullfile(testCase.Folder, 'duplicate'));
            testCase.verifyError(@() testCase.read(testCase.Epochs), 'ULMO:seismicgldas:DuplicateMonth');
        end
        function incompatibleUnitsFail(testCase)
            ncwriteatt(testCase.Files(1), 'SWE_inst', 'units', 'm');
            testCase.verifyError(@() testCase.read(testCase.Epochs), 'ULMO:seismicgldas:Unit');
        end
        function changingMissingMaskFails(testCase)
            ncwrite(testCase.Files(2), 'SWE_inst', single(-9999), [1, 1, 1]);
            testCase.verifyError(@() testCase.read(testCase.Epochs), 'ULMO:seismicgldas:MissingValues');
        end
        function nativeMaskHasZeroContribution(testCase)
            original = testCase.read(testCase.Epochs);
            for iMonth = 1:2
                for iVariable = 1:numel(testCase.Variables)
                    ncwrite(testCase.Files(iMonth), testCase.Variables{iVariable}, single(-9999), [1, 1, 1]);
                end
            end
            actual = testCase.read(testCase.Epochs);
            testCase.verifyTrue(all(isfinite(actual), 'all'));
            testCase.verifyNotEqual(actual, original);
            testCase.verifyEqual(actual(2, :, 3:4), 3 * actual(1, :, 3:4), AbsTol=1e-12);
        end
        function missingValueAttributeIsHonoured(testCase)
            for iMonth = 1:2
                for iVariable = 1:numel(testCase.Variables)
                    name = testCase.Variables{iVariable};
                    ncwriteatt(testCase.Files(iMonth), name, 'missing_value', single(-7777));
                    ncwrite(testCase.Files(iMonth), name, single(-7777), [1, 1, 1]);
                end
            end
            p = testCase.read(testCase.Epochs);
            testCase.verifyLessThan(max(abs(p(:, :, 3:4)), [], 'all'), 10);
        end
        function cosineAndSineOrderingMatchesLongitude(testCase)
            longitude = deg2rad((-179.5:179.5).');
            for iVariable = 1:numel(testCase.Variables)
                ncwrite(testCase.Files(1), testCase.Variables{iVariable}, ...
                    single(repmat(cos(longitude) + 2 * sin(longitude), 1, 150)));
            end
            p = testCase.read(testCase.Epochs(1, :));
            testCase.verifyGreaterThan(p(1, 3, 3), 0);
            testCase.verifyEqual(p(1, 3, 4), 2 * p(1, 3, 3), AbsTol=1e-7);
            % Exact integral of cos(latitude)^2 between -60 and 90 degrees.
            south = -pi / 3;
            integral = pi / 4 - south / 2 - sin(2 * south) / 4;
            expected = 6 * sqrt(3) * integral / 4 * sinc(1 / 360);
            testCase.verifyEqual(p(1, 3, 3), expected, AbsTol=1e-7);
        end
    end
    methods
        function p = read(~, epochs)
            p = seismicgldas(epochs, 2, 'SD', 'ISSM', 6378136.3);
        end
    end
end
