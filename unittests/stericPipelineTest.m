%% STERICPIPELINETEST - Tests shared density, integration, and aggregation stages.
% Uses direct GSW density calculations and analytic density-ratio integrals
% to check scientific results independently of native-product readers.
%
% Syntax
%   results = runtests(fullfile('unittests', 'stericPipelineTest.m'));
%   assertSuccess(results);
%
% See also
%   computeStericDensity, computeStericDensityVar, computeStericSeaLevel,
%   computeStericClimatology, aggregateStericSeaLevel, convertStericSourceTSTest
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

classdef stericPipelineTest < matlab.unittest.TestCase

    properties (Access = private)
        Folder
        MonthPath
        ClimatologyPath
        ExpectedDensity
    end

    methods (TestClassSetup)
        % Use this checkout's shared helpers and restore the path after the suite.
        function addSourcePath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'aux')));
        end

    end

    methods (TestMethodSetup)
        % Create a three-latitude reference with density computed directly by GSW.
        function createFixture(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            tc.Folder = folder.Folder;
            lat = [-70; 0; 70]; lon = [20, 120]; depth = [0; 1200; 2400];
            pres = gsw_p_from_z(repmat(-depth', numel(lat), 1), lat);
            salinity = 35 * ones(3, 2, 3);
            consTemp = repmat(reshape([15, 7, 2], 1, 1, []), 3, 2, 1);
            date = datetime(2020, 1, 16);
            path = fullfile(folder.Folder, 'month.mat');
            save(path, 'salinity', 'consTemp', 'date', 'lat', 'lon', 'pres', 'depth');
            densityClim = nan(size(salinity));

            for k = 1:numel(depth)
                densityClim(:, :, k) = gsw_rho(salinity(:, :, k), consTemp(:, :, k), pres(:, k));
            end

            salinityClim = salinity; consTempClim = consTemp;
            clim = fullfile(folder.Folder, 'clim.mat');
            save(clim, 'salinityClim', 'consTempClim', 'densityClim');
            tc.MonthPath = path;
            tc.ClimatologyPath = clim;
            tc.ExpectedDensity = densityClim;
        end

    end

    methods (Test)
        % Check monthly density against the independent GSW reference.
        function testDensityAgainstGsw(tc)
            computeStericDensity(tc.MonthPath, BeQuiet = true);
            data = load(tc.MonthPath, 'density');
            tc.verifyEqual(double(data.density), tc.ExpectedDensity, AbsTol = 1e-4);
        end

        % With identical source/reference T/S, both component densities equal density.
        function testComponentDensitiesAgainstGsw(tc)
            computeStericDensityVar(tc.MonthPath, tc.ClimatologyPath, BeQuiet = true);
            data = load(tc.MonthPath);
            tc.verifyEqual(double(data.haloDensity), tc.ExpectedDensity, AbsTol = 1e-4);
            tc.verifyEqual(double(data.thermoDensity), tc.ExpectedDensity, AbsTol = 1e-4);
        end

        % A 1% density ratio gives 26 m total, 20 m shallow, and 6 m deep height.
        function testAnalyticLayerIntegralsAndMissingColumns(tc)
            tc.setAnalyticDensity();
            computeStericSeaLevel(tc.MonthPath, tc.ClimatologyPath, ...
                Bottom = 2600, BeQuiet = true, ForceNew = true);
            data = load(tc.MonthPath);
            tc.verifyEqual(data.stericSl(:, 1), 26 * ones(3, 1), AbsTol = 1e-9);
            tc.verifyEqual(data.shallowStericSl(:, 1), 20 * ones(3, 1), AbsTol = 1e-9);
            tc.verifyEqual(data.deepStericSl(:, 1), 6 * ones(3, 1), AbsTol = 1e-9);
            tc.verifyTrue(all(isnan(data.stericSl(:, 2))));
            computeStericSeaLevel(tc.MonthPath, tc.ClimatologyPath, ...
                Bottom = 2000, BeQuiet = true, ForceNew = true);
            data = load(tc.MonthPath);
            tc.verifyTrue(all(isnan(data.deepStericSl), 'all'));
        end

        % A one-month climatology must reproduce source values and missing columns.
        function testClimatologyPreservesFieldsAndMissingColumns(tc)
            tc.setAnalyticDensity();
            target = fullfile(tc.Folder, 'mean.mat');
            computeStericClimatology([datetime(2020, 1, 1), datetime(2020, 12, 31)], ...
                tc.Folder, {'month.mat'}, target, BeQuiet = true, ForceNew = true);
            source = load(tc.MonthPath);
            result = load(target);
            tc.verifyEqual(result.salinityClim, source.salinity);
            tc.verifyEqual(result.consTempClim, source.consTemp);
            tc.verifyEqual(result.densityClim, source.density);
        end

        % Reverse the input file order to check date sorting and NaN preservation.
        function testAggregationSortsMonthsAndPreservesMissingColumns(tc)
            tc.setAnalyticDensity();
            computeStericSeaLevel(tc.MonthPath, tc.ClimatologyPath, ...
                Bottom = 2000, BeQuiet = true, ForceNew = true);
            second = fullfile(tc.Folder, 'month2.mat');
            copyfile(tc.MonthPath, second);
            date = datetime(2020, 2, 15); save(second, 'date', '-append');
            target = fullfile(tc.Folder, 'aggregate.mat');
            aggregateStericSeaLevel(tc.Folder, {'month2.mat', 'month.mat'}, ...
                target, BeQuiet = true);
            result = load(target);
            tc.verifySize(result.stericSl, [3, 2, 2]);
            tc.verifyEqual(result.dates(:), [datetime(2020, 1, 16); date]);
            tc.verifyTrue(all(isnan(result.stericSl(:, 2, :)), 'all'));
        end

    end

    methods (Access = private)
        % Set rho=1000 and reference rho=1010, with one entirely missing column.
        function setAnalyticDensity(tc)
            density = 1000 * ones(3, 2, 3);
            density(:, 2, :) = NaN;
            haloDensity = density; thermoDensity = density;
            densityClim = 1010 * ones(3, 2, 3);
            save(tc.MonthPath, 'density', 'haloDensity', 'thermoDensity', '-append');
            save(tc.ClimatologyPath, 'densityClim', '-append');
        end

    end

end
