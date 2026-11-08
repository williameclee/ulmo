%% STERICFIELDSCHEMATEST - Test shared steric metadata and caller requirements.
%
% Syntax
%   results = runtests('unittests/stericFieldSchemaTest.m')
%
% Notes
%   Uses temporary MAT files under an isolated IFILES folder. Tests distinguish
%   absent, empty and zero-valued grids, residual support, invalid layer splits,
%   component-only products and estimator-specific numeric validation.
%   Production data is not read or changed.
%
% Created by
%   2026/11/05, OpenAI Codex
%
% Last modified
%   2026/11/08, OpenAI Codex

classdef stericFieldSchemaTest < matlab.unittest.TestCase
    properties
        nativePath
        fixtureData
    end

    methods (TestMethodSetup)
        function prepareNativeFile(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            testCase.nativePath = stericfilepath('EN4c14');
            mkdir(fileparts(testCase.nativePath));
            testCase.fixtureData = struct('lon', [0, 120, 240], 'lat', [-30; 30], ...
                'dates', datetime(2020, 1, (1:3)'), 'stericSl', 6 * ones(2, 3, 3), ...
                'thermostericSl', 2 * ones(2, 3, 3), 'halostericSl', ones(2, 3, 3));
        end
    end

    methods (Test)
        function namesWithoutMetadata(testCase)
            schema = stericfieldschema();
            testCase.verifyEqual(schema.types, ...
                ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"]);
            testCase.verifyEqual(schema.depths, ["full", "shallow", "deep"]);
            testCase.verifySize(schema.gridNames, [3, 3]);
            testCase.verifySize(schema.meanNames, [5, 3]);
            testCase.verifyEqual(schema.gridNames(3, 2), "shallowHalostericSl");
            testCase.verifyEqual(schema.meanNames(5, 3), "deepNonHalostericGmsl");
            testCase.verifyFalse(any(schema.available, 'all'));
            testCase.verifyFalse(any(schema.supported, 'all'));
            testCase.verifyFalse(schema.incompleteSplit || schema.incompatibleSplits);
        end

        function missingEmptyAndZeroGrids(testCase)
            data = rmfield(testCase.fixtureData, 'thermostericSl');
            data.halostericSl = [];
            data.shallowStericSl = data.stericSl;
            data.deepStericSl = zeros(size(data.stericSl));
            data.halostericGmsl = ones(3, 1);
            schema = testCase.writeSchema(data);
            testCase.verifyEqual(schema.available, logical([1, 1, 1; 0, 0, 0; 0, 0, 0]));
            testCase.verifyEqual(schema.supported, logical([1, 1, 1; zeros(4, 3)]));
            testCase.verifyFalse(schema.incompleteSplit || schema.incompatibleSplits);
            info = stericproductinfo('EN4c14');
            testCase.verifyTrue(info.hasDeep && info.hasShallowField && info.hasDeepField);
            testCase.verifyFalse(info.hasHalosteric);
        end

        function completeSplitsSupportAllMeans(testCase)
            data = testCase.fixtureData;
            components = {'Steric', 'Thermosteric', 'Halosteric'};
            for iComponent = 1:numel(components)
                component = components{iComponent};
                fullName = [lower(component(1)), component(2:end), 'Sl'];
                data.(['shallow', component, 'Sl']) = data.(fullName) / 4;
                data.(['deep', component, 'Sl']) = 3 * data.(fullName) / 4;
            end
            schema = testCase.writeSchema(data);
            testCase.verifyTrue(all(schema.available, 'all'));
            testCase.verifyTrue(all(schema.supported, 'all'));
            testCase.verifyTrue(all(schema.split));
            testCase.verifyFalse(schema.incompleteSplit || schema.incompatibleSplits);
        end

        function residualsRequireTotal(testCase)
            schema = testCase.writeSchema(rmfield(testCase.fixtureData, 'stericSl'));
            testCase.verifyEqual(schema.supported(:, 1), logical([0; 1; 1; 0; 0]));
            schema = testCase.writeSchema(rmfield(testCase.fixtureData, 'halostericSl'));
            testCase.verifyEqual(schema.supported(:, 1), logical([1; 1; 0; 1; 0]));
            schema = testCase.writeSchema(rmfield(testCase.fixtureData, 'thermostericSl'));
            testCase.verifyEqual(schema.supported(:, 1), logical([1; 0; 1; 0; 1]));
        end

        function partialSplitsRetainCallerErrors(testCase)
            components = {'Steric', 'Thermosteric', 'Halosteric'};
            for iComponent = 1:numel(components)
                data = testCase.fixtureData;
                data.(['shallow', components{iComponent}, 'Sl']) = ones(2, 3, 3);
                schema = testCase.writeSchema(data);
                testCase.verifyTrue(schema.incompleteSplit);
                testCase.verifyCallerSchemaErrors;
            end
        end

        function incompatibleSplitsRetainCallerErrors(testCase)
            components = {'Steric', 'Halosteric'};
            for iComponent = 1:numel(components)
                data = testCase.fixtureData;
                data.(['shallow', components{iComponent}, 'Sl']) = ones(2, 3, 3);
                data.(['deep', components{iComponent}, 'Sl']) = ones(2, 3, 3);
                schema = testCase.writeSchema(data);
                testCase.verifyFalse(schema.incompleteSplit);
                testCase.verifyTrue(schema.incompatibleSplits);
                testCase.verifyCallerSchemaErrors;
            end
        end

        function componentOnlyCallerRequirements(testCase)
            data = rmfield(testCase.fixtureData, {'stericSl', 'halostericSl'});
            testCase.writeSchema(data);
            testCase.verifyError(@() stericproductinfo('EN4c14'), 'ULMO:LoadData:VariableNotFound');
            gmsl = steric2gmsl('EN4c14', 'thermosteric', SaveData=false, BeQuiet=true);
            testCase.verifyEqual(gmsl, 2 * ones(3, 1), AbsTol=1e-12);
            grid = steric2lonlatt('EN4c14', Type='thermosteric', SaveData=false, BeQuiet=true);
            testCase.verifyEqual(grid, data.thermostericSl);
        end

        function numericValidationRemainsWithEstimator(testCase)
            grids = {uint8(ones(2, 3)), complex(ones(2, 3), ones(2, 3)), sparse(ones(2, 3))};
            for iGrid = 1:numel(grids)
                data = testCase.fixtureData;
                data.dates = data.dates(1);
                data.stericSl = grids{iGrid};
                data.thermostericSl = data.thermostericSl(:, :, 1);
                data.halostericSl = data.halostericSl(:, :, 1);
                schema = testCase.writeSchema(data);
                testCase.verifyTrue(schema.available(1, 1));
                testCase.verifyWarningFree(@() stericproductinfo('EN4c14'));
                testCase.verifyError(@() steric2gmsl('EN4c14', 'total', ...
                    SaveData=false, BeQuiet=true), 'ULMO:steric2gmsl:InvalidSchema');
            end
        end
    end

    methods (Access=private)
        function schema = writeSchema(testCase, data)
            save(testCase.nativePath, '-struct', 'data', '-v7.3');
            schema = stericfieldschema(whos('-file', testCase.nativePath));
        end

        function verifyCallerSchemaErrors(testCase)
            testCase.verifyError(@() stericproductinfo('EN4c14'), 'ULMO:stericproductinfo:InvalidSchema');
            testCase.verifyError(@() steric2gmsl('EN4c14', 'total', SaveData=false, ...
                BeQuiet=true), 'ULMO:steric2gmsl:InvalidSchema');
            testCase.verifyError(@() steric2lonlatt('EN4c14', ForceNew=true, SaveData=false, ...
                BeQuiet=true), 'ULMO:steric2lonlatt:InvalidSchema');
        end
    end
end
