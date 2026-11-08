classdef stericUnsupportedFieldsTest < matlab.unittest.TestCase
    % Verify loader support using prepared aggregates, without generating data.
    methods (Test)
        function emptyLayersAndCaches(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            productFolder = fullfile(folder.Folder, 'STERIC', 'SIO'); mkdir(productFolder);
            native = fullfile(productFolder, 'SIO-StericSeaLevel.mat');
            lon = [0, 120, 240]; lat = [-30; 30];
            dates = datetime(2005, 1, 16) + calmonths((0:3)');
            stericSl = 15 * ones(2, 3, numel(dates));
            shallowStericSl = stericSl; halostericSl = stericSl;
            shallowHalostericSl = halostericSl;
            thermostericSl = stericSl; shallowThermostericSl = stericSl;
            deepStericSl = []; deepHalostericSl = []; deepThermostericSl = [];
            save(native, 'dates', 'lon', 'lat', 'stericSl', 'shallowStericSl', ...
                'deepStericSl', 'halostericSl', 'shallowHalostericSl', ...
                'deepHalostericSl', 'thermostericSl', 'shallowThermostericSl', ...
                'deepThermostericSl', '-v7.3');
            testCase.verifyFalse(stericproductinfo('SIO').hasDeep);
            types = {'total', 'thermosteric', 'halosteric', 'no halosteric drift', ...
                'non-thermosteric', 'non-halosteric'};
            for iType = 1:numel(types)
                testCase.verifyError(@() steric2lonlatt('SIO', [], [], [], ...
                    Depth='deep', Type=types{iType}, BeQuiet=true), ...
                    'ULMO:steric2lonlatt:UnsupportedField');
            end
            shallow = steric2lonlatt('SIO', [], [], [], Depth='shallow', BeQuiet=true);
            testCase.verifyEqual(shallow, stericSl);
            % Interpolated field caches preserve empty arrays.
            full = steric2lonlatt('SIO', days(20), [], [], BeQuiet=true);
            testCase.verifyEqual(full, 15 * ones(size(full)));
            cache = fullfile(productFolder, 'SIO-StericSeaLevel-T20days-interp-linear.mat');
            cached = load(cache);
            testCase.verifyEmpty(cached.deepStericSl);
            % Native declarations override stale interpolated placeholders.
            deepStericSl = nan(size(cached.stericSl));
            save(cache, 'deepStericSl', '-append');
            testCase.verifyError(@() steric2lonlatt('SIO', days(20), [], [], ...
                Depth='deep', BeQuiet=true), 'ULMO:steric2lonlatt:UnsupportedField');
            residual = steric2lonlatt('SIO', days(20), [], [], ...
                Depth='shallow', Type='non-halosteric', BeQuiet=true);
            testCase.verifyEqual(residual, zeros(size(residual)));
            % Epoch caches are selected-component arrays, following main.
            monthly = steric2lonlatt('SIO', 'midmonth', [], [], BeQuiet=true);
            testCase.verifyEqual(monthly, 15 * ones(size(monthly)));
            epochCache = fullfile(productFolder, 'SIO-StericSeaLevel-Tmidmonth-linear-deep-total.mat');
            steric = monthly;
            save(epochCache, 'lon', 'lat', 'dates', 'steric');
            testCase.verifyError(@() steric2lonlatt('SIO', 'midmonth', [], [], ...
                Depth='deep', BeQuiet=true), 'ULMO:steric2lonlatt:UnsupportedField');
        end

        function newProductIdentifiers(testCase)
            for product = ["JAMSTEC", "ISAS", "IPRC", "IAPv4"]
                testCase.verifyEqual(string(StericProduct.(product)), product);
            end
        end
    end
end
