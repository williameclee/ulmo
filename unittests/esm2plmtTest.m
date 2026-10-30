%% ESM2PLMTTEST - Tests raw/monthly caching and ESM coefficient conventions.
% Synthetic fields exercise boundary selection, missing-data rejection and
% physical unit conversion independently of the downloaded ESM samples.
%
% Created by
%   2026/10/29
%
% Last modified
%   2026/10/30

classdef esm2plmtTest < matlab.unittest.TestCase
    properties
        Scratch
        RawFolder
        Start = datetime(2003, 1, 1)
    end
    methods (TestMethodSetup)
        function setup(testCase)
            testCase.Scratch = tempname;
            mkdir(testCase.Scratch);
            testCase.addTeardown(@() rmdir(testCase.Scratch, 's'));
            original = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', original));
            setenv('IFILES', testCase.Scratch);
            testCase.RawFolder = fullfile(testCase.Scratch, 'ESM', 'ESAESM2', 'raw', 'H');
            mkdir(testCase.RawFolder);
            mkdir(fullfile(testCase.Scratch, 'EARTHMODELS', 'CONSTANTS'));
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(testCase.Scratch, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
            mkdir(fullfile(testCase.Scratch, 'LOVENUMS'));
            loadinggravitationalpotential = zeros(181, 1);
            save(fullfile(testCase.Scratch, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential');
        end
    end
    methods (Test)
        function usesConfiguredDataRoot(testCase)
            % A duplicate outside the fixed raw folder must not be discovered.
            file = testCase.fixture(0, 1);
            copyfile(file, testCase.Scratch);
            range = [testCase.Start, testCase.Start + hours(6)];
            p = esm2plmt('H', [], range, 'Unit', 'GRAV');
            testCase.verifyEqual(p(1, 1, 3), 1e-10);
        end
        function componentFirstDefaultsAndCacheHitFormatting(testCase)
            testCase.fixture(0, 2);
            p = esm2plmt('H', 'Unit', 'GRAV');
            testCase.verifyEqual(p(1, 1, 3), 2e-10);
            rmdir(testCase.RawFolder, 's');
            [cached, ~, dates] = esm2plmt('H', [], [], 'Unit', 'GRAV', ...
                'OutputFormat', 'traditional', 'TimeFormat', 'datetime', 'Loutput', 1);
            testCase.verifyEqual(cached, squeeze(p(:, 1:3, :)));
            testCase.verifyEqual(dates, testCase.Start);
        end
        function oceanComponentDoesNotMatchOutputFormat(testCase)
            testCase.fullmonth(0, 2);
            esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            folder = fileparts(testCase.cachepath());
            copyfile(testCase.cachepath('midmonth'), ...
                fullfile(folder, 'ESAESM2_O_GRAV-Tmidmonth.mat'));
            [p, ~, dates] = esm2plmt('O', 'midmonth', 'Unit', 'GRAV');
            testCase.verifyEqual(p(1, 1, 3), 2e-10, 'RelTol', 2e-14);
            testCase.verifyEqual(dates, datenum(datetime(2003, 1, 16, 12, 0, 0)));
        end
        function rawRefreshInvalidatesAllEpochVariants(testCase)
            testCase.fixture(0, 1);
            folder = fileparts(testCase.cachepath());
            variants = ["ESAESM2_H_GRAV-Tmidmonth.mat", ...
                "ESAESM2_H_GRAV-Tmidmonth-mean.mat", "ESAESM2_H_GRAV-Tfuture.mat"];
            sentinel = 1;
            for variant = variants
                save(fullfile(folder, variant), 'sentinel');
            end
            unrelated = fullfile(folder, 'ESAESM2_O_GRAV-Tmidmonth.mat');
            save(unrelated, 'sentinel');
            esm2plmt('H', [], [], 'Unit', 'GRAV', 'ForceNew', true);
            for variant = variants
                testCase.verifyFalse(isfile(fullfile(folder, variant)));
            end
            testCase.verifyTrue(isfile(unrelated));
        end
        function requiresConfiguredDataRoot(testCase)
            range = [testCase.Start, testCase.Start + hours(6)];
            setenv('IFILES', '');
            testCase.verifyError(@() esm2plmt('O', [], range), 'ULMO:esm2plmt:MissingDataRoot');
            setenv('IFILES', testCase.Scratch);
            testCase.verifyError(@() esm2plmt('O', [], range), ...
                'ULMO:esm2plmt:MissingDataFolder');
        end
        function rawTimeRangeIsInclusive(testCase)
            for k = 0:4
                testCase.fixture(k, k + 1);
            end
            range = testCase.Start + hours([6, 18]);
            [p, s, dates] = esm2plmt('H', [], datenum(range), ...
                'Unit', 'GRAV');
            testCase.verifySize(p, [3, 6, 4]);
            testCase.verifyEqual(p(:, 1, 3), [2; 3; 4] * 1e-10, 'RelTol', 1e-14);
            testCase.verifyEqual(s(:, :, 1:2), p(:, :, 1:2));
            testCase.verifyTrue(all(isnan(s(:, :, 3:4)), 'all'));
            testCase.verifyEqual(dates, datenum(testCase.Start + hours([6; 12; 18])));
            testCase.verifyEqual(nargout('esm2plmt'), 3);
            cached = load(testCase.cachepath());
            testCase.verifySize(cached.esmPlmt, [5, 6, 4]);
        end
        function monthlyCachePrecedesSelectionAndTruncation(testCase)
            % Unequal month lengths and a leap day exercise true calendar bins.
            testCase.Start = datetime(2004, 1, 1);
            testCase.fullmonth(0, 2);
            testCase.fullmonth(1, 7);
            range = [datetime(2004, 2, 15), datetime(2004, 2, 16)];
            [p, s, dates] = esm2plmt('H', 'midmonth', range, ...
                'Unit', 'GRAV', 'Loutput', 1, 'TimeFormat', 'datetime', ...
                'OutputFormat', 'traditional');
            testCase.verifySize(p, [3, 4]);
            testCase.verifyEqual(p(:, 3), 7e-10 * (1:3)', 'RelTol', 2e-14);
            testCase.verifyEqual(dates, datetime(2004, 2, 15, 12, 0, 0));
            testCase.verifyTrue(all(isnan(s(:, 3:4)), 'all'));
            cached = load(testCase.cachepath('midmonth'));
            testCase.verifySize(cached.esmPlmt, [2, 6, 4]);
            testCase.verifyEqual(cached.esmPlmt(:, 1, 3), [2; 7] * 1e-10, 'RelTol', 2e-14);
            testCase.verifyEqual(cached.dates, ...
                [datetime(2004, 1, 16, 12, 0, 0); datetime(2004, 2, 15, 12, 0, 0)]);
            testCase.verifyEqual(sort(string(fieldnames(cached))), ...
                sort(["dates"; "esmPlmt"; "esmStdPlmt"]));
            raw = load(testCase.cachepath());
            testCase.verifySize(raw.esmPlmt, [240, 6, 4]);
            % The monthly cache remains usable without either underlying input.
            delete(testCase.cachepath());
            rmdir(testCase.RawFolder, 's');
            [allFields, ~, d] = esm2plmt('H', "midmonth", [], 'Unit', 'GRAV', ...
                'TimeFormat', 'datetime');
            testCase.verifyEqual(allFields, cached.esmPlmt);
            testCase.verifyEqual(d, cached.dates);
        end
        function rawCacheAvoidsSourceReads(testCase)
            testCase.fixture(0, 1);
            testCase.fixture(2, 3); % Raw data may have gaps.
            p = esm2plmt('H', [], [], 'Unit', 'GRAV', 'Loutput', 1);
            rmdir(testCase.RawFolder, 's');
            full = esm2plmt('H', [], [], 'Unit', 'GRAV');
            testCase.verifySize(full, [2, 6, 4]);
            testCase.verifyEqual(p, full(:, 1:3, :));
        end
        function forceNewRefreshesBothCaches(testCase)
            testCase.fullmonth(0, 1);
            first = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.fixture(0, 125);
            cached = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.verifyEqual(cached, first);
            updated = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV', 'ForceNew', true);
            testCase.verifyEqual(updated(1, 1, 3), 2e-10, 'RelTol', 2e-14);
            testCase.fixture(0, 249);
            esm2plmt('H', [], [], 'Unit', 'GRAV', 'ForceNew', true);
            testCase.verifyFalse(isfile(testCase.cachepath('midmonth')));
            updated = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.verifyEqual(updated(1, 1, 3), 3e-10, 'RelTol', 2e-14);
        end
        function saveDataFalseWritesNeitherCache(testCase)
            testCase.fullmonth(0, 1);
            esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV', 'SaveData', false);
            testCase.verifyFalse(isfile(testCase.cachepath()));
            testCase.verifyFalse(isfile(testCase.cachepath('midmonth')));
        end
        function monthlyUsesRawCache(testCase)
            testCase.fullmonth(0, 3);
            esm2plmt('H', [], [], 'Unit', 'GRAV');
            rmdir(testCase.RawFolder, 's');
            p = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.verifyEqual(p(1, 1, 3), 3e-10, 'RelTol', 2e-14);
        end
        function incompleteMonthsDoNotBecomeBiasedMeans(testCase)
            testCase.fixture(0, 100);
            testCase.verifyWarning(@() esm2plmt('H', 'midmonth', [], ...
                'Unit', 'GRAV'), 'ULMO:esm2plmt:IncompleteEpoch');
            p = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.verifyTrue(all(isnan(p(:, :, 3:4)), 'all'));
            testCase.verifyEqual(squeeze(p(1, :, 1:2)), [0 0; 1 0; 1 1; 2 0; 2 1; 2 2]);
        end
        function emptySelectionRetainsOutputShape(testCase)
            testCase.fixture(0, 1);
            range = testCase.Start + days([2, 3]);
            [p, s, d] = esm2plmt('H', [], range, 'Unit', 'GRAV', 'Loutput', 3);
            testCase.verifySize(p, [0, 10, 4]);
            testCase.verifySize(s, [0, 10, 4]);
            testCase.verifySize(d, [0, 1]);
        end
        function forceWithoutSavingPreservesCaches(testCase)
            testCase.fullmonth(0, 1);
            original = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.fixture(0, 125);
            fresh = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV', ...
                'ForceNew', true, 'SaveData', false);
            testCase.verifyEqual(fresh(1, 1, 3), 2e-10, 'RelTol', 2e-14);
            cached = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV');
            testCase.verifyEqual(cached, original);
            raw = load(testCase.cachepath());
            testCase.verifyEqual(raw.esmPlmt(1, 1, 3), 1e-10);
        end
        function incompleteCacheIsRebuilt(testCase)
            testCase.fixture(0, 2);
            dates = testCase.Start;
            save(testCase.cachepath(), 'dates');
            p = esm2plmt('H', [], [], 'Unit', 'GRAV');
            testCase.verifyEqual(p(1, 1, 3), 2e-10);
            testCase.verifyEqual(numel(who('-file', testCase.cachepath())), 3);
        end
        function missingInteriorMonthHasNaNs(testCase)
            testCase.fullmonth(0, 2);
            testCase.fullmonth(2, 7);
            testCase.verifyWarning(@() esm2plmt('H', 'midmonth', [], ...
                'Unit', 'GRAV'), 'ULMO:esm2plmt:IncompleteEpoch');
            [p, ~, d] = esm2plmt('H', 'midmonth', [], 'Unit', 'GRAV', ...
                'TimeFormat', 'datetime');
            testCase.verifySize(p, [3, 6, 4]);
            testCase.verifyTrue(all(isnan(p(2, :, 3:4)), 'all'));
            testCase.verifyEqual(p([1, 3], 1, 3), [2; 7] * 1e-10, 'RelTol', 2e-14);
            testCase.verifyEqual(d(2), datetime(2003, 2, 15));
        end
        function densityConversionPreservesLowDegrees(testCase)
            % Explicit physical factors include the CF k1 = 0.021 term.
            testCase.fixture(0, 2);
            [p, ~] = esm2plmt('H', [], ...
                [testCase.Start, testCase.Start + hours(6)]);
            degree = [0; 1; 1; 2; 2; 2];
            love = [0; 0.021; 0.021; 0; 0; 0];
            factor = 6378137 * 5517 / 3 .* (2 * degree + 1) ./ (1 + love);
            testCase.verifyEqual(reshape(p(1, :, 3), [], 1), 2e-10 * (1:6)' .* factor, ...
                'RelTol', 2e-14);
        end
        function formatsTruncationAndPadding(testCase)
            testCase.fixture(0, 1);
            range = [testCase.Start, testCase.Start + hours(6)];
            [p, s, d] = esm2plmt('H', [], range, ...
                'Unit', 'POT', 'Loutput', 1, 'OutputFormat', 'traditional', ...
                'TimeFormat', 'datetime');
            testCase.verifySize(p, [3, 4]);
            testCase.verifyEqual(p(:, 3), 6378137e-10 * (1:3)', 'RelTol', 1e-14);
            testCase.verifyEqual(s(:, 1:2), p(:, 1:2));
            testCase.verifyEqual(d, testCase.Start);
            p = esm2plmt('H', [], range, ...
                'Unit', 'GRAV', 'Loutput', 3);
            testCase.verifySize(p, [1, 10, 4]);
            testCase.verifyTrue(all(isnan(p(1, 7:10, 3:4)), 'all'));
            testCase.verifyEqual(squeeze(p(1, 7:10, 1:2)), [3, 0; 3, 1; 3, 2; 3, 3]);
        end
        function rejectDuplicateSamples(testCase)
            file = testCase.fixture(0, 1);
            mkdir(fullfile(testCase.RawFolder, 'duplicate'));
            copyfile(file, fullfile(testCase.RawFolder, 'duplicate'));
            testCase.verifyError(@() esm2plmt('H', [], [], 'Unit', 'GRAV'), ...
                'ULMO:esm2plmt:InvalidSampling');
        end
        function ignoresHeaderMetadata(testCase)
            % Header descriptions are not validated; radius still sets units.
            testCase.fixture(0, 1);
            file = testCase.fixture(1, 2);
            testCase.rewrite(file, '6378137', '6378136');
            testCase.rewrite(file, 'fully_normalized', 'unused-description');
            testCase.rewrite(file, 'max_degree 2', 'max_degree 99');
            testCase.rewrite(file, 'mtmshc_H_20030101_06', 'different-model-name');
            testCase.rewrite(file, 'earth_gravity_constant 3.986005D14', 'unused_key ignored');
            range = [testCase.Start, testCase.Start + hours(12)];
            p = esm2plmt('H', [], range, 'Unit', 'GRAV');
            testCase.verifySize(p, [2, 6, 4]);
            testCase.verifyEqual(p(:, 1, 3), [1; 2] * 1e-10, 'RelTol', 1e-14);
            p = esm2plmt('H', [], range, 'Unit', 'POT');
            testCase.verifyEqual(p(:, 1, 3), [6378137; 2 * 6378136] * 1e-10, ...
                'RelTol', 1e-14);
        end
        function rejectBadHarmonicsAndInputs(testCase)
            file = testCase.fixture(0, 1);
            range = [testCase.Start, testCase.Start + hours(6)];
            testCase.rewrite(file, 'gfc 2 2', 'gfc 2 1');
            testCase.verifyError(@() esm2plmt('H', [], range, ...
                'Unit', 'GRAV'), 'ULMO:esm2plmt:InvalidCoefficients');
            testCase.verifyError(@() esm2plmt('H', hours(6)), 'ULMO:esm2plmt:InvalidTimeStep');
            testCase.verifyError(@() esm2plmt('H', [], [range; range]), 'ULMO:notValidTimeRange');
        end
    end
    methods (Access = private)
        function path = cachepath(testCase, timestep)
            if nargin < 2 || isempty(timestep)
                name = 'ESAESM2_H_GRAV.mat';
            else
                name = 'ESAESM2_H_GRAV-Tmidmonth.mat';
            end
            path = fullfile(testCase.Scratch, 'ESM', 'ESAESM2', name);
        end
        function fullmonth(testCase, offset, amplitude)
            first = testCase.Start + calmonths(offset);
            last = first + calmonths(1);
            steps = hours(first - testCase.Start) / 6:hours(last - testCase.Start) / 6 - 1;
            for step = steps
                testCase.fixture(step, amplitude);
            end
        end
        function file = fixture(testCase, step, amplitude)
            timestamp = testCase.Start + hours(6 * step);
            timestamp.Format = 'yyyyMMdd_HH';
            name = ['mtmshc_H_', char(timestamp)];
            file = fullfile(testCase.RawFolder, [name, '.gfc']);
            fid = fopen(file, 'w');
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, ['product_type gravity_field\nmodelname %s\n', ...
                'earth_gravity_constant 3.986005D14\nradius 6378137\n', ...
                'max_degree 2\nerror no\nnorm fully_normalized\nend_of_head\n'], name);
            % Deliberately reverse the coefficient rows to check harmonic sorting.
            labels = [0, 0; 1, 0; 1, 1; 2, 0; 2, 1; 2, 2];
            for k = 6:-1:1
                fprintf(fid, 'gfc %d %d %.15e %.15e\n', labels(k, :), ...
                    amplitude * k * 1e-10, amplitude * labels(k, 2) * 1e-10);
            end
        end
        function rewrite(~, file, old, new)
            content = strrep(fileread(file), old, new);
            fid = fopen(file, 'w');
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s', content);
        end
    end
end
