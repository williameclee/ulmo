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
            testCase.RawFolder = fullfile(testCase.Scratch, 'external-source', 'H');
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
        function importerUsesExternalFolderAndNativeMat(testCase)
            testCase.fixture(0, 1);
            testCase.fixture(1, 2);
            testCase.prepare();
            data = load(testCase.inputpath());
            testCase.verifyEqual(sort(string(fieldnames(data))), sort(["esmPlmt"; "dates"; "equatorRadius"]));
            testCase.verifySize(data.esmPlmt, [2, 6, 4]);
            testCase.verifyEqual(data.esmPlmt(:, 1, 3), [1; 2] * 1e-10, 'RelTol', 1e-14);
            testCase.verifyEqual(data.equatorRadius, [6378137; 6378137]);
            rmdir(testCase.RawFolder, 's');
            [p, s, d] = esm2plmt('H', 'Unit', 'GRAV');
            testCase.verifyEqual(p, data.esmPlmt);
            testCase.verifyTrue(all(isnan(s(:, :, 3:4)), 'all'));
            testCase.verifyEqual(d, datenum(data.dates));
        end
        function readerNeverImportsGfcFiles(testCase)
            testCase.fixture(0, 2);
            testCase.verifyError(@() esm2plmt('H', 'Unit', 'GRAV'), 'ULMO:esm2plmt:MissingPreparedData');
            testCase.verifyFalse(isfile(testCase.inputpath()));
        end
        function requiresConfiguredRoot(testCase)
            setenv('IFILES', '');
            testCase.verifyError(@() esm2plmt('O'), 'ULMO:esm2plmt:MissingDataRoot');
            testCase.verifyError(@() prepareesm(testCase.RawFolder, 'H'), 'ULMO:prepareesm:MissingDataRoot');
        end
        function defaultOceanAndCharStringOptions(testCase)
            testCase.fixture(0, 2);
            testCase.prepare();
            copyfile(testCase.inputpath(), fullfile(fileparts(testCase.inputpath()), 'ESAESM2_O.mat'));
            default = esm2plmt();
            explicit = esm2plmt("O", [], "", "Unit", "SD", "LoveNumSource", "ISSM");
            testCase.verifyEqual(default, explicit);
            [p, ~, d] = esm2plmt('O', [], '', 'Unit', 'GRAV', ...
                'OutputFormat', 'traditional', 'TimeFormat', 'datetime');
            testCase.verifyEqual(p(1, 3), 2e-10);
            testCase.verifyEqual(d, testCase.Start);
        end
        function inclusiveRawSelectionAndPadding(testCase)
            for k = 0:4
                testCase.fixture(k, k + 1);
            end
            testCase.prepare();
            range = testCase.Start + hours([6, 18]);
            [p, s, dates] = esm2plmt('H', datenum(range), '', 'Unit', 'GRAV', 'Loutput', 3);
            testCase.verifySize(p, [3, 10, 4]);
            testCase.verifyEqual(p(:, 1, 3), [2; 3; 4] * 1e-10, 'RelTol', 1e-14);
            testCase.verifyEqual(s(:, :, 1:2), p(:, :, 1:2));
            testCase.verifyTrue(all(isnan(p(:, 7:10, 3:4)), 'all'));
            testCase.verifyEqual(dates, datenum(testCase.Start + hours([6; 12; 18])));
            data = load(testCase.inputpath());
            testCase.verifySize(data.esmPlmt, [5, 6, 4]);
        end
        function monthlyCachesAllDatesAndNativeDegreesBeforeSelection(testCase)
            testCase.Start = datetime(2004, 1, 1);
            testCase.fullmonth(0, 2);
            testCase.fullmonth(1, 7);
            testCase.prepare();
            range = [datetime(2004, 2, 15), datetime(2004, 2, 16)];
            [p, s, dates] = esm2plmt('H', range, 'midmonth', 'Unit', 'GRAV', ...
                'Loutput', 1, 'TimeFormat', 'datetime', 'OutputFormat', 'traditional');
            testCase.verifySize(p, [3, 4]);
            testCase.verifyEqual(p(:, 3), 7e-10 * (1:3)', 'RelTol', 2e-14);
            testCase.verifyEqual(dates, datetime(2004, 2, 15, 12, 0, 0));
            testCase.verifyEqual(s(:, 3:4), zeros(3, 2), 'AbsTol', 1e-23);
            cached = load(testCase.cachepath());
            testCase.verifySize(cached.esmPlmt, [2, 6, 4]);
            testCase.verifyEqual(cached.esmPlmt(:, 1, 3), [2; 7] * 1e-10, 'RelTol', 2e-14);
            testCase.verifyEqual(sort(string(fieldnames(cached))), sort(["dates"; "esmPlmt"; "esmStdPlmt"]));
            delete(testCase.inputpath());
            rmdir(testCase.RawFolder, 's');
            [p, s, d] = esm2plmt('H', [], "midmonth", 'Unit', 'GRAV', 'TimeFormat', 'datetime');
            testCase.verifyEqual(p, cached.esmPlmt);
            testCase.verifyEqual(s, cached.esmStdPlmt);
            testCase.verifyEqual(d, cached.dates);
        end
        function monthlyEstimatesSpreadAndCachesIt(testCase)
            testCase.fullmonth(0, 1);
            testCase.fixture(0, 125);
            testCase.prepare();
            rmdir(testCase.RawFolder, 's');
            [p, s] = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV');
            testCase.verifyEqual(p(1, 1, 3), 2e-10, 'RelTol', 2e-14);
            testCase.verifyEqual(s(1, 1, 3), 1e-10, 'RelTol', 2e-14);
            data = load(testCase.cachepath());
            testCase.verifyEqual(s, data.esmStdPlmt);
        end
        function monthlyPropagatesSuppliedFormalErrors(testCase)
            testCase.fullmonth(0, 2);
            testCase.prepare();
            source = matfile(testCase.inputpath(), 'Writable', true);
            sigma = source.esmPlmt;
            sigma(:, :, 3:4) = 3;
            source.esmStdPlmt = sigma;
            [~, s] = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV');
            testCase.verifyEqual(s(:, :, 3:4), repmat(3 / sqrt(124), 1, 6, 2), 'RelTol', 1e-14);
        end
        function forceNewReadsPreparedDataOnly(testCase)
            testCase.fullmonth(0, 1);
            testCase.prepare();
            first = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV');
            source = matfile(testCase.inputpath(), 'Writable', true);
            source.esmPlmt(1, 1, 3) = 125e-10;
            rmdir(testCase.RawFolder, 's');
            cached = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV');
            testCase.verifyEqual(cached, first);
            fresh = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV', 'ForceNew', true);
            testCase.verifyEqual(fresh(1, 1, 3), 2e-10, 'RelTol', 2e-14);
        end
        function saveDataFalsePreservesInputAndCache(testCase)
            testCase.fullmonth(0, 1);
            testCase.prepare();
            original = load(testCase.inputpath());
            esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV', 'SaveData', false);
            testCase.verifyFalse(isfile(testCase.cachepath()));
            esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV');
            cached = load(testCase.cachepath());
            source = matfile(testCase.inputpath(), 'Writable', true);
            source.esmPlmt(1, 1, 3) = 125e-10;
            fresh = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV', 'ForceNew', true, 'SaveData', false);
            testCase.verifyEqual(fresh(1, 1, 3), 2e-10, 'RelTol', 2e-14);
            testCase.verifyEqual(load(testCase.cachepath()), cached);
            % Raw force/read options never rewrite the prepared input.
            before = load(testCase.inputpath());
            esm2plmt('H', [], '', 'Unit', 'SD', 'ForceNew', true);
            testCase.verifyEqual(load(testCase.inputpath()), before);
            testCase.verifyEqual(before.dates, original.dates);
            testCase.verifyFalse(isfile(fullfile(fileparts(testCase.inputpath()), 'ESAESM2_H_SD_ISSM.mat')));
        end
        function reimportInvalidatesEveryUnitAndEpochVariant(testCase)
            testCase.fixture(0, 1);
            testCase.prepare();
            folder = fileparts(testCase.inputpath());
            variants = ["ESAESM2_H_GRAV-Tmidmonth.mat", "ESAESM2_H_SD_ISSM-Tfuture.mat", ...
                "ESAESM2_H_POT-Tmidmonth-mean.mat"];
            sentinel = 1;
            for name = variants
                save(fullfile(folder, name), 'sentinel');
            end
            other = fullfile(folder, 'ESAESM2_O_GRAV-Tmidmonth.mat');
            save(other, 'sentinel');
            testCase.fixture(0, 2);
            testCase.prepare('ForceNew', true);
            for name = variants
                testCase.verifyFalse(isfile(fullfile(folder, name)));
            end
            testCase.verifyTrue(isfile(other));
            data = load(testCase.inputpath());
            testCase.verifyEqual(data.esmPlmt(1, 1, 3), 2e-10);
        end
        function failedImportPreservesPreparedDataAndCaches(testCase)
            testCase.fixture(0, 1);
            testCase.prepare();
            original = load(testCase.inputpath());
            sentinel = 1;
            save(testCase.cachepath(), 'sentinel');
            testCase.verifyError(@() testCase.prepare(), 'ULMO:prepareesm:FileExists');
            file = testCase.fixture(1, 2);
            testCase.rewrite(file, 'gfc 2 2', 'gfc 2 1');
            testCase.verifyError(@() testCase.prepare('ForceNew', true), 'ULMO:prepareesm:InvalidCoefficients');
            testCase.verifyEqual(load(testCase.inputpath()), original);
            testCase.verifyEqual(load(testCase.cachepath()), struct('sentinel', 1));
            testCase.verifyEqual(numel(dir(fullfile(fileparts(testCase.inputpath()), '*.mat'))), 2);
        end
        function incompleteMonthAndInteriorGapHaveNaNs(testCase)
            testCase.fullmonth(0, 2);
            testCase.fullmonth(2, 7);
            testCase.prepare();
            testCase.verifyWarning(@() esm2plmt('H', [], 'midmonth', ...
                'Unit', 'GRAV'), 'ULMO:averageepochs:IncompleteEpoch');
            [p, s, d] = esm2plmt('H', [], 'midmonth', 'Unit', 'GRAV', 'TimeFormat', 'datetime');
            testCase.verifySize(p, [3, 6, 4]);
            testCase.verifyTrue(all(isnan(p(2, :, 3:4)), 'all'));
            testCase.verifyTrue(all(isnan(s(2, :, 3:4)), 'all'));
            testCase.verifyEqual(d(2), datetime(2003, 2, 15));
            source = matfile(testCase.inputpath(), 'Writable', true);
            dates = source.dates;
            dates(1) = dates(1) - hours(6);
            source.dates = dates;
            testCase.verifyWarning(@() esm2plmt('H', [], 'midmonth', ...
                'Unit', 'GRAV', 'ForceNew', true), 'ULMO:averageepochs:IncompleteEpoch');
            data = load(testCase.cachepath());
            testCase.verifyTrue(all(isnan(data.esmPlmt(2, :, 3:4)), 'all'));
        end
        function emptySelectionRetainsShape(testCase)
            testCase.fixture(0, 1);
            testCase.prepare();
            [p, s, d] = esm2plmt('H', testCase.Start + days([2, 3]), '', 'Unit', 'GRAV', 'Loutput', 3);
            testCase.verifySize(p, [0, 10, 4]);
            testCase.verifySize(s, [0, 10, 4]);
            testCase.verifySize(d, [0, 1]);
        end
        function densityAndPotentialConversionPreserveLowDegrees(testCase)
            testCase.fixture(0, 2);
            file = testCase.fixture(1, 3);
            testCase.rewrite(file, '6378137', '6378136');
            testCase.prepare();
            p = esm2plmt('H');
            degree = [0; 1; 1; 2; 2; 2];
            love = [0; 0.021; 0.021; 0; 0; 0];
            factor = 6378137 * 5517 / 3 .* (2 * degree + 1) ./ (1 + love);
            testCase.verifyEqual(reshape(p(1, :, 3), [], 1), 2e-10 * (1:6)' .* factor, 'RelTol', 2e-14);
            [p, ~, d] = esm2plmt('H', [], '', 'Unit', 'POT', 'Loutput', 1, ...
                'OutputFormat', 'traditional', 'TimeFormat', 'datetime');
            testCase.verifySize(p, [3, 4, 2]);
            testCase.verifyEqual(squeeze(p(1, 3, :)), [2 * 6378137; 3 * 6378136] * 1e-10, 'RelTol', 1e-14);
            testCase.verifyEqual(d, testCase.Start + hours([0; 6]));
        end
        function invalidPreparedInputNeverFallsBackToGfc(testCase)
            testCase.fixture(0, 1);
            testCase.prepare();
            dates = testCase.Start;
            save(testCase.inputpath(), 'dates', '-v7.3');
            testCase.verifyError(@() esm2plmt('H'), 'ULMO:esm2plmt:InvalidPreparedData');
        end
        function rejectDuplicateAndOffGridSources(testCase)
            file = testCase.fixture(0, 1);
            mkdir(fullfile(testCase.RawFolder, 'duplicate'));
            copyfile(file, fullfile(testCase.RawFolder, 'duplicate'));
            testCase.verifyError(@() testCase.prepare(), 'ULMO:prepareesm:InvalidSampling');
            rmdir(fullfile(testCase.RawFolder, 'duplicate'), 's');
            movefile(file, strrep(file, '_00.gfc', '_01.gfc'));
            testCase.verifyError(@() testCase.prepare(), 'ULMO:prepareesm:InvalidSampling');
        end
        function ignoresDescriptiveMetadataButRequiresRadius(testCase)
            file = testCase.fixture(0, 1);
            testCase.rewrite(file, 'fully_normalized', 'unused-description');
            testCase.rewrite(file, 'max_degree 2', 'max_degree 99');
            testCase.prepare();
            testCase.rewrite(file, 'radius 6378137', 'unused-radius 6378137');
            testCase.verifyError(@() testCase.prepare('ForceNew', true), 'ULMO:prepareesm:InvalidRadius');
        end
        function rejectsInvalidTimeArguments(testCase)
            testCase.verifyError(@() esm2plmt('H', [], 'GRACE'), 'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() esm2plmt('H', [], []), 'MATLAB:validators:mustBeTextScalar');
            testCase.verifyError(@() esm2plmt('H', [testCase.Start, NaT]), 'ULMO:notValidTimeRange');
            testCase.verifyError(@() esm2plmt('H', repmat(testCase.Start, 2)), 'ULMO:notValidTimeRange');
        end
        function demoPropagatesFormalErrorsFromPreparedData(testCase)
            for k = 0:3
                testCase.fixture(k, k + 1);
            end
            testCase.prepare();
            source = matfile(testCase.inputpath(), 'Writable', true);
            sigma = source.esmPlmt;
            sigma(:, :, 3) = 3;
            sigma(:, :, 4) = 4;
            source.esmStdPlmt = sigma;
            folder = fileparts(testCase.inputpath());
            for component = ["I", "O"]
                copyfile(testCase.inputpath(), fullfile(folder, "ESAESM2_" + component + ".mat"));
            end
            graceFile = fullfile(testCase.Scratch, 'GSM-fixture');
            fid = fopen(graceFile, 'w');
            for row = 1:6
                fprintf(fid, 'GRCOF2 %d %d 1 0 0 0 20030101.0100 20030101.1300 nnnn\n', ...
                    sigma(1, row, 1:2));
            end
            fclose(fid);
            products = esm2plmt_demo(graceFile, fullfile(testCase.Scratch, 'products'), ...
                'Loutput', 1, 'MakeFigure', false);
            h = load(products.H.nativeGRAV);
            hi = load(products.HI.nativeSD);
            hio = load(products.HIO.nativeSD);
            testCase.verifyEqual(h.esmStdPlmt(1, 1, 3), 3 / sqrt(2), 'RelTol', 1e-14);
            factor = 6378137 * 5517 / 3;
            testCase.verifyEqual(hi.esmStdPlmt(1, 1, 3), 3 * factor, 'RelTol', 1e-14);
            testCase.verifyEqual(hio.esmStdPlmt(1, 1, 3), 3 * sqrt(3 / 2) * factor, 'RelTol', 1e-14);
            testCase.verifyEqual(h.dates, datenum(datetime(2003, 1, 1, 7, 0, 0)));
            % Perfectly correlated component variability must add before the
            % empirical standard error, rather than combine in quadrature.
            for component = ["H", "I", "O"]
                path = fullfile(folder, "ESAESM2_" + component + ".mat");
                data = rmfield(load(path), 'esmStdPlmt');
                save(path, '-struct', 'data', '-v7.3');
            end
            products = esm2plmt_demo(graceFile, fullfile(testCase.Scratch, 'empirical'), ...
                'Loutput', 1, 'MakeFigure', false);
            hio = load(products.HIO.nativeSD);
            testCase.verifyEqual(hio.esmStdPlmt(1, 1, 3), 1.5e-10 * factor, 'RelTol', 2e-14);
        end
    end
    methods (Access = private)
        function path = inputpath(testCase)
            path = fullfile(testCase.Scratch, 'ESM', 'ESAESM2', 'ESAESM2_H.mat');
        end
        function path = cachepath(testCase)
            path = fullfile(testCase.Scratch, 'ESM', 'ESAESM2', 'ESAESM2_H_GRAV-Tmidmonth.mat');
        end
        function prepare(testCase, varargin)
            prepareesm(testCase.RawFolder, 'H', 'BeQuiet', true, varargin{:});
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
