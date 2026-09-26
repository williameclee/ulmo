%% graceInputFormatsTest - Synthetic input/correction fixtures.
% Created by
%   2026/09/26, En-Chi Lee (williameclee@azirona.edu)

classdef graceInputFormatsTest < matlab.unittest.TestCase

    properties
        Scratch
    end

    methods (TestMethodSetup)

        function setup(testCase)
            testCase.Scratch = tempname;
            mkdir(testCase.Scratch);
            testCase.addTeardown(@() rmdir(testCase.Scratch, 's'));

            for name = {'GRACEDATA', 'ORIGINALGRACEDATA', 'IFILES'}
                original = getenv(name{1});
                testCase.addTeardown(@() setenv(name{1}, original));
                setenv(name{1}, testCase.Scratch);
            end

            visible = get(0, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(0, 'DefaultFigureVisible', visible));
            set(0, 'DefaultFigureVisible', 'off');
            constants = fullfile(testCase.Scratch, 'EARTHMODELS', 'CONSTANTS');
            mkdir(constants);
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(constants, 'Earth.mat'), 'Earth');
            mkdir(fullfile(testCase.Scratch, 'LOVENUMS'));
            loadinggravitationalpotential = zeros(721, 1);
            save(fullfile(testCase.Scratch, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
            'loadinggravitationalpotential');
        end

    end

    methods (Test)

        function legacyEarthParameters(testCase)
            % Verify parsing and processing of legacy GRACE source file formats (e.g., GFZ RL05),
            % ensuring Earth parameters (GM, radius), dates, and POT unit scaling are correctly applied.
            f = testCase.fixture('GFZ', 'RL05', 90, 1, 'legacy', 6378136.46);
            [g, s, d, dr, gm, radius] = parsegracesourcefile(f);
            testCase.verifySize(g, [4186, 4]);
            testCase.verifySize(s, size(g));
            testCase.verifyEqual(gm, 3.986004415e14);
            testCase.verifyEqual(radius, 6378136.46);
            testCase.verifyEqual(d, mean(dr));
            [p, ~] = testCase.process('GFZ', 'RL05', 90, 'POT');
            testCase.verifyEqual(p(1, 5, 3), g(5, 3) * radius, 'RelTol', 1e-12);
        end

        function yamlMetadataAndValidation(testCase)
            % Verify YAML header parsing for RL06 metadata and confirm that invalid/corrupted
            % header keys raise the expected parsing error.
            f = testCase.fixture('CSR', 'RL06', 60, 1, 'yaml', 6378136.3);
            [~, ~, ~, ~, gm, radius] = parsegracesourcefile(f);
            testCase.verifyEqual(gm, 3.986004415e14);
            testCase.verifyEqual(radius, 6378136.3);
            text = fileread(f);
            % Do not accidentally use the next attribute's value as GM.
            text = strrep(text, 'value: 3.986004415e14', 'missing: 3.986004415e14');
            fid = fopen(f, 'w'); fprintf(fid, '%s', text); fclose(fid);
            testCase.verifyError(@() sixOutputs(f), 'ULMO:parsegracesourcefile:InvalidHeader');
        end

        function mixedJplDegrees(testCase)
            % Verify ingestion of monthly GRACE solutions with differing maximum degrees (e.g., JPL RL05
            % mixed 60/90), ensuring proper NaN padding in standard format and truncation in traditional format.
            testCase.fixture('JPL', 'RL05', 60, 1, 'legacy', 6378136.3);
            testCase.fixture('JPL', 'RL05', 90, 2, 'legacy', 6378136.3);
            [p, s] = testCase.process('JPL', 'RL05', 90, 'SD');
            testCase.verifySize(p, [2, 4186, 4]);
            testCase.verifyEqual(p(:, :, 1:2), s(:, :, 1:2));
            testCase.verifyTrue(all(isfinite(p(:, :, 1:2)), 'all'));
            testCase.verifyTrue(all(isnan(p(1, 1892:end, 3:4)), 'all'));
            testCase.verifyTrue(all(isnan(s(1, 1892:end, 3:4)), 'all'));
            testCase.verifyTrue(all(isfinite(p(2, :, :)), 'all'));
            [p, s] = testCase.process('JPL', 'RL05', 90, 'SD', 'Loutput', 60, ...
                'OutputFormat', 'traditional');
            testCase.verifySize(p, [1891, 4, 2]);
            testCase.verifyTrue(all(isfinite(p), 'all'));
            testCase.verifyTrue(all(isfinite(s), 'all'));
        end

        function csrBandwidthAndOutputPadding(testCase)
            % Verify bandwidth limits and output zero/NaN padding for CSR RL05, checking that exceeding
            % input data degree issues a warning and specifying Loutput pads higher degrees appropriately.
            testCase.fixture('CSR', 'RL05', 60, 1, 'legacy', 6378136.3);
            testCase.verifyWarning(@() testCase.process('CSR', 'RL05', 96, 'SD'), ...
            'ULMO:grace2plmt:InvalidLdata');
            [p, s] = testCase.process('CSR', 'RL05', 60, 'SD', 'Loutput', 96);
            testCase.verifySize(p, [1, 4753, 4]);
            testCase.verifySize(s, size(p));
            testCase.verifyTrue(all(isfinite(p(:, :, 1:2)), 'all'));
            testCase.verifyEqual(p(:, :, 1:2), s(:, :, 1:2));
            testCase.verifyTrue(all(isnan(p(:, 1892:end, 3:4)), 'all'));
            testCase.verifyTrue(all(isnan(s(:, 1892:end, 3:4)), 'all'));
        end

        function rl05CorrectionContract(testCase)
            % Verify the correction contract for RL05 data, ensuring default/disabled corrections behave
            % identically and requesting unsupported corrections (Degree 1, C20, C30) raises an error.
            testCase.fixture('CSR', 'RL05', 60, 1, 'legacy', 6378136.3);
            expected = testCase.process('CSR', 'RL05', 60, 'GRAV');
            actual = testCase.process('CSR', 'RL05', 60, 'GRAV', ...
                'Deg1Correction', false, 'C20Correction', false, 'C30Correction', false);
            testCase.verifyEqual(actual, expected);
            actual = testCase.process('CSR', 'RL05', 60, 'GRAV', ...
                'Deg1Correction', [], 'C20Correction', [], 'C30Correction', []);
            testCase.verifyEqual(actual, expected);

            for name = {'Deg1Correction', 'C20Correction', 'C30Correction'}
                testCase.verifyError(@() testCase.process('CSR', 'RL05', 60, 'GRAV', name{1}, true), ...
                'ULMO:grace2plmt:UnsupportedCorrection');
            end

        end

        function masconReferenceAndCorrections(testCase)
            % Verify processing of CSR RL06 mascon (GSU) solutions with and without external Degree 1
            % (TN-13) and C20/C30 (TN-14) corrections applied to spherical harmonic coefficients and errors.
            f = testCase.fixture('CSR mascon', 'RL06', 720, 1, 'yaml', 6378136.3);
            raw = parsegracesourcefile(f);
            [p, s] = testCase.process('CSR mascon', 'RL06', 720, 'GRAV', ...
                'Deg1Correction', false, 'C20Correction', false, 'C30Correction', false);
            raw(1, 3) = 0;
            testCase.verifyEqual(squeeze(p), raw);
            testCase.corrections();
            [p, correctedStd] = testCase.process('CSR mascon', 'RL06', 720, 'GRAV');
            testCase.verifyEqual(p(1, 4, 3), 2e-10, 'AbsTol', 1e-19);
            testCase.verifyEqual(p(1, 7, 3), 3e-10, 'AbsTol', 1e-19);
            testCase.verifyEqual(p(1, 11, 3), raw(11, 3));
            testCase.verifyEqual(squeeze(p(1, 2:3, 3:4)), [1e-10 0; 2e-10 3e-10]);
            testCase.verifyEqual(correctedStd(1, 4, 3), 4e-12, 'AbsTol', 1e-25);
            testCase.verifyEqual(correctedStd(1, 7, 3), 5e-12, 'AbsTol', 1e-25);
            testCase.verifyEqual(correctedStd(1, 11, :), s(1, 11, :));
        end

        function rl06GsmKeepsFullFieldConvention(testCase)
            % Verify that RL06 GSM spherical harmonic fields maintain the full-field gravity convention
            % when applying SLR C20/C30 replacements by including the background static field component.
            testCase.fixture('CSR', 'RL06', 60, 1, 'yaml', 6378136.3);
            testCase.corrections();
            p = testCase.process('CSR', 'RL06', 60, 'GRAV');
            testCase.verifyEqual(p(1, 4, 3), ...
                -4.841694573200e-4 +2e-10 +0.108262982131e-2 / sqrt(5), 'AbsTol', 1e-19);
            testCase.verifyEqual(p(1, 7, 3), 9.571647583412e-7 +3e-10);
        end

        function bypassOldCache(testCase)
            % Verify that legacy cache files are bypassed, new cache files (*_inputV2.mat) are created,
            % and data can be loaded directly from the new cache in the absence of source files.
            testCase.fixture('CSR', 'RL05', 60, 1, 'legacy', 6378136.3);
            gracePlmt = -999; graceStdPlmt = -999; dates = NaT; %#ok<NASGU>
            oldCache = fullfile(testCase.Scratch, 'CSR_RL05_alldata_nDeg1_nC20_nC30_60_GRAV.mat');
            save(oldCache, 'gracePlmt', 'graceStdPlmt', 'dates');
            p = grace2plmt_new('CSR', 'RL05', 60, 'GRAV', 'ForceNew', false, ...
                'SaveData', true, 'BeQuiet', true);
            testCase.verifySize(p, [1, 1891, 4]);
            cache = dir(fullfile(testCase.Scratch, '*_inputV2.mat'));
            testCase.verifyNumElements(cache, 1);
            % Prove that the new cache can be read without the source files.
            rmdir(fullfile(testCase.Scratch, 'RL05'), 's');
            loaded = grace2plmt_new('CSR', 'RL05', 60, 'GRAV', ...
                'ForceNew', false, 'BeQuiet', true);
            testCase.verifyEqual(loaded, p);
        end

    end

    methods

        function varargout = process(~, center, release, L, unit, varargin)
            [varargout{1:nargout}] = grace2plmt_new(center, release, L, unit, ...
                'ForceNew', true, 'SaveData', false, 'BeQuiet', true, varargin{:});
        end

        function f = fixture(testCase, center, release, L, month, format, radius)
            folder = fullfile(testCase.Scratch, release, center);
            if ~isfolder(folder), mkdir(folder); end

            if strcmp(release, 'RL05')

                switch center
                    case 'CSR', suffix = 'UTCSR_0060_0005.txt';
                    case 'GFZ', suffix = 'EIGEN_G---_005a';
                    case 'JPL', suffix = 'JPLEM_0000_0005.txt';
                end

                prefix = 'GSM';
            elseif strcmp(center, 'CSR mascon')
                prefix = 'GSU'; suffix = 'UTCSR_B---_0600.txt';
            else
                prefix = 'GSM'; suffix = 'UTCSR_BA01_0600.txt';
            end

            f = fullfile(folder, sprintf('%s-2_2020%03d-2020%03d_TEST_%s', ...
                prefix, month * 30, month * 30 + 29, suffix));
            fid = fopen(f, 'w'); cleanup = onCleanup(@() fclose(fid));

            if strcmp(format, 'legacy')
                fprintf(fid, 'FIRST synthetic GRACE fixture\nEARTH 3.986004415D+14 %.10E\nSHM %d %d\n', radius, L, L);
            else
                fprintf(fid, ['header:\n  non-standard_attributes:\n' ...
                                  '    earth_gravity_param:\n      value: 3.986004415e14\n' ...
                              '    mean_equator_radius:\n      value: %.10E\n# End of YAML header\n'], radius);
            end

            [order, degree] = addmon(L);
            first = 1;
            if strcmp(center, 'JPL'), first = 4; end
            startDate = datetime(2020, month, 1, 'Format', 'yyyyMMdd');
            endDate = startDate + calmonths(1);
            records = [degree(first:end), order(first:end), ...
                           1e-10 * (degree(first:end) + 1), 1e-11 * order(first:end), ...
                           repmat([1e-12 2e-12], numel(order) - first + 1, 1)]';
            fmt = sprintf('GRCOF2 %%d %%d %%.14e %%.14e %%.14e %%.14e %s.0000 %s.0000 nnnn\n', ...
                string(startDate), string(endDate));
            fprintf(fid, fmt, records);
        end

        function corrections(testCase)
            mkdir(fullfile(testCase.Scratch, 'Degree1'));
            mkdir(fullfile(testCase.Scratch, 'Degree2'));
            fid = fopen(fullfile(testCase.Scratch, 'Degree1', 'TN-13_GEOC_CSR_RL06.txt'), 'w');
            fprintf(fid, 'end of header ===============================================================================\n');

            for month = 1:4
                first = string(datetime(2020, month, 1, 'Format', 'yyyyMMdd'));
                last = string(datetime(2020, month + 1, 1, 'Format', 'yyyyMMdd'));
                fprintf(fid, 'GRCOF2 1 0 1e-10 0 1e-12 0 %s.0000 %s.0000 0\n', first, last);
                fprintf(fid, 'GRCOF2 1 1 2e-10 3e-10 2e-12 3e-12 %s.0000 %s.0000 0\n', first, last);
            end

            fclose(fid);
            fid = fopen(fullfile(testCase.Scratch, 'Degree2', 'TN-14_C30_C20_GSFC_SLR.txt'), 'w');
            fprintf(fid, 'Product:\n');

            for month = 1:4
                first = convertTo(datetime(2020, month, 1), 'modifiedjuliandate');
                last = convertTo(datetime(2020, month + 1, 1), 'modifiedjuliandate');
                fprintf(fid, '%.0f 0 %.16e 0 0.04 %.16e 0 0.05 %.0f 0\n', ...
                    first, -4.841694573200e-4 +2e-10, 9.571647583412e-7 +3e-10, last);
            end

            fclose(fid);
        end

    end

end

function sixOutputs(f)
    [~, ~, ~, ~, ~, ~] = parsegracesourcefile(f);
end
