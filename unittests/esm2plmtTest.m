%% ESM2PLMTTEST - Tests interval averaging and ESM coefficient conventions.
% Synthetic fields exercise boundary selection, missing-data rejection and
% physical unit conversion independently of the downloaded ESM samples.
%
% Created by
%   2026/10/29

classdef esm2plmtTest < matlab.unittest.TestCase
    properties
        Scratch
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
        function averageAndExclusiveBoundary(testCase)
            % A large sample at the end belongs only to the next interval.
            for k = 0:4
                testCase.fixture(k, k + 1);
            end
            ranges = [testCase.Start, testCase.Start + hours(12); ...
                testCase.Start + hours(12), testCase.Start + hours(24)];
            [p, s, dates, meta] = esm2plmt(testCase.Scratch, ranges, ...
                'Component', 'H', 'Unit', 'GRAV');
            testCase.verifySize(p, [2, 6, 4]);
            testCase.verifyEqual(p(:, 1, 3), [1.5; 3.5] * 1e-10, 'RelTol', 1e-14);
            testCase.verifyEqual(s(:, :, 1:2), p(:, :, 1:2));
            testCase.verifyTrue(all(isnan(s(:, :, 3:4)), 'all'));
            testCase.verifyEqual(dates, datenum(testCase.Start + hours([6; 18])));
            testCase.verifyEqual(meta.sampleCounts, [2; 2]);
            testCase.verifyEqual(meta.gravityParam, 3.986005e14);
        end
        function densityConversionPreservesLowDegrees(testCase)
            % Explicit physical factors include the CF k1 = 0.021 term.
            testCase.fixture(0, 2);
            [p, ~] = esm2plmt(testCase.Scratch, ...
                [testCase.Start, testCase.Start + hours(6)], 'Component', 'H');
            degree = [0; 1; 1; 2; 2; 2];
            love = [0; 0.021; 0.021; 0; 0; 0];
            factor = 6378137 * 5517 / 3 .* (2 * degree + 1) ./ (1 + love);
            testCase.verifyEqual(reshape(p(1, :, 3), [], 1), 2e-10 * (1:6)' .* factor, ...
                'RelTol', 2e-14);
        end
        function formatsTruncationAndPadding(testCase)
            testCase.fixture(0, 1);
            range = [testCase.Start, testCase.Start + hours(6)];
            [p, s, d] = esm2plmt(testCase.Scratch, range, 'Component', 'H', ...
                'Unit', 'POT', 'Loutput', 1, 'OutputFormat', 'traditional', ...
                'TimeFormat', 'datetime');
            testCase.verifySize(p, [3, 4]);
            testCase.verifyEqual(p(:, 3), 6378137e-10 * (1:3)', 'RelTol', 1e-14);
            testCase.verifyEqual(s(:, 1:2), p(:, 1:2));
            testCase.verifyEqual(d, testCase.Start + hours(3));
            p = esm2plmt(testCase.Scratch, range, 'Component', 'H', ...
                'Unit', 'GRAV', 'Loutput', 3);
            testCase.verifySize(p, [1, 10, 4]);
            testCase.verifyTrue(all(isnan(p(1, 7:10, 3:4)), 'all'));
            testCase.verifyEqual(squeeze(p(1, 7:10, 1:2)), [3, 0; 3, 1; 3, 2; 3, 3]);
        end
        function rejectMissingAndDuplicateSamples(testCase)
            file = testCase.fixture(0, 1);
            range = [testCase.Start, testCase.Start + hours(12)];
            call = @() esm2plmt(testCase.Scratch, range, 'Component', 'H', 'Unit', 'GRAV');
            testCase.verifyError(call, 'ULMO:esm2plmt:IncompleteSampling');
            testCase.fixture(1, 2);
            mkdir(fullfile(testCase.Scratch, 'duplicate'));
            copyfile(file, fullfile(testCase.Scratch, 'duplicate'));
            testCase.verifyError(call, 'ULMO:esm2plmt:IncompleteSampling');
        end
        function rejectInconsistentHeaders(testCase)
            testCase.fixture(0, 1);
            file = testCase.fixture(1, 2);
            testCase.rewrite(file, '6378137', '6378136');
            range = [testCase.Start, testCase.Start + hours(12)];
            testCase.verifyError(@() esm2plmt(testCase.Scratch, range, ...
                'Component', 'H', 'Unit', 'GRAV'), 'ULMO:esm2plmt:MixedMetadata');
            testCase.rewrite(file, 'fully_normalized', 'unnormalized');
            testCase.verifyError(@() esm2plmt(testCase.Scratch, range, ...
                'Component', 'H', 'Unit', 'GRAV'), 'ULMO:esm2plmt:InvalidHeader');
        end
        function rejectBadHarmonicsAndIntervals(testCase)
            file = testCase.fixture(0, 1);
            range = [testCase.Start, testCase.Start + hours(6)];
            testCase.rewrite(file, 'gfc 2 2', 'gfc 2 1');
            testCase.verifyError(@() esm2plmt(testCase.Scratch, range, ...
                'Component', 'H', 'Unit', 'GRAV'), 'ULMO:esm2plmt:InvalidCoefficients');
            range(2) = range(2) - minutes(1);
            testCase.verifyError(@() esm2plmt(testCase.Scratch, range), ...
                'ULMO:esm2plmt:InvalidIntervals');
        end
    end
    methods (Access = private)
        function file = fixture(testCase, step, amplitude)
            timestamp = testCase.Start + hours(6 * step);
            timestamp.Format = 'yyyyMMdd_HH';
            name = ['mtmshc_H_', char(timestamp)];
            file = fullfile(testCase.Scratch, [name, '.gfc']);
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
