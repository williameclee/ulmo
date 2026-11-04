classdef stericHaloCorrectionTest < matlab.unittest.TestCase
    % Verify native support, depth offsets and loader temporal/cache behavior.
    methods (Test)
        function nativeAndFixedStepCorrections(~)
            folder = tempname; mkdir(folder); previous = getenv('IFILES');
            cleanup = onCleanup(@() restore(folder, previous));
            setenv('IFILES', folder);
            productFolder = fullfile(folder, 'STERIC', 'EN4c14'); mkdir(productFolder);
            native = fullfile(productFolder, 'EN4c14-StericSeaLevel.mat');
            lon = [0, 120, 240]; lat = [-60; 30];
            dates = datetime(2005, 1, 15) + calmonths((0:9)');
            t = reshape(0:9, 1, 1, []);
            shallowHalostericSl = [1, 2, 3; 4, 5, 6] .* (.01 + .001 * t.^2);
            deepHalostericSl = [0, 2, 1; 0, 1, 3] .* (.02 + .002 * t);
            halostericSl = shallowHalostericSl + deepHalostericSl;
            % Excluded by the full-column 80% rule despite finite values later.
            halostericSl(1, 3, 1:3) = 0;
            shallowHalostericSl(1, 3, 1:3) = 0;
            deepHalostericSl = halostericSl - shallowHalostericSl;
            shallowStericSl = shallowHalostericSl + .2;
            deepStericSl = deepHalostericSl + .1;
            stericSl = shallowStericSl + deepStericSl;
            save(native, 'lon', 'lat', 'dates', 'stericSl', 'shallowStericSl', ...
                'deepStericSl', 'halostericSl', 'shallowHalostericSl', 'deepHalostericSl', '-v7.3');
            h = sterichalocorrection(native, [0; (1:10)'], true, true);
            assert(isnan(h.total(1)) && ~h.supportMask(1, 3));
            weights = cosd(lat) .* ones(2, 3); weights(1, 3) = 0;
            expected = sum(halostericSl .* weights, [1, 2]) / sum(weights, 'all');
            assert(max(abs(h.total(2:end) / 1000 - expected(:))) < 1e-12);
            assert(max(abs(h.total(2:end) - h.shallow(2:end) - h.deep(2:end))) < 1e-12);
            disabled = sterichalocorrection(native, (1:10)', true, false);
            assert(all(disabled.total == 0) && isequal(disabled.estimatedTotal, h.total(2:end)));

            % Native output and units use the same global correction.
            corrected = steric2lonlatt('EN4c14', [], [], [], BeQuiet=true);
            assert(max(abs(corrected - (stericSl - expected)), [], 'all') < 1e-12);
            shallow = steric2lonlatt('EN4c14', [], [], [], Depth='shallow', BeQuiet=true);
            deep = steric2lonlatt('EN4c14', [], [], [], Depth='deep', BeQuiet=true);
            assert(max(abs(corrected - shallow - deep), [], 'all') < 1e-12);
            mm = steric2lonlatt('EN4c14', [], [], [], Unit='mm', OutputFormat='ndgrid', BeQuiet=true);
            assert(max(abs(mm - permute(corrected * 1000, [2, 1, 3])), [], 'all') < 1e-9);

            % Native correction must be invariant to the spatial output mesh. Exercise
            % both processing and cached paths, time truncation and nonlinear time interpolation.
            limits = dates([2, 8]);
            for method = ["linear", "pchip"]
                mesh = 60;
                if method == "pchip", mesh = []; end
                layers = cell(1, 3);
                names = {'full', 'shallow', 'deep'};
                [offsets, target] = interptemporal(dates, reshape([h.total(2:end), h.shallow(2:end)]', 2, 1, []) / 1000, ...
                    days(20), method, TimeRange=limits, BeQuiet=true);
                expectedOffsets = {offsets(1, 1, :), offsets(2, 1, :), offsets(1, 1, :) - offsets(2, 1, :)};
                for iDepth = 1:3
                    args = {'Depth', names{iDepth}, 'Interpolation', method, 'TruncationBeforeInterpolation', true, 'BeQuiet', true};
                    [fresh, ~, actualDates] = steric2lonlatt('EN4c14', days(20), mesh, limits, args{:}, ForceNew=true);
                    cached = steric2lonlatt('EN4c14', days(20), mesh, limits, args{:});
                    raw = steric2lonlatt('EN4c14', days(20), mesh, limits, args{:}, Type='total');
                    assert(isequaln(fresh, cached) && isequal(actualDates(:), target(:)));
                    delta = raw - fresh - expectedOffsets{iDepth};
                    assert(any(isfinite(delta), 'all') && max(abs(delta), [], 'all', 'omitmissing') < 1e-12);
                    layers{iDepth} = raw - fresh;
                end
                closure = layers{1} - layers{2} - layers{3};
                assert(max(abs(closure), [], 'all', 'omitmissing') < 1e-12);
            end
            % Shallow may also be derived from full minus deep in the shared helper.
            data = load(native); data = rmfield(data, 'shallowHalostericSl');
            save(native, '-struct', 'data', '-v7.3');
            derived = sterichalocorrection(native, (1:10)', true, true);
            assert(max(abs(derived.shallow - h.shallow(2:end))) < 1e-10);
            derivedShallow = steric2lonlatt('EN4c14', [], [], [], Depth='shallow', BeQuiet=true);
            assert(max(abs(derivedShallow - shallow), [], 'all') < 1e-12);
            % SIO defaults still work without halosteric fields or deep coverage.
            sioFolder = fullfile(folder, 'STERIC', 'SIO'); mkdir(sioFolder);
            save(fullfile(sioFolder, 'SIO-StericSeaLevel.mat'), 'lon', 'lat', 'dates', 'stericSl');
            sio = steric2lonlatt('SIO', [], [], [], BeQuiet=true);
            assert(isequal(sio, stericSl));
        end

        function epochMeansUseNativeCorrectionAndFieldSupport(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            productFolder = fullfile(folder.Folder, 'STERIC', 'EN4c14');
            mkdir(productFolder);
            native = fullfile(productFolder, 'EN4c14-StericSeaLevel.mat');
            dates = (datetime(2020, 1, 1):days(3):datetime(2020, 2, 29))';
            lon = [0, 120, 240]; lat = [-60; 30];
            t = reshape(1:numel(dates), 1, 1, []);
            halostericSl = repmat(.001 * t, 2, 3, 1);
            stericSl = repmat(.01 * t.^2, 2, 3, 1);
            stericSl(1, 1, 1:3) = NaN;
            save(native, 'lon', 'lat', 'dates', 'stericSl', 'halostericSl', '-v7.3');
            epochs = temporalmeanepochs(dates, 'midmonth', []);
            [expected, target] = averageepochs(dates, stericSl - halostericSl, epochs);
            [fresh, ~, labels] = steric2lonlatt('EN4c14', 'midmonth', [], [], BeQuiet=true, ForceNew=true);
            cached = steric2lonlatt('EN4c14', 'midmonth', [], [], BeQuiet=true);
            testCase.verifyEqual(labels, target);
            testCase.verifyEqual(fresh, expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(cached, fresh);
            raw = steric2lonlatt('EN4c14', 'midmonth', [], [], Type='total', BeQuiet=true);
            effective = raw - fresh;
            testCase.verifyGreaterThan(abs(effective(1, 1, 1) - effective(2, 1, 1)), 1e-4);
            % Native halo means must not depend on the output mesh.
            correctedGrid = steric2lonlatt('EN4c14', 'midmonth', 60, [], BeQuiet=true, ForceNew=true);
            rawGrid = steric2lonlatt('EN4c14', [], 60, [], Type='total', BeQuiet=true, ForceNew=true);
            testCase.verifyEqual(correctedGrid, averageepochs(dates, rawGrid - .001 * t, epochs), 'AbsTol', 1e-12);
            selected = steric2lonlatt('EN4c14', 'midmonth', [], [target(1), target(1) + days(1)], ...
                TruncationBeforeInterpolation=true, BeQuiet=true);
            testCase.verifyEqual(selected, expected(:, :, 1), 'AbsTol', 1e-12);

            % GRACE uses the same correction-before-averaging operation.
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('ORIGINALGRACEDATA', folder.Folder));
            refFolder = fullfile(folder.Folder, 'RL06', 'CSR'); mkdir(refFolder);
            writelines(["# End of YAML header"; ...
                "GRCOF2 0 0 1 0 0 0 20200101.0000 20200201.0000 0"], ...
                fullfile(refFolder, 'GSM-2_test_BA01_0600'));
            graceEpochs = temporalmeanepochs(dates, 'GRACE', []);
            expectedGrace = averageepochs(dates, stericSl - halostericSl, graceEpochs);
            actualGrace = steric2lonlatt('EN4c14', 'GRACE', [], [], BeQuiet=true);
            testCase.verifyEqual(actualGrace, expectedGrace, 'AbsTol', 1e-12);
            testCase.verifyEqual(steric2lonlatt('EN4c14', 'GRACE', [], [], BeQuiet=true), actualGrace);
            % Completed epoch caches retain main's standalone-read behavior.
            delete(native);
            testCase.verifyEqual(steric2lonlatt('EN4c14', 'midmonth', [], [], BeQuiet=true), fresh);
        end
    end
end

function restore(folder, previous)
    setenv('IFILES', previous); rmdir(folder, 's');
end
