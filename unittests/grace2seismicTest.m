%% GRACE2SEISMICTEST - Checks independent fits, event identities and NetCDF caches.
% Uses synthetic GRACE harmonics and real source-window parsing, cap bases,
% temporal fits and NetCDF I/O. A loader fixture records refresh/save options.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/13, En-Chi Lee (williameclee@arizona.edu)

classdef grace2seismicTest < matlab.unittest.TestCase
    properties
        Root
        Scratch
        Epochs
        Dates
        Events = [datetime(2004, 12, 26); datetime(2010, 2, 27)]
        Product = {'CSR', 'RL06', 60}
        L = 4
    end
    methods (TestMethodSetup)
        function setup(testCase)
            testCase.Root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testCase.Root, IncludingSubfolders=true));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            testCase.Scratch = fixture.Folder;
            for name = ["IFILES", "GRACEDATA", "ORIGINALGRACEDATA"]
                old = getenv(name);
                testCase.addTeardown(@() setenv(name, old));
            end
            setenv('IFILES', testCase.Scratch);
            setenv('GRACEDATA', fullfile(testCase.Scratch, 'GRACE'));
            setenv('ORIGINALGRACEDATA', fullfile(testCase.Scratch, 'raw'));
            mkdir(fullfile(testCase.Scratch, 'DLMB'));
            raw = fullfile(getenv('ORIGINALGRACEDATA'), 'RL06', 'CSR');
            mkdir(raw);
            starts = datetime(2002, 1, 1) + calmonths((0:143).');
            starts(80) = []; % The missing month must never become an observation.
            testCase.Epochs = [starts, starts + calmonths(1)];
            testCase.Dates = mean(testCase.Epochs, 2);
            [~, ~, ~, template, ~, ~, ~, ~, ~, mapping] = addmon(testCase.L);
            nHarmonic = (testCase.L + 1) ^ 2;
            time = years(testCase.Dates - testCase.Dates(1));
            values = 1 + time * (1:nHarmonic) / 30 + ...
                sin(time * (1:nHarmonic) * 3.17) / 40 + ...
                (testCase.Dates >= testCase.Events(1)) * cos(1:nHarmonic) + ...
                (testCase.Dates >= testCase.Events(2)) * sin(1:nHarmonic);
            coefficients = zeros(numel(starts), 2 * size(template, 1));
            coefficients(:, mapping) = values;
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            ULMO_SEISMIC_FIXTURE = struct('calls', 0, 'dates', testCase.Dates);
            ULMO_SEISMIC_FIXTURE.plmt = cat(3, ...
                repmat(reshape(template(:, 1:2), 1, [], 2), numel(starts), 1, 1), ...
                reshape(coefficients, numel(starts), [], 2));
            ULMO_SEISMIC_FIXTURE.sigma = ULMO_SEISMIC_FIXTURE.plmt;
            ULMO_SEISMIC_FIXTURE.sigma(:, :, 3:4) = 0.01;
            testCase.addTeardown(@() clearfixture());
            copyfile(fullfile(testCase.Root, 'unittests', 'fixtures', 'seismic', ...
                'grace2plmt_new.m.txt'), fullfile(testCase.Scratch, 'grace2plmt_new.m'));
            clear grace2plmt_new
            testCase.addTeardown(@() clearloader());
            for iDate = 1:numel(starts)
                path = fullfile(raw, sprintf('GSM-2_%03d_BA01_0600', iDate));
                fid = fopen(path, 'w');
                fprintf(fid, 'EARTH 3.986004415e14 6378136.3\n');
                for iRow = 1:size(template, 1)
                    fprintf(fid, 'GRCOF2 %d %d 0 0 1e-12 1e-12 %s %s 0\n', ...
                        template(iRow, 1:2), ...
                        string(testCase.Epochs(iDate, 1), 'yyyyMMdd.HHmm'), ...
                        string(testCase.Epochs(iDate, 2), 'yyyyMMdd.HHmm'));
                end
                fclose(fid);
            end
        end
    end
    methods (Test)
        function excludedDegreesCannotAffectFitOrUncertainty(testCase)
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            [a, sa, ~, ea] = testCase.call('Sumatra', 60, testCase.Events(1), SaveData=false);
            ULMO_SEISMIC_FIXTURE.plmt(:, 1:3, 3:4) = NaN;
            [b, sb, ~, eb] = testCase.call('Sumatra', 60, testCase.Events(1), SaveData=false);
            testCase.verifyEqual(b, a);
            testCase.verifyEqual(sb, sa);
            testCase.verifyEqual(eb, ea);
            testCase.verifyEqual(a(:, 1:3, 3:4), zeros(numel(testCase.Dates), 3, 2));
            testCase.verifyEqual(sa(:, 1:3, 3:4), zeros(numel(testCase.Dates), 3, 2));
            testCase.verifyEqual(ea(:, 1:3, 3:4), zeros(1, 3, 2));
        end
        function maximumDegreeMustBeAtLeastTwo(testCase)
            for degree = [0, 1]
                testCase.verifyError(@() grace2seismic(testCase.Product, 'Sumatra', 60, ...
                    testCase.Events(1), degree), 'MATLAB:validators:mustBeGreaterThanOrEqual');
            end
            [p, s] = grace2seismic(testCase.Product, 'Sumatra', 60, ...
                testCase.Events(1), 2, SaveData=false, BeQuiet=true);
            testCase.verifySize(p, [numel(testCase.Dates), 6, 4]);
            testCase.verifyEqual(s(:, 1:3, 3:4), zeros(numel(testCase.Dates), 3, 2));
        end
        function matchingProductDegreeOmitsCacheSuffix(testCase)
            % Seed a degree-60 cache to test naming/reuse without a large cap solve.
            [~, ~, ~, field] = addmon(60);
            field(:, 3:4) = 0;
            field(end, 3) = 2;
            sigma = field;
            sigma(:, 3:4) = 0;
            folder = fullfile(getenv('GRACEDATA'), 'Seismic');
            mkdir(folder);
            path = fullfile(folder, 'grace2seismic-CSR0660-Sumatra_60-SD.nc');
            writeshnc(path, field, sigma, testCase.Events(1), 'SD');
            [p, s, ~, event] = grace2seismic(testCase.Product, 'Sumatra', 60, ...
                testCase.Events(1), 60, BeQuiet=true);
            testCase.verifyEqual(squeeze(event), field);
            testCase.verifyEqual(p(end, :, 3:4), reshape(field(:, 3:4), 1, [], 2));
            testCase.verifyEqual(s(:, :, 3:4), zeros(numel(testCase.Dates), size(field, 1), 2));
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.calls, 0);
            testCase.verifyNumElements(dir(fullfile(folder, '*.nc')), 1);
        end
        function uncachedFitsRemainIndependent(testCase)
            [combined, sigma, ~, events] = grace2seismic(testCase.Product, {'Sumatra', 'Chile'}, ...
                60, testCase.Events, testCase.L, ForceNew=true, SaveData=false, BeQuiet=true);
            [a, sa, ~, ea] = testCase.call('Sumatra', 60, testCase.Events(1), ForceNew=true, SaveData=false);
            [b, sb, ~, eb] = testCase.call('Chile', 60, testCase.Events(2), ForceNew=true, SaveData=false);
            testCase.verifyEqual(combined(:, :, 3:4), a(:, :, 3:4) + b(:, :, 3:4), AbsTol=1e-13);
            testCase.verifyEqual(sigma(:, :, 3:4).^2, sa(:, :, 3:4).^2 + sb(:, :, 3:4).^2, AbsTol=1e-13);
            testCase.verifyEqual(events, cat(1, ea, eb), AbsTol=1e-13);
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.calls, 3);
        end
        function unitsAreForwardedAndHaveSeparateCaches(testCase)
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            for unit = ["SD", "POT", "GRAV"]
                [~, ~, ~, field] = testCase.call('Sumatra', 60, testCase.Events(1), Unit=unit);
                path = testCase.cache("Sumatra_60-L4-" + unit);
                testCase.verifyTrue(isfile(path));
                testCase.verifyEqual(string(ncreadatt(path, '/', 'coefficientUnit')), unit);
                args = ULMO_SEISMIC_FIXTURE.arguments;
                testCase.verifyEqual(string(args{2}), unit);
                testCase.verifyTrue(all(isfinite(field), 'all'));
            end
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.calls, 3);
        end
        function realLoaderConsumesSyntheticSourceFiles(testCase)
            delete(fullfile(testCase.Scratch, 'grace2plmt_new.m'));
            clear grace2plmt_new
            mkdir(getenv('GRACEDATA'));
            constants = fullfile(testCase.Scratch, 'EARTHMODELS', 'CONSTANTS');
            mkdir(constants);
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(constants, 'Earth.mat'), 'Earth');
            [p, sigma, dates] = testCase.call('Sumatra', 60, testCase.Events(1), ...
                Unit='GRAV', Deg1Correction=false, C20Correction=false, C30Correction=false, ...
                ForceNew=true, SaveData=false, TimeFormat='datetime');
            testCase.verifyEqual(dates, testCase.Dates);
            testCase.verifyLessThan(max(abs(p(:, :, 3:4)), [], 'all'), 1e-14);
            testCase.verifyTrue(all(isfinite(sigma), 'all'));
            testCase.verifyFalse(isfolder(fullfile(getenv('GRACEDATA'), 'Seismic')));
        end
        function gldasIsSubtractedBeforeFittingAndCachedSeparately(testCase)
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            ULMO_SEISMIC_FIXTURE.storage = ULMO_SEISMIC_FIXTURE.plmt;
            ULMO_SEISMIC_FIXTURE.storage(:, :, 3:4) = 0.25 * ULMO_SEISMIC_FIXTURE.plmt(:, :, 3:4);
            ULMO_SEISMIC_FIXTURE.storage(:, 1:3, 3:4) = NaN;
            ULMO_SEISMIC_FIXTURE.gldasCalls = 0;
            copyfile(fullfile(testCase.Root, 'unittests', 'fixtures', 'seismic', ...
                'seismicgldas.m.txt'), fullfile(testCase.Scratch, 'seismicgldas.m'));
            clear seismicgldas
            testCase.addTeardown(@() cleargldas());
            [original, originalStd, ~, originalEvent] = testCase.call('Sumatra', 60, testCase.Events(1));
            [corrected, correctedStd, ~, correctedEvent] = testCase.call('Sumatra', 60, ...
                testCase.Events(1), RemoveGLDAS=true);
            testCase.verifyEqual(corrected(:, :, 3:4), 0.75 * original(:, :, 3:4), AbsTol=1e-12);
            testCase.verifyEqual(correctedStd(:, :, 3:4), 0.75 * originalStd(:, :, 3:4), AbsTol=1e-12);
            testCase.verifyEqual(correctedEvent(:, :, 3:4), 0.75 * originalEvent(:, :, 3:4), AbsTol=1e-12);
            testCase.verifyTrue(isfile(testCase.cache('Sumatra_60-L4-SD-GLDAS')));
            testCase.verifyTrue(contains(ncreadatt(testCase.cache('Sumatra_60-L4-SD-GLDAS'), ...
                '/', 'gldasProduct'), 'GLDAS_NOAH10_M.2.1'));
            testCase.call('Sumatra', 60, testCase.Events(1), RemoveGLDAS=true);
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.gldasCalls, 1);
        end
        function mixedEventsAreIndependentAndAdditive(testCase)
            [combined, uncertainty, dates, events] = grace2seismic(testCase.Product, ...
                {'Sumatra', [-72.733, -35.909]}, [60, 70], testCase.Events, testCase.L, BeQuiet=true);
            [first, firstStd, ~, firstEvent] = testCase.call('Sumatra', 60, testCase.Events(1));
            [second, secondStd, ~, secondEvent] = testCase.call('Chile', 70, testCase.Events(2));
            testCase.verifyEqual(combined(:, :, 3:4), first(:, :, 3:4) + second(:, :, 3:4), AbsTol=1e-13);
            testCase.verifyEqual(uncertainty(:, :, 3:4) .^ 2, ...
                firstStd(:, :, 3:4) .^ 2 + secondStd(:, :, 3:4) .^ 2, AbsTol=1e-13);
            testCase.verifyEqual(events, cat(1, firstEvent, secondEvent));
            testCase.verifyEqual(dates, datenum(testCase.Dates)); %#ok<DATNM>
            [reversed, ~, ~, reversedEvents] = grace2seismic(testCase.Product, ...
                {'Chile', 'Sumatra'}, [70, 60], flip(testCase.Events), testCase.L, BeQuiet=true);
            testCase.verifyEqual(reversed, combined, AbsTol=1e-13);
            testCase.verifyEqual(reversedEvents, flip(events, 1));
        end
        function defaultFitMatchesDirectBasisCalculation(testCase)
            [~, ~, ~, event] = testCase.call('Sumatra', 60, testCase.Events(1));
            GVN = seismiccapbasis([95.982, 3.295], 60, testCase.L, [], SaveData=false);
            G = GVN{1};
            [~, ~, ~, template, ~, ~, ~, ~, ~, mapping] = addmon(testCase.L);
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            values = reshape(ULMO_SEISMIC_FIXTURE.plmt(:, :, 3:4), numel(testCase.Dates), []);
            [step, sigma] = slept2seismic(values(:, mapping) * G, testCase.Dates, testCase.Events(1));
            actual = reshape(event(1, :, 3:4), [], 2);
            testCase.verifyEqual(actual(mapping), G * step.', AbsTol=1e-12);
            [cached, cachedSigma, date] = readshnc(testCase.cache('Sumatra_60-L4-SD'));
            s = cachedSigma(:, 3:4);
            testCase.verifyEqual(s(mapping), sqrt((G .^ 2) * sigma.'.^2), AbsTol=1e-12);
            testCase.verifyEqual(cached(:, 1:2), template(:, 1:2));
            testCase.verifyEqual(date, testCase.Events(1));
            testCase.verifyEqual(ncreadatt(testCase.cache('Sumatra_60-L4-SD'), '/', 'retainedModeCount'), size(G, 2));
        end
        function outputSelectionAndLayoutsReuseStaticCache(testCase)
            [all, allStd, ~, allEvents] = testCase.call('Sumatra', 60, testCase.Events(1));
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            calls = ULMO_SEISMIC_FIXTURE.calls;
            selected = [20, 30];
            [part, partStd, dates, events] = testCase.call('Sumatra', 60, testCase.Events(1), ...
                TimeRange=testCase.Dates(selected).', OutputFormat='traditional', TimeFormat='datetime');
            testCase.verifyEqual(part, permute(all(20:30, :, :), [2, 3, 1]));
            testCase.verifyEqual(partStd, permute(allStd(20:30, :, :), [2, 3, 1]));
            testCase.verifyEqual(events, permute(allEvents, [2, 3, 1]));
            testCase.verifyEqual(dates, testCase.Dates(20:30));
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.calls, calls);
            [empty, emptyStd, dates] = testCase.call('Sumatra', 60, testCase.Events(1), ...
                TimeRange=[datetime(1990, 1, 1), datetime(1991, 1, 1)]);
            testCase.verifySize(empty, [0, 15, 4]);
            testCase.verifySize(emptyStd, [0, 15, 4]);
            testCase.verifyEmpty(dates);
        end
        function namesAliasesAndCoordinatesShareCache(testCase)
            expected = testCase.call('Sumatra', 60, testCase.Events(1));
            actual = testCase.call([95.982, 3.295], 60, testCase.Events(1));
            alias = testCase.call('Indian Ocean', 60, datenum(testCase.Events(1))); %#ok<DATNM>
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(alias, expected);
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.calls, 1);
            testCase.verifyNumElements(dir(fullfile(getenv('GRACEDATA'), 'Seismic', '*.nc')), 1);
        end
        function namedOverrideAndCoordinateFilenames(testCase)
            testCase.call([121.48, -30.72], 20, datetime(2009, 10, 27, 12, 0, 0));
            testCase.verifyTrue(isfile(testCase.cache('121p48_n30p72_20-T20091028-L4-SD')));
            testCase.call('Sumatra', 20, datetime(2004, 12, 26, 11, 59, 59));
            testCase.verifyTrue(isfile(testCase.cache('Sumatra_20-L4-SD')));
            testCase.call('Sumatra', 20, datetime(2004, 12, 26, 12, 0, 0));
            testCase.verifyTrue(isfile(testCase.cache('Sumatra_20-T20041227-L4-SD')));
        end
        function intervalWeightsIncludeExactBoundaries(testCase)
            event = datetime(2004, 12, 26);
            [p, s, ~, a] = testCase.call('Sumatra', 60, event);
            index = find(testCase.Epochs(:, 1) == datetime(2004, 12, 1));
            weight = 6 / 31;
            testCase.verifyEqual(p(index, :, 3:4), weight * a(:, :, 3:4), AbsTol=1e-13);
            testCase.verifyEqual(p(index - 1, :, 3:4), zeros(1, 15, 2));
            testCase.verifyEqual(p(index + 1, :, 3:4), a(:, :, 3:4));
            testCase.verifyEqual(s(index, :, 3:4), weight * s(index + 1, :, 3:4), AbsTol=1e-13);
            [boundary, ~, ~, jump] = testCase.call([121, 20], 60, datetime(2005, 1, 1));
            testCase.verifyEqual(boundary(index, :, 3:4), zeros(1, 15, 2));
            testCase.verifyEqual(boundary(index + 1, :, 3:4), jump(:, :, 3:4));
        end
        function incompleteCacheRegeneratesInPlace(testCase)
            [~, ~, ~, event] = testCase.call('Sumatra', 60, testCase.Events(1));
            path = testCase.cache('Sumatra_60-L4-SD');
            writeshnc(path, squeeze(event), [], testCase.Events(1), 'SD');
            testCase.call('Sumatra', 60, testCase.Events(1));
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            testCase.verifyEqual(ULMO_SEISMIC_FIXTURE.calls, 2);
            info = ncinfo(path);
            testCase.verifyTrue(ismember('uncertainty', {info.Variables.Name}));
            testCase.verifyNumElements(dir(fullfile(getenv('GRACEDATA'), 'Seismic', '*.nc')), 1);
        end
        function forcedRefreshAndSaveDisabled(testCase)
            original = testCase.call('Sumatra', 60, testCase.Events(1));
            path = testCase.cache('Sumatra_60-L4-SD');
            contents = readbytes(path);
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            ULMO_SEISMIC_FIXTURE.plmt(:, :, 3:4) = 2 * ULMO_SEISMIC_FIXTURE.plmt(:, :, 3:4);
            fresh = testCase.call('Sumatra', 60, testCase.Events(1), ForceNew=true, SaveData=false);
            testCase.verifyEqual(fresh(:, :, 3:4), 2 * original(:, :, 3:4), AbsTol=1e-12);
            testCase.verifyEqual(readbytes(path), contents);
            args = ULMO_SEISMIC_FIXTURE.arguments;
            testCase.verifyFalse(args{find(strcmp(args(1:2:end), 'SaveData')) * 2});
            testCase.verifyTrue(args{find(strcmp(args(1:2:end), 'ForceNew')) * 2});
            testCase.verifyEqual(testCase.call('Sumatra', 60, testCase.Events(1)), original);
            fresh = testCase.call('Sumatra', 60, testCase.Events(1), ForceNew=true);
            testCase.verifyEqual(testCase.call('Sumatra', 60, testCase.Events(1)), fresh);
        end
        function noSaveDoesNotCreateEventOrBasisCache(testCase)
            testCase.call('Sumatra', 60, testCase.Events(1), SaveData=false);
            testCase.verifyFalse(isfolder(fullfile(getenv('GRACEDATA'), 'Seismic')));
            testCase.verifyEmpty(dir(fullfile(testCase.Scratch, 'GLMALPHA', '*.mat')));
            testCase.verifyFalse(isfolder(fullfile(testCase.Scratch, 'SDWCAP')));
        end
        function nondefaultSettingsRequireUncachedMode(testCase)
            alternatives = {{'polynomialDegree', 1}, {'periods', []}, {'truncation', 3}, ...
                {'FitRange', testCase.Dates([10, end]).'}, {'Deg1Correction', false}, ...
                {'LoveNumSource', 'Wahr'}};
            for iOption = 1:numel(alternatives)
                args = alternatives{iOption};
                testCase.verifyError(@() testCase.call('Sumatra', 60, testCase.Events(1), args{:}), ...
                    'ULMO:grace2seismic:NondefaultCache');
            end
            actual = testCase.call('Sumatra', 60, testCase.Events(1), ...
                polynomialDegree=1, periods=[], truncation=3, ForceNew=true, SaveData=false);
            testCase.verifyTrue(all(isfinite(actual), 'all'));
            testCase.verifyFalse(isfolder(fullfile(getenv('GRACEDATA'), 'Seismic')));
        end
        function eventPairingAndInvalidCells(testCase)
            testCase.verifyError(@() testCase.call({'Sumatra', [1, 2]}, 60, testCase.Events(1)), ...
                'ULMO:grace2seismic:EventCount');
            testCase.verifyError(@() testCase.call('Sumatra', [10, 20], testCase.Events(1)), ...
                'ULMO:grace2seismic:EventCount');
            testCase.verifyError(@() testCase.call({[1; 2]}, 60, testCase.Events(1)), ...
                'ULMO:resolveearthquakes:InvalidInput');
            testCase.verifyError(@() testCase.call([1, 100], 60, testCase.Events(1)), ...
                'ULMO:resolveearthquakes:InvalidInput');
            testCase.verifyError(@() testCase.call('unknown', 60, testCase.Events(1)), ...
                'ULMO:resolveearthquakes:UnknownName');
        end
        function fallbackRootAndProductNormalisation(testCase)
            setenv('GRACEDATA', '');
            grace2seismic({'csr', 6, 60}, 'Sumatra', 60, testCase.Events(1), 4, BeQuiet=true);
            testCase.verifyTrue(isfile(fullfile(testCase.Scratch, 'GRACE', 'Seismic', ...
                'grace2seismic-CSR0660-Sumatra_60-L4-SD.nc')));
        end
        function coefficientDatesMustMatchActualWindows(testCase)
            global ULMO_SEISMIC_FIXTURE %#ok<GVMIS>
            ULMO_SEISMIC_FIXTURE.dates(1) = ULMO_SEISMIC_FIXTURE.dates(1) + days(1);
            testCase.verifyError(@() testCase.call('Sumatra', 60, testCase.Events(1)), ...
                'ULMO:grace2seismic:Dates');
        end
        function missingGldasCannotReuseUncorrectedCache(testCase)
            testCase.call('Sumatra', 60, testCase.Events(1));
            old = getenv('GLDASDATA');
            testCase.addTeardown(@() setenv('GLDASDATA', old));
            setenv('GLDASDATA', fullfile(testCase.Scratch, 'missing'));
            testCase.verifyError(@() testCase.call('Sumatra', 60, testCase.Events(1), RemoveGLDAS=true), ...
                'ULMO:seismicgldas:MissingData');
        end
    end
    methods
        function varargout = call(testCase, earthquakes, radius, dates, varargin)
            [varargout{1:nargout}] = grace2seismic(testCase.Product, earthquakes, radius, dates, ...
                testCase.L, 'BeQuiet', true, varargin{:});
        end
        function path = cache(~, stem)
            path = fullfile(getenv('GRACEDATA'), 'Seismic', "grace2seismic-CSR0660-" + stem + ".nc");
        end
    end
end

function clearfixture()
    clear global ULMO_SEISMIC_FIXTURE
end
function clearloader()
    clear grace2plmt_new
end
function cleargldas()
    clear seismicgldas
end
function data = readbytes(path)
    fid = fopen(path, 'r');
    cleanup = onCleanup(@() fclose(fid));
    data = fread(fid, Inf, '*uint8');
end
