%% SEISMICHELPERSTEST - Checks readable identities and circular-cap geometry.
% Exercises coordinate canonicalisation, precise filename numbers, polar
% cap dispatch, retained-mode ordering and rotation to the epicentre.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

classdef seismicHelpersTest < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function setup(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root, IncludingSubfolders=true));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            old = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', old));
            setenv('IFILES', fixture.Folder);
            mkdir(fullfile(fixture.Folder, 'DLMB'));
        end
    end
    methods (Test)
        function numbersHaveReadableRoundTripRepresentation(testCase)
            testCase.verifyEqual(num2filename(121.48), '121p48');
            testCase.verifyEqual(num2filename(-30.72), 'n30p72');
            testCase.verifyEqual(num2filename(-0), '0');
            values = [pi, -1e-20, eps, realmin, realmax, 1.000000000000001, 20, 180];
            for value = values
                text = num2filename(value);
                testCase.verifyEqual(str2double(strrep(strrep(text, 'p', '.'), 'n', '-')), value);
                testCase.verifyEmpty(regexp(text, '[.+-]', 'once'));
            end
            testCase.verifyNotEqual(num2filename(121.480001), num2filename(121.480002));
        end
        function aliasesAndMixedCellsResolveInOrder(testCase)
            [xy, names, dates] = resolveearthquakes({'Maule', [142.373, 38.297], 'Indian Ocean'});
            testCase.verifyEqual(xy, [-72.733, -35.909; 142.373, 38.297; 95.982, 3.295]);
            testCase.verifyEqual(names, ["Chile"; "Tohoku"; "Sumatra"]);
            testCase.verifyEqual(dates, [datetime(2010, 2, 27); datetime(2011, 3, 11); datetime(2004, 12, 26)]);
            [~, names] = resolveearthquakes(["Tōhoku", "SUMATRA-2004"]);
            testCase.verifyEqual(names, ["Tohoku"; "Sumatra"]);
        end
        function wrappedLongitudesAndPolesHaveCanonicalNames(testCase)
            [xy, names] = resolveearthquakes([287.267, -35.909; 180, 0; -180, 0; 77, 90; -42, 90]);
            testCase.verifyEqual(names(1), "Chile");
            testCase.verifyEqual(names(2), names(3));
            testCase.verifyEqual(names(4), names(5));
            testCase.verifyEqual(xy(4:5, 1), [0; 0]);
        end
        function invalidEntriesAreRejected(testCase)
            cases = {{[]}, {NaN, 'Sumatra'}, {1}, {struct()}, {[0, 91]}, {""}, {missing}};
            for iCase = 1:numel(cases)
                testCase.verifyError(@() resolveearthquakes(cases{iCase}), ...
                    errorid(cases{iCase}));
            end
        end
        function polarBasisHasExpectedConcentrationAndOrthonormality(testCase)
            [G, V, ~, ~, N] = glmalpha_new(60, 4, SaveData=false, BeQuiet=true);
            testCase.verifyEqual(N, 25 * (1 - cosd(60)) / 2, AbsTol=1e-12);
            testCase.verifyEqual(G.' * G, eye(25), AbsTol=1e-12);
            testCase.verifyEqual(sum(V), N, AbsTol=1e-12);
            basis = seismiccapbasis([0, 90], 60, 4, [], SaveData=false);
            testCase.verifySize(basis{1}, [25, round(N)]);
            testCase.verifyGreaterThanOrEqual(-diff(basis{2}), zeros(round(N)-1, 1));
        end
        function capRotatesToEpicentreAndPreservesModeSpectrum(testCase)
            L = 4;
            centre = [121.48, -30.72];
            polar = seismiccapbasis([0, 90], 60, L, 4, SaveData=false);
            rotated = seismiccapbasis(centre, 60, L, 4, SaveData=false);
            testCase.verifyEqual(rotated{1}.' * rotated{1}, eye(4), AbsTol=1e-12);
            testCase.verifyEqual(rotated{2}, polar{2}, AbsTol=1e-12);
            orders = addmout(L);
            scale = sqrt(4*pi) * (-1).^orders(:);
            Y = ylm([0, L], [], deg2rad(90 - [centre(2); -centre(2)]), ...
                deg2rad([centre(1); centre(1) + 180]), [], [], 0, 1);
            amplitudes = (Y .* scale).' * rotated{1}(:, 1);
            testCase.verifyGreaterThan(abs(amplitudes(1)), 5 * abs(amplitudes(2)));
        end
        function fractionalRadiiDoNotSharePolarCache(testCase)
            [~, ~, ~, ~, first] = glmalpha_new(20.1, 4, SaveData=true, BeQuiet=true);
            [~, ~, ~, ~, second] = glmalpha_new(20.2, 4, SaveData=true, BeQuiet=true);
            testCase.verifyNotEqual(first, second);
            folder = fullfile(getenv('IFILES'), 'GLMALPHA');
            testCase.verifyTrue(isfile(fullfile(folder, 'glmalpha-20p1-4-1-0.mat')));
            testCase.verifyTrue(isfile(fullfile(folder, 'glmalpha-20p2-4-1-0.mat')));
        end
        function zeroDegreeBasisIsConstant(testCase)
            basis = seismiccapbasis([10, -20], 60, 0, [], SaveData=false);
            testCase.verifyEqual(abs(basis{1}), 1, AbsTol=1e-12);
            testCase.verifyEqual(basis{3}, 0.25, AbsTol=1e-12);
        end
    end
end

function id = errorid(input)
    if isequal(input, {""})
        id = 'ULMO:resolveearthquakes:UnknownName';
    else
        id = 'ULMO:resolveearthquakes:InvalidInput';
    end
end
