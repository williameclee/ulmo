%% SEISMICHELPERSTEST - Checks readable identities and circular-cap geometry.
% Exercises coordinate canonicalisation, precise filename numbers, polar
% cap dispatch, retained-mode ordering and rotation to the epicentre.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/13, En-Chi Lee (williameclee@arizona.edu)

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
            basis = seismiccapbasis([0, 90], 60, 4, [], SaveData=false);
            N = 21 * (1 - cosd(60)) / 2;
            testCase.verifyEqual(basis{3}, N, AbsTol=1e-12);
            testCase.verifySize(basis{1}, [25, round(N)]);
            testCase.verifyEqual(basis{1}.' * basis{1}, eye(round(N)), AbsTol=1e-12);
            testCase.verifyEqual(basis{1}(1:4, :), zeros(4, round(N)));
            testCase.verifyGreaterThanOrEqual(-diff(basis{2}), zeros(round(N)-1, 1));
            % Integrate independently to verify the restricted concentration eigenproblem.
            [weights, nodes] = gausslegendrecof(20);
            z = (1 + cosd(60))/2 + (1 - cosd(60))/2 * nodes(:);
            weights = weights(:) * (1 - cosd(60))/2;
            longitude = (0:19).' * 2*pi/20;
            [theta, phi] = ndgrid(acos(z), longitude);
            orders = addmout(4);
            Y = ylm([0, 4], [], theta(:), phi(:), [], [], 0, 1);
            Y = Y .* ((-1).^orders(:));
            D = (Y .* repmat(weights, 20, 1).') * Y.' * 2*pi/20;
            G = basis{1}(5:end, :);
            testCase.verifyEqual(D(5:end, 5:end) * G, G .* basis{2}.', AbsTol=1e-11);
        end
        function capRotatesToEpicentreAndPreservesModeSpectrum(testCase)
            L = 4;
            centre = [121.48, -30.72];
            polar = seismiccapbasis([0, 90], 60, L, 4, SaveData=false);
            rotated = seismiccapbasis(centre, 60, L, 4, SaveData=false);
            testCase.verifyEqual(rotated{1}.' * rotated{1}, eye(4), AbsTol=1e-12);
            testCase.verifyEqual(rotated{2}, polar{2}, AbsTol=1e-12);
            testCase.verifyEqual(rotated{1}(1:4, :), zeros(4, 4));
            % Integrate the rotated modes over the actual epicentral cap.
            [weights, z] = gausslegendrecof(2*L, [], [cosd(60), 1]);
            azimuth = (0:19) * 2*pi/20;
            east = [-sind(centre(1)); cosd(centre(1)); 0];
            north = [-sind(centre(2))*cosd(centre(1)); ...
                -sind(centre(2))*sind(centre(1)); cosd(centre(2))];
            axis = [cosd(centre(2))*cosd(centre(1)); ...
                cosd(centre(2))*sind(centre(1)); sind(centre(2))];
            [z, azimuth] = ndgrid(z, azimuth);
            points = axis*z(:).' + east*(sqrt(1-z(:).^2).*cos(azimuth(:))).' + ...
                north*(sqrt(1-z(:).^2).*sin(azimuth(:))).';
            Y = ylm([0, L], [], acos(points(3, :)), atan2(points(2, :), points(1, :)), ...
                [], [], 0, 1);
            orders = addmout(L);
            values = (Y .* ((-1).^orders(:))).' * rotated{1};
            concentration = values.' * (values .* repmat(weights(:), 20, 1)) * 2*pi/20;
            testCase.verifyEqual(concentration, diag(rotated{2}), AbsTol=1e-11);
        end
        function fractionalRadiiDoNotSharePolarCache(testCase)
            [~, ~, ~, ~, first] = glmalpha_new(20.1, 4, SaveData=true, BeQuiet=true);
            [~, ~, ~, ~, second] = glmalpha_new(20.2, 4, SaveData=true, BeQuiet=true);
            testCase.verifyNotEqual(first, second);
            folder = fullfile(getenv('IFILES'), 'GLMALPHA');
            testCase.verifyTrue(isfile(fullfile(folder, 'glmalpha-20p1-4-1-0.mat')));
            testCase.verifyTrue(isfile(fullfile(folder, 'glmalpha-20p2-4-1-0.mat')));
        end
        function bandpassWithZeroMinimumMatchesLowpass(testCase)
            [G, V, ~, ~, N] = glmalpha_new(60, [0, 4], SaveData=false, BeQuiet=true);
            [lowpass, lowpassV, ~, ~, lowpassN] = glmalpha_new(60, 4, SaveData=false, BeQuiet=true);
            testCase.verifyEqual(G.' * G, eye(25), AbsTol=1e-12);
            testCase.verifyEqual(sort(V), sort(lowpassV), AbsTol=1e-12);
            testCase.verifyEqual(N, lowpassN, AbsTol=1e-12);
            testCase.verifyEqual(G * diag(V) * G.', ...
                lowpass * diag(lowpassV) * lowpass.', AbsTol=1e-12);
            testCase.verifyFalse(isfolder(fullfile(getenv('IFILES'), 'SDWCAP')));
        end
        function bandpassCacheReusesRestrictedBasis(testCase)
            first = seismiccapbasis([0, 90], 60, 4, [], SaveData=true);
            second = seismiccapbasis([0, 90], 60, 4, [], SaveData=false);
            testCase.verifyEqual(second, first);
            testCase.verifyNumElements(dir(fullfile(getenv('IFILES'), 'GLMALPHA', '*.mat')), 1);
            testCase.verifyFalse(isfolder(fullfile(getenv('IFILES'), 'SDWCAP')));
        end
        function minimumDegreeAndMaximumTruncation(testCase)
            basis = seismiccapbasis([10, -20], 60, 2, 5, SaveData=false);
            testCase.verifySize(basis{1}, [9, 5]);
            testCase.verifyEqual(basis{1}(1:4, :), zeros(4, 5));
            testCase.verifyEqual(basis{1}.' * basis{1}, eye(5), AbsTol=1e-12);
            testCase.verifyEqual(basis{3}, 1.25, AbsTol=1e-12);
            testCase.verifyError(@() seismiccapbasis([0, 90], 60, 2, 6, SaveData=false), ...
                'ULMO:prepareslepianbasis:Truncation');
            for L = [0, 1]
                testCase.verifyError(@() seismiccapbasis([0, 90], 60, L, [], SaveData=false), ...
                    'MATLAB:validators:mustBeGreaterThanOrEqual');
            end
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
