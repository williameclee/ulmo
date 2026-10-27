%% jointLowDegreeTest - Checks harmonic ordering and cached joint-solution unit conversion.
%
% Created by
%   2026/10/26, En-Chi Lee (williameclee@gmail.com)

classdef jointLowDegreeTest < matlab.unittest.TestCase

    methods (Test)

        function kernelOrdering(testCase)

            for c20 = [false, true]

                for c30 = [false, true]

                    for c21 = [false, true]
                        [locs, ids, degrees] = lowdegreecoefficients(c20, c30, c21);
                        % Independently index a labelled lmcosi array in kernel order.
                        labels = reshape(1:20, 10, 2);
                        packed = labels(kernelorder(3));
                        linear = sub2ind(size(labels), locs(:, 1), locs(:, 2));
                        testCase.verifyEqual(packed(ids(:)), labels(linear));
                        testCase.verifyEqual(numel(degrees), 3 + c20 + c30 + 2 * c21);
                    end

                end

            end

        end

        function cachedJointUnitConversion(testCase)
            folder = tempname;
            mkdir(folder);
            original = getenv('GRACEDATA');
            cleanup = onCleanup(@() restoreFolder(original, folder));
            setenv('GRACEDATA', folder);
            previousData = getenv('IFILES');
            testCase.addTeardown(@() setenv('IFILES', previousData));
            setenv('IFILES', folder);
            mkdir(fullfile(folder, 'EARTHMODELS', 'CONSTANTS'));
            Earth = struct('GM_EGM96', 3.986004415e14, 'a_EGM96', 6378136.3);
            save(fullfile(folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
            mkdir(fullfile(folder, 'LOVENUMS'));
            loadinggravitationalpotential = [0; .021; -.3];
            save(fullfile(folder, 'LOVENUMS', 'lovenumbers-ISSM.mat'), ...
                'loadinggravitationalpotential');
            product = {'CSR', 'RL06', 60};
            domain = GeoDomain('alloceans', 'Buffer', 0.5);
            cacheFolder = fullfile(folder, 'Degree1', 'new');
            mkdir(cacheFolder);
            % Seed a known cache fixture to exercise the real unit conversion.
            cacheName = sprintf('CSRRL0660_ICE6GD-WC21S21-Ld60_Ls96-%s-SD.mat', domain.Id);
            coeffs = [1:5; 6:10];
            coeffStds = ones(2, 5);
            dates = datetime(2003, [1; 2], 15);
            save(fullfile(cacheFolder, cacheName), ...
                'coeffs', 'coeffStds', 'dates');
            [actual, uncertainty, actualDates] = solvedegree1(product, ...
                'Solver', 'iteration', 'OceanDomain', domain, 'IncludeC21S21', true, 'Unit', 'GRAV', 'BeQuiet', true);
            expected = [convertgravity(coeffs(:, 1:3), 'SD', 'GRAV', 'InputFormat', 'L', 'L', 1), ...
                            convertgravity(coeffs(:, 4:5), 'SD', 'GRAV', 'InputFormat', 'L', 'L', 2)];
            testCase.verifyEqual(actual, expected, 'RelTol', 1e-14);
            testCase.verifyEqual(actualDates, dates);
            testCase.verifyGreaterThan(uncertainty, zeros(2, 5));
        end

    end

end

function restoreFolder(original, folder)
    setenv('GRACEDATA', original);
    rmdir(folder, 's');
end
