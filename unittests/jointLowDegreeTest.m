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

        function separateCacheAndConversion(testCase)
            folder = tempname;
            mkdir(folder);
            original = getenv('GRACEDATA');
            cleanup = onCleanup(@() restoreFolder(original, folder));
            setenv('GRACEDATA', folder);
            product = {'CSR', 'RL06', 60};
            baseId = solvedegree1(product, 'OnlyId', true);
            jointId = solvedegree1(product, 'IncludeC21S21', true, 'OnlyId', true);
            testCase.verifyNotEqual(baseId, jointId);
            coeffs = [1:5; 6:10];
            coeffStds = ones(2, 5);
            dates = datetime(2003, [1; 2], 15);
            save(fullfile(folder, 'Degree1', 'new', ['CSRRL0660_', jointId, '-SD.mat']), ...
                'coeffs', 'coeffStds', 'dates');
            [actual, uncertainty, actualDates] = solvedegree1(product, ...
                'IncludeC21S21', true, 'Unit', 'GRAV', 'BeQuiet', true);
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
