%% graceSourceEpochTest - Check GRACE source epochs across header formats.
%
% Created by
%   2026/11/01

classdef graceSourceEpochTest < matlab.unittest.TestCase

    properties (TestParameter)
        invalidYaml = struct( ...
            'missingGM', ["  earth_gravity_param:"; "    units: m^3/s^2"; ...
                          "  mean_equator_radius:"; "    value: 6378136.3"], ...
            'missingRadius', ["  earth_gravity_param:"; "    value: 3.986004415e14"; ...
                              "  mean_equator_radius:"; "  other_attribute:"; "    value: 42"], ...
            'invalidGM', ["  earth_gravity_param:"; "    value: NaN"; ...
                          "  mean_equator_radius:"; "    value: 6378136.3"])
        header = struct( ...
            'legacy', ["FIRST GRACE SHM"; "EARTH 3.986004415D+14 6.3781363D+06"; "SHM 0 0"], ...
            'yaml', ["header:"; "  earth_gravity_param:"; "    value: 3.986004415e14"; ...
                     "  mean_equator_radius:"; "    value: 6378136.3"; "# End of YAML header"])
    end

    methods (Test)
        function rejectsInvalidYamlReferenceParameters(testCase, invalidYaml)
            % Each requested parameter needs a finite positive value in its own block.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            sourceFile = fullfile(folder.Folder, 'GSM.txt');
            writelines(["header:"; invalidYaml; "# End of YAML header"; ...
                        "GRCOF2 0 0 1 0 0 0 20020404.0000 20020501.0000 nnnn"], sourceFile);
            testCase.verifyError(@() readReferenceParameters(sourceFile), ...
                'ULMO:parsegracesourcefile:InvalidHeader');
        end

        function padsJplFullField(testCase)
            % A JPL file starting at degree 2 needs a unit C00 and zero uncertainties below degree 2.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            sourceFile = fullfile(folder.Folder, 'GSM.txt');
            records = compose("GRCOF2 2 %d 0.5 0 0.25 0 20020404.0000 20020501.0000 nnnn", (0:2)');
            writelines(["# End of YAML header"; records], sourceFile);

            [gravity, uncertainty] = parsegracesourcefile(sourceFile);
            testCase.verifyEqual(gravity(1:3, :), [0 0 1 0; 1 0 0 0; 1 1 0 0]);
            testCase.verifyEqual(uncertainty(1:3, :), [0 0 0 0; 1 0 0 0; 1 1 0 0]);
        end

        function readsEpochAndReferenceParameters(testCase, header)
            % Exercise the four-output epoch call and the six-output gravity call.
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            sourceFile = fullfile(folder.Folder, 'GSM.txt');
            record = "GRCOF2 0 0 1 0 0.25 0 20020404.1100 20020501.0000 nnnn";
            writelines([header; record], sourceFile);

            [~, ~, midpoint, bounds] = parsegracesourcefile(sourceFile);
            expected = [datetime(2002, 4, 4, 11, 0, 0), datetime(2002, 5, 1)];
            testCase.verifyEqual(bounds, expected);
            testCase.verifyEqual(midpoint, mean(expected));

            [gravity, uncertainty, ~, ~, gm, radius] = parsegracesourcefile(sourceFile);
            testCase.verifyEqual(gravity, [0, 0, 1, 0]);
            testCase.verifyEqual(uncertainty, [0, 0, 0.25, 0]);
            testCase.verifyEqual(gm, 3.986004415e14);
            testCase.verifyEqual(radius, 6378136.3);
        end
    end
end

function readReferenceParameters(sourceFile)
    [~, ~, ~, ~, ~, ~] = parsegracesourcefile(sourceFile);
end
