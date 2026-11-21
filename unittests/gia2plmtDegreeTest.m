%% GIA2PLMTDEGREETEST - Verifies degree selection for native GIA model files.
% Uses temporary surface-density models with and without uncertainty bounds.
%
% Created by
%   2026/11/09, En-Chi Lee (williameclee@gmail.com)

classdef gia2plmtDegreeTest < matlab.unittest.TestCase

    methods (Test)

        function resizesModelAndBoundsIndependently(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            path = fullfile(folder.Folder, 'model.mat');
            lmcosiM = single([0 0 0 0; 1 0 2 0; 1 1 3 4]);
            degreeTwo = single([2 0 0 0; 2 1 0 0; 2 2 0 0]);
            lmcosiU = [lmcosiM; degreeTwo; single([3 0 5 0; 3 1 6 7; 3 2 8 9; 3 3 10 11])];
            lmcosiL = lmcosiM;
            lmcosiU(4:6, 3) = [5; 6; 7];
            lmcosiL(2:3, 3) = [-2; -3];
            save(path, 'lmcosiM', 'lmcosiU', 'lmcosiL');

            [model, upperBound, lowerBound] = testCase.verifyWarning(@() gia2plmt(path, 2), ...
                'ULMO:LoadData:InsufficientResolution');
            testCase.verifyEqual(model, [lmcosiM; degreeTwo]);
            testCase.verifyEqual(upperBound, lmcosiU(1:6, :));
            testCase.verifyEqual(lowerBound, [lmcosiL; degreeTwo]);
        end

        function truncatesModelWithoutBounds(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            path = fullfile(folder.Folder, 'model.mat');
            lmcosiM = [0 0 0 0; 1 0 2 0; 1 1 3 4; 2 0 5 0; 2 1 6 7; 2 2 8 9];
            save(path, 'lmcosiM');
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:load:variableNotFound'));

            [model, upperBound, lowerBound] = gia2plmt(path, 1, 'BeQuiet', true);
            testCase.verifyEqual(model, lmcosiM(1:3, :));
            testCase.verifyEmpty(upperBound);
            testCase.verifyEmpty(lowerBound);
        end

    end

end
