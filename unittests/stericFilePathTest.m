%% STERICFILEPATHTEST - Checks native and derived steric cache filenames.
%
% Syntax
%   results = runtests('unittests/stericFilePathTest.m')
%
% Notes
%   Checks exact native, temporal and spatial
%   filenames. A temporary IFILES root is used, without creating product
%   directories or reading native data.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/20, En-Chi Lee (williameclee@gmail.com)

classdef stericFilePathTest < matlab.unittest.TestCase

    methods (Test)

        function nativeAndInterpolationPaths(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            prefix = fullfile(folder.Folder, 'STERIC', 'EN4c14');
            testCase.verifyEqual(stericfilepath('EN4c14'), fullfile(prefix, 'EN4c14-StericSeaLevel.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', days(10), [], [], 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-T10days-interp-linear.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 10, [], [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-T10-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 'GRACE', 1, 180, 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-Tgrace-H1O180-linear.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 'midmonth', [], [], 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-Tmidmonth-linear.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', [], 0.5, [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-H0.5-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', [], [], 180, 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-O180-linear.mat'));
            testCase.verifyFalse(isfolder(prefix));
        end

    end

end
