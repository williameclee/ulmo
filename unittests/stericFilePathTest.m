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
%   2026/11/22, En-Chi Lee (williameclee@gmail.com)

classdef stericFilePathTest < matlab.unittest.TestCase

    methods (Test)

        function nativeAndInterpolationPaths(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
            prefix = fullfile(folder.Folder, 'STERIC', 'EN4c14');
            testCase.verifyEqual(stericfilepath('EN4c14'), fullfile(prefix, 'EN4c14-StericSeaLevel.mat'));
            fixedPath = fullfile(prefix, 'EN4c14-StericSeaLevel-T864000seconds-linear.mat');
            testCase.verifyEqual(stericfilepath('EN4c14', days(10)), fixedPath);
            testCase.verifyEqual(stericfilepath('EN4c14', 10, [], [], 'LINEAR'), fixedPath);
            testCase.verifyEqual(stericfilepath('EN4c14', hours(240), [], [], 'linear'), fixedPath);
            testCase.verifyEqual(stericfilepath('EN4c14', 10, [], [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-T864000seconds-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 10, 1, 180, 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-T864000seconds-nearest-H1O180.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 'GRACE', 1, 180, 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-grace-H1O180-linear.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 'midmonth', [], [], 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-midmonth.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', 'MIDMONTH', [], [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-midmonth.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', [], 0.5, [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-H0.5-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', duration.empty, 0.5, [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-H0.5-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', '', 0.5, [], 'nearest'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-H0.5-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14', [], [], 180, 'linear'), ...
                fullfile(prefix, 'EN4c14-StericSeaLevel-O180-linear.mat'));
            testCase.verifyFalse(isfolder(prefix));
        end

    end

end
