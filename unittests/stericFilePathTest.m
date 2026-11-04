%% STERICFILEPATHTEST - Preserve existing native and derived cache filenames.
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

classdef stericFilePathTest < matlab.unittest.TestCase
    methods (Test)
        function nativeAndInterpolationPaths(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES',folder.Folder));
            prefix = fullfile(folder.Folder,'STERIC','EN4c14');
            testCase.verifyEqual(stericfilepath('EN4c14'), fullfile(prefix,'EN4c14-StericSeaLevel.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14',days(10),[],[],'linear'), ...
                fullfile(prefix,'EN4c14-StericSeaLevel-T10days-interp-linear.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14',10,[],[],'nearest'), ...
                fullfile(prefix,'EN4c14-StericSeaLevel-T10-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14','GRACE',1,180,'linear',false,[],'-deep-total'), ...
                fullfile(prefix,'EN4c14-StericSeaLevel-Tgrace-H1O180-linear-deep-total.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14','midmonth',[],[],'linear'), ...
                fullfile(prefix,'EN4c14-StericSeaLevel-Tmidmonth-linear.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14',[],0.5,[],'nearest'), ...
                fullfile(prefix,'EN4c14-StericSeaLevel-H0.5-nearest.mat'));
            testCase.verifyEqual(stericfilepath('EN4c14',[],[],180,'linear'), ...
                fullfile(prefix,'EN4c14-StericSeaLevel-O180-linear.mat'));
            testCase.verifyFalse(isfolder(prefix));
        end

        function truncatedPaths(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES',folder.Folder));
            prefix = fullfile(folder.Folder,'STERIC','SIO');
            limits = datetime(2020,1,[1,31]);
            expected = fullfile(prefix,'SIO-StericSeaLevel-Tmidmonth-linear-20200101_20200131-full-total.mat');
            testCase.verifyEqual(stericfilepath('SIO','midmonth',[],[],'linear',true,limits,'-full-total'),expected);
            testCase.verifyEqual(stericfilepath('SIO','midmonth',[],[],'linear',true,datenum(limits),'-full-total'),expected); %#ok<DATNM>
            limits = limits + hours(12);
            expected = fullfile(prefix,'SIO-StericSeaLevel-Tmidmonth-linear-20200101T120000000000000_20200131T120000000000000.mat');
            testCase.verifyEqual(stericfilepath('SIO','midmonth',[],[],'linear',true,limits),expected);
        end
    end
end
