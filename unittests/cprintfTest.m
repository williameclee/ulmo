classdef cprintfTest < matlab.unittest.TestCase
    % Unit tests for the cprintf function

    methods (Test)
        function testHotLinksEnabled(testCase)
            % Test when HotLinks are enabled
            txt = 'txt <a href="matlab: open(''/path/to/file.m'')">file</a>'; %#ok<NASGU>
            expectedOutput = 'txt <a href="matlab: open(''/path/to/file.m'')">file</a>';
            testCase.verifyEqual(evalc('feature(''HotLinks'', true); cprintf(''%s'', txt)'), expectedOutput);
        end

        function testHotLinksDisabledScriptFile(testCase)
            % Test when HotLinks are disabled for .m files
            txt = 'txt <a href="matlab: open(''/path/to/file.m'')">file</a>'; %#ok<NASGU>
            expectedOutput = 'txt file';
            testCase.verifyEqual(evalc('feature(''HotLinks'', false); cprintf(''%s'', txt)'), expectedOutput);
        end

        function testHotLinksDisabledDataFile(testCase)
            % Test when HotLinks are disabled for data files
            txt = 'txt <a href="matlab: fprintf(''data/file.mat\n'');open(''data/file.mat'')">data</a>'; %#ok<NASGU>
            expectedOutput = 'txt data (data/file.mat)';
            testCase.verifyEqual(evalc('feature(''HotLinks'', false); cprintf(''%s'', txt)'), expectedOutput);
        end

        function testHotLinksDisabledExternalLink(testCase)
            % Test when HotLinks are disabled for external links
            txt = 'txt <a href="https://example.com">example</a>'; %#ok<NASGU>
            expectedOutput = 'txt example (https://example.com)';
            testCase.verifyEqual(evalc('feature(''HotLinks'', false); cprintf(''%s'', txt)'), expectedOutput);
        end

        function testHotLinksDisabledRemainingTags(testCase)
            % Test when HotLinks are disabled for remaining tags
            txt = 'txt <a href="#">link</a>'; %#ok<NASGU>
            expectedOutput = 'txt link';
            testCase.verifyEqual(evalc('feature(''HotLinks'', false); cprintf(''%s'', txt)'), expectedOutput);
        end
    end
end