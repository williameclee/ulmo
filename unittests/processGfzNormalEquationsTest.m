%% PROCESSGFZNORMALEQUATIONSTEST - Checks the multi-epoch GFZ covariance factors.
% Synthetic SINEX solutions use a known covariance with cross-degree coupling.
%
% Created by
%   2026/11/10
classdef processGfzNormalEquationsTest < matlab.unittest.TestCase
    properties
        Source
        Root
        Early = 'GSM-2_2018002-2018031_GRFO_GFZOP_BB01_0603.snx'
        Late = 'GSM-2_2018001-2018090_GRAC_GFZOP_BB01_0600.snx'
    end
    methods (TestMethodSetup)
        function isolateFiles(testCase)
            import matlab.unittest.fixtures.*
            temp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.Root = fullfile(temp.Folder, 'local');
            testCase.Source = fullfile(temp.Folder, 'source');
            mkdir(testCase.Source);
            testCase.applyFixture(EnvironmentVariableFixture('GRACEDATA', testCase.Root));
        end
    end
    methods (Test)
        function reducedAndFullDegreeSeries(testCase)
            early = fullfile(testCase.Source, testCase.Early);
            late = fullfile(testCase.Source, testCase.Late);
            [earlyCov, labels] = writeSolution(early, 15, 1);
            lateCov = writeSolution(late, 75, 2);
            gzip(late);
            delete(late);
            writeText([early '.gz.part'], 'partial download');
            before = posixtime(datetime('now','TimeZone','UTC'));
            paths = cell(2, 1);
            paths{1} = processGfzNormalEquations(char(testCase.Source), 2);
            paths{2} = processGfzNormalEquations(string(testCase.Source), 3);
            after = posixtime(datetime('now','TimeZone','UTC'));
            rmdir(testCase.Source, 's');

            for L = 2:3
                path = paths{L-1};
                keep = labels(:,1) <= L;
                n = nnz(keep);
                testCase.verifyEqual(path, fullfile(testCase.Root, 'GFZ', 'NEQ', ...
                    sprintf('GFZ_NEQ_L%d.nc', L)));
                testCase.verifyEqual(ncread(path,'degree'), uint16(labels(keep,1)));
                testCase.verifyEqual(ncread(path,'order'), uint16(labels(keep,2)));
                testCase.verifyEqual(ncread(path,'cs'), uint8(labels(keep,3)));
                testCase.verifyEqual(ncread(path,'sourceFilename'), [string(testCase.Early); string(testCase.Late)+".gz"]);
                testCase.verifyEqual(ncread(path,'varianceFactor'), [1;2], 'AbsTol', 1e-12);
                testCase.verifyEqual(ncread(path,'sourceMaxDegree'), uint16([3;3]));
                testCase.verifyEqual(ncread(path,'referenceRadius'), [6378136;6378137]);
                testCase.verifyEqual(ncreadatt(path,'/','maxDegree'), L);
                testCase.verifyEqual(ncreadatt(path,'/','normalisation'), 'fully_normalized');
                testCase.verifyEqual(ncreadatt(path,'/','tideSystem'), 'tide_free');
                start = seconds(datetime(2018,1,1)-datetime(2000,1,1));
                testCase.verifyEqual(ncread(path,'time'), start + [14;74]*86400);
                testCase.verifyEqual(ncread(path,'epochStart'), start + [1;0]*86400);
                testCase.verifyEqual(ncread(path,'epochEnd'), start + [31;90]*86400);
                created = ncinfo(path, 'creationDate');
                testCase.verifyEmpty(created.Dimensions);
                testCase.verifyGreaterThanOrEqual(ncread(path, 'creationDate'), before);
                testCase.verifyLessThanOrEqual(ncread(path, 'creationDate'), after);
                testCase.verifyEqual(ncreadatt(path, 'creationDate','units'), 'seconds since 1970-01-01 00:00:00 UTC');
                info = ncinfo(path, 'noiseFactor');
                testCase.verifyEqual({info.Dimensions.Name}, {'packed','time'});
                testCase.verifyEqual(info.Size, [n*(n+1)/2, 2]);
                testCase.verifyEqual(info.ChunkSize(2), 1);
                covs = {earlyCov, lateCov};
                for t = 1:2
                    expected = covs{t}(keep, keep);
                    lower = zeros(n);
                    lower(tril(true(n))) = ncread(path, 'noiseFactor', [1,t], [n*(n+1)/2,1]);
                    testCase.verifyEqual(lower * lower.', expected, 'AbsTol', 1e-11);
                    testCase.verifyEqual(ncread(path, 'formalStd', [1,t], [n,1]), ...
                        sqrt(diag(expected)), 'AbsTol', 1e-12);
                end
            end
        end

        function rejectsNonscalarFolder(testCase)
            folders = repmat(string(testCase.Source), 1, 2);
            testCase.verifyError(@() processGfzNormalEquations(folders), ...
                'MATLAB:validators:mustBeTextScalar');
        end

        function rejectsDuplicateEpoch(testCase)
            writeSolution(fullfile(testCase.Source, testCase.Early), 15, 1);
            writeSolution(fullfile(testCase.Source, testCase.Late), 15, 1);
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source, 2), ...
                'ULMO:gfzneq:DuplicateEpoch');
        end

        function rejectsEmptyFolderAndInvalidDegrees(testCase)
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source), 'ULMO:gfzneq:NoFiles');
            writeSolution(fullfile(testCase.Source,testCase.Early), 15, 1);
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source), 'ULMO:gfzneq:Degree');
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source, 1), ...
                'MATLAB:validators:mustBeGreaterThanOrEqual');
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source, 2.5), ...
                'MATLAB:validators:mustBeInteger');
        end

        function detectsInconsistentFormalScale(testCase)
            path = fullfile(testCase.Source, testCase.Early);
            writeSolution(path, 15, 1);
            % The final printed estimate is native index 11, a degree-2 term.
            text = regexprep(fileread(path), '[^ \n]+\n-SOLUTION/ESTIMATE', '1.0\n-SOLUTION/ESTIMATE');
            writeText(path, text);
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source, 2), ...
                'ULMO:gfzneq:Scale');
        end

        function preservesPriorSeriesOnTruncatedSecondMonth(testCase)
            writeSolution(fullfile(testCase.Source,testCase.Early), 15, 1);
            out = processGfzNormalEquations(testCase.Source, 2);
            original = ncread(out,'noiseFactor');
            created = ncread(out,'creationDate');
            late = fullfile(testCase.Source,testCase.Late);
            writeSolution(late, 75, 2);
            body = extractBefore(string(fileread(late)), '-SOLUTION/NORMAL_EQUATION_MATRIX L');
            writeText(late, body);
            testCase.verifyError(@() processGfzNormalEquations(testCase.Source, 2), ...
                'ULMO:gfzneq:InvalidSource');
            testCase.verifyEqual(ncread(out,'noiseFactor'), original);
            testCase.verifyEqual(ncread(out,'creationDate'), created);
            files = dir(fullfile(fileparts(out), '*.nc'));
            testCase.verifyEqual({files.name}, {'GFZ_NEQ_L2.nc'});
        end
    end
end

function [covariance, labels] = writeSolution(path, day, scale)
    labels = [2 0 0; 2 1 0; 2 1 1; 2 2 0; 2 2 1; ...
        3 0 0; 3 1 0; 3 1 1; 3 2 0; 3 2 1; 3 3 0; 3 3 1];
    n = size(labels,1);
    weights = (1:n)'/n;
    covariance = diag(1 + scale*weights) + scale*(weights*weights.');
    % Reverse the native order and print estimates in a third order.
    native = n:-1:1;
    normal = scale * (covariance(native,native) \ eye(n));
    sigma = sqrt(diag(covariance(native,native)));
    file = fopen(path, 'w');
    closeFile = onCleanup(@() fclose(file));
    if day <= 31
        interval = '18:002:00000 18:032:00000';
    else
        interval = '18:001:00000 18:091:00000';
    end
    fprintf(file, '%%=SNX 2.02 GFZ 18:100:00000 GFZ %s C %d 2\n', interval, n);
    fprintf(file, '+FILE/COMMENT\n max_degree 3\n radius %d\n', 6378135+scale);
    fprintf(file, ' norm fully_normalized\n tide_system tide_free\n-FILE/COMMENT\n');
    fprintf(file, '+SOLUTION/ESTIMATE\n*INDEX TYPE CODE PT SOLN EPOCH UNIT C EST STD\n');
    kinds = ["CN","SN"];
    for idx = [2:2:n, 1:2:n]
        label = labels(native(idx), :);
        fprintf(file, ' %d %s %d -- %d 18:%03d:00000 ---- 2 0 %.16e\n', ...
            idx, kinds(label(3)+1), label(1), label(2), day, sigma(idx));
    end
    fprintf(file, '-SOLUTION/ESTIMATE\n+SOLUTION/NORMAL_EQUATION_MATRIX L\n*PARA1 PARA2 VALUES\n');
    for row = 1:n
        for col = 1:3:row
            fprintf(file, ' %d %d', row, col);
            fprintf(file, ' %.16e', normal(row, col:min(col+2,row)));
            fprintf(file, '\n');
        end
    end
    fprintf(file, '-SOLUTION/NORMAL_EQUATION_MATRIX L\n%%ENDSNX\n');
end

function writeText(path, content)
    file = fopen(path, 'w');
    closeFile = onCleanup(@() fclose(file));
    fprintf(file, '%s', content);
end
