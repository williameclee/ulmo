%% SHNCTEST - Tests compact storage and observes actual NetCDF subset reads.
% A transparent NCREAD wrapper records requested slabs while delegating all
% I/O to MATLAB's original function. Tests use isolated temporary files.
%
% Created by
%   2026/10/30

classdef shncTest < matlab.unittest.TestCase
    properties
        Scratch
        Folder
        Labels = [0 0; 1 0; 1 1; 2 0; 2 1; 2 2]
    end
    methods (TestMethodSetup)
        function setup(testCase)
            testCase.Scratch = tempname;
            testCase.Folder = fullfile(testCase.Scratch, 'ESM', 'ESAESM2');
            mkdir(testCase.Folder);
            testCase.addTeardown(@() rmdir(testCase.Scratch, 's'));
            oldRoot = getenv('IFILES');
            setenv('IFILES', testCase.Scratch);
            testCase.addTeardown(@() setenv('IFILES', oldRoot));
        end
    end
    methods (Test)
        function nativeReadsTimeFirstAndOnlyRequestedSlab(testCase)
            dates = datetime(2003, 1, 1) + hours((0:9)' * 6);
            path = fullfile(testCase.Folder, 'ESAESM2_O.nc');
            c = reshape(1:120, 6, 2, 10);
            createshnc(path, testCase.Labels, dates, 'GRAV', 'EquatorRadius', 6378137);
            ncwrite(path, 'coefficients', c);
            testCase.recordReads();
            [p, s, d] = esm2plmt('O', dates([3, 4]), '', 'Unit', 'GRAV', ...
                'Loutput', 1, 'TimeFormat', 'datetime');
            testCase.verifyEqual(p(:, :, 3:4), permute(c(1:3, :, 3:4), [3, 1, 2]));
            testCase.verifyEqual(d, dates(3:4));
            testCase.verifyTrue(all(isnan(s(:, :, 3:4)), 'all'));
            testCase.verifySlab('coefficients', [1, 1, 3], [3, 2, 2]);
        end
        function cacheHitReadsOnlySelectedEpochAndDegree(testCase)
            dates = [datetime(2003, 1, 16, 12, 0, 0); datetime(2003, 2, 15); datetime(2003, 3, 16, 12, 0, 0)];
            c = reshape(1:36, 6, 2, 3);
            path = fullfile(testCase.Folder, 'ESAESM2_O_GRAV-Tmidmonth.nc');
            writeshnc(path, c, c / 10, testCase.Labels, dates, 'GRAV');
            testCase.recordReads();
            [p, s] = esm2plmt('O', [datetime(2003, 2, 1), datetime(2003, 2, 28)], ...
                'midmonth', 'Unit', 'GRAV', 'Loutput', 1);
            testCase.verifyEqual(p(:, :, 3:4), permute(c(1:3, :, 2), [3, 1, 2]));
            testCase.verifyEqual(s(:, :, 3:4), permute(c(1:3, :, 2) / 10, [3, 1, 2]));
            testCase.verifySlab('coefficients', [1, 1, 2], [3, 2, 1]);
            testCase.verifySlab('uncertainty', [1, 1, 2], [3, 2, 1]);
        end
        function emptySelectionDoesNotReadCoefficientVariables(testCase)
            path = fullfile(testCase.Folder, 'ESAESM2_O.nc');
            date = datetime(2003, 1, 1);
            writeshnc(path, ones(6, 2), ones(6, 2), testCase.Labels, date, 'GRAV');
            testCase.recordReads();
            [p, s, d] = esm2plmt('O', date + days([1, 2]), '', 'Unit', 'GRAV');
            testCase.verifySize(p, [0, 6, 4]);
            testCase.verifySize(s, [0, 6, 4]);
            testCase.verifyEmpty(d);
            global ULMO_NC_READ_LOG
            names = cellfun(@(args) string(args{2}), ULMO_NC_READ_LOG);
            testCase.verifyEqual(names(1), "time");
            testCase.verifyFalse(any(ismember(names, ["coefficients", "uncertainty"])));
        end
        function singletonAndMultipleEpochsRoundTripExactly(testCase)
            for n = [1, 3]
                path = fullfile(testCase.Folder, "roundtrip" + n + ".nc");
                c = reshape((1:12*n) * pi, 6, 2, n);
                c(3, 2, 1) = NaN;
                dates = datetime(1999, 12, 31, 12, 0, 0) + seconds((0:n-1)' * 43201.5);
                writeshnc(path, c, c .^ 2, testCase.Labels, dates, 'POT');
                [actual, sigma, actualDates, labels] = readshnc(path);
                testCase.verifyEqual(actual, c);
                testCase.verifyEqual(sigma, c .^ 2);
                testCase.verifyEqual(actualDates, dates);
                testCase.verifyEqual(labels, testCase.Labels);
                info = ncinfo(path, 'coefficients');
                testCase.verifyEqual({info.Dimensions.Name}, {'lm', 'cs', 'time'});
                testCase.verifyEqual(info.Size, [6, 2, n]);
                testCase.verifyEqual(info.Datatype, 'double');
                testCase.verifyClass(ncread(path, 'degree'), 'uint16');
                testCase.verifyEqual(ncread(path, 'cs'), uint8([0; 1]));
            end
        end
        function uncertaintyAxesMustMatchEvenWhenSizesAgree(testCase)
            path = fullfile(testCase.Folder, 'wrong-axes.nc');
            dates = datetime(2003, 1, 1) + days((0:5)');
            createshnc(path, testCase.Labels, dates, 'GRAV');
            nccreate(path, 'uncertainty', 'Dimensions', {'time', 6, 'cs', 2, 'lm', 6});
            testCase.verifyError(@() readshnc(path), 'ULMO:readshnc:InvalidDimensions');
        end
        function exclusiveUpperBoundOmitsNextEpochSample(testCase)
            dates = datetime(2003, 1, 1) + hours([0; 6; 12]);
            path = fullfile(testCase.Folder, 'bounds.nc');
            writeshnc(path, ones(6, 2, 3), [], testCase.Labels, dates, 'GRAV');
            [c, ~, d] = readshnc(path, dates([1, 3]), [], 'EndExclusive', true);
            testCase.verifySize(c, [6, 2, 2]);
            testCase.verifyEqual(d, dates(1:2));
        end
    end
    methods (Access = private)
        function recordReads(testCase)
            global ULMO_NC_READ_LOG ULMO_ORIGINAL_NCREAD
            ULMO_NC_READ_LOG = {};
            ULMO_ORIGINAL_NCREAD = @ncread;
            shadow = fullfile(testCase.Scratch, 'spy');
            mkdir(shadow);
            fid = fopen(fullfile(shadow, 'ncread.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function result = ncread(varargin)', ...
                'global ULMO_NC_READ_LOG ULMO_ORIGINAL_NCREAD', ...
                'ULMO_NC_READ_LOG{end + 1} = varargin;', ...
                'result = ULMO_ORIGINAL_NCREAD(varargin{:});', ...
                'end');
            fclose(fid);
            oldPath = path;
            addpath(shadow, '-begin');
            testCase.addTeardown(@() path(oldPath));
            testCase.addTeardown(@() clearReadLog());
        end
        function verifySlab(testCase, variable, start, count)
            global ULMO_NC_READ_LOG
            names = cellfun(@(args) string(args{2}), ULMO_NC_READ_LOG);
            testCase.verifyEqual(names(1), "time");
            calls = ULMO_NC_READ_LOG(names == string(variable));
            testCase.assertNumElements(calls, 1);
            testCase.verifyEqual(calls{1}{3}, start);
            testCase.verifyEqual(calls{1}{4}, count);
        end
    end
end

function clearReadLog()
    clear global ULMO_NC_READ_LOG ULMO_ORIGINAL_NCREAD
end
