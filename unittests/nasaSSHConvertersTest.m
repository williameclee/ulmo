%% NASASSHCONVERTERSTEST - Checks native and monthly NASA-SSH conversion.
% Each test uses temporary NetCDF inputs and a private IFILES directory.
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

classdef nasaSSHConvertersTest < matlab.unittest.TestCase

    properties (TestParameter)
        filenameType = {'char', 'string'}
    end

    methods (TestMethodSetup)

        function createFolders(testCase)
            % Exercise filename-only outputs from the temporary working folder.
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', pwd));
            mkdir('input');
        end

    end

    methods (Test)

        function singletonRoundTrip(testCase, filenameType)
            % One map is a complete input for both converters and a partial month.
            date = datetime(2000, 2, 7);
            values = [1.125, -9999; 3, 11];
            writeMap(date, values);
            weekly = feval(filenameType, 'weekly.mat');
            monthly = feval(filenameType, 'monthly.mat');
            testCase.verifyEqual(processSSHDataNasa('input', weekly), weekly);
            testCase.verifyEqual(processSSHMonthlyNasa(weekly, monthly), monthly);
            w = load(weekly);
            m = load(monthly);
            expected = single([1.125, NaN; 3, NaN]);
            testCase.verifyEqual(w.sshs, expected);
            testCase.verifyEqual(w.dates, date);
            testCase.verifyEqual(m.sshs, expected);
            testCase.verifyEqual(m.validMapCounts, uint8(isfinite(expected)));
            testCase.verifyEqual(m.mapCounts, uint8(1));
            testCase.verifyEqual(m.expectedMapCounts, uint8(4));
            testCase.verifyTrue(m.partialMonths);
            testCase.verifyTrue(all(m.incompleteCoverage, 'all'));
            testCase.verifyEqual(m.dates, datetime(2000, 2, 15, 12, 0, 0));
        end

        function multipleMapsInOneMonth(testCase)
            % A 3-D weekly file may produce just one 2-D monthly map.
            writeMap(datetime(2000, 2, 7), [1 2; 3 4]);
            writeMap(datetime(2000, 2, 14), [5 -9999; 7 8]);
            processSSHDataNasa('input', BeQuiet = true);
            m = load(processSSHMonthlyNasa(BeQuiet = true));
            testCase.verifyEqual(m.sshs, single([3 2; 5 6]));
            testCase.verifyEqual(m.validMapCounts, uint8([2 1; 2 2]));
            testCase.verifyEqual(m.metadata.sourceIndices, {[1; 2]});
        end

        function completeMonthCoverageAndPrecision(testCase)
            % Leap February has four Mondays. One cell is missing its last sample.
            dates = datetime(2000, 2, [7, 14, 21, 28]);
            values = [8, 2 ^ -22, -8, 0];

            for iDate = 1:numel(dates)
                field = [values(iDate), 0; 1, -9999];

                if iDate == numel(dates)
                    field(2, 1) = -9999;
                end

                writeMap(dates(iDate), field);
            end

            weekly = processSSHDataNasa('input', BeQuiet = true);
            m = load(processSSHMonthlyNasa(weekly, BeQuiet = true));
            testCase.verifyFalse(m.partialMonths);
            testCase.verifyEqual(m.sshs, single([2 ^ -24, 0; 1, NaN]));
            testCase.verifyEqual(m.validMapCounts, uint8([4 4; 3 0]));
            testCase.verifyEqual(m.incompleteCoverage, [false false; true true]);
        end

        function defaultPathsAndReader(testCase)
            % Exercise default logging, partial edges, empty maps and a missing month.
            dates = [datetime(2000, 1, 31); datetime(2000, 2, 7); ...
                         datetime(2000, 2, 14); datetime(2000, 4, 3)];
            values = cat(3, [2 3; 4 5], [1 2; 3 4], -9999 * ones(2), [5 6; 7 8]);

            for iDate = 1:numel(dates)
                writeMap(dates(iDate), values(:, :, iDate));
            end

            weekly = processSSHDataNasa('input');
            monthly = processSSHMonthlyNasa();
            w = load(weekly);
            m = load(monthly);
            testCase.verifyEqual(weekly, fullfile(pwd, 'SSH', 'NASASSH.mat'));
            testCase.verifyEqual(monthly, fullfile(pwd, 'SSH', 'NASASSHMonthly.mat'));
            testCase.verifyEqual(w.dates, dates);
            testCase.verifyEqual(w.lon, [0; 180]);
            testCase.verifyEqual(w.lat, [-90; 90]);
            testCase.verifyEmpty(w.sshErrors);
            testCase.verifyEqual(w.creationDate.TimeZone, 'UTC');
            testCase.verifyEqual(m.creationDate.TimeZone, 'UTC');
            testCase.verifyEqual(w.metadata.coverageStart, string(dates - days(5), 'yyyy-MM-dd'));
            testCase.verifyEqual(w.metadata.coverageEnd, string(dates + days(5), 'yyyy-MM-dd'));
            testCase.verifyFalse(any(isfield(w.metadata, {'inputFolder', 'sourceFiles', ...
                                                              'sourceBytes', 'sourceModifiedDatenum', 'finiteCounts', 'created'})));
            testCase.verifyFalse(any(isfield(m.metadata, {'sourcePath', 'created'})));
            testCase.verifyEqual(m.metadata.sourceMetadata, w.metadata);
            testCase.verifyEqual(m.metadata.sourceDates, dates);
            testCase.verifyEqual(m.metadata.sourceIndices, {1; [2; 3]; zeros(0, 1); 4});
            testCase.verifyEqual(m.mapCounts, uint8([1; 2; 0; 1]));
            testCase.verifyEqual(m.expectedMapCounts, uint8([5; 4; 4; 4]));
            testCase.verifyEqual(m.partialMonths, [true; false; false; true]);
            testCase.verifyEqual(m.validMapCounts, repmat(reshape(uint8([1 1 0 1]), 1, 1, []), 2, 2));
            testCase.verifyTrue(all(m.incompleteCoverage, 'all'));
            values(:, :, 3) = NaN;
            [ssh, sigma, t] = ssh2lonlatt('NASASSH', OutputFormat = 'ndgrid', ...
                Unit = 'mm', SaveData = false, BeQuiet = true);
            testCase.verifyEmpty(sigma);
            testCase.verifyEqual(t, dates);
            testCase.verifyEqual(ssh, single(1000 * values));
            [ssh, sigma, t] = ssh2lonlatt('NASASSHMonthly', OutputFormat = 'ndgrid', ...
                SaveData = false, BeQuiet = true);
            testCase.verifyEmpty(sigma);
            testCase.verifyEqual(t, m.monthStarts + (m.monthEnds - m.monthStarts) / 2);
            testCase.verifyEqual(ssh, single(values));
        end

        function preserveExistingDestinations(testCase)
            % Reject existing files and directories through both public converters.
            writeMap(datetime(2000, 1, 3), ones(2));
            weekly = processSSHDataNasa('input');
            mkdir('existing');
            writelines("original", 'existing.mat');

            for target = ["existing", "existing.mat"]
                testCase.verifyError(@() processSSHDataNasa('input', target), 'ULMO:OutputExists');
                testCase.verifyError(@() processSSHMonthlyNasa(weekly, target), 'ULMO:OutputExists');
            end

            testCase.verifyEqual(readlines('existing.mat', EmptyLineRule = 'skip'), "original");
            testCase.verifyEmpty(dir(fullfile('existing', '*.mat')));
        end

        function destinationAppearsBeforePublication(testCase)
            % Stage both candidates before publishing either to exercise the race.
            writelines("first", 'first.tmp');
            writelines("second", 'second.tmp');
            publishnewfile('first.tmp', 'product.mat');
            deletefileifpresent('first.tmp');
            testCase.verifyError(@() publishnewfile('second.tmp', 'product.mat'), 'ULMO:OutputExists');
            testCase.verifyEqual(readlines('product.mat', EmptyLineRule = 'skip'), "first");
        end

        function discardPartialOutput(testCase)
            % Abort on the second map after the first has been written to disk.
            writeMap(datetime(2000, 1, 3), ones(2));
            bad = writeMap(datetime(2000, 1, 10), ones(2));
            ncwriteatt(bad, 'ssha', 'units', 'km');
            target = fullfile('products', 'weekly.mat');
            testCase.verifyError(@() processSSHDataNasa('input', target), 'ULMO:InvalidSSHUnits');
            testCase.verifyEmpty(dir(fullfile('products', '*.mat')));
            ncwriteatt(bad, 'ssha', 'units', 'm');
            testCase.verifyEqual(processSSHDataNasa('input', target), target);
            testCase.verifyEqual({dir(fullfile('products', '*.mat')).name}, {'weekly.mat'});
        end

    end

end

function path = writeMap(date, values)
    % Write a minimal reference-mission granule with the native NetCDF attributes.
    path = fullfile('input', ['NASA-SSH_alt_ref_simple_grid_v1_', char(date, 'yyyyMMdd'), '.nc']);
    nccreate(path, 'longitude', 'Dimensions', {'longitude', 2});
    ncwrite(path, 'longitude', [0; 180]);
    nccreate(path, 'latitude', 'Dimensions', {'latitude', 2});
    ncwrite(path, 'latitude', [-90; 90]);
    nccreate(path, 'time');
    ncwrite(path, 'time', seconds(date - datetime(1990, 1, 1)));
    ncwriteatt(path, 'time', 'units', 'seconds since 1990-01-01');
    nccreate(path, 'ssha', 'Dimensions', {'longitude', 2, 'latitude', 2}, 'FillValue', -9999);
    ncwrite(path, 'ssha', values);
    ncwriteatt(path, 'ssha', 'units', 'm');
    ncwriteatt(path, 'ssha', 'valid_min', -10);
    ncwriteatt(path, 'ssha', 'valid_max', 10);
    ncwriteatt(path, '/', 'product_short_name', 'NASA_SSH_REF_SIMPLE_GRID_V1');
    ncwriteatt(path, '/', 'mean_sea_surface', 'DTU21');
    ncwriteatt(path, '/', 'gridding_method', 'synthetic fixture');
    ncwriteatt(path, '/', 'time_coverage_start', char(date - days(5), 'yyyy-MM-dd'));
    ncwriteatt(path, '/', 'time_coverage_end', char(date + days(5), 'yyyy-MM-dd'));
end
