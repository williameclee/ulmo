%% PROCESSSSHMONTHLYNASA - Averages weekly NASA-SSH maps by calendar month.
% Give each finite value equal weight and assign maps by their central epoch.
% Accumulate in double precision. Retain empty months and partial edge months.
%
% Syntax
%   outPath = processSSHMonthlyNasa()
%   outPath = processSSHMonthlyNasa(inPath, outPath)
%   outPath = processSSHMonthlyNasa(__, "BeQuiet", true)
%
% Input arguments
%   inPath (optional) - Weekly MAT file produced by processSSHDataNasa.
%       Default path: fullfile(getenv('IFILES'), 'SSH', 'NASASSH.mat')
%   outPath (optional) - Destination MAT file.
%       Default path: fullfile(getenv('IFILES'), 'SSH', 'NASASSHMonthly.mat')
%   BeQuiet (name-value) - Whether to suppress progress and completion messages.
%       Default value: false
%
% Input file content
%   dates - Increasing central epochs on a common seven-day sampling grid.
%       Size: (nWeek, 1).
%   lon, lat - Native longitude and latitude vectors in degrees.
%       Size: (nLon, 1) and (nLat, 1), respectively.
%   sshs - Weekly anomalies in metres, with NaN for missing values.
%       Size: (nLon, nLat, nWeek).
%   metadata - Source metadata with product set to NASA_SSH_REF_SIMPLE_GRID_V1.
%
% Output arguments
%   outPath - Path to the monthly MAT file.
%
% Output file content
%   sshs - Single-precision monthly means, with NaN for cells with no samples.
%       Size: (nLon, nLat, nMonth).
%       Unit: metres.
%   lon, lat - Longitude and latitude vectors copied from the weekly file.
%   dates - Midpoints between consecutive calendar-month boundaries.
%       Size: (nMonth, 1).
%   monthStarts, monthEnds - Inclusive start and exclusive end of each month.
%       Size: (nMonth, 1).
%   mapCounts, expectedMapCounts - Numbers of available and scheduled weekly
%       epochs per month. Available epochs include entirely missing maps.
%       Size: (nMonth, 1).
%   validMapCounts - Number of finite weekly values per cell and month.
%       Size: (nLon, nLat, nMonth).
%   incompleteCoverage - True where validMapCounts is below expectedMapCounts.
%       Size: (nLon, nLat, nMonth).
%   partialMonths - True where scheduled epochs fall outside the input span.
%       Size: (nMonth, 1).
%   sshErrors - Empty array. Overlapping source windows share observations.
%   metadata - Averaging method, coverage policy, weekly dates,
%       indices of contributing epochs and weekly product metadata.
%   creationDate - UTC datetime when the converter wrote the file.
%
% See also
%   PROCESSSSHDATANASA, SSH2LONLATT
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function outPath = processSSHMonthlyNasa(inPath, outPath, options)

    arguments (Input)
        inPath {mustBeTextScalar} = fullfile(getenv('IFILES'), 'SSH', 'NASASSH.mat')
        outPath {mustBeTextScalar} = fullfile(getenv('IFILES'), 'SSH', 'NASASSHMonthly.mat')
        options.BeQuiet (1, 1) logical = false
    end

    assert(~isfile(outPath) && ~isfolder(outPath), 'ULMO:OutputExists', 'Output file already exists: %s', outPath);
    d = load(inPath, 'dates', 'lon', 'lat', 'metadata');
    assert(strcmp(d.metadata.product, 'NASA_SSH_REF_SIMPLE_GRID_V1'));
    assert(~isempty(d.dates) && all(diff(d.dates) > days(0)));
    assert(all(mod(days(d.dates - d.dates(1)), 7) == 0), 'Input is not on a weekly grid.');
    src = matfile(inPath);
    gridSize = size(src, 'sshs');
    gridSize(end + 1:3) = 1;
    assert(isequal(gridSize, [numel(d.lon), numel(d.lat), numel(d.dates)]));
    monthStarts = (dateshift(d.dates(1), 'start', 'month'):calmonths(1): ...
        dateshift(d.dates(end), 'start', 'month'))';
    monthEnds = monthStarts + calmonths(1);
    dates = monthStarts + (monthEnds - monthStarts) / 2;
    lon = d.lon; lat = d.lat; sshErrors = [];
    nMonths = numel(dates);
    mapCounts = zeros(nMonths, 1, 'uint8');
    expectedMapCounts = mapCounts;
    partialMonths = false(nMonths, 1); srcIdxs = cell(nMonths, 1);
    % Extend the observed weekly phase across both edge months.
    weeklyDates = (d.dates(1) - days(35):days(7):d.dates(end) + days(35))';
    folder = fileparts(char(outPath));

    if isempty(folder)
        folder = pwd;
    end

    if ~isfolder(folder)
        mkdir(folder);
    end

    tmp = [tempname(folder), '.mat'];
    cleanup = onCleanup(@() deletefileifpresent(tmp));
    save(tmp, 'dates', 'lon', 'lat', 'sshErrors', 'monthStarts', 'monthEnds', '-v7.3');
    out = matfile(tmp, 'Writable', true);
    out.sshs(numel(lon), numel(lat), nMonths) = single(NaN);
    out.validMapCounts(numel(lon), numel(lat), nMonths) = uint8(0);
    out.incompleteCoverage(numel(lon), numel(lat), nMonths) = false;

    for iMonth = 1:nMonths
        ix = find(d.dates >= monthStarts(iMonth) & d.dates < monthEnds(iMonth));
        srcIdxs{iMonth} = ix; mapCounts(iMonth) = numel(ix);
        expected = weeklyDates(weeklyDates >= monthStarts(iMonth) & weeklyDates < monthEnds(iMonth));
        expectedMapCounts(iMonth) = numel(expected);
        partialMonths(iMonth) = any(expected < d.dates(1) | expected > d.dates(end));

        if isempty(ix)
            counts = zeros(numel(lon), numel(lat), 'uint8');
            avg = nan(numel(lon), numel(lat), 'single');
        else

            if isscalar(d.dates)
                block = double(src.sshs);
            else
                block = double(src.sshs(:, :, ix));
            end

            valid = isfinite(block);
            counts = uint8(sum(valid, 3));
            block(~valid) = 0;
            avg = single(sum(block, 3) ./ double(counts));
            avg(counts == 0) = NaN;
        end

        % MAT files omit the trailing dimension for a single month.
        page = {':', ':'};

        if nMonths > 1
            page{3} = iMonth;
        end

        out.sshs(page{:}) = avg;
        out.validMapCounts(page{:}) = counts;
        out.incompleteCoverage(page{:}) = counts < expectedMapCounts(iMonth);

        if ~options.BeQuiet && (mod(iMonth, 60) == 0 || iMonth == nMonths)
            fprintf('Monthly NASA-SSH: %d/%d\n', iMonth, nMonths);
        end

    end

    metadata = struct('product', 'NASASSHMonthly', 'sourceProduct', d.metadata.product, ...
        'sourceMetadata', d.metadata, 'sourceDates', d.dates, ...
        'sourceIndices', {srcIdxs}, ...
        'units', 'm', 'meanSeaSurface', 'DTU21', ...
        'method', 'Equal-weight finite-map calendar-month average by central epoch; double accumulation, single output.', ...
        'coveragePolicy', 'Retain all months; NaN if zero valid maps; incompleteCoverage if fewer than expected weekly maps; no threshold applied.', ...
        'caveat', 'Overlapping 10-day source windows; approximate monthly means; counts are not independent sample counts or uncertainty.');
    out.mapCounts = mapCounts;
    out.expectedMapCounts = expectedMapCounts;
    out.partialMonths = partialMonths;
    out.metadata = metadata;
    out.creationDate = datetime('now', 'TimeZone', 'UTC');
    clear out

    publishnewfile(tmp, outPath);

    if ~options.BeQuiet
        fprintf('[ULMO>%s] Saved %s (%s to %s)\n', ...
            mfilename, filehref(outPath, 'monthly averaged NASA-SSH data'), string(dates(1)), string(dates(end)));
    end

end
