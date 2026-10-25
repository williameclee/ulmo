% PROCESSSSHMONTHLYNASA Calendar-month averages of native NASA-SSH maps.
% Equal weights for maps with central timestamps in each calendar month.
% Finite values are averaged in double precision, saved as single metres.
% No valid values => NaN. Overlapping 10-day windows are not deconvolved.
% validMapCounts and expectedMapCounts describe sampling, not uncertainty.
% incompleteCoverage flags cells missing any expected weekly sample.
% partialMonths flags months whose expected epochs extend beyond input dates.
% Retains partial edge months; does not impose an arbitrary coverage cutoff.

function outputPath = processSSHMonthlyNasa(inputPath, outputPath)

    arguments (Input)
        inputPath (1, :) char = fullfile(getenv('IFILES'), 'SSH', 'NASASSH.mat')
        outputPath (1, :) char = fullfile(getenv('IFILES'), 'SSH', 'NASASSHMonthly.mat')
    end

    assert(~isfile(outputPath), 'ULMO:OutputExists', 'Output exists: %s', outputPath);
    d = load(inputPath, 'dates', 'lon', 'lat', 'metadata');
    assert(strcmp(d.metadata.product, 'NASA_SSH_REF_SIMPLE_GRID_V1'));
    assert(numel(d.dates) > 1 && all(diff(d.dates) > days(0)));
    assert(all(mod(days(d.dates - d.dates(1)), 7) == 0), 'Input is not on a weekly grid.');
    src = matfile(inputPath);
    assert(isequal(size(src, 'sshs'), [numel(d.lon), numel(d.lat), numel(d.dates)]));
    monthStarts = (dateshift(d.dates(1), 'start', 'month'):calmonths(1): ...
        dateshift(d.dates(end), 'start', 'month'))';
    monthEnds = monthStarts + calmonths(1);
    dates = monthStarts + (monthEnds - monthStarts) / 2;
    lon = d.lon; lat = d.lat; sshErrors = [];
    n = numel(dates); mapCounts = zeros(n, 1, 'uint8'); expectedMapCounts = mapCounts;
    partialMonths = false(n, 1); sourceIndices = cell(n, 1);
    % Extend the observed weekly phase across both edge months.
    weeklyDates = (d.dates(1) - days(35):days(7):d.dates(end) + days(35))';
    folder = fileparts(outputPath);

    if isempty(folder)
        folder = pwd;
    end

    if ~isfolder(folder)
        mkdir(folder);
    end

    tmp = [tempname(folder), '.mat'];
    cleanup = onCleanup(@() removePartial(tmp));
    save(tmp, 'dates', 'lon', 'lat', 'sshErrors', 'monthStarts', 'monthEnds', '-v7.3');
    out = matfile(tmp, 'Writable', true);
    out.sshs(numel(lon), numel(lat), n) = single(NaN);
    out.validMapCounts(numel(lon), numel(lat), n) = uint8(0);
    out.incompleteCoverage(numel(lon), numel(lat), n) = false;

    for k = 1:n
        ix = find(d.dates >= monthStarts(k) & d.dates < monthEnds(k));
        sourceIndices{k} = ix; mapCounts(k) = numel(ix);
        expected = weeklyDates(weeklyDates >= monthStarts(k) & weeklyDates < monthEnds(k));
        expectedMapCounts(k) = numel(expected);
        partialMonths(k) = any(expected < d.dates(1) | expected > d.dates(end));

        if isempty(ix)
            counts = zeros(numel(lon), numel(lat), 'uint8');
            avg = nan(numel(lon), numel(lat), 'single');
        else
            block = double(src.sshs(:, :, ix));
            valid = isfinite(block); counts = uint8(sum(valid, 3));
            block(~valid) = 0;
            avg = single(sum(block, 3) ./ double(counts));
            avg(counts == 0) = NaN;
        end

        out.sshs(:, :, k) = avg; out.validMapCounts(:, :, k) = counts;
        out.incompleteCoverage(:, :, k) = counts < expectedMapCounts(k);

        if mod(k, 60) == 0 || k == n
            fprintf('Monthly NASA-SSH: %d/%d\n', k, n);
        end

    end

    metadata = struct('product', 'NASASSHMonthly', 'sourceProduct', d.metadata.product, ...
        'sourcePath', inputPath, 'sourceMetadata', d.metadata, 'sourceDates', d.dates, ...
        'sourceIndices', {sourceIndices}, 'created', datetime('now'), ...
        'units', 'm', 'meanSeaSurface', 'DTU21', ...
        'method', 'Equal-weight finite-map calendar-month average by central epoch; double accumulation, single output.', ...
        'coveragePolicy', 'Retain all months; NaN if zero valid maps; incompleteCoverage if fewer than expected weekly maps; no threshold applied.', ...
        'caveat', 'Overlapping 10-day source windows; approximate monthly means; counts are not independent sample counts or uncertainty.');
    out.mapCounts = mapCounts; out.expectedMapCounts = expectedMapCounts;
    out.partialMonths = partialMonths; out.metadata = metadata;
    clear out
    assert(~isfile(outputPath), 'Output appeared during processing.');
    movefile(tmp, outputPath);
    fprintf('Saved %s\n', outputPath);
end

function removePartial(path)

    if isfile(path)
        delete(path);
    end

end
