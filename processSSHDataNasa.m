% PROCESSSSHDATANASA Aggregate NASA_SSH_REF_SIMPLE_GRID_V1 for ssh2lonlatt.
% Preserves native weekly epochs, DTU21 anomalies (metres), and missing data.
% Output is single precision (longitude, latitude, time), MATLAB v7.3.
% No uncertainty is supplied by this product; sshErrors is empty.
% Existing outputs are never overwritten. Partial output is not published.
function outputPath = processSSHDataNasa(inputFolder, outputPath)

    arguments (Input)
        inputFolder (1, :) char
        outputPath (1, :) char = fullfile(getenv('IFILES'), 'SSH', 'NASASSH.mat')
    end

    assert(~isfile(outputPath), 'ULMO:OutputExists', 'Output already exists: %s', outputPath);
    files = dir(fullfile(inputFolder, 'NASA-SSH_alt_ref_simple_grid_v1_*.nc'));
    assert(~isempty(files), 'ULMO:InputDataNotFound', 'No NASA-SSH files in %s', inputFolder);
    [~, order] = sort({files.name});
    files = files(order);
    first = fullfile(inputFolder, files(1).name);
    lon = ncread(first, 'longitude'); lat = ncread(first, 'latitude');
    lon = lon(:); lat = lat(:);
    assert(all(diff(lon) > 0) && all(diff(lat) > 0), 'Coordinates must be increasing.');
    n = numel(files); dates = NaT(n, 1); finiteCounts = zeros(n, 1);
    coverageStart = strings(n, 1); coverageEnd = strings(n, 1);
    folder = fileparts(outputPath);

    if ~isfolder(folder)
        mkdir(folder);
    end

    tmpPath = [tempname(folder), '.mat'];
    cleanup = onCleanup(@() removePartial(tmpPath));
    sshErrors = [];
    save(tmpPath, 'lon', 'lat', 'sshErrors', '-v7.3');
    m = matfile(tmpPath, 'Writable', true);
    m.sshs(numel(lon), numel(lat), n) = single(NaN);

    for k = 1:n
        p = fullfile(inputFolder, files(k).name);
        assert(strcmp(ncreadatt(p, '/', 'product_short_name'), 'NASA_SSH_REF_SIMPLE_GRID_V1'));
        assert(strcmp(ncreadatt(p, '/', 'mean_sea_surface'), 'DTU21'));
        assert(strcmp(ncreadatt(p, 'ssha', 'units'), 'm'));
        assert(strcmp(ncreadatt(p, 'time', 'units'), 'seconds since 1990-01-01'));
        assert(isequal(lon, ncread(p, 'longitude')) && isequal(lat, ncread(p, 'latitude')));
        info = ncinfo(p, 'ssha');
        assert(isequal({info.Dimensions.Name}, {'longitude', 'latitude'}));
        field = ncread(p, 'ssha');
        field(field < ncreadatt(p, 'ssha', 'valid_min') | ...
            field > ncreadatt(p, 'ssha', 'valid_max') | ~isfinite(field)) = NaN;
        assert(isequal(size(field), [numel(lon), numel(lat)]));
        dates(k) = datetime(1990, 1, 1) + seconds(ncread(p, 'time'));
        stamp = regexp(files(k).name, '(\d{8})\.nc$', 'tokens', 'once');
        assert(dates(k) == datetime(stamp{1}, 'InputFormat', 'yyyyMMdd'));

        if k > 1
            assert(dates(k) > dates(k - 1), 'Duplicate or unordered epochs.');
        end

        finiteCounts(k) = nnz(isfinite(field));
        % Entirely missing maps are legitimate source gaps; retain their epochs.
        coverageStart(k) = string(ncreadatt(p, '/', 'time_coverage_start'));
        coverageEnd(k) = string(ncreadatt(p, '/', 'time_coverage_end'));
        m.sshs(:, :, k) = single(field);

        if mod(k, 100) == 0 || k == n
            fprintf('NASA-SSH: %d/%d maps\n', k, n);
        end

    end

    metadata = struct('product', 'NASA_SSH_REF_SIMPLE_GRID_V1', ...
        'doi', '10.5067/NSREF-SG0V1', 'meanSeaSurface', 'DTU21', ...
        'units', 'm', 'ordering', 'longitude,latitude,time', ...
        'inputFolder', inputFolder, 'sourceFiles', {{files.name}}, ...
        'sourceBytes', [files.bytes], 'sourceModifiedDatenum', [files.datenum], ...
        'finiteCounts', finiteCounts, 'coverageStart', coverageStart, ...
        'coverageEnd', coverageEnd, 'created', datetime('now'), ...
        'griddingMethod', ncreadatt(first, '/', 'gridding_method'), ...
        'processing', 'Native ssha; invalid values to NaN; single precision; no additional corrections or rebaselining.');
    m.dates = dates; m.metadata = metadata;
    clear m
    assert(~isfile(outputPath), 'Output appeared while processing.');
    movefile(tmpPath, outputPath);
    fprintf('Saved %s (%s to %s)\n', outputPath, string(dates(1)), string(dates(end)));
end

function removePartial(path)

    if isfile(path)
        delete(path);
    end

end
