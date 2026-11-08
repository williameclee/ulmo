%% PROCESSSSHDATANASA - Converts NASA-SSH weekly anomaly maps to a MAT file.
% Preserve the native grid, central epochs and DTU21 reference surface.
% Replace invalid cells with NaN and retain epochs with entirely missing maps.
%
% Syntax
%   outPath = processSSHDataNasa(inFolder)
%   outPath = processSSHDataNasa(inFolder, outPath)
%   outPath = processSSHDataNasa(__, "BeQuiet", true)
%
% Input arguments
%   inFolder - Folder containing NASA-SSH reference-mission NetCDF maps
%       named NASA-SSH_alt_ref_simple_grid_v1_YYYYMMDD.nc.
%   outPath (optional) - Destination MAT file.
%       Default path: fullfile(getenv('IFILES'), 'SSH', 'NASASSH.mat')
%   BeQuiet (name-value) - Whether to suppress progress and completion messages.
%       Default value: false
%
% Output arguments
%   outPath - Path to the weekly MAT file.
%
% Output file content
%   sshs - Single-precision sea level anomalies relative to DTU21.
%       Size: (nLon, nLat, nWeek).
%       Unit: metres.
%   lon, lat - Increasing longitude and latitude column vectors.
%       Size: (nLon, 1) and (nLat, 1), respectively.
%       Unit: degrees.
%   dates - Central epoch of each source map, in chronological order.
%       Size: (nWeek, 1).
%   sshErrors - Empty array. The source supplies no uncertainty.
%   metadata - Product identity, coverage intervals, gridding method,
%       units, dimension order and conversion details.
%   creationDate - UTC datetime when the converter wrote the file.
%
% Data sources
%   NASA-SSH reference-mission simple grid, version 1. https://doi.org/10.5067/NSREF-SG0V1
%
% See also
%   PROCESSSSHMONTHLYNASA, SSH2LONLATT
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function outPath = processSSHDataNasa(inFolder, outPath, options)

    arguments (Input)
        inFolder {mustBeTextScalar, mustBeFolder}
        outPath {mustBeTextScalar} = fullfile(getenv('IFILES'), 'SSH', 'NASASSH.mat')
        options.BeQuiet (1, 1) logical = false
    end

    assert(~isfile(outPath) && ~isfolder(outPath), 'ULMO:OutputExists', 'Output already exists: %s', outPath);
    files = dir(fullfile(inFolder, 'NASA-SSH_alt_ref_simple_grid_v1_*.nc'));
    assert(~isempty(files), 'ULMO:InputDataNotFound', 'No NASA-SSH files in %s', inFolder);
    [~, order] = sort({files.name});
    files = files(order);
    first = fullfile(inFolder, files(1).name);
    lon = ncread(first, 'longitude');
    lat = ncread(first, 'latitude');
    lon = lon(:);
    lat = lat(:);
    assert(all(diff(lon) > 0) && all(diff(lat) > 0), 'Coordinates must be increasing.');
    nFiles = numel(files);
    dates = NaT(nFiles, 1);
    covStart = strings(nFiles, 1); covEnd = strings(nFiles, 1);
    folder = fileparts(char(outPath));

    if isempty(folder)
        folder = pwd;
    elseif ~isfolder(folder)
        mkdir(folder);
    end

    tmpPath = [tempname(folder), '.mat'];
    cleanup = onCleanup(@() deletefileifpresent(tmpPath));
    sshErrors = [];
    save(tmpPath, 'lon', 'lat', 'sshErrors', '-v7.3');
    m = matfile(tmpPath, 'Writable', true);
    m.sshs(numel(lon), numel(lat), nFiles) = single(NaN);

    for iFile = 1:nFiles
        p = fullfile(inFolder, files(iFile).name);
        assert(strcmp(ncreadatt(p, '/', 'product_short_name'), 'NASA_SSH_REF_SIMPLE_GRID_V1'));
        assert(strcmp(ncreadatt(p, '/', 'mean_sea_surface'), 'DTU21'));
        assert(strcmp(ncreadatt(p, 'ssha', 'units'), 'm'), ...
            'ULMO:InvalidSSHUnits', 'Expected ssha in metres in %s.', p);
        assert(strcmp(ncreadatt(p, 'time', 'units'), 'seconds since 1990-01-01'));
        assert(isequal(lon, ncread(p, 'longitude')) && isequal(lat, ncread(p, 'latitude')));
        info = ncinfo(p, 'ssha');
        assert(isequal({info.Dimensions.Name}, {'longitude', 'latitude'}));
        field = ncread(p, 'ssha');
        field(field < ncreadatt(p, 'ssha', 'valid_min') | ...
            field > ncreadatt(p, 'ssha', 'valid_max') | ~isfinite(field)) = NaN;
        assert(isequal(size(field), [numel(lon), numel(lat)]));
        dates(iFile) = datetime(1990, 1, 1) + seconds(ncread(p, 'time'));
        stamp = regexp(files(iFile).name, '(\d{8})\.nc$', 'tokens', 'once');
        assert(dates(iFile) == datetime(stamp{1}, 'InputFormat', 'yyyyMMdd'));

        if iFile > 1
            assert(dates(iFile) > dates(iFile - 1), 'Duplicate or unordered epochs.');
        end

        % Entirely missing maps are legitimate source gaps; retain their epochs.
        covStart(iFile) = string(ncreadatt(p, '/', 'time_coverage_start'));
        covEnd(iFile) = string(ncreadatt(p, '/', 'time_coverage_end'));

        if nFiles == 1
            % A MAT file stores a singleton time dimension as a 2-D map.
            m.sshs = single(field);
        else
            m.sshs(:, :, iFile) = single(field);
        end

        if ~options.BeQuiet && (mod(iFile, 100) == 0 || iFile == nFiles)
            fprintf('NASA-SSH: %d/%d maps\n', iFile, nFiles);
        end

    end

    metadata = struct('product', 'NASA_SSH_REF_SIMPLE_GRID_V1', ...
        'doi', '10.5067/NSREF-SG0V1', 'meanSeaSurface', 'DTU21', ...
        'units', 'm', 'ordering', 'longitude,latitude,time', ...
        'coverageStart', covStart, ...
        'coverageEnd', covEnd, ...
        'griddingMethod', ncreadatt(first, '/', 'gridding_method'), ...
        'processing', 'Native ssha; invalid values to NaN; single precision; no additional corrections or rebaselining.');

    m.dates = dates;
    m.metadata = metadata;
    m.creationDate = datetime('now', 'TimeZone', 'UTC');
    clear m

    publishnewfile(tmpPath, outPath);

    if ~options.BeQuiet
        fprintf('[ULMO>%s] Saved %s (%s to %s)\n', ...
            mfilename, filehref(outPath, 'NASA-SSH data'), string(dates(1)), string(dates(end)));
    end

end
