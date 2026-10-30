%% PREPAREESM - Imports ESAESM2 GFC files into a compact harmonic NetCDF file.
% Reads a source folder, including subdirectories, separately from ESM2PLMT.
% The source may reside on an external drive and is never modified.
%
% Syntax
%   path = prepareesm(folder, component)
%   path = prepareesm(folder, component, "ForceNew", true)
%
% Input arguments
%   folder - Folder containing mtmshc_COMPONENT_YYYYMMDD_HH.gfc files.
%   component - 'A', 'O', 'H', 'I', 'S' or 'AOHIS'.
%   ForceNew (name-value) - Whether to replace an existing prepared NetCDF file.
%       Default option: false (an existing destination raises an error).
%   BeQuiet (name-value) - Whether to suppress the completion message.
%       Default option: false.
%
% Output arguments
%   path - Prepared NetCDF file at $IFILES/ESM/ESAESM2/ESAESM2_COMPONENT.nc.
%
% Output file content
%   degree, order - uint16 vectors on the shared lm dimension. Size: (M, 1).
%   coefficients - Fully normalized Stokes values, double. Size: (M, 2, N).
%       Dimensions: lm, cs (cosine then sine), time.
%   time - Double seconds since 2000-01-01 00:00:00 UTC. Size: (N, 1).
%   equatorRadius (global attribute) - Constant reference radius, metres.
%
% Notes
%   Requires IFILES. Imports at most 32 samples per batch to a temporary NetCDF-4
%   file, then replaces the destination after successful validation/import.
%   Duplicate/off-grid timestamps and missing harmonic rows are rejected.
%   Formal errors are absent from ESAESM2; no redundant NaN array is stored.
%   Refreshing a component removes its derived -T*.nc caches in all units.
%   Reimport the entire desired source collection when adding data; this
%   replaces the prepared record rather than appending to it.
%
% See also
%   ESM2PLMT
%
% Created by
%   2026/10/30

function path = prepareesm(folder, component, options)
    arguments (Input)
        folder {mustBeTextScalar, mustBeFolder}
        component {mustBeTextScalar, mustBeMember(component, ...
            ["A", "O", "H", "I", "S", "AOHIS"])}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
    end
    arguments (Output)
        path (1, 1) string
    end
    root = getenv('IFILES');
    assert(~isempty(root), 'ULMO:prepareesm:MissingDataRoot', 'Set IFILES before importing ESM.');
    destination = fullfile(root, 'ESM', 'ESAESM2');
    path = fullfile(destination, "ESAESM2_" + string(component) + ".nc");
    assert(options.ForceNew || ~isfile(path), 'ULMO:prepareesm:FileExists', ...
        'Prepared data already exist. Use ForceNew=true to replace %s.', path);
    files = dir(fullfile(folder, '**', "mtmshc_" + string(component) + "_*.gfc"));
    assert(~isempty(files), 'ULMO:prepareesm:MissingData', 'No ESM files found in %s.', folder);
    dates = NaT(numel(files), 1);
    for k = 1:numel(files)
        token = regexp(files(k).name, '^mtmshc_' + string(component) + ...
            '_(\d{8}_\d{2})\.gfc$', 'tokens', 'once');
        assert(~isempty(token), 'ULMO:prepareesm:InvalidName', ...
            'Invalid ESM filename: %s', files(k).name);
        dates(k) = datetime(token{1}, 'InputFormat', 'yyyyMMdd_HH');
    end
    [dates, order] = sort(dates);
    files = files(order);
    assert(numel(unique(dates)) == numel(dates) && all(mod(hour(dates), 6) == 0), ...
        'ULMO:prepareesm:InvalidSampling', 'Duplicate or off-grid %s samples.', component);
    [first, equatorRadius] = readesmfile(fullfile(files(1).folder, files(1).name));
    labels = first(:, 1:2);
    if ~isfolder(destination)
        mkdir(destination);
    end
    temporary = string(tempname(destination)) + ".nc";
    cleanup = onCleanup(@() removetemporary(temporary));
    createshnc(temporary, labels, dates, 'GRAV', 'Attributes', {'equatorRadius', equatorRadius});
    for firstIndex = 1:32:numel(files)
        ids = firstIndex:min(firstIndex + 31, numel(files));
        batch = nan(size(labels, 1), 2, numel(ids));
        for j = 1:numel(ids)
            k = ids(j);
            [plm, radius] = readesmfile(fullfile(files(k).folder, files(k).name));
            assert(isequal(plm(:, 1:2), labels), 'ULMO:prepareesm:InvalidCoefficients', ...
                'Harmonic rows differ in %s.', files(k).name);
            assert(radius == equatorRadius, 'ULMO:prepareesm:InconsistentRadius', ...
                'Reference radius differs in %s.', files(k).name);
            batch(:, :, j) = plm(:, 3:4);
        end
        ncwrite(temporary, 'coefficients', batch, [1, 1, ids(1)]);
    end
    movefile(temporary, path, 'f');
    caches = dir(fullfile(destination, "ESAESM2_" + string(component) + "_*-T*.nc"));
    for k = 1:numel(caches)
        delete(fullfile(caches(k).folder, caches(k).name));
    end
    if ~options.BeQuiet
        fprintf('[ULMO>prepareesm] Saved %s (%d samples).\n', path, numel(dates));
    end
end

function removetemporary(path)
    if isfile(path)
        delete(path);
    end
end

function [plm, equatorRadius] = readesmfile(path)
    contents = fileread(path);
    endHeader = regexp(contents, '(?m)^end_of_head[^\r\n]*\r?\n', 'end', 'once');
    assert(~isempty(endHeader), 'ULMO:prepareesm:InvalidHeader', ...
        'Missing end_of_head in %s', path);
    radiusText = regexp(contents(1:endHeader), '(?m)^radius[ \t]+(\S+)', 'tokens', 'once');
    equatorRadius = NaN;
    if ~isempty(radiusText)
        equatorRadius = str2double(regexprep(radiusText{1}, '[dD]', 'E'));
    end
    assert(isscalar(equatorRadius) && isfinite(equatorRadius) && equatorRadius > 0, ...
        'ULMO:prepareesm:InvalidRadius', 'Missing or invalid radius in %s', path);
    data = regexprep(contents(endHeader + 1:end), '[dD]', 'E');
    columns = textscan(data, '%s%f%f%f%f', 'ReturnOnError', false);
    plm = [columns{2:5}];
    plm = sortrows(plm, [1, 2]);
    assert(~isempty(plm) && all(strcmp(columns{1}, 'gfc')) && all(isfinite(plm), 'all') && ...
        all(plm(:, 1:2) >= 0 & plm(:, 1:2) == fix(plm(:, 1:2)), 'all') && ...
        isequal(plm(:, 1:2), harmoniclabels(plm(end, 1))), ...
        'ULMO:prepareesm:InvalidCoefficients', 'Incomplete or invalid harmonics in %s', path);
end

function labels = harmoniclabels(degree)
    labels = zeros((degree + 1) * (degree + 2) / 2, 2);
    for l = 0:degree
        rows = l * (l + 1) / 2 + (1:l + 1);
        labels(rows, :) = [repmat(l, l + 1, 1), (0:l)'];
    end
end
