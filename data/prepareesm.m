%% PREPAREESM - Imports ESAESM2 GFC files into a compact harmonic NetCDF file.
% Reads a source folder, including subdirectories, separately from ESM2PLMT.
% The source may reside on an external drive and is never modified.
%
% Syntax
%   outPath = prepareesm(inFolder, comp)
%   outPath = prepareesm(inFolder, comp, "ForceNew", true)
%
% Input arguments
%   inFolder - Folder containing mtmshc_COMPONENT_YYYYMMDD_HH.gfc files.
%   comp - 'A', 'O', 'H', 'I', 'S' or 'AOHIS'.
%   ForceNew (name-value) - Whether to replace an existing prepared NetCDF file.
%       Default option: false (an existing destination raises an error).
%   BeQuiet (name-value) - Whether to suppress the completion message.
%       Default option: false.
%
% Output arguments
%   outPath - Prepared NetCDF file at $IFILES/ESM/ESAESM2/ESAESM2_COMPONENT.nc.
%
% Output file content
%   degree, order - uint16 vectors on the shared lm dimension. Size: (M, 1).
%   coefficients - Fully normalized Stokes values, double. Size: (M, 2, N).
%       Dimensions: lm, cs (cosine then sine), time.
%   time - Double seconds since 2000-01-01 00:00:00 UTC. Size: (N, 1).
%   equatorRadius (global attribute) - Constant reference radius, metres.
%
% Notes
%   Imports at most 32 samples per batch to a temporary NetCDF-4
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
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)

function outPath = prepareesm(inFolder, comp, options)

    arguments (Input)
        inFolder {mustBeTextScalar, mustBeFolder}
        comp ...
            {mustBeTextScalar, mustBeMember(comp, ["A", "O", "H", "I", "S", "AOHIS"])}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
    end

    arguments (Output)
        outPath (1, 1) string
    end

    % Discover component files and resolve the configured prepared-data destination.
    [inFiles, outPath] = getIOpaths(inFolder, comp);

    % Parse timestamps from filenames before ordering the coefficient samples.
    dates = NaT(numel(inFiles), 1);

    for k = 1:numel(inFiles)
        token = regexp(inFiles(k).name, '^mtmshc_' + string(comp) + ...
            '_(\d{8}_\d{2})\.gfc$', 'tokens', 'once');
        assert(~isempty(token), 'ULMO:prepareesm:InvalidName', ...
            'Invalid ESM filename: %s', inFiles(k).name);
        dates(k) = datetime(token{1}, 'InputFormat', 'yyyyMMdd_HH');
    end

    % Keep files and timestamps in the same order; reject duplicate/off-grid samples.
    [dates, order] = sort(dates);
    inFiles = inFiles(order);
    assert(numel(unique(dates)) == numel(dates) && all(mod(hour(dates), 6) == 0), ...
        'ULMO:prepareesm:InvalidSampling', 'Duplicate or off-grid %s samples.', comp);
    % Use the first sample to establish the shared harmonic labels and radius.
    [first, eqRad] = readesmfile(fullfile(inFiles(1).folder, inFiles(1).name));
    labels = first(:, 1:2);

    % Build a temporary output so failed imports do not replace prepared data.
    tmpPath = string(tempname(outFolder)) + ".nc";
    cleanup = onCleanup(@() removetemporary(tmpPath));
    createshnc(tmpPath, labels, dates, 'GRAV', 'Attributes', {'equatorRadius', eqRad});

    % Read and write file in batches
    for firstIdx = 1:32:numel(inFiles)
        ids = firstIdx:min(firstIdx + 31, numel(inFiles));
        batch = nan(size(labels, 1), 2, numel(ids));

        % Require every sample to use the same harmonic rows and reference radius.
        for j = 1:numel(ids)
            k = ids(j);
            [plm, fileEqRad] = readesmfile(fullfile(inFiles(k).folder, inFiles(k).name));
            assert(isequal(plm(:, 1:2), labels), 'ULMO:prepareesm:InvalidCoefficients', ...
                'Harmonic rows differ in %s.', inFiles(k).name);
            assert(fileEqRad == eqRad, 'ULMO:prepareesm:InconsistentRadius', ...
                'Reference radius differs in %s.', inFiles(k).name);
            batch(:, :, j) = plm(:, 3:4);
        end

        % Write the batch directly into its time slab instead of retaining the record.
        ncwrite(tmpPath, 'coefficients', batch, [1, 1, ids(1)]);
    end

    % Publish the completed input before removing derived caches that are now stale.
    movefile(tmpPath, outPath, 'f');
    caches = dir(fullfile(outFolder, "ESAESM2_" + string(comp) + "-*-T*.nc"));

    for k = 1:numel(caches)
        delete(fullfile(caches(k).folder, caches(k).name));
    end

    if ~options.BeQuiet
        fprintf('[ULMO>%s] Saved %s (%d samples).\n', mfilename, outPath, numel(dates));
    end

end

% Delete partial output on failure; do nothing after a successful move.
function removetemporary(path)

    if isfile(path)
        delete(path);
    end

end

% Parse one GFC header and coefficient table without changing its normalization.
function [plm, eqRad] = readesmfile(path)
    contents = fileread(path);
    endHeader = regexp(contents, '(?m)^end_of_head[^\r\n]*\r?\n', 'end', 'once');
    assert(~isempty(endHeader), ...
        'ULMO:prepareesm:InvalidHeader', 'Missing end_of_head in %s', path);
    radTxt = regexp(contents(1:endHeader), '(?m)^radius[ \t]+(\S+)', 'tokens', 'once');
    eqRad = NaN;

    if ~isempty(radTxt)
        eqRad = str2double(regexprep(radTxt{1}, '[dD]', 'E'));
    end

    assert(isscalar(eqRad) && isfinite(eqRad) && eqRad > 0, ...
        'ULMO:prepareesm:InvalidRadius', 'Missing or invalid radius in %s', path);
    % Accept Fortran D exponents, then sort coefficients into degree-major order.
    data = regexprep(contents(endHeader + 1:end), '[dD]', 'E');
    cols = textscan(data, '%s%f%f%f%f', 'ReturnOnError', false);
    plm = [cols{2:5}];
    plm = sortrows(plm, [1, 2]);
    % Compare against the full triangular label set to detect missing or duplicate rows.
    assert(~isempty(plm) && all(strcmp(cols{1}, 'gfc')) && all(isfinite(plm), 'all') && ...
        all(plm(:, 1:2) >= 0 & plm(:, 1:2) == fix(plm(:, 1:2)), 'all') && ...
        isequal(plm(:, 1:2), harmoniclabels(plm(end, 1))), ...
        'ULMO:prepareesm:InvalidCoefficients', 'Incomplete or invalid harmonics in %s', path);
end

% Resolve the output under IFILES and recursively discover source component files.
function [inFiles, outPath] = getIOpaths(inFolder, comp)
    dataRoot = getenv('IFILES');
    assert(~isempty(dataRoot), 'ULMO:prepareesm:MissingDataRoot', ...
    'Set IFILES before importing ESM fields.');
    outFolder = fullfile(dataRoot, 'ESM', 'ESAESM2');
    outPath = fullfile(outFolder, "ESAESM2_" + string(comp) + ".nc");
    assert(options.ForceNew || ~isfile(outPath), ...
        'ULMO:prepareesm:FileExists', ...
        'Prepared data already exist. Use ForceNew=true to replace %s.', outPath);
    inFiles = dir(fullfile(inFolder, '**', "mtmshc_" + string(comp) + "_*.gfc"));
    assert(~isempty(inFiles), ...
        'ULMO:prepareesm:MissingData', 'No ESM files found in %s.', inFolder);

    if ~isfolder(outFolder)
        mkdir(outFolder);
    end

end

% Construct the complete (degree, order) sequence expected in each GFC file.
function labels = harmoniclabels(degree)
    labels = zeros((degree + 1) * (degree + 2) / 2, 2);

    for l = 0:degree
        rows = l * (l + 1) / 2 + (1:l + 1);
        labels(rows, :) = [repmat(l, l + 1, 1), (0:l)'];
    end

end
