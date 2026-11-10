%% PROCESSGFZNORMALEQUATIONS - Converts monthly GFZ normal equations to a NetCDF series.
% Stores the scaled marginal covariance factor for each available GSM solution.
%
% Syntax
%   ncPath = processGfzNormalEquations(sourceFolder)
%   ncPath = processGfzNormalEquations(sourceFolder, "Loutput", 60)
%
% Input arguments
%   sourceFolder - Folder containing GFZ GRACE RL06 and/or GRACE-FO RL06.3
%       degree-96 GSM SINEX files (.snx or .snx.gz).
%   Loutput (name-value) - Highest retained spherical harmonic degree.
%       Default degree: 60. Supported range: 2 through the source maximum.
%
% Output arguments
%   ncPath - Completed GFZ/NEQ/GFZ_NEQ_L<degree>.nc under GRACEDATA,
%       or under IFILES/GRACE when GRACEDATA is unset. Rebuilds the series
%       from the available epochs on each call.
%
% Output file content
%   noiseFactor - Lower Cholesky factor of each monthly marginal covariance
%       in dimensionless Stokes units. Packed by columns from the diagonal
%       to the last row. Size: (n*(n+1)/2, nTime), n = (Loutput+1)^2-4.
%   degree, order, cs - Shared parameter labels; cs is 0 for C and 1 for S.
%   formalStd - Published dimensionless standard deviations. Size: (n, nTime).
%   time, epochStart, epochEnd - SINEX epochs in seconds since 2000-01-01 UTC.
%   sourceFilename, sourceMaxDegree, referenceRadius, varianceFactor - Monthly
%       source basename, maximum degree, radius in metres and variance scale.
%   creationDate - Scalar completion time in seconds since 1970-01-01 UTC.
%
% Data sources
%   GFZ GRACE RL06 and GRACE-FO RL06.3 GSM normal equations, GFZ ISDC.
%
% Created by
%   2026/11/10
function ncPath = processGfzNormalEquations(sourceFolder, options)
    arguments (Input)
        sourceFolder (1, 1) string {mustBeFolder}
        options.Loutput (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.Loutput, 2)} = 60
    end

    files = dir(fullfile(sourceFolder, 'GSM-2_*.snx*'));
    pattern = '^GSM-2_\d{7}-\d{7}_(GRAC_GFZOP_BB01_0600|GRFO_GFZOP_BB01_0603)\.snx(\.gz)?$';
    files = files(~[files.isdir] & ~cellfun(@isempty, regexp({files.name}, pattern, 'once')));
    assert(~isempty(files), 'ULMO:gfzneq:NoFiles', 'No GFZ normal equations in %s.', sourceFolder);
    for iFile = 1:numel(files)
        path = fullfile(sourceFolder, files(iFile).name);
        months(iFile) = readHeader(path, options.Loutput); %#ok<AGROW>
    end
    [times, order] = sort([months.time]);
    months = months(order);
    duplicate = find(diff(times) == 0, 1);
    if ~isempty(duplicate)
        error('ULMO:gfzneq:DuplicateEpoch', 'Duplicate epoch in %s and %s.', ...
            months(duplicate).name, months(duplicate+1).name);
    end
    assert(all(strcmp({months.norm}, months(1).norm)) && all(strcmp({months.tide}, months(1).tide)), ...
        'ULMO:gfzneq:Conventions', 'Source normalisations or tide systems differ.');

    root = getenv('GRACEDATA');
    if isempty(root)
        root = getenv('IFILES');
        assert(~isempty(root), 'ULMO:gfzneq:NoRoot', 'Set GRACEDATA or IFILES.');
        root = fullfile(root, 'GRACE');
    end
    folder = fullfile(root, 'GFZ', 'NEQ');
    if ~isfolder(folder)
        mkdir(folder);
    end
    ncPath = fullfile(folder, sprintf('GFZ_NEQ_L%d.nc', options.Loutput));
    tmp = [tempname(folder), '.nc'];
    cleanup = onCleanup(@() deletefileifpresent(tmp));
    createSeries(tmp, months, options.Loutput);
    for iMonth = 1:numel(months)
        month = months(iMonth);
        [packed, scale] = monthlyFactor(fullfile(sourceFolder, month.name), month);
        ncwrite(tmp, 'noiseFactor', packed, [1 iMonth]);
        ncwrite(tmp, 'formalStd', month.std, [1 iMonth]);
        ncwrite(tmp, 'varianceFactor', scale, iMonth);
    end
    ncwrite(tmp, 'creationDate', posixtime(datetime('now', 'TimeZone', 'UTC')));
    movefile(tmp, ncPath, 'f');
end

function month = readHeader(path, L)
    % Small header reads keep preflight independent of matrix size.
    stream = openSource(path);
    cleanup = onCleanup(@() stream.close());
    header = '';
    marker = '(?m)^\+SOLUTION/NORMAL_EQUATION_MATRIX[^\r\n]*\r?\n';
    last = [];
    while isempty(last)
        block = char(stream.readNBytes(128*1024)).';
        assert(~isempty(block), 'ULMO:gfzneq:InvalidSource', 'Missing normal matrix in %s.', path);
        header = [header block]; %#ok<AGROW>
        [first, last] = regexp(header, marker, 'once');
        % A damaged opening marker must fail before buffering the dense matrix.
        if isempty(last) && ~isempty(regexp(header, '(?m)^ +\d+ +\d+ +[+.\d-]', 'once'))
            error('ULMO:gfzneq:InvalidSource', 'Matrix entries before the normal matrix header in %s.', path);
        end
    end
    assert(strcmp(strtrim(header(first:last)), '+SOLUTION/NORMAL_EQUATION_MATRIX L'), ...
        'ULMO:gfzneq:InvalidSource', 'Expected a lower normal matrix in %s.', path);
    month.offset = last;
    [~, stem, ext] = fileparts(path);
    month.name = char(stem + ext);
    firstLine = split(strtrim(extractBefore(string(header), newline)));
    assert(numel(firstLine) >= 9 && firstLine(1) == "%=SNX", ...
        'ULMO:gfzneq:InvalidSource', 'Invalid SINEX header in %s.', path);
    month.start = epochSeconds(firstLine(6));
    month.end = epochSeconds(firstLine(7));
    month.count = str2double(firstLine(9));
    month.maxDegree = str2double(headerField(header, 'max_degree', path));
    month.radius = str2double(headerField(header, 'radius', path));
    month.norm = headerField(header, 'norm', path);
    month.tide = headerField(header, 'tide_system', path);
    assert(month.maxDegree == fix(month.maxDegree) && month.count == (month.maxDegree+1)^2-4 ...
        && month.radius > 0 && isfinite(month.radius), 'ULMO:gfzneq:InvalidSource', ...
        'Invalid source degree, parameter count or radius in %s.', path);
    assert(L <= month.maxDegree, 'ULMO:gfzneq:Degree', 'Requested degree exceeds degree %d in %s.', ...
        month.maxDegree, path);

    block = regexp(header, '(?ms)^\+SOLUTION/ESTIMATE\s*?\r?\n(.*?)^-SOLUTION/ESTIMATE', 'tokens', 'once');
    assert(~isempty(block), 'ULMO:gfzneq:InvalidSource', 'Missing estimates in %s.', path);
    est = textscan(block{1}, '%f%s%f%s%f%s%s%f%f%f', 'CommentStyle', '*', 'ReturnOnError', false);
    idx = est{1};
    cs = strcmp(est{2}, 'SN');
    degree = est{3};
    order = est{5};
    [labels, perm] = sortrows([degree order cs]);
    valid = degree >= 2 & degree <= month.maxDegree & order >= 0 & order <= degree & ...
        degree == fix(degree) & order == fix(order) & (cs | strcmp(est{2}, 'CN')) & ~(order == 0 & cs);
    assert(isequal(sort(idx), (1:month.count)') && all(valid) && ...
        ~any(all(diff(labels) == 0, 2)) && all(isfinite(est{10}) & est{10} > 0) && ...
        all(strcmp(est{6}, est{6}{1})), 'ULMO:gfzneq:InvalidSource', ...
        'Invalid coefficient map, uncertainties or epochs in %s.', path);
    month.time = epochSeconds(est{6}{1});
    assert(month.start <= month.time && month.time <= month.end, ...
        'ULMO:gfzneq:InvalidSource', 'Reference epoch outside solution interval in %s.', path);
    keep = labels(:, 1) <= L;
    month.labels = labels(keep, :);
    month.idx = idx(perm(keep));
    month.std = est{10}(perm(keep));
end

function value = headerField(header, name, path)
    token = regexp(header, ['(?m)^[ \t]+' name '[ \t]+(\S+)'], 'tokens', 'once');
    assert(~isempty(token), 'ULMO:gfzneq:InvalidSource', 'Missing %s in %s.', name, path);
    value = token{1};
end

function value = epochSeconds(token)
    assert(~isempty(regexp(token, '^\d{2}:\d{3}:\d{5}$', 'once')), ...
        'ULMO:gfzneq:InvalidEpoch', 'Invalid SINEX epoch %s.', token);
    parts = sscanf(token, '%d:%d:%d');
    year = 2000 + parts(1) - 100*(parts(1) >= 50);
    start = datetime(year, 1, 1);
    assert(parts(2) >= 1 && parts(2) <= days(datetime(year+1,1,1)-start) && parts(3) < 86400, ...
        'ULMO:gfzneq:InvalidEpoch', 'Invalid SINEX epoch %s.', token);
    value = seconds(start + days(parts(2)-1) - datetime(2000,1,1)) + parts(3);
end

function stream = openSource(path)
    stream = java.io.FileInputStream(char(path));
    if endsWith(path, '.gz')
        try
            stream = java.util.zip.GZIPInputStream(stream, 65536);
        catch exception
            stream.close();
            rethrow(exception);
        end
    end
end

function [packed, scale] = monthlyFactor(path, month)
    upper = readNormalTriangle(path, month);
    R = chol(upper, 'upper');
    clear upper
    n = numel(month.idx);
    E = full(sparse(month.idx, 1:n, 1, month.count, n));
    % Q is the selected block of the inverse of the full normal matrix.
    Y = linsolve(R, E, struct('UT', true, 'TRANSA', true));
    clear R E
    Q = Y.' * Y;
    clear Y
    ratios = month.std.^2 ./ diag(Q);
    scale = median(ratios);
    % GFZ prints six significant digits. Two rounded variances can differ
    % by up to 2e-5 in relative terms even when their true scales agree.
    assert(all(abs(ratios / scale - 1) <= 2e-5), 'ULMO:gfzneq:Scale', ...
        'Formal standard deviations disagree with the normal matrix in %s.', path);
    lower = chol(scale * Q, 'lower');
    packed = lower(tril(true(n)));
end

function upper = readNormalTriangle(path, month)
    stream = openSource(path);
    cleanup = onCleanup(@() stream.close());
    stream.readNBytes(month.offset);
    n = month.count;
    upper = nan(n);
    carry = '';
    count = 0;
    while true
        bytes = stream.readNBytes(2^21);
        assert(~isempty(bytes), 'ULMO:gfzneq:InvalidSource', 'Incomplete normal matrix in %s.', path);
        block = [carry char(bytes).'];
        tail = regexp(block, '(?m)^-SOLUTION/NORMAL_EQUATION_MATRIX L[ \t]*\r?\n', 'once');
        if isempty(tail)
            last = find(block == newline, 1, 'last');
            carry = block(last+1:end);
        else
            last = tail-1;
        end
        data = textscan(block(1:last), '%f%f%f%f%f', 'Delimiter', ' ', ...
            'MultipleDelimsAsOne', true, 'CollectOutput', true, ...
            'CommentStyle', '*', 'ReturnOnError', false);
        rows = data{1};
        row = rows(:, 1);
        col = rows(:, 2);
        assert(all(isfinite(rows(:,3))) && all(row <= n & row == fix(row) & ...
            col >= 1 & col == fix(col)), 'ULMO:gfzneq:InvalidSource', 'Invalid matrix record in %s.', path);
        for offset = 0:2
            values = rows(:, 3+offset);
            used = ~isnan(values);
            assert(all(isfinite(values(used))) && all(col(used)+offset <= row(used)), ...
                'ULMO:gfzneq:InvalidSource', 'Invalid lower triangle in %s.', path);
            upper(col(used)+offset + n*(row(used)-1)) = values(used);
            count = count + nnz(used);
        end
        if ~isempty(tail)
            % Finish the stream to verify the SINEX terminator and gzip CRC.
            footer = [block(tail:end) char(stream.readAllBytes()).'];
            assert(endsWith(strtrim(footer), '%ENDSNX'), 'ULMO:gfzneq:InvalidSource', ...
                'Missing SINEX terminator in %s.', path);
            break
        end
    end
    assert(count == n*(n+1)/2, 'ULMO:gfzneq:InvalidSource', 'Wrong matrix entry count in %s.', path);
    for col = 1:n
        assert(all(isfinite(upper(1:col,col))), 'ULMO:gfzneq:InvalidSource', ...
            'Missing normal matrix entry in %s.', path);
    end
end

function createSeries(path, months, L)
    labels = months(1).labels;
    n = size(labels, 1);
    nPacked = n*(n+1)/2;
    timeDim = {'time', numel(months)};
    paramDim = {'parameter', n};
    epochs = [[months.time]; [months.start]; [months.end]];
    names = ["time", "epochStart", "epochEnd"];
    for col = 1:3
        nccreate(path, names(col), 'Dimensions', timeDim, 'Format', 'netcdf4');
        ncwrite(path, names(col), epochs(col,:).');
        ncwriteatt(path, names(col), 'units', 'seconds since 2000-01-01 00:00:00 UTC');
        ncwriteatt(path, names(col), 'calendar', 'proleptic_gregorian');
    end
    nccreate(path, 'noiseFactor', 'Dimensions', [{'packed', nPacked} timeDim], ...
        'Format', 'netcdf4', 'Datatype', 'double', 'ChunkSize', [min(nPacked, 131072), 1]);
    ncwriteatt(path, 'noiseFactor', 'packing', 'Lower triangle by columns, diagonal to last row.');
    ncwriteatt(path, 'noiseFactor', 'units', '1');
    nccreate(path, 'formalStd', 'Dimensions', [paramDim timeDim], 'Datatype', 'double');
    ncwriteatt(path, 'formalStd', 'units', '1');
    names = ["degree", "order", "cs"];
    types = ["uint16", "uint16", "uint8"];
    for col = 1:3
        nccreate(path, names(col), 'Dimensions', paramDim, 'Datatype', types(col));
        ncwrite(path, names(col), labels(:, col));
    end
    ncwriteatt(path, 'cs', 'description', '0=cosine, 1=sine');
    nccreate(path, 'sourceFilename', 'Dimensions', timeDim, 'Datatype', 'string');
    ncwrite(path, 'sourceFilename', string({months.name}).');
    nccreate(path, 'sourceMaxDegree', 'Dimensions', timeDim, 'Datatype', 'uint16');
    ncwrite(path, 'sourceMaxDegree', [months.maxDegree].');
    nccreate(path, 'referenceRadius', 'Dimensions', timeDim);
    ncwrite(path, 'referenceRadius', [months.radius].');
    ncwriteatt(path, 'referenceRadius', 'units', 'm');
    nccreate(path, 'varianceFactor', 'Dimensions', timeDim);
    nccreate(path, 'creationDate', 'Datatype', 'double');
    ncwriteatt(path, 'creationDate', 'units', 'seconds since 1970-01-01 00:00:00 UTC');
    ncwriteatt(path, '/', 'maxDegree', L);
    ncwriteatt(path, '/', 'normalisation', months(1).norm);
    ncwriteatt(path, '/', 'tideSystem', months(1).tide);
end
