%% PREPARELOWDEGREEMAP - Builds the linear map from known to unknown loads.
%
% Syntax
%   map = preparelowdegreemap(kernels, unknownIds)
%   map = preparelowdegreemap(kernels, unknownIds, Method="uniform")
%
% Input arguments
%   kernels - Output of PREPARESLEKERNELS. For Method="uniform", only O
%       is required; no Earth response or SLE solve is used.
%   unknownIds - Distinct kernel-order indices of the estimated coefficients.
%       Size: (1, U). Order is retained in the output.
%   Method (name-value) - Ocean constraint, fingerprint or uniform.
%       Default method: fingerprint.
%   RelativeTolerance (name-value) - Relative residual tolerance for SLE solves.
%       Default tolerance: 1e-12.
%   MaxIter (name-value) - Maximum iterations in each SLE system solve.
%       Default limit: 100.
%
%   CacheFolder (name-value) - Directory for numeric response-map caches.
%       Default folder: "" (disabled). Identity hashes the actual matrices,
%       coefficient selection, method, tolerance, iteration cap and version.
%   SaveData (name-value) - Whether to save a newly prepared map.
%       Default option: true. Existing matching maps can still be loaded.
%   ForceNew (name-value) - Whether to bypass a matching map cache.
%       Default option: false.
%
% Output arguments
%   map - Numeric response matrices and index metadata.
%       D: selected ocean response to a full load, size (U,M).
%       Dunknown, Dknown: the small-system blocks, sizes (U,U), (U,M-U).
%       T: maps known loads to estimated loads, size (U,M-U).
%       unknownIds, knownIds: row and column index order for T.
%       reciprocalCondition: rcond(Dunknown); diagnostics: SLE solve accuracy.
%       cacheHit: whether a stored map was reused; cacheKey: numeric identity.
%
% Notes
%   For fingerprint loading, D=(O(unknownIds,:)/A)*KF*(I-O).
%   Only U transposed systems are solved. All epochs then use T*knownLoads.
%   Independent input variances propagate as (T.^2)*knownVariances.
%   Reuse this map only with the same ocean matrix, Earth response, bandwidth,
%   coefficient selection and numerical accuracy. It is independent of epochs.
%   Finite-bandwidth O is not assumed idempotent.
%
% See also
%   PREPARESLEKERNELS, SOLVESLESYSTEM, SOLVEDEGREE1
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function map = preparelowdegreemap(kernels, unknownIds, options)
    arguments (Input)
        kernels (1, 1) struct
        unknownIds (1, :) double {mustBeInteger, mustBePositive}
        options.Method (1, 1) string {mustBeMember(options.Method, ["fingerprint", "uniform"])} = "fingerprint"
        options.RelativeTolerance (1, 1) double {mustBePositive, mustBeFinite, mustBeLessThan(options.RelativeTolerance, 1)} = 1e-12
        options.MaxIter (1, 1) double {mustBeInteger, mustBePositive} = 100
        options.CacheFolder (1, 1) string = ""
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
    end
    O = kernels.O;
    M = size(O, 1);
    if isempty(unknownIds) || numel(unique(unknownIds)) ~= numel(unknownIds) || any(unknownIds > M)
        error('ULMO:preparelowdegreemap:InvalidIndices', 'Unknown indices must be distinct and within the bandwidth.');
    end
    if size(O, 2) ~= M || ~isfinite(O(1, 1)) || O(1, 1) <= 0
        error('ULMO:preparelowdegreemap:InvalidOceanKernel', 'O must be square with positive finite area.');
    end
    knownIds = setdiff(1:M, unknownIds);
    cacheKey = "";
    cacheFile = "";
    if options.CacheFolder ~= "" && (options.SaveData || isfolder(options.CacheFolder))
        cacheKey = mapidentity(kernels, unknownIds, options);
        cacheFile = fullfile(options.CacheFolder, "sle-map-v1-" + cacheKey + ".mat");
        if ~options.ForceNew && isfile(cacheFile)
            saved = load(cacheFile, 'map');
            if isfield(saved, 'map') && isfield(saved.map, 'cacheKey') && ...
                    saved.map.cacheKey == cacheKey && isequal(saved.map.unknownIds, unknownIds) && ...
                    isequal(saved.map.knownIds, knownIds) && ...
                    isequal(size(saved.map.T), [numel(unknownIds), numel(knownIds)]) && ...
                    all(isfinite(saved.map.T), 'all')
                map = saved.map;
                map.cacheHit = true;
                return
            end
        end
    end
    Orows = O(unknownIds, :);
    diagnostics = struct;
    if options.Method == "uniform"
        landC00 = -O(1, :);
        landC00(1) = landC00(1) + 1;
        D = -(O(unknownIds, 1) / O(1, 1)) * landC00;
    else
        [Ht, diagnostics] = solveslesystem(kernels, Orows.', Transpose=true, ...
            RelativeTolerance=options.RelativeTolerance, MaxIter=options.MaxIter);
        H = Ht.' * kernels.KF;
        D = H - H * O;
    end
    Dunknown = D(:, unknownIds) - Orows(:, unknownIds);
    Dknown = Orows(:, knownIds) - D(:, knownIds);
    reciprocalCondition = rcond(Dunknown);
    if ~isfinite(reciprocalCondition) || reciprocalCondition <= eps
        error('ULMO:preparelowdegreemap:SingularSystem', 'The selected ocean constraint is singular to working precision.');
    elseif reciprocalCondition < 1e-10
        warning('ULMO:preparelowdegreemap:IllConditionedSystem', ...
            'The selected ocean constraint has reciprocal condition %.6g.', reciprocalCondition);
    end
    T = Dunknown \ Dknown;
    if any(~isfinite(T), 'all')
        error('ULMO:preparelowdegreemap:NonFiniteMap', 'The response map is nonfinite.');
    end
    map = struct('D', D, 'Dunknown', Dunknown, 'Dknown', Dknown, 'T', T, ...
        'unknownIds', unknownIds, 'knownIds', knownIds, ...
        'reciprocalCondition', reciprocalCondition, 'diagnostics', diagnostics, ...
        'cacheHit', false, 'cacheKey', cacheKey);
    if cacheFile ~= "" && options.SaveData
        if ~isfolder(options.CacheFolder)
            mkdir(options.CacheFolder);
        end
        temporary = string(tempname(options.CacheFolder)) + ".mat";
        cleanup = onCleanup(@() removetemporary(temporary));
        save(temporary, 'map');
        movefile(temporary, cacheFile, 'f');
    end
end

function key = mapidentity(kernels, unknownIds, options)
    digest = java.security.MessageDigest.getInstance('SHA-256');
    metadata = struct('Version', 1, 'UnknownIds', unknownIds, ...
        'Method', options.Method, 'RelativeTolerance', options.RelativeTolerance, ...
        'MaxIter', options.MaxIter, 'Size', size(kernels.O));
    digest.update(typecast(uint8(jsonencode(metadata)), 'int8'));
    % Stream the large matrix, keeping temporary byte arrays bounded.
    for first = 1:64:size(kernels.O, 2)
        block = full(kernels.O(:, first:min(first+63, size(kernels.O, 2))));
        digest.update(typecast(block(:), 'int8'));
    end
    if options.Method == "fingerprint"
        for field = {'response', 'KF'}
            matrix = kernels.(field{1});
            header = struct('Field', field{1}, 'Size', size(matrix), 'Count', nnz(matrix));
            digest.update(typecast(uint8(jsonencode(header)), 'int8'));
            [i,j,v] = find(matrix);
            values = [i,j,v];
            digest.update(typecast(values(:), 'int8'));
        end
    end
    bytes = typecast(digest.digest(), 'uint8');
    key = string(lower(reshape(dec2hex(bytes, 2).', 1, [])));
end

function removetemporary(file)
    if isfile(file)
        delete(file);
    end
end
