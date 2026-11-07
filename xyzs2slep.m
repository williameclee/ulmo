%% XYZS2SLEP - Fits scalar fields on paired coordinates to a regional Slepian basis.
%
% Syntax
%   [C, V, N, fit, design] = xyzs2slep(vals, lon, lat, domain, L, ...)
%   [C, V, N, fit] = xyzs2slep(vals, design, ...)
%
% Input arguments
%   vals - Real nPoint-by-nField matrix in field units. 
%       A vector matching the coordinate count is one field. Each column uses the same grid.
%   lon, lat - Paired coordinate vectors in degrees. 
%       Longitude wraps at 360.
%   domain - GeoDomain, region name/cell, closed [lon,lat] polygon, or
%       scalar north-polar cap radius in degrees. 
%       All supplied points are used.
%   L - Nonnegative integer maximum harmonic degree.
%   design - Fifth output of an earlier call on the same ordered grid.
%       Reuses the unweighted matrix. Supply new weights/errors on each call.
%   GVN (name-value) - Optional {G,V,N} for the requested domain and L.
%       G has (L+1)^2 rows in ADDMOUT order. Columns retain their supplied
%       order and signs. Supply an already-retained basis to fit exactly it.
%   truncation (name-value) - Number of leading basis columns. Defaults to
%       max(1,round(N)), bounded by basis size, for a computed or supplied
%       full basis. For a supplied subset, defaults to all its columns.
%       A reused design already fixes the truncation.
%   missingPolicy (name-value) - "error" (default) or "omit". Omit NaN
%       values per field and NaN coordinates across fields. Inf is invalid.
%   weight (name-value) - Finite nonnegative vector, one per grid point,
%       shared across fields. Minimize sum(weight .* residual.^2). Zero
%       excludes a point. Uniform centered lat/lon grids can use cosd(lat).
%       Default is equal weighting. Exclusive with all uncertainty options.
%   dataStd, dataCovariance (name-value) - Optional absolute observation
%       standard deviations (positive vector) or positive-definite covariance
%       (nPoint-by-nPoint). Shared across fields and finite on the full grid.
%   estimateNoiseVariance (name-value) - Estimate iid residual variance
%       separately per field. Requires more observations than coefficients.
%       Default false. Select at most one uncertainty source.
%   uncertaintyOutput (name-value) - "none" (default), "std", or "covariance".
%       With a noise model, the last choice returns both covariance and std.
%   rankTolerance (name-value) - Relative SVD cutoff in [0,1).
%       Default max(nObs,J)*eps per mask. Rank deficiency raises an error.
%   blockSize (name-value) - Maximum points per basis evaluation block.
%       Default 1024. Only rows used by at least one field are sampled
%       unless the fifth output requests a reusable design.
%
% Output arguments
%   C - J-by-nField coefficients in field units, using 4*pi real harmonics.
%   V, N - Retained concentration eigenvalues and full Shannon number.
%   fit - nObs, rank, condNum, residualRms (one entry per field). RMS is
%       in field units, area-weighted when weight is supplied and ordinary
%       otherwise. uncSrc identifies the shared error model. coeffStds
%       (J-by-nField) and coeffCov (J-by-J-by-nField) are empty unless
%       requested with an error model. Covariances are marginal per field.
%   design - A, lon, lat, and GVN. A retains original row order, including
%       points missing from any field. Undefined coordinate rows are NaN.
%
% Notes
%   Fields with identical usable rows share one factorisation. A field with
%   insufficient coverage fails the call with its column index. Propagated
%   errors are conditional on the chosen basis and observation error model.
%
% See also
%   GLMALPHA_NEW, PLM2SLEP_NEW, SLEP2PLM_NEW
%
% Created by
%   2026/10/06, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/07, En-Chi Lee (williameclee@arizona.edu)

function [C, V, N, fit, design] = xyzs2slep(vals, grid, lat, domain, L, options)

    arguments (Input)
        vals (:, :) {mustBeNumeric, mustBeReal}
        grid
        lat = []
        domain = []
        L = []
        options.GVN cell = {}
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.rankTolerance double ...
            {mustBeScalarOrEmpty, mustBeFinite, mustBeNonnegative, mustBeLessThan(options.rankTolerance, 1)} = []
        options.missingPolicy (1, 1) string {mustBeMember(options.missingPolicy, ["error", "omit"])} = "error"
        options.blockSize (1, 1) double {mustBeInteger, mustBeFinite, mustBePositive} = 1024
        options.weight {mustBeNumeric, mustBeReal} = []
        options.dataStd {mustBeNumeric, mustBeReal} = []
        options.dataCovariance {mustBeNumeric, mustBeReal} = []
        options.estimateNoiseVariance (1, 1) logical = false
        options.uncertaintyOutput (1, 1) string ...
            {mustBeMember(options.uncertaintyOutput, ["none", "std", "covariance"])} = "none"
    end

    reuseDesign = isstruct(grid);

    if reuseDesign
        design = grid;

        if ~isempty(lat) || ~isempty(domain) || ~isempty(L)
            error('ULMO:xyzs2slep:InvalidDesign', 'Do not supply geometry when reusing a design.');
        end

        if ~isscalar(design) || ~all(isfield(design, {'A', 'GVN', 'lon', 'lat'})) || ...
                ~isempty(options.GVN) || ~isempty(options.truncation)
            error('ULMO:xyzs2slep:InvalidDesign', 'Supply a returned design without GVN or truncation.');
        end

        [lon, lat] = coordinates(design.lon, design.lat);
    else

        if isempty(L)
            error('ULMO:xyzs2slep:MissingGeometry', 'Supply longitude, latitude, domain, and L.');
        end

        design = [];
        [lon, lat] = coordinates(grid, lat);
    end

    nPts = numel(lon);

    if isvector(vals) && numel(vals) == nPts
        vals = vals(:);
    end

    if isempty(vals) || size(vals, 1) ~= nPts
        error('ULMO:xyzs2slep:ObservationSize', 'Supply one row of values per coordinate pair.');
    end

    if any(isinf(vals), 'all')
        error('ULMO:xyzs2slep:InvalidObservation', 'Infinite observations are invalid.');
    end

    wght = double(options.weight(:));

    if ~isempty(wght) && (~isvector(options.weight) || numel(wght) ~= nPts || ...
            any(~isfinite(wght) | wght < 0))
        error('ULMO:xyzs2slep:InvalidWeight', 'Supply one finite nonnegative weight per point.');
    end

    if sum([~isempty(wght), ~isempty(options.dataStd), ...
                ~isempty(options.dataCovariance), options.estimateNoiseVariance]) > 1
        error('ULMO:xyzs2slep:ConflictingUncertainty', ...
        'Select only one weight or uncertainty model.');
    end

    valid = ~isnan(vals) & ~isnan(lon) & ~isnan(lat);

    if ~isempty(wght)
        actv = wght > 0;
    else
        actv = true(nPts, 1);
    end

    if strcmpi(options.missingPolicy, "error") && any(~valid & actv, 'all')
        error('ULMO:xyzs2slep:MissingObservation', ...
        'NaN values or coordinates require missingPolicy="omit".');
    end

    valid = valid & actv;
    nObs = sum(valid, 1);

    if ~isempty(find(nObs == 0, 1))
        emptyField = find(nObs == 0, 1);
        error('ULMO:xyzs2slep:ObservationSize', ...
            'Field %d has no usable observations.', emptyField);
    end

    sampleRows = true(nPts, 1);

    if reuseDesign
        validateattributes(design.A, {'double'}, {'real', '2d', 'nrows', nPts, 'nonempty'});
        validateBasis(design.GVN);

        if size(design.A, 2) ~= numel(design.GVN{2})
            error('ULMO:xyzs2slep:InvalidDesign', ...
            'Design columns must match the retained basis.');
        end

        basis = design.GVN;
        sampled = design.A;
    else
        basis = selectBasis(domain, L, options);
    end

    [G, V, N] = basis{:};
    J = numel(V);
    field = find(nObs < J, 1);

    if ~isempty(field)
        error('ULMO:xyzs2slep:InvalidTruncation', ...
            'Field %d has %d usable points for %d coefficients.', field, nObs(field), J);
    end

    if options.estimateNoiseVariance
        field = find(nObs == J, 1);

        if ~isempty(field)
            error('ULMO:xyzs2slep:NoiseDegreesOfFreedom', ...
                'Field %d needs positive residual degrees of freedom.', field);
        end

    end

    [sigma, cov, fullFactor, source] = prepareUncertainty(options.dataStd, ...
        options.dataCovariance, options.estimateNoiseVariance, nPts);

    if ~reuseDesign

        if nargout < 5
            sampleRows = any(valid, 2);
        end

        sampled = sampleBasis(lon(sampleRows), lat(sampleRows), G, L, options.blockSize);

        if nargout >= 5
            design = struct('A', sampled, 'lon', lon, 'lat', lat, 'GVN', {basis});
        end

    end

    nField = size(vals, 2);
    C = zeros(J, nField);
    fit = struct('nObs', nObs, 'rank', zeros(1, nField), 'condNum', zeros(1, nField), ...
        'residualRms', zeros(1, nField), 'uncSrc', source, 'coeffStds', [], 'coeffCov', []);
    wantStd = nargout >= 4 && options.uncertaintyOutput ~= "none" && source ~= "none";
    wantCov = wantStd && options.uncertaintyOutput == "covariance";

    if wantStd
        fit.coeffStds = zeros(J, nField);
    end

    if wantCov
        fit.coeffCov = zeros(J, J, nField);
    end

    [masks, ~, group] = unique(valid.', 'rows');

    for iMask = 1:size(masks, 1)
        rows = masks(iMask, :).';
        fields = find(group == iMask);
        count = nObs(fields(1));
        % rows indexes the original observations. Only A may be compact.
        A = sampled(rows(sampleRows), :);
        data = double(vals(rows, fields));
        B = A;
        rhs = data;

        if ~isempty(wght)
            rowScale = sqrt(wght(rows));
            B = A .* rowScale;
            rhs = data .* rowScale;
        elseif ~isempty(sigma)
            B = A ./ sigma(rows);
            rhs = data ./ sigma(rows);
        elseif ~isempty(fullFactor)
            factor = fullFactor;

            if count ~= nPts
                factor = chol(cov(rows, rows), 'lower');
            end

            B = factor \ A;
            rhs = factor \ data;
        end

        if any(~isfinite(B), 'all') || any(~isfinite(rhs), 'all')
            error('ULMO:xyzs2slep:WhiteningOverflow', 'Nonfinite scaled values in field %d.', fields(1));
        end

        [C(:, fields), rnk, condition, invFactor] = solveMask(B, rhs, ...
            fields(1), options.rankTolerance, wantStd);

        if nargout < 4
            continue
        end

        residual = data - A * C(:, fields);

        if isempty(wght)
            rms = vecnorm(residual, 2, 1) / sqrt(count);
        else
            rms = vecnorm(rowScale .* residual, 2, 1) / norm(rowScale);
        end

        fit.residualRms(fields) = rms;
        fit.rank(fields) = rnk;
        fit.condNum(fields) = condition;

        if wantStd
            noiseStd = ones(1, numel(fields));

            if options.estimateNoiseVariance
                noiseStd = vecnorm(residual, 2, 1) / sqrt(count - J);
            end

            fit.coeffStds(:, fields) = vecnorm(invFactor, 2, 2) .* noiseStd;

            if wantCov
                fit.coeffCov(:, :, fields) = (invFactor * invFactor.') .* reshape(noiseStd .^ 2, 1, 1, []);
            end

            if any(~isfinite(fit.coeffStds(:, fields)), 'all') || ...
                    (wantCov && any(~isfinite(fit.coeffCov(:, :, fields)), 'all'))
                error('ULMO:xyzs2slep:CovarianceOverflow', 'Coefficient uncertainty overflowed.');
            end

        end

    end

end

function [coeffs, rnk, cond, invFactor] = solveMask(B, rhs, field, tol, wantStd)
    [count, J] = size(B);
    [U, S, Q] = svd(B, 'econ');
    s = diag(S);

    if isempty(tol)
        tol = max(count, J) * eps;
    end

    rnk = sum(s > tol * s(1));

    if rnk < J
        error('ULMO:xyzs2slep:RankDeficient', ...
            'Field %d has sampled basis rank %d of %d.', field, rnk, J);
    end

    cond = s(1) / s(end);

    if cond > 1 / sqrt(eps)
        warning('ULMO:xyzs2slep:IllConditioned', ...
            'Sampled basis condition number %.3g for field %d.', cond, field);
    end

    coeffs = Q * ((U.' * rhs) ./ s);
    invFactor = [];

    if wantStd
        invFactor = Q ./ s.';
    end

end

function [lon, lat] = coordinates(lon, lat)
    validateattributes(lon, {'numeric'}, {'real', 'vector'});
    validateattributes(lat, {'numeric'}, {'real', 'vector'});

    if numel(lat) ~= numel(lon)
        error('ULMO:xyzs2slep:ObservationSize', 'Supply paired longitude and latitude vectors.');
    end

    if any(isinf(lon)) || any(abs(lat) > 90)
        error('ULMO:xyzs2slep:InvalidObservation', 'Invalid longitude or latitude.');
    end

    lon = mod(double(lon(:)), 360);
    lat = double(lat(:));
end

function validateBasis(GVN)

    if ~iscell(GVN) || numel(GVN) ~= 3
        error('ULMO:xyzs2slep:InvalidBasis', 'Supply {G, V, N}.');
    end

    [G, V, N] = GVN{:};
    validateattributes(G, {'double', 'single'}, {'real', 'finite', '2d', 'nonempty'});
    validateattributes(V, {'numeric'}, {'real', 'finite', 'vector', 'numel', size(G, 2)});
    validateattributes(N, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});

    if mod(sqrt(size(G, 1)), 1) ~= 0 || size(G, 2) > size(G, 1)
        error('ULMO:xyzs2slep:InvalidBasis', 'G must have (L+1)^2 rows and at most that many columns.');
    end

end

function basis = selectBasis(domain, L, options)
    validateattributes(L, {'numeric'}, {'real', 'scalar', 'integer', 'nonnegative', 'finite'});
    J = options.truncation;

    if ~isempty(options.GVN)
        validateBasis(options.GVN);
        [G, V, N] = options.GVN{:};

        if size(G, 1) ~= (L + 1) ^ 2
            error('ULMO:xyzs2slep:InvalidBasis', 'Supplied basis degree differs from L.');
        end

        order = 1:size(G, 2);
    else

        if isnumeric(domain) && isscalar(domain)
            validateattributes(domain, {'numeric'}, {'real', 'finite', '>', 0, '<=', 180});
            [G, V, ~, ~, N] = glmalpha(domain, L);
        elseif (isa(domain, 'GeoDomain') && isscalar(domain)) || ischar(domain) || ...
                (isstring(domain) && isscalar(domain)) || iscell(domain) || ...
                (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2)
            [G, V, ~, ~, N] = glmalpha_new(domain, L, 'BeQuiet', true);
        else
            error('ULMO:xyzs2slep:InvalidDomain', 'Supply a region, polygon or cap radius.');
        end

        [V, order] = sort(V(:), 'descend');
    end

    if isempty(J)
        J = size(G, 2);

        if J == (L + 1) ^ 2
            J = min(J, max(1, round(N)));
        end

    end

    if J > size(G, 2)
        error('ULMO:xyzs2slep:InvalidTruncation', 'Truncation exceeds the available basis columns.');
    end

    G = double(G(:, order(1:J)));
    V = V(:);
    basis = {G, V(1:J), N};
end

function A = sampleBasis(lon, lat, G, L, blockSize)
    A = nan(numel(lon), size(G, 2));
    points = find(~isnan(lon) & ~isnan(lat));
    orders = addmout(L);
    scale = sqrt(4 * pi) * (-1) .^ orders(:);

    for first = 1:blockSize:numel(points)
        rows = points(first:min(first + blockSize - 1, end));

        if L == 0
            A(rows, :) = G;
        else
            Y = ylm([0 L], [], deg2rad(90 - lat(rows)), deg2rad(lon(rows)), [], [], 0, 1);
            A(rows, :) = (Y .* scale).' * G;
        end

    end

end

function [sigma, cov, fac, src] = prepareUncertainty(sigma, cov, estimate, cnt)
    fac = [];
    src = "none";

    if ~isempty(sigma)
        src = "dataStd";

        if ~isvector(sigma) || numel(sigma) ~= cnt
            error('ULMO:xyzs2slep:InvalidDataStd', ...
                ['Size of dataStd must match the number of input values. ', ...
             'But got %d values and %d STDs.'], ...
                cnt, numel(sigma));
        elseif any(~isfinite(sigma(:))) || any(sigma(:) <= 0)
            error('ULMO:xyzs2slep:InvalidDataStd', ...
            'dataStd must be finite and positive.');
        end

        sigma = full(double(sigma(:)));
    elseif ~isempty(cov)
        src = "dataCovariance";

        if ~ismatrix(cov) || ~isequal(size(cov), [cnt cnt])
            error('ULMO:xyzs2slep:InvalidDataCovariance', ...
            'dataCovariance must be a square matrix for the original observations.');
        elseif any(~isfinite(cov(:)))
            error('ULMO:xyzs2slep:InvalidDataCovariance', ...
            'dataCovariance must be finite.');
        end

        cov = full(double(cov));
        maxCov = max(abs(cov(:)));

        if maxCov == 0
            error('ULMO:xyzs2slep:NonPositiveCovariance', ...
            'dataCovariance must be positive definite.');
        end

        % Make sure the covariance matrix is symmetric
        normedCov = cov / maxCov;

        if norm(normedCov - normedCov.', 'fro') > 100 * eps * norm(normedCov, 'fro')
            error('ULMO:xyzs2slep:AsymmetricCovariance', ...
            'dataCovariance must be symmetric within 100*eps relative Frobenius tolerance.');
        end

        cov = cov / 2 + cov.' / 2;
        [fac, flag] = chol(cov, 'lower');

        if flag ~= 0
            error('ULMO:xyzs2slep:NonPositiveCovariance', ...
            'dataCovariance must be positive definite; no jitter is added.');
        end

    elseif estimate
        src = "estimatedIid";
    end

end
