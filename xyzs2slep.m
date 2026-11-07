%% XYZS2SLEP - Fits scalar fields on paired coordinates to a regional Slepian basis.
% Each field is fitted by least squares on its usable observations. Fields
% share coordinates, a retained basis, and any weights or observation errors.
%
% Syntax
%   [coeffs, V, N, fit, design] = xyzs2slep(vals, lon, lat, domain, L)
%   [coeffs, V, N, fit] = xyzs2slep(vals, design = d)
%   [__] = xyzs2slep(__, "Name", Value)
%
% Input arguments
%   vals - Observations to fit, with one field per column.
%       Rows follow the supplied coordinate order. A vector matching the
%       coordinate count is treated as one field, in either orientation.
%       Size: (nPts, nField), or a vector of length nPts for one field.
%       Unit: field units.
%   lon (optional) - Paired longitudes, wrapped to [0, 360).
%       Supply lon, lat, domain, and L together to build a design.
%       Size: (nPts, 1) or (1, nPts).
%       Unit: degrees east.
%       Default value: [].
%   lat (optional) - Paired latitudes in [-90, 90].
%       Required with lon; omit when supplying design.
%       Size: (nPts, 1) or (1, nPts).
%       Unit: degrees north.
%       Default value: [].
%   domain (optional) - Region used to construct the Slepian basis.
%       Accepts a scalar GeoDomain, region name/cell, closed [lon, lat]
%       polygon, or scalar north-polar cap radius in (0, 180]. Required
%       with lon; omit when supplying design.
%       Unit: degrees for polygon coordinates and cap radius.
%       Default value: [].
%   L (optional) - Nonnegative integer maximum harmonic degree.
%       Required with lon; omit when supplying design.
%       Size: scalar.
%       Unit: dimensionless.
%       Default value: [].
%   design (name-value) - Fifth output of an earlier call on the same grid.
%       Pass design = d to reuse its coordinates and retained basis. Omit
%       positional geometry, GVN, and truncation. Supply weights or errors
%       for the current observations. See the output fields below.
%       Size: scalar structure.
%       Default design: omitted (evaluated from lon, lat, domain, and L).
%   GVN (name-value) - Supplied {G, V, N} for the requested domain and L.
%       G contains real harmonic coefficients in ADDMOUT order. V contains
%       their concentration eigenvalues; N is the full Shannon number.
%       Supplied columns retain their order and signs. A supplied basis
%       bypasses domain construction; its geographical origin is not checked.
%       Size: three-element cell; G is ((L+1)^2, nBasis), V has nBasis
%       elements, and N is scalar.
%       Unit: dimensionless.
%       Default basis: {}, computed from domain and L.
%   truncation (name-value) - Positive integer number J of leading modes.
%       Cannot exceed the available basis columns. A reused design fixes J.
%       Size: scalar or empty.
%       Default count: [], using min(nBasis, max(1, round(N))) for a full
%       basis with (L+1)^2 columns, or all columns of a supplied subset.
%   missingPolicy (name-value) - Treatment of NaN observations/coordinates.
%       "error" rejects NaNs at points with nonzero weight. "omit" removes
%       NaN observations per field and NaN coordinates across all fields.
%       Infinite observations or coordinates are always invalid.
%       Size: text scalar.
%       Default policy: "error".
%   weight (name-value) - Finite nonnegative weights shared across fields.
%       Minimises sum(weight .* residual.^2); zero excludes a point.
%       Uniform centred latitude-longitude grids can use cosd(lat) for
%       approximate area weighting. Cannot be combined with an error model.
%       Size: (nPts, 1), (1, nPts), or empty.
%       Unit: relative weights; a common positive scale does not affect C.
%       Default weights: [], equal weighting.
%   dataStd (name-value) - Absolute observation standard deviations.
%       Must be finite and positive on the full grid, including omitted
%       points. Shared across fields; assumes uncorrelated point errors.
%       Size: (nPts, 1), (1, nPts), or empty.
%       Unit: field units.
%       Default value: [], no supplied standard deviations.
%   dataCovariance (name-value) - Absolute observation error covariance.
%       Must be finite, symmetric, and positive definite on the full grid.
%       Shared across fields; rows and columns follow the coordinate order.
%       Size: (nPts, nPts) or empty.
%       Unit: squared field units.
%       Default value: [], no supplied covariance.
%   estimateNoiseVariance (name-value) - Whether to estimate residual variance.
%       Assumes independent, identically distributed errors within each
%       field. Uses sum(residual.^2)/(nObs-J) separately for each field;
%       each field must have more usable observations than coefficients.
%       Size: logical scalar.
%       Default value: false.
%   uncertaintyOutput (name-value) - Coefficient uncertainty to return in fit.
%       "none" returns neither array; "std" returns standard deviations;
%       "covariance" returns both standard deviations and covariance.
%       Requires an error model and the fourth output to produce arrays.
%       Size: text scalar.
%       Default output: "none".
%   rankTolerance (name-value) - Relative singular-value cutoff in [0, 1).
%       Singular values must exceed this fraction of the largest value.
%       Applied after weighting or whitening; rank deficiency is an error.
%       Size: scalar or empty.
%       Unit: dimensionless.
%       Default tolerance: [], using max(nObs, J)*eps for each usable mask.
%   blockSize (name-value) - Maximum points per basis evaluation block.
%       Must be a finite positive integer. Limits harmonic evaluation
%       storage; the sampled Slepian matrix is still retained for the fit.
%       Size: scalar.
%       Default count: 1024.
%
% Output arguments
%   coeffs - Coefficients in the retained Slepian basis.
%       Uses real spherical harmonics with squared norm 4*pi on the sphere.
%       Size: (J, nField).
%       Unit: field units.
%   V - Concentration eigenvalues of the retained modes.
%       Computed modes are sorted in descending concentration. Supplied
%       modes retain the order in GVN.
%       Size: (J, 1).
%       Unit: dimensionless.
%   N - Full Shannon number, before truncation.
%       Size: scalar.
%       Unit: dimensionless.
%   fit - Diagnostics and optional coefficient uncertainties.
%       nObs, rank, and condNum: (1, nField) usable counts, numerical ranks,
%           and 2-norm condition numbers of the weighted/whitened matrices.
%       residualRms: (1, nField), in field units. Uses the supplied weight,
%           or ordinary RMS even when observation uncertainties are supplied.
%       uncSrc: "none", "dataStd", "dataCovariance", or "estimatedIid".
%       coeffStds: (J, nField), in field units, or [] when not requested.
%       coeffCov: (J, J, nField), in squared field units, or [] when not
%           requested. Each slice is a marginal covariance for one field.
%       Both uncertainty arrays are [] when no error model is selected.
%       Size: scalar structure.
%   design - Unweighted sampled basis for reuse on the same ordered grid.
%       A: (nPts, J) dimensionless matrix, including points unused by this
%           fit. Rows with NaN coordinates are NaN.
%       lon, lat: (nPts, 1) coordinate vectors in degrees.
%       GVN: {G, V, N} for the retained basis.
%       Size: scalar structure.
%
% Notes
%   Select at most one of weight, dataStd, dataCovariance, or
%   estimateNoiseVariance. Supplied uncertainties are absolute and are not
%   rescaled by residual variance. Errors are conditional on the retained
%   basis and the selected observation error model; cross-field covariance
%   is not computed.
%   Fields with identical usable rows share one SVD. Insufficient coverage
%   or rank deficiency fails the call and identifies an affected column.
%   Without the fifth output, only rows used by at least one field are
%   sampled. Requesting design samples every defined coordinate so later
%   calls can use observations or positive weights at previously unused rows.
%
% See also
%   GLMALPHA_NEW, PLM2SLEP_NEW, SLEP2PLM_NEW
%
% Created by
%   2026/10/06, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/07, En-Chi Lee (williameclee@arizona.edu)

function [coeffs, V, N, fit, design] = xyzs2slep(vals, lon, lat, domain, L, options)

    arguments (Input)
        vals (:, :) {mustBeNumeric, mustBeReal}
        lon = []
        lat = []
        domain = []
        L = []
        options.design (1, 1) struct
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

    %% Resolve coordinates or reuse a sampled design
    % With an arguments block, nargin counts only supplied positional inputs.
    reuseDesign = isfield(options, 'design');

    if reuseDesign

        if nargin ~= 1 || ~isempty(options.GVN) || ~isempty(options.truncation)
            error('ULMO:xyzs2slep:InvalidDesign', ...
            'A saved design fixes the coordinates, basis, and truncation.');
        end

        design = options.design;

        if ~all(isfield(design, {'A', 'GVN', 'lon', 'lat'}))
            error('ULMO:xyzs2slep:InvalidDesign', ...
            'Supply the design returned by xyzs2slep. It must contain the fields ''A'', ''GVN'', ''lon'', and ''lat''.');
        end

        lon = design.lon;
        lat = design.lat;
    elseif any([isempty(lon), isempty(lat), isempty(domain), isempty(L)])
        error('ULMO:xyzs2slep:MissingGeometry', ...
        'When a fit design is not provided for reuse, lon, lat, domain, and L must all be provided.');
    end

    [lon, lat] = preparecoords(lon, lat);

    %% Validating observations and construct each field's usable mask
    nPts = numel(lon);

    % Validate values
    if isvector(vals) && numel(vals) == nPts
        vals = vals(:);
    end

    if isempty(vals) || size(vals, 1) ~= nPts
        error('ULMO:xyzs2slep:ObservationSize', 'Supply one row of values per coordinate pair.');
    end

    if any(isinf(vals), 'all')
        error('ULMO:xyzs2slep:InvalidObservation', 'Infinite observations are invalid.');
    end

    % Validating weights
    wghts = double(options.weight(:));

    if ~isempty(wghts) && (~isvector(options.weight) || numel(wghts) ~= nPts || ...
            any(~isfinite(wghts) | wghts < 0))
        error('ULMO:xyzs2slep:InvalidWeight', 'Supply one finite nonnegative weight per point.');
    end

    % Using prescribed weights and defining weights from data uncertainties are mutually exclusive
    if sum([~isempty(wghts), ~isempty(options.dataStd), ~isempty(options.dataCovariance), options.estimateNoiseVariance]) > 1
        error('ULMO:xyzs2slep:ConflictingUncertainty', ...
        'Select only one weight or uncertainty model.');
    end

    % Coordinate validity is shared; missing values can differ by field.
    valid = ~isnan(vals) & ~isnan(lon) & ~isnan(lat);

    if ~isempty(wghts)
        actv = wghts > 0;
    else
        actv = true(nPts, 1);
    end

    % Zero-weight rows are excluded even under the strict missing policy.
    if strcmpi(options.missingPolicy, "error") && any(~valid & actv, 'all')
        error('ULMO:xyzs2slep:MissingObservation', ...
            ['NaN are detected in the values or coordinates, ', ...
             'which are prohibited under the selected missing policy ''error''. ', ...
         'Use missingPolicy=''omit'' instead.']);
    end

    valid = valid & actv;
    nObs = sum(valid, 1);

    if ~isempty(find(nObs == 0, 1))
        emptyField = find(nObs == 0, 1);
        error('ULMO:xyzs2slep:ObservationSize', ...
            'Field number %d has no usable observations.', emptyField);
    end

    %% Selecting the retained basis and check coverage before sampling
    sampleRows = true(nPts, 1);

    if reuseDesign
        validateattributes(design.A, {'double'}, {'real', '2d', 'nrows', nPts, 'nonempty'});
        mustBeGVN(design.GVN);

        if size(design.A, 2) ~= numel(design.GVN{2})
            error('ULMO:xyzs2slep:InvalidDesign', ...
            'Design columns must match the retained basis.');
        end

        GVN = design.GVN;
        sampled = design.A;
    else
        GVN = prepareGVN(domain, L, options);
    end

    [G, V, N] = GVN{:};
    J = numel(V);

    invalidFieldId = find(nObs < J, 1);

    if ~isempty(invalidFieldId)
        error('ULMO:xyzs2slep:InvalidTruncation', ...
            'Field number %d has only %d usable points, which cannot be used to fit %d coefficients.', ...
            invalidFieldId, nObs(invalidFieldId), J);
    end

    % Residual variance needs positive degrees of freedom; supplied absolute
    % errors can also be propagated for an exactly determined fit.
    if options.estimateNoiseVariance
        invalidFieldId = find(nObs == J, 1);

        if ~isempty(invalidFieldId)
            error('ULMO:xyzs2slep:NoiseDegreesOfFreedom', ...
                ['Field number %d needs positive residual degrees of freedom, ', ...
             'but the number of observations is equal to the number of coefficients.'], ...
                invalidFieldId);
        end

    end

    % Validate the full error model before selecting any field's submatrix.
    [sigma, cov, fsqrtCov, source] = prepareUncertainty(options.dataStd, ...
        options.dataCovariance, options.estimateNoiseVariance, nPts);

    %% Sampling the basis and allocate requested outputs
    if ~reuseDesign

        % A reusable design must include rows that later fields may activate.
        % Otherwise, sampling only the union of usable rows avoids extra work.
        if nargout < 5
            sampleRows = any(valid, 2);
        end

        sampled = sampleBasis(lon(sampleRows), lat(sampleRows), G, L, options.blockSize);

        if nargout >= 5
            design = struct('A', sampled, 'lon', lon, 'lat', lat, 'GVN', {GVN});
        end

    end

    nFields = size(vals, 2);
    coeffs = zeros(J, nFields);
    fit = struct('nObs', nObs, 'rank', zeros(1, nFields), 'condNum', zeros(1, nFields), ...
        'residualRms', zeros(1, nFields), 'uncSrc', source, 'coeffStds', [], 'coeffCov', []);
    wantStd = nargout >= 4 && ~strcmpi(options.uncertaintyOutput, "none") && ~strcmpi(source, "none");
    wantCov = wantStd && strcmpi(options.uncertaintyOutput, "covariance");

    if wantStd
        fit.coeffStds = zeros(J, nFields);
    end

    if wantCov
        fit.coeffCov = zeros(J, J, nFields);
    end

    %% Fitting fields with identical masks together
    % Shared weights/errors make one factorisation sufficient for each mask.
    [masks, ~, group] = unique(valid.', 'rows');

    for iMask = 1:size(masks, 1)
        rows = masks(iMask, :).';
        fields = find(group == iMask);
        count = nObs(fields(1));
        % rows refers to the original grid; sampleRows maps it to sampled.
        % Observation arrays retain their original rows even when sampled
        % contains only the union of rows needed by the current fields.
        A = sampled(rows(sampleRows), :);
        data = double(vals(rows, fields));
        B = A;
        rhs = data;

        % Apply the same transformation to A and data. Standard deviations
        % or a Cholesky factor whiten errors to unit covariance; geometric
        % weights instead define the least-squares objective without errors.
        if ~isempty(wghts)
            rowScale = sqrt(wghts(rows));
            B = A .* rowScale;
            rhs = data .* rowScale;
        elseif ~isempty(sigma)
            B = A ./ sigma(rows);
            rhs = data ./ sigma(rows);
        elseif ~isempty(fsqrtCov)
            factor = fsqrtCov;

            % A subset of a Cholesky factor is generally not the factor of
            % the corresponding covariance submatrix, so refactor the subset.
            if count ~= nPts
                factor = chol(cov(rows, rows), 'lower');
            end

            B = factor \ A;
            rhs = factor \ data;
        end

        if any(~isfinite(B), 'all') || any(~isfinite(rhs), 'all')
            error('ULMO:xyzs2slep:WhiteningOverflow', 'Nonfinite scaled values in field %d.', fields(1));
        end

        [coeffs(:, fields), rnk, condition, invFactor] = solveMask(B, rhs, ...
            fields(1), options.rankTolerance, wantStd);

        if nargout < 4
            continue
        end

        % Report residuals in field units, before whitening. Only geometric
        % weights enter the reported RMS; observation errors do not.
        res = data - A * coeffs(:, fields);

        if isempty(wghts)
            rms = vecnorm(res, 2, 1) / sqrt(count);
        else
            rms = vecnorm(rowScale .* res, 2, 1) / norm(rowScale);
        end

        fit.residualRms(fields) = rms;
        fit.rank(fields) = rnk;
        fit.condNum(fields) = condition;

        if wantStd
            % Absolute error scales are already included in invFactor.
            % Only estimated iid noise adds a field-specific residual scale.
            noiseStd = ones(1, numel(fields));

            if options.estimateNoiseVariance
                noiseStd = vecnorm(res, 2, 1) / sqrt(count - J);
            end

            % Row norms give sqrt(diag(invFactor*invFactor.')) without
            % forming a covariance matrix when only standard errors are needed.
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

%% SOLVEMASK - Solves all fields sharing one transformed design matrix.
% Solves B * coeffs = rhs.
%
% Input arguments
%   B - Weighted or whitened design. Size: (nObs, J).
%   rhs - Matching transformed observations. Size: (nObs, nField).
%   field - Original field index used in errors and warnings. Size: scalar.
%   tol - Relative singular-value cutoff, or [] for max(nObs, J)*eps.
%   wantStd - Whether to construct the coefficient covariance factor.
%
% Output arguments
%   coeffs - Fitted coefficients in field units. Size: (J, nField).
%   rnk, cond - Numerical rank and 2-norm condition number of B.
%   invFactor - Q*diag(1./s), where Q contains right singular vectors.
%       invFactor*invFactor.' is the coefficient covariance for unit
%       transformed observation variance. Size: (J, J), or [] if unused.
function [coeffs, rnk, cond, invFactor] = solveMask(B, rhs, field, tol, wantStd)

    [count, J] = size(B);
    % Solve with the SVD directly rather than forming B.'*B, which would
    % square the condition number and obscure the rank check.
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

%% PREPARECOORDS - Validates paired coordinates and returns column vectors.
%
% Input arguments
%   lon, lat - Coordinate vectors of equal length, in degrees.
%       Either orientation is accepted. NaNs are left for missingPolicy
%       to handle; Inf and latitudes outside [-90, 90] are invalid.
%
% Output arguments
%   lon, lat - Double column vectors in degrees. Size: (nPts, 1).
%       Longitude is wrapped to [0, 360); point order is preserved.
function [lon, lat] = preparecoords(lon, lat)

    arguments (Input)
        lon {mustBeNumeric, mustBeReal, mustBeVector}
        lat {mustBeNumeric, mustBeReal, mustBeVector}
    end

    arguments (Output)
        lon {mustBeNumeric, mustBeReal, mustBeColumn}
        lat {mustBeNumeric, mustBeReal, mustBeColumn}
    end

    if numel(lat) ~= numel(lon)
        error('ULMO:xyzs2slep:ObservationSize', ...
            'Length of longitude and latitude vectors must agree, but got lengths %d and %d, respectively.', ...
            numel(lon), numel(lat));
    end

    if any(isinf(lon)) || any(abs(lat) > 90)
        error('ULMO:xyzs2slep:InvalidObservation', 'Invalid longitude or latitude.');
    end

    lon = mod(double(lon(:)), 360);
    lat = double(lat(:));
end

%% PREPAREGVN - Constructs or accepts a basis and retains its leading modes.
%
% Input arguments
%   domain, L - Region and maximum degree accepted by XYZS2SLEP.
%   options - Validated options structure; uses GVN and truncation.
%
% Output arguments
%   GVN - Cell {G, V, N} for J retained modes. G is ((L+1)^2, J),
%       V is (J, 1), and N is the full scalar Shannon number.
function GVN = prepareGVN(domain, L, options)

    validateattributes(L, {'numeric'}, {'real', 'scalar', 'integer', 'nonnegative', 'finite'});
    J = options.truncation;

    if ~isempty(options.GVN)
        mustBeGVN(options.GVN);
        [G, V, N] = options.GVN{:};

        if size(G, 1) ~= (L + 1) ^ 2
            error('ULMO:xyzs2slep:InvalidBasis', ...
                'Supplied basis degree differs from L. For L = %d, the basis must have %d rows, but the received array has %d.', ...
                L, (L + 1) ^ 2, size(G, 1));
        end

        order = 1:size(G, 2);
    else

        if isnumeric(domain) && isscalar(domain)
            % Use the slepian_alpha routine for polar cap
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

    % A full basis defaults to the rounded Shannon count. A supplied subset
    % is treated as an intentional selection and retained in full by default.
    if isempty(J)
        J = size(G, 2);

        if J == (L + 1) ^ 2
            J = min(J, max(1, round(N)));
        end

    end

    if J > size(G, 2)
        error('ULMO:xyzs2slep:InvalidTruncation', ...
            'Truncation exceeds the available basis columns. The provided basis only contains %d functions, but %d are requested.', ...
            size(G, 2), J);
    end

    G = double(G(:, order(1:J)));
    V = V(:);
    GVN = {G, V(1:J), N};
end

%% SAMPLEBASIS - Evaluates retained Slepian functions at paired coordinates.
%
% Input arguments
%   lon, lat - Paired column vectors in degrees.
%       Size: (nPts, 1).
%   G - Retained real harmonic coefficients.
%       Size: ((L+1)^2, J).
%   L - Nonnegative maximum harmonic degree.
%   blockSize - Positive integer limit on points per harmonic evaluation.
%
% Output arguments
%   A - Dimensionless sampled basis.
%       Size: (nPts, J).
function A = sampleBasis(lon, lat, G, L, blockSize)

    A = nan(numel(lon), size(G, 2));
    pts = find(~isnan(lon) & ~isnan(lat));
    orders = addmout(L);
    % Convert YLM's unit-normalised harmonics and order-dependent phase to
    % the 4*pi real-harmonic convention used by the Slepian coefficients.
    scale = sqrt(4 * pi) * (-1) .^ orders(:);

    for first = 1:blockSize:numel(pts)
        rows = pts(first:min(first + blockSize - 1, end));

        % The degree-zero harmonic is one in this convention.
        if L == 0
            A(rows, :) = G;
        else
            Y = ylm([0, L], [], deg2rad(90 - lat(rows)), deg2rad(lon(rows)), [], [], 0, 1);
            A(rows, :) = (Y .* scale).' * G;
        end

    end

end

%% PREPAREUNCERTAINTY - Validates the observation error model on the full grid.
%
% Input arguments
%   sigma - Absolute standard-deviation vector in field units, or [].
%   cov - Absolute covariance in squared field units, or [].
%       Size: (cnt, cnt) when supplied.
%   estimate - Whether to estimate iid residual variance for each field.
%   cnt - Number of original coordinate pairs. Size: scalar.
%
% Output arguments
%   sigma - Validated full double column vector. Size: (cnt, 1), or [].
%   cov - Validated full double covariance, with roundoff asymmetry removed.
%   fac - Lower Cholesky factor satisfying cov = fac*fac.', or [].
%   src - Error source: "none", "dataStd", "dataCovariance", or "estimatedIid".
%
% Notes
%   The caller checks that sources are mutually exclusive. Validation
%   includes points omitted from the fit; missing data do not excuse
%   invalid uncertainties. No jitter or residual rescaling is applied.
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

        % Scale before the Frobenius-norm comparison to avoid overflow while
        % applying the same relative symmetry tolerance at every error scale.
        normedCov = cov / maxCov;

        if norm(normedCov - normedCov.', 'fro') > 100 * eps * norm(normedCov, 'fro')
            error('ULMO:xyzs2slep:AsymmetricCovariance', ...
            'dataCovariance must be symmetric within 100*eps relative Frobenius tolerance.');
        end

        % Average roundoff asymmetry without overflowing the sum of entries.
        % Cholesky then checks positive definiteness without regularisation.
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
