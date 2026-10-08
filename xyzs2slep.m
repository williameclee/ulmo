%% XYZS2SLEP - Fits scattered scalar observations to regional Slepian functions.
% Uses ordinary or generalised least squares at paired observation locations.
%
% Syntax
%   [falpha, V, N] = xyzs2slep(vals, lon, lat, domain, L)
%   [falpha, V, N, J, nObs, validIdxs, rnk, SVs, condNum, fitVals, resids, residDof] = xyzs2slep(__)
%   [__] = xyzs2slep(__, truncation=J, missingPolicy="omit")
%   [__, coefficientCovariance, coefficientStd, uncertaintySource] = ...
%       xyzs2slep(__, dataStd=sigma)
%   [__] = xyzs2slep(__, dataCovariance=Cdata)
%   [__] = xyzs2slep(__, estimateNoiseVariance=true)
%
% Input arguments
%   vals - Real vector of scalar observations, in field units.
%   lon, lat - Paired longitude and latitude vectors, in degrees.
%       Same length as values. Longitudes are wrapped to [0, 360).
%       Latitudes must lie in [-90, 90]. Row vectors are accepted.
%   domain - Concentration region accepted by glmalpha_new.
%       GeoDomain, region name, buffered-region cell, or closed [lon, lat]
%       polygon in degrees. A scalar is a north-polar cap radius (0, 180].
%       All supplied observations are used, even outside the domain.
%   L - Nonnegative integer maximum spherical-harmonic degree.
%   truncation (name-value) - Positive integer number of retained functions.
%       Default: max(1, round(N)), limited to the available basis dimension.
%       An explicit value is never silently reduced.
%   rankTolerance (name-value) - Relative singular-value threshold in [0, 1).
%       Default: max(nObs, J)*eps. Singular values must exceed the threshold
%       times the largest singular value; deficient fits raise an error.
%   missingPolicy (name-value) - "error" (default) or "omit".
%       "omit" removes NaN values/coordinates together. Inf is always invalid.
%   blockSize (name-value) - Positive integer maximum points per harmonic
%       evaluation block. Default: 1024. The full nObs-by-J matrix is stored.
%
%   dataStd (name-value) - Optional positive finite standard-deviation vector.
%       Same length and units as the original vals; mutually exclusive with
%       dataCovariance. Omission subsets it using validIdxs.
%   dataCovariance (name-value) - Optional positive-definite real covariance.
%       Original-observation-by-original-observation matrix in squared field
%       units. All entries must be finite, including omitted observations.
%       Relative Frobenius asymmetry up to 100*eps is averaged; larger
%       asymmetry is rejected. The full supplied matrix must be positive
%       definite. Omission subsets both axes using validIdxs.
%   estimateNoiseVariance (name-value) - Estimate an iid residual variance.
%       Default: false. Requires nObs>J and no supplied uncertainty. Uses
%       sum(resids.^2)/(nObs-J); this can include model mismatch.
%
% Output arguments
%   falpha - J-by-1 fitted coefficients, compatible with plm2slep_new and
%       slep2plm_new (4*pi-normalized real harmonics), in field units.
%   V - J-by-1 concentration eigenvalues in descending order.
%   N - Shannon number of the full concentration problem.
%   J - Number of retained basis functions.
%   nObs - Number of retained observations.
%   validIdxs - nObs-by-1 indices into the original observation vectors.
%   rnk - Numerical rank of the (whitened, when weighted) design matrix.
%   SVs - J-by-1 singular values of that same matrix, descending.
%   condNum - Ratio of largest to smallest singular value.
%   fitVals - nObs-by-1 predictions in field units.
%   resids - nObs-by-1 observed minus fitted values, in field units.
%   residDof - Residual degrees of freedom, nObs-J.
%   coefficientCovariance - J-by-J propagated covariance in squared field
%       units. Empty if no uncertainty is supplied or estimated. Supplied
%       absolute uncertainty is never rescaled using residual variance.
%   coefficientStd - J-by-1 square roots of the covariance diagonal, or [].
%   uncertaintySource - "dataStd", "dataCovariance", "estimatedIid", or "none".
%
% Notes
%   Without uncertainty, each observation has equal weight. Supplied errors
%   whiten the observations and design matrix before SVD. No area weighting
%   or regularisation is performed. Covariance is conditional on the fixed
%   basis/error model and excludes leakage, truncation bias, and coordinate
%   uncertainty. Singular covariance/exact constraints are not supported.
%   Rank-deficient fits are rejected; full-rank fits with condition number > 1/sqrt(eps) warn.
%   The full regional basis is requested before sorting, to avoid truncating
%   block-ordered cap eigenfunctions prematurely.
%
% See also
%   GLMALPHA_NEW, PLM2SLEP_NEW, SLEP2PLM_NEW, YLM, XYZ2SLEP
%
% Created by
%   2026/10/06, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/08, En-Chi Lee (williameclee@arizona.edu)

function [falpha, V, N, J, nObs, validIdxs, rnk, SVs, condNum, fitVals, resids, residDof, ...
        coefficientCovariance, coefficientStd, uncertaintySource] = ...
        xyzs2slep(vals, lon, lat, domain, L, options)

    arguments (Input)
        vals {mustBeNumeric, mustBeReal, mustBeVector}
        lon {mustBeNumeric, mustBeReal, mustBeVector}
        lat {mustBeNumeric, mustBeReal, mustBeVector}
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.rankTolerance double {mustBeScalarOrEmpty, mustBeFinite, mustBeNonnegative, mustBeLessThan(options.rankTolerance, 1)} = []
        options.missingPolicy (1, 1) string {mustBeMember(options.missingPolicy, ["error", "omit"])} = "error"
        options.blockSize (1, 1) double {mustBeFinite, mustBeInteger, mustBePositive} = 1024
        options.dataStd {mustBeNumeric, mustBeReal} = []
        options.dataCovariance {mustBeNumeric, mustBeReal} = []
        options.estimateNoiseVariance (1, 1) logical = false
    end

    originalCount = numel(vals);
    [vals, lon, lat, nObs, validIdxs, J, Jmax, isCap] = ...
        preprocessInputs(vals, lon, lat, L, domain, options.truncation, options.missingPolicy);

    [sigma, covarianceFactor, uncertaintySource] = prepareUncertainty( ...
        options.dataStd, options.dataCovariance, options.estimateNoiseVariance, ...
        originalCount, validIdxs);

    %% Computing the Slepian functions evaluated at each data point
    if isCap
        validateattributes(domain, {'numeric'}, {'real', 'finite', '>', 0, '<=', 180});
        % ULMO's extracted axisymmetric helper is incomplete; use Alpha.
        [G, V, ~, ~, N] = glmalpha(domain, L);
    else
        [G, V, ~, ~, N] = glmalpha_new(domain, L, 'BeQuiet', true);
    end

    [V, order] = sort(V(:), 'descend');

    if isempty(J)
        J = min(Jmax, max(1, round(N)));
    end

    if J > nObs
        error('ULMO:xyzs2slep:InvalidTruncation', ...
            ['The default truncation exceeds the observation count. ', ...
             'The Shannon number of the domain is %d, but there are only %d data points available.'
         'Specify a smaller truncation.'], ...
            J, nObs);
    end

    residDof = nObs - J;
    if options.estimateNoiseVariance && residDof <= 0
        error('ULMO:xyzs2slep:NoiseDegreesOfFreedom', ...
            'Estimating noise variance requires more observations than coefficients.');
    end

    G = G(:, order(1:J));
    V = V(1:J);

    % ylm uses unit normalisation and the Condon-Shortley phase. ULMO's
    % coefficient transforms use 4*pi normalisation and omit that phase.
    orders = addmout(L);
    scale = sqrt(4 * pi) * (-1) .^ orders(:);
    slepVals = zeros(nObs, J); % Slepian functions evaluated at each data point

    for first = 1:options.blockSize:nObs
        points = first:min(first + options.blockSize - 1, nObs);

        if L == 0
            slepVals(points, :) = repmat(G, numel(points), 1);
        else
            Y = ylm([0 L], [], deg2rad(90 - lat(points)), ...
                deg2rad(lon(points)), [], [], 0, 1);
            slepVals(points, :) = (Y .* scale).' * G;
        end

    end

    %% Inverting for the Slepian coefficients
    % The svd(A, 0) syntax in datafit (slepian_bravo) is not recommended per MathWork's documentation
    design = slepVals;
    rhs = vals;
    if ~isempty(sigma)
        design = slepVals ./ sigma;
        rhs = vals ./ sigma;
    elseif ~isempty(covarianceFactor)
        design = covarianceFactor \ slepVals;
        rhs = covarianceFactor \ vals;
    end
    if any(~isfinite(design(:))) || any(~isfinite(rhs))
        error('ULMO:xyzs2slep:WhiteningOverflow', ...
            'Whitening produced nonfinite values; rescale the field and uncertainty units.');
    end
    [U, S, Q] = svd(design, 'econ');
    SVs = diag(S);
    tol = options.rankTolerance;

    if isempty(tol)
        tol = max(nObs, J) * eps;
    end

    rnk = sum(SVs > tol * SVs(1));

    if rnk < J
        error('ULMO:xyzs2slep:RankDeficient', ...
            ['Sampled basis has rank %d of %d. ', ...
         'Reduce truncation or improve spatial coverage.'], ...
            rnk, J);
    end

    condNum = SVs(1) / SVs(end);

    if condNum > 1 / sqrt(eps)
        warning('ULMO:xyzs2slep:IllConditioned', ...
            ['Sampled basis condition number is large (%.3g); ', ...
         'coefficients may be poorly constrained.'] ...
            , condNum);
    end

    falpha = Q * ((U.' * rhs) ./ SVs);
    fitVals = [];
    resids = [];

    if nargout >= 10 || options.estimateNoiseVariance
        fitVals = slepVals * falpha;
    end

    if nargout >= 11 || options.estimateNoiseVariance
        resids = vals - fitVals;
    end

    coefficientCovariance = [];
    coefficientStd = [];
    if nargout >= 13 && uncertaintySource ~= "none"
        % Factor form avoids normal equations and an explicit matrix inverse.
        inverseFactor = Q ./ SVs.';
        if options.estimateNoiseVariance
            noiseStd = norm(resids) / sqrt(residDof);
            inverseFactor = inverseFactor * noiseStd;
        end
        coefficientCovariance = inverseFactor * inverseFactor.';
        if any(~isfinite(coefficientCovariance(:)))
            error('ULMO:xyzs2slep:CovarianceOverflow', ...
                'Coefficient covariance overflowed; rescale the field and uncertainty units.');
        end
        if nargout >= 14
            coefficientStd = sqrt(diag(coefficientCovariance));
        end
    end
end

%% Subfunctions

function [sigma, factor, source] = prepareUncertainty(sigma, covariance, estimate, count, validIdxs)
    factor = [];
    source = "none";
    if (~isempty(sigma) && ~isempty(covariance)) || ...
            (estimate && (~isempty(sigma) || ~isempty(covariance)))
        error('ULMO:xyzs2slep:ConflictingUncertainty', ...
            'Choose only one of dataStd, dataCovariance, or estimateNoiseVariance.');
    end

    if ~isempty(sigma)
        if ~isvector(sigma) || numel(sigma) ~= count || ...
                any(~isfinite(sigma(:))) || any(sigma(:) <= 0)
            error('ULMO:xyzs2slep:InvalidDataStd', ...
                'dataStd must have one finite positive value per original observation.');
        end
        sigma = full(double(sigma(:)));
        sigma = sigma(validIdxs);
        source = "dataStd";
    elseif ~isempty(covariance)
        if ~ismatrix(covariance) || ~isequal(size(covariance), [count count]) || ...
                any(~isfinite(covariance(:)))
            error('ULMO:xyzs2slep:InvalidDataCovariance', ...
                'dataCovariance must be a finite square matrix for the original observations.');
        end
        covariance = full(double(covariance));
        magnitude = max(abs(covariance(:)));
        if magnitude == 0
            error('ULMO:xyzs2slep:NonPositiveCovariance', ...
                'dataCovariance must be positive definite.');
        end
        normalised = covariance / magnitude;
        if norm(normalised - normalised.', 'fro') > 100 * eps * norm(normalised, 'fro')
            error('ULMO:xyzs2slep:AsymmetricCovariance', ...
                'dataCovariance must be symmetric within 100*eps relative Frobenius tolerance.');
        end
        covariance = covariance / 2 + covariance.' / 2;
        [factor, flag] = chol(covariance, 'lower');
        if flag ~= 0
            error('ULMO:xyzs2slep:NonPositiveCovariance', ...
                'dataCovariance must be positive definite; no jitter is added.');
        end
        if numel(validIdxs) ~= count
            % Subset the covariance, not its Cholesky factor: these operations
            % do not commute for correlated observations.
            factor = chol(covariance(validIdxs, validIdxs), 'lower');
        end
        source = "dataCovariance";
    elseif estimate
        source = "estimatedIid";
    end
end

function [vals, lon, lat, nObs, validIdxs, J, Jmax, domainIsCap] = ...
        preprocessInputs(vals, lon, lat, L, domain, J, missingPolicy)
    vals = double(vals(:));
    lon = double(lon(:));
    lat = double(lat(:));

    if isempty(vals)
        error('ULMO:xyzs2slep:ObservationSize', ...
        'Provide nonempty values.');
    elseif numel(vals) ~= numel(lon) || numel(vals) ~= numel(lat)
        error('ULMO:xyzs2slep:ObservationSize', ...
            ['The data points and their coordinate arrays have different sizes. ', ...
             'There are %d values, %d longitudes, and %d latitudes. ', ...
         'Provide same number for all three arrays.'], ...
            numel(vals), numel(lon), numel(lat));
    end

    if any(isinf(vals))
        error('ULMO:xyzs2slep:InvalidObservation', ...
        'Infinite observations/coordinates are invalid.');
    elseif any(isinf(lon) | isinf(lat))
        error('ULMO:xyzs2slep:InvalidObservation', ...
        'Infinite coordinates and latitudes outside [-90, 90] are invalid.');
    elseif any(abs(lat) > 90)
        error('ULMO:xyzs2slep:InvalidObservation', ...
        'Latitudes outside [-90, 90] are invalid.');
    end

    isValid = ~(isnan(vals) | isnan(lon) | isnan(lat));

    if ~all(isValid) && missingPolicy == "error"
        error('ULMO:xyzs2slep:MissingObservation', ...
        'NaN observations/coordinates require missingPolicy="omit".');
    end

    validIdxs = find(isValid);
    vals = vals(validIdxs);
    lon = mod(lon(validIdxs), 360);
    lat = lat(validIdxs);
    nObs = numel(vals);

    if nObs == 0
        error('ULMO:xyzs2slep:ObservationSize', 'No usable observations remain.');
    end

    Jmax = (L + 1) ^ 2;

    if ~isempty(J)

        if J > Jmax
            error('ULMO:xyzs2slep:InvalidTruncation', ...
                ['Truncation must not exceed (L+1)^2. ', ...
             'Maximum functions at degree %d is %d, but %d functions are requested.'], ...
                L, Jmax, J);
        elseif J > nObs
            error('ULMO:xyzs2slep:InvalidTruncation', ...
                ['Truncation must not exceed the usable observation count. ', ...
             'There are only %d valid data points, but %d functions are requested.'], ...
                nObs, J);
        end

    end

    domainIsCap = isnumeric(domain) && isscalar(domain);

    if ~domainIsCap && ~((isa(domain, 'GeoDomain') && isscalar(domain)) || ...
            ischar(domain) || (isstring(domain) && isscalar(domain)) || ...
            iscell(domain) || (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2))
        error('ULMO:xyzs2slep:InvalidDomain', ...
            ['Use a GeoDomain, region name/cell, [lon, lat] polygon, or scalar cap radius.', ...
         'The input domain format of type %s is not supported.'], class(domain));
    end

end
