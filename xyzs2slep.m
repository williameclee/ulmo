%% XYZS2SLEP - Fits scattered scalar observations to regional Slepian functions.
% Uses ordinary least squares at paired observation locations.
%
% Syntax
%   [falpha, V, N] = xyzs2slep(vals, lon, lat, domain, L)
%   [falpha, V, N, J, nObs, validIdxs, rnk, SVs, condNum, fitVals, resids, residDof] = xyzs2slep(__)
%   [__] = xyzs2slep(__, truncation=J, missingPolicy="omit")
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
% Output arguments
%   falpha - J-by-1 fitted coefficients, compatible with plm2slep_new and
%       slep2plm_new (4*pi-normalized real harmonics), in field units.
%   V - J-by-1 concentration eigenvalues in descending order.
%   N - Shannon number of the full concentration problem.
%   J - Number of retained basis functions.
%   nObs - Number of retained observations.
%   validIdx - nObs-by-1 indices into the original observation vectors.
%   rnk - Numerical rank of the sampled design matrix.
%   SVs - J-by-1 singular values of that matrix, descending.
%   condNum - Ratio of largest to smallest singular value.
%   fitVals - nObs-by-1 predictions in field units.
%   resids - nObs-by-1 observed minus fitted values, in field units.
%   residDof - Residual degrees of freedom, nObs-J.
%
% Notes
%   Each observation has equal weight, so densely sampled areas contribute more.
%   No error propagation, area weighting, or regularisation is performed.
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

function [falpha, V, N, J, nObs, validIdxs, rnk, SVs, condNum, fitVals, resids, residDof] = ...
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
    end

    [vals, lon, lat, nObs, validIdxs, J, Jmax] = ...
        preprocessInputs(vals, lon, lat, L, options.truncation, options.missingPolicy);

    isCap = isnumeric(domain) && isscalar(domain);

    if isCap
        validateattributes(domain, {'numeric'}, {'real', 'finite', '>', 0, '<=', 180});
        % ULMO's extracted axisymmetric helper is incomplete; use Alpha.
        [G, V, ~, ~, N] = glmalpha(domain, L);
    else

        if ~((isa(domain, 'GeoDomain') && isscalar(domain)) || ...
                ischar(domain) || (isstring(domain) && isscalar(domain)) || ...
                iscell(domain) || (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2))
            error('ULMO:xyzs2slep:InvalidDomain', ...
                ['Use a GeoDomain, region name/cell, [lon, lat] polygon, or scalar cap radius.', ...
             'The input domain format of type %s is not supported.'], class(domain));
        end

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

    G = G(:, order(1:J));
    V = V(1:J);

    % ylm uses unit normalisation and the Condon-Shortley phase. ULMO's
    % coefficient transforms use 4*pi normalisation and omit that phase.
    orders = addmout(L);
    scale = sqrt(4 * pi) * (-1) .^ orders(:);
    A = zeros(nObs, J);

    for first = 1:options.blockSize:nObs
        points = first:min(first + options.blockSize - 1, nObs);

        if L == 0
            A(points, :) = repmat(G, numel(points), 1);
        else
            Y = ylm([0 L], [], deg2rad(90 - lat(points)), ...
                deg2rad(lon(points)), [], [], 0, 1);
            A(points, :) = (Y .* scale).' * G;
        end

    end

    % The svd(A, 0) syntax in datafit (slepian_bravo) is not recommended per MathWork's documentation
    [U, S, Q] = svd(A, 'econ');
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

    falpha = Q * ((U.' * vals) ./ SVs);
    fitVals = [];
    resids = [];

    if nargout >= 10
        fitVals = A * falpha;
    end

    if nargout >= 11
        resids = vals - fitVals;
    end

    residDof = nObs - J;
end

%% Subfunctions

function [vals, lon, lat, nObs, validIdxs, J, Jmax] = ...
        preprocessInputs(vals, lon, lat, L, J, missingPolicy)
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

end
