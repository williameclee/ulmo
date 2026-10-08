%% XYZS2SLEP - Fits scattered scalar observations to regional Slepian functions.
% Uses ordinary least squares at paired observation locations without gridding.
%
% Syntax
%   [falpha, V, N] = xyzs2slep(values, lon, lat, domain, L)
%   [falpha, V, N, J, nObs, usedIndices, numericalRank, singularValues, ...
%       conditionNumber, fittedValues, residuals, residualDof] = xyzs2slep(__)
%   [__] = xyzs2slep(__, truncation=J, missingPolicy="omit")
%
% Input arguments
%   values - Real vector of scalar observations, in field units.
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
%   usedIndices - nObs-by-1 indices into the original observation vectors.
%   numericalRank - Numerical rank of the sampled design matrix.
%   singularValues - J-by-1 singular values of that matrix, descending.
%   conditionNumber - Ratio of largest to smallest singular value.
%   fittedValues - nObs-by-1 predictions in field units.
%   residuals - nObs-by-1 observed minus fitted values, in field units.
%   residualDof - Residual degrees of freedom, nObs-J.
%
% Notes
%   s means scattered, not MATLAB sparse storage. Each observation has equal
%   weight, so densely sampled areas contribute more. No error propagation,
%   area weighting, or regularization is performed. Rank-deficient fits are
%   rejected; full-rank fits with condition number > 1/sqrt(eps) warn.
%   The full regional basis is requested before sorting, to avoid truncating
%   block-ordered cap eigenfunctions prematurely. Basis construction uses
%   existing Slepian caches under IFILES. Scalar caps use Alpha's glmalpha;
%   geographic domains use ULMO's glmalpha_new. The latter's domain topology
%   and geometric integration limitations apply unchanged.
%   Bravo's xyz2slep remains unchanged and still resolves under that name.
%
% See also
%   GLMALPHA_NEW, PLM2SLEP_NEW, SLEP2PLM_NEW, YLM, XYZ2SLEP
%
% Created by
%   2026/10/06, En-Chi Lee (williameclee@arizona.edu)

function [falpha, V, N, J, nObs, usedIndices, numericalRank, singularValues, ...
        conditionNumber, fittedValues, residuals, residualDof] = ...
        xyzs2slep(values, lon, lat, domain, L, options)
    arguments (Input)
        values {mustBeNumeric, mustBeReal, mustBeVector}
        lon {mustBeNumeric, mustBeReal, mustBeVector}
        lat {mustBeNumeric, mustBeReal, mustBeVector}
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        options.truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        options.rankTolerance double {mustBeScalarOrEmpty, mustBeFinite, mustBeNonnegative, mustBeLessThan(options.rankTolerance, 1)} = []
        options.missingPolicy (1, 1) string {mustBeMember(options.missingPolicy, ["error", "omit"])} = "error"
        options.blockSize (1, 1) double {mustBeFinite, mustBeInteger, mustBePositive} = 1024
    end

    values = double(values(:));
    lon = double(lon(:));
    lat = double(lat(:));
    if isempty(values) || numel(values) ~= numel(lon) || numel(values) ~= numel(lat)
        error('ULMO:xyzs2slep:ObservationSize', ...
            'Provide nonempty values, lon, and lat vectors of equal length.');
    end
    if any(isinf(values) | isinf(lon) | isinf(lat)) || any(abs(lat) > 90)
        error('ULMO:xyzs2slep:InvalidObservation', ...
            'Infinite observations/coordinates and latitudes outside [-90, 90] are invalid.');
    end
    missing = isnan(values) | isnan(lon) | isnan(lat);
    if any(missing) && options.missingPolicy == "error"
        error('ULMO:xyzs2slep:MissingObservation', ...
            'NaN observations/coordinates require missingPolicy="omit".');
    end
    usedIndices = find(~missing);
    values = values(usedIndices);
    lon = mod(lon(usedIndices), 360);
    lat = lat(usedIndices);
    nObs = numel(values);
    if nObs == 0
        error('ULMO:xyzs2slep:ObservationSize', 'No usable observations remain.');
    end
    dimension = (L + 1)^2;
    J = options.truncation;
    if ~isempty(J) && (J > dimension || J > nObs)
        error('ULMO:xyzs2slep:InvalidTruncation', ...
            'truncation must not exceed (L+1)^2 or the usable observation count.');
    end

    isCap = isnumeric(domain) && isscalar(domain);
    if isCap
        validateattributes(domain, {'numeric'}, {'real', 'finite', '>', 0, '<=', 180});
        if L == 0
            % Degree zero has one constant basis function, including at 180 deg.
            G = 1;
            V = (1 - cosd(domain))/2;
            N = V;
        else
            % ULMO's extracted axisymmetric helper is incomplete; use Alpha.
            [G, V, ~, ~, N] = glmalpha(domain, L);
        end
    else
        if ~((isa(domain, 'GeoDomain') && isscalar(domain)) || ...
                ischar(domain) || (isstring(domain) && isscalar(domain)) || ...
                iscell(domain) || (isnumeric(domain) && ismatrix(domain) && size(domain, 2) == 2))
            error('ULMO:xyzs2slep:InvalidDomain', ...
                'Use a GeoDomain, region name/cell, [lon, lat] polygon, or scalar cap radius.');
        end
        [G, V, ~, ~, N] = glmalpha_new(domain, L, 'BeQuiet', true);
    end
    [V, order] = sort(V(:), 'descend');
    if isempty(J)
        J = min(dimension, max(1, round(N)));
    end
    if J > nObs
        error('ULMO:xyzs2slep:InvalidTruncation', ...
            'The default truncation exceeds the observation count; specify a smaller truncation.');
    end
    G = G(:, order(1:J));
    V = V(1:J);

    % ylm uses unit normalization and the Condon-Shortley phase. ULMO's
    % coefficient transforms use 4*pi normalization and omit that phase.
    orders = addmout(L);
    scale = sqrt(4*pi) * (-1).^orders(:);
    A = zeros(nObs, J);
    for first = 1:options.blockSize:nObs
        points = first:min(first + options.blockSize - 1, nObs);
        if L == 0
            A(points, :) = repmat(G, numel(points), 1);
        else
            Y = ylm([0 L], [], deg2rad(90-lat(points)), ...
                deg2rad(lon(points)), [], [], 0, 1);
            A(points, :) = (Y .* scale).' * G;
        end
    end

    [U, S, Q] = svd(A, 'econ');
    singularValues = diag(S);
    tolerance = options.rankTolerance;
    if isempty(tolerance)
        tolerance = max(nObs, J)*eps;
    end
    numericalRank = sum(singularValues > tolerance*singularValues(1));
    if numericalRank < J
        error('ULMO:xyzs2slep:RankDeficient', ...
            'Sampled basis has rank %d of %d. Reduce truncation or improve spatial coverage.', ...
            numericalRank, J);
    end
    conditionNumber = singularValues(1)/singularValues(end);
    if conditionNumber > 1/sqrt(eps)
        warning('ULMO:xyzs2slep:IllConditioned', ...
            'Sampled basis condition number is %.3g; coefficients may be poorly constrained.', conditionNumber);
    end
    falpha = Q * ((U.' * values)./singularValues);
    fittedValues = [];
    residuals = [];
    if nargout >= 10
        fittedValues = A*falpha;
    end
    if nargout >= 11
        residuals = values-fittedValues;
    end
    residualDof = nObs-J;
end
