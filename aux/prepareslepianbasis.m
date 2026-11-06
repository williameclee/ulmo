%% PREPARESLEPIANBASIS - Prepare retained harmonic Slepian basis columns.
%
% Syntax
%   basis = PREPARESLEPIANBASIS(domain, L)
%   basis = PREPARESLEPIANBASIS(domain, L, truncation, supplied, SaveData=false)
%
% Input arguments
%   domain - Geographic domain accepted by GLMALPHA_NEW.
%   L - Nonnegative integer maximum harmonic degree.
%   truncation (optional) - Positive retained column count, or [].
%       Default count: bounded round of N for full bases; all columns for
%       an already truncated supplied basis.
%   supplied (optional) - Cell array {G,V,N}, or {} to calculate the basis.
%       G has (L+1)^2 harmonic rows; V has one value per column; N is the
%       full Shannon number. Supplied column order and signs are preserved.
%       Default basis: {}.
%   SaveData (name-value) - Whether computed basis/kernel caches may be saved.
%       Default flag: true.
%   RequireOrthonormal (name-value) - Whether retained columns must be
%       orthonormal, using the precision of the supplied harmonic matrix.
%       Default flag: false (the XYZS2SLEP validation convention).
%
% Output arguments
%   basis - Retained {G,V,N}; G is double, V is a column, N is unchanged.
%       Size: G is ((L+1)^2, J), V is (J, 1), N is scalar.
%
% Notes
%   Computed bases are sorted by descending concentration. This helper does
%   not read or write projection caches. Callers own supplied-basis bypass.
%
% See also
%   glmalpha_new, xyzs2slep, fitgriddedslepian
%
% Created by
%   2026/11/06, En-Chi Lee

function basis = prepareslepianbasis(domain, L, truncation, supplied, options)
    arguments (Input)
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        truncation double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        supplied cell = {}
        options.SaveData (1, 1) logical = true
        options.RequireOrthonormal (1, 1) logical = false
    end
    if isempty(supplied)
        [G, V, ~, ~, N] = glmalpha_new(domain, L, SaveData=options.SaveData, BeQuiet=true);
        [V, order] = sort(V(:), 'descend');
        G = G(:, order);
    else
        assert(numel(supplied) == 3, 'ULMO:prepareslepianbasis:Basis', 'Supply {G,V,N}.');
        [G, V, N] = supplied{:};
    end
    validateattributes(G, {'double', 'single'}, {'real', 'finite', '2d', 'nonempty', 'nrows', (L+1)^2});
    validateattributes(V, {'numeric'}, {'real', 'finite', 'vector', 'numel', size(G, 2)});
    validateattributes(N, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});
    assert(size(G, 2) <= size(G, 1), 'ULMO:prepareslepianbasis:Basis', 'Too many basis columns.');
    J = truncation;
    if isempty(J)
        J = size(G, 2);
        if J == (L+1)^2, J = min(J, max(1, round(N))); end
    end
    assert(J <= size(G, 2), 'ULMO:prepareslepianbasis:Truncation', 'Truncation exceeds basis size.');
    precision = double(eps(class(G)));
    G = double(G(:, 1:J));
    if options.RequireOrthonormal
        assert(norm(G.'*G - eye(J), 'fro') <= 100*precision*max(size(G)), ...
            'ULMO:prepareslepianbasis:Basis', 'The retained harmonic basis must be orthonormal.');
    end
    V = V(:);
    basis = {G, V(1:J), N};
end
