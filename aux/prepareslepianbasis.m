%% PREPARESLEPIANBASIS - Selects harmonic Slepian basis columns for a fit.
%
% Syntax
%   GVN = PREPARESLEPIANBASIS(domain, L)
%   GVN = PREPARESLEPIANBASIS(domain, L, J, GVN, SaveData=false)
%
% Input arguments
%   domain - Geographic domain accepted by GLMALPHA_NEW.
%   L - Maximum harmonic degree.
%   J (optional) - Number of leading functions to retain.
%       With [], a full basis uses round(N) bounded by 1 and its column
%       count. An already truncated basis keeps all its columns.
%       Default count: [].
%   GVN (optional) - Supplied basis {G, V, N}, or {} to compute the basis.
%       G is the harmonic coefficient matrix, V contains one concentration
%       per column, and N is the full Shannon number.
%       Default basis: {}.
%   SaveData (name-value) - Whether to save computed basis and kernel caches.
%       Default flag: true.
%   RequireOrthonormal (name-value) - Whether to check the retained
%       columns for orthonormality using the input precision of G.
%       Default flag: false.
%
% Output arguments
%   GVN - Retained basis {G, V, N}. G is double and N is unchanged.
%       Computed columns are sorted by descending concentration. Supplied
%       columns retain their input order and signs.
%       Size: G is ((L+1)^2, J), V is (J, 1), and N is scalar.
%
% See also
%   glmalpha_new, mustBeGVN, xyzs2slep, fitgriddedslepian
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function GVN = prepareslepianbasis(domain, L, J, GVN, options)

    arguments (Input)
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        J double {mustBeScalarOrEmpty, mustBeFinite, mustBeInteger, mustBePositive} = []
        GVN cell = {}
        options.SaveData (1, 1) logical = true
        options.RequireOrthonormal (1, 1) logical = false
    end

    if isempty(GVN)
        [G, V, ~, ~, N] = glmalpha_new(domain, L, SaveData = options.SaveData, BeQuiet = true);
        [V, order] = sort(V(:), 'descend');
        G = G(:, order);
    else
        mustBeGVN(GVN, L);
        [G, V, N] = GVN{:};
    end

    if isempty(J)
        J = size(G, 2);

        if J == (L + 1) ^ 2
            J = min(J, max(1, round(N)));
        end

    end

    assert(J <= size(G, 2), 'ULMO:prepareslepianbasis:Truncation', 'Truncation exceeds basis size.');
    precision = double(eps(class(G)));
    G = double(G(:, 1:J));

    if options.RequireOrthonormal
        % Each inner product includes input rounding and double-precision summation.
        tol = 100 * max(precision, size(G, 1) * eps('double'));
        % For an orthogonal matrix, G^T * G is the identity
        assert(max(abs(G.' * G - eye(J)), [], 'all') <= tol, ...
            'ULMO:prepareslepianbasis:Basis', 'The retained harmonic basis must be orthonormal.');
    end

    V = V(:);
    GVN = {G, V(1:J), N};
end
