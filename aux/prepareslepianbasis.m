%% PREPARESLEPIANBASIS - Prepares retained harmonic Slepian basis columns.
%
% Syntax
%   GVN = PREPARESLEPIANBASIS(domain, L)
%   GVN = PREPARESLEPIANBASIS(domain, L, J, GVN, SaveData=false)
%
% Input arguments
%   domain - Geographic domain.
%   L - Maximum harmonic degree for the basis.
%   J (optional) - Positive retained column count.
%       Default count: bounded round of N for full bases; all columns for
%       an already truncated supplied basis.
%   GVN (optional) - Cell array {G,V,N}, or {} to calculate the basis.
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
%   2026/11/07, En-Chi Lee (williameclee@gmail.com)

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
        assert(norm(G.' * G - eye(J), 'fro') <= 100 * precision * max(size(G)), ...
            'ULMO:prepareslepianbasis:Basis', 'The retained harmonic basis must be orthonormal.');
    end

    V = V(:);
    GVN = {G, V(1:J), N};
end
