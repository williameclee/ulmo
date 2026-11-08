%% MUSTBEGVN - Checks the structure and dimensions of a supplied Slepian basis.
%
% Syntax
%   MUSTBEGVN(GVN)
%   MUSTBEGVN(GVN, L)
%
% Input arguments
%   GVN - Three-element cell array {G, V, N}.
%       G is a nonempty, finite, real single- or double-precision matrix
%       with a square number of rows and at most that many columns.
%       V contains one finite real eigenvalue per column. N is a finite,
%       nonnegative scalar Shannon number. All entries are dimensionless.
%   L (optional) - Maximum harmonic degree. When supplied, G must have
%       exactly (L+1)^2 rows.
%       Default degree: [] (infer from the row count).
%
% Created by
%   2026/11/07, En-Chi Lee (williameclee@arizona.edu)

function mustBeGVN(GVN, L)

    arguments (Input)
        GVN cell
        L double = []
    end

    if numel(GVN) ~= 3
        error('ULMO:xyzs2slep:InvalidBasis', ...
            'GVN must be a 3-element cell, but only received %d elements.', ...
            numel(GVN));
    end

    [G, V, N] = GVN{:};
    validateattributes(G, {'double', 'single'}, {'real', 'finite', '2d', 'nonempty'});
    validateattributes(V, {'numeric'}, {'real', 'finite', 'vector', 'numel', size(G, 2)});
    validateattributes(N, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});

    if ~isempty(L) && size(G, 1) ~= (L + 1) ^ 2
        error('ULMO:xyzs2slep:InvalidBasis', ...
            'G must have %d rows for Lmax = %d, but the received array only has %d rows.', ...
            (L + 1) ^ 2, L, size(G, 1));
    elseif isempty(L) && mod(sqrt(size(G, 1)), 1) ~= 0
        error('ULMO:xyzs2slep:InvalidBasis', ...
            'G must have (L+1)^2 rows, but the received array has %d rows, which is not a square number.', ...
            size(G, 1));
    end

    if size(G, 2) > size(G, 1)
        error('ULMO:xyzs2slep:InvalidBasis', ...
            'G must not have more columns than rows, but received of size (%d, %d).', size(G, 1), size(G, 2));
    end

end
