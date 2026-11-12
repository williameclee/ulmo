%% SEISMICCAPBASIS - Constructs a retained Slepian basis for an epicentral cap.
%
% Syntax
%   GVN = seismiccapbasis(epicentre, radius, L, truncation, SaveData=false)
%
% Input arguments
%   epicentre - [longitude, latitude] in degrees.
%   radius - Cap radius in degrees, including both land and ocean.
%   L - Maximum harmonic degree, at least 2. The basis uses degrees 2:L.
%   truncation - Positive mode count, or [] for max(1, round(N)).
%   SaveData (name-value) - Whether to save the polar basis where supported.
%       Default flag: true. PLM2ROT's dependency DLMB has no save control.
%
% Output arguments
%   GVN - {G, V, N}, with G in ADDMOUT order and descending concentrations.
%       Only retained columns are rotated; their order and signs are kept.
%
% See also
%   prepareslepianbasis, glmalpha_new, plm2rot
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/13, En-Chi Lee (williameclee@arizona.edu)

function GVN = seismiccapbasis(epicentre, radius, L, truncation, options)
    arguments (Input)
        epicentre (1, 2) double {mustBeFinite, mustBeReal}
        radius (1, 1) double {mustBePositive, mustBeLessThanOrEqual(radius, 180)}
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeGreaterThanOrEqual(L, 2)}
        truncation = []
        options.SaveData (1, 1) logical = true
    end
    [G, V, ~, ~, N] = glmalpha_new(radius, [2, L], SaveData = options.SaveData, BeQuiet = true);
    [V, order] = sort(V(:), 'descend');
    if isempty(truncation), truncation = min(size(G, 2), max(1, round(N))); end
    GVN = prepareslepianbasis(radius, L, truncation, {G(:, order), V, N});
    G = GVN{1};
    [~, ~, ~, template, ~, ~, ~, ~, ~, mapping] = addmon(L);
    for iMode = 1:size(G, 2)
        coefficients = zeros(size(template, 1), 2);
        coefficients(mapping) = G(:, iMode);
        rotated = plm2rot([template(:, 1:2), coefficients], 0, epicentre(2) - 90, -epicentre(1));
        coefficients = rotated(:, 3:4);
        G(:, iMode) = coefficients(mapping);
    end
    GVN{1} = G;
end
