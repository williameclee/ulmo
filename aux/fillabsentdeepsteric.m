%% FILLABSENTDEEPSTERIC - Fills absent deep-ocean contributions with zero.
%
% Syntax
%   deep = FILLABSENTDEEPSTERIC(full, shallow, deep, oceanMask, deepOceanMask)
%
% Input arguments
%   full, shallow, deep - Steric grids or same-epoch grid slices.
%   oceanMask, deepOceanMask - Static logical grids of product coverage.
%
% Output arguments
%   deep - Deep-ocean field with structural NaNs filled at shared valid cells.
%       The full-depth and upper-ocean values must both be finite at the epoch.
%       Cells outside ocean coverage are missing.
%
% See also
%   stericcoverage, steric2lonlatt, steric2gmsl
%
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function deep = fillabsentdeepsteric(full, shallow, deep, oceanMask, deepOceanMask)
    shallowWater = oceanMask & ~deepOceanMask;
    outsideOcean = ~oceanMask;

    if ~any(shallowWater | outsideOcean, 'all')
        return;
    end

    for iEpoch = 1:size(deep, 3)
        deepSlice = deep(:, :, iEpoch);
        fill = shallowWater & isfinite(full(:, :, iEpoch)) & ...
            isfinite(shallow(:, :, iEpoch)) & isnan(deepSlice);
        deepSlice(fill) = 0;
        deepSlice(outsideOcean) = NaN;
        deep(:, :, iEpoch) = deepSlice;
    end

end
