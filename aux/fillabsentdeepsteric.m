%% FILLABSENTDEEPSTERIC - Fills absent deep-ocean contributions with zero.
%
% Syntax
%   deep = FILLABSENTDEEPSTERIC(fullDepth, shallow, deep, oceanMask, deepOceanMask)
%
% Input arguments
%   fullDepth, shallow, deep - Native steric grids or same-epoch grid slices.
%   oceanMask, deepOceanMask - Static logical grids of native product coverage.
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

function deep = fillabsentdeepsteric(fullDepth, shallow, deep, oceanMask, deepOceanMask)
    shallowWater = oceanMask & ~deepOceanMask;
    outsideOcean = ~oceanMask;

    if ~any(shallowWater | outsideOcean, 'all')
        return;
    end

    for iDate = 1:size(deep, 3)
        field = deep(:, :, iDate);
        fill = shallowWater & isfinite(fullDepth(:, :, iDate)) & ...
            isfinite(shallow(:, :, iDate)) & isnan(field);
        field(fill) = 0;
        field(outsideOcean) = NaN;
        deep(:, :, iDate) = field;
    end
end
