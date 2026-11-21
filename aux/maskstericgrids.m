%% MASKSTERICGRIDS - Applies shared ocean validity rules to steric grids.
%
% Syntax
%   [data, masks] = MASKSTERICGRIDS(data, schema, nEpochs)
%
% Input arguments
%   data - Structure containing the supplied steric grids.
%   schema - Field schema returned by STERICFIELDSCHEMA.
%   nEpochs - Number of input timestamps.
%
% Output arguments
%   data - Steric grids with absent deep-ocean contributions filled.
%   masks - Static ocean and, for split products, deep-ocean coverage masks.
%
% See also
%   stericcoverage, fillabsentdeepsteric, steric2lonlatt, steric2slept
%
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function [data, masks] = maskstericgrids(data, schema, nEpochs)
    masks = stericcoverage(data, schema.gridNames, schema.available, nEpochs);

    for iComp = find(schema.split).'
        names = schema.gridNames(iComp, :);
        data.(names(3)) = fillabsentdeepsteric( ...
            data.(names(1)), data.(names(2)), data.(names(3)), masks.oceanMask, masks.deepOceanMask);
    end

end
