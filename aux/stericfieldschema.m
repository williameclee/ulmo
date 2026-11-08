%% STERICFIELDSCHEMA - Describe supplied steric components and depth coverage.
% Resolves native grid names, global-mean names and depth consistency from
% MAT-file metadata without reading grids or inferring coverage from products.
%
% Syntax
%   schema = STERICFIELDSCHEMA()
%   schema = STERICFIELDSCHEMA(metadata)
%
% Input arguments
%   metadata (optional) - Structure array returned by whos('-file', path).
%       Only the name and size fields are inspected. Missing and empty grids
%       are unavailable; nonempty zero-valued grids remain available.
%       Size: structure array with one entry per stored variable.
%       Default metadata: empty structure array (no supplied grids).
%
% Output arguments
%   schema - Scalar structure with the following fields.
%   schema.types - Component labels: total, thermosteric, halosteric,
%       non-thermosteric and non-halosteric, in row order.
%       Size: (1, 5) string.
%   schema.depths - Depth labels: full, shallow and deep, in column order.
%       Size: (1, 3) string.
%   schema.gridNames - Native grid names for the first three component rows.
%       Size: (3, 3) string.
%   schema.meanNames - Cached global-mean names for all component rows.
%       Size: (5, 3) string.
%   schema.available - Whether each native grid is present and nonempty.
%       Size: (3, 3) logical.
%   schema.supported - Whether each mean has its required native grids.
%       Residual rows require total and the subtracted component at that depth.
%       Size: (5, 3) logical.
%   schema.split - Whether either shallow or deep is supplied per component.
%       Size: (3, 1) logical.
%   schema.incompleteSplit - Whether any split lacks full, shallow or deep.
%       Size: scalar logical.
%   schema.incompatibleSplits - Whether supplied full total and halosteric
%       components disagree on the presence of a depth split.
%       Size: scalar logical.
%
% Notes
%   Availability and support describe stored fields even for invalid splits.
%   Callers must check the split flags before using native schemas; the flags
%   allow callers to retain their own error identifiers and validation order.
%   No component is required here. Coordinate requirements, numeric types,
%   grid geometry, finite support and numerical layer closure belong to callers.
%   Residuals have mean names but no native grid names. They are total minus
%   the named component, not aliases for another physical component.
%
% See also
%   steric2lonlatt, steric2gmsl, stericproductinfo
%
% Created by
%   2026/11/05, OpenAI Codex

function schema = stericfieldschema(metadata)
    arguments (Input)
        metadata struct = struct('name', {}, 'size', {})
    end

    schema.types = ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"];
    schema.depths = ["full", "shallow", "deep"];
    schema.gridNames = ["stericSl", "shallowStericSl", "deepStericSl"; ...
                        "thermostericSl", "shallowThermostericSl", "deepThermostericSl"; ...
                        "halostericSl", "shallowHalostericSl", "deepHalostericSl"];
    schema.meanNames = ["stericGmsl", "shallowStericGmsl", "deepStericGmsl"; ...
                        "thermostericGmsl", "shallowThermostericGmsl", "deepThermostericGmsl"; ...
                        "halostericGmsl", "shallowHalostericGmsl", "deepHalostericGmsl"; ...
                        "nonThermostericGmsl", "shallowNonThermostericGmsl", "deepNonThermostericGmsl"; ...
                        "nonHalostericGmsl", "shallowNonHalostericGmsl", "deepNonHalostericGmsl"];
    available = false(size(schema.gridNames));
    for iField = 1:numel(schema.gridNames)
        entry = metadata(strcmp({metadata.name}, schema.gridNames(iField)));
        available(iField) = ~isempty(entry) && all(entry.size > 0);
    end
    schema.available = available;
    schema.supported = [available; available(1, :) & available(2, :); ...
                       available(1, :) & available(3, :)];
    schema.split = any(available(:, 2:3), 2);
    schema.incompleteSplit = any(schema.split & ~all(available, 2));
    schema.incompatibleSplits = available(1, 1) && available(3, 1) && ...
        schema.split(1) ~= schema.split(3);
end
