%% STERICFIELDSCHEMA - Describe supplied steric components and depth coverage.
% Resolves native grid names, global-mean names and depth consistency from
% MAT-file metadata and the declared native depth range without reading grids.
%
% Syntax
%   schema = STERICFIELDSCHEMA()
%   schema = STERICFIELDSCHEMA(metadata)
%   schema = STERICFIELDSCHEMA(metadata, nativeDepth)
%
% Input arguments
%   metadata (optional) - Structure array returned by whos('-file', path).
%       Only missing and empty grids are unavailable; zero-valued ones remain available.
%       Default metadata: empty structure array (no supplied grids).
%   nativeDepth (optional) - Depth range represented by unprefixed native fields.
%       full-depth or shallow. With shallow coverage, unprefixed fields map
%       to the shallow column.
%       Default depth range: full-depth.
%
% Output arguments
%   schema - Scalar structure with the following fields.
%   schema.components - Component labels: total, thermosteric, halosteric,
%       non-thermosteric and non-halosteric, in row order.
%       Size: (1, 5) string.
%   schema.depths - Depth labels: full-depth, shallow and deep, in column order.
%       Full-depth means integration from the sea surface to the seafloor.
%       Size: (1, 3) string.
%   schema.gridNames - Native grid names for the first three component rows.
%       For shallow native coverage, the shallow column uses unprefixed names.
%       Size: (3, 3) string.
%   schema.meanNames - Global-mean names for all component rows.
%       Only the first three rows are cached; difference components are derived.
%       Size: (5, 3) string.
%   schema.cacheGridNames - Processed grid names for the three supplied components.
%       Depth names are explicit even for products with shallow native coverage.
%       Size: (3, 3) string.
%   schema.available - Whether each native grid is present and nonempty.
%       Size: (3, 3) logical.
%   schema.supported - Whether each component/depth has its required native grids.
%       Residual rows require total and the subtracted component at that depth.
%       Size: (5, 3) logical.
%   schema.split - Whether a full-depth native component supplies either layer.
%       Size: (3, 1) logical.
%   schema.incompleteSplit - Whether any split lacks full-depth, shallow or deep.
%       Size: scalar logical.
%   schema.incompatibleSplits - Whether supplied full-depth total steric and halosteric
%       components disagree on the presence of a depth split.
%       Size: scalar logical.
%
% Notes
%   Availability and support describe stored fields even for invalid splits.
%   Callers must check the split flags before using native schemas; the flags
%   allow callers to retain their own error identifiers and validation order.
%   Residuals have mean names but no native grid names.
%
% See also
%   steric2lonlatt, steric2gmsl, stericproductinfo
%
% Created by
%   2026/11/05, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/20, En-Chi Lee (williameclee@gmail.com)

function schema = stericfieldschema(metadata, nativeDepth)

    arguments (Input)
        metadata struct = struct('name', {}, 'size', {})
        nativeDepth (1, 1) string {mustBeMember(nativeDepth, ["full-depth", "shallow"])} = "full-depth"
    end

    arguments (Output)
        schema (1, 1) struct
    end

    schema.components = ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"];
    schema.depths = ["full-depth", "shallow", "deep"];
    schema.gridNames = ...
        ["stericSl", "shallowStericSl", "deepStericSl"; ...
         "thermostericSl", "shallowThermostericSl", "deepThermostericSl"; ...
         "halostericSl", "shallowHalostericSl", "deepHalostericSl"];
    schema.meanNames = ...
        ["stericGmsl", "shallowStericGmsl", "deepStericGmsl"; ...
         "thermostericGmsl", "shallowThermostericGmsl", "deepThermostericGmsl"; ...
         "halostericGmsl", "shallowHalostericGmsl", "deepHalostericGmsl"; ...
         "nonThermostericGmsl", "shallowNonThermostericGmsl", "deepNonThermostericGmsl"; ...
         "nonHalostericGmsl", "shallowNonHalostericGmsl", "deepNonHalostericGmsl"];
    schema.cacheGridNames = replace(schema.meanNames(1:3, :), "Gmsl", "Sl");

    if nativeDepth == "shallow"
        schema.gridNames(:, 2) = schema.gridNames(:, 1);
    end

    avail = false(size(schema.gridNames));

    for iField = 1:numel(schema.gridNames)
        entry = metadata(strcmp({metadata.name}, schema.gridNames(iField)));
        avail(iField) = ~isempty(entry) && all(entry.size > 0);
    end

    if nativeDepth == "shallow"
        avail(:, [1, 3]) = false;
    end

    schema.available = avail;
    schema.supported = ...
        [avail; avail(1, :) & avail(2, :); ... % non-thermosteric
         avail(1, :) & avail(3, :)]; % non-halosteric
    schema.split = nativeDepth == "full-depth" & any(avail(:, 2:3), 2);
    schema.incompleteSplit = any(schema.split & ~all(avail, 2));
    schema.incompatibleSplits = avail(1, 1) && avail(3, 1) && ...
        schema.split(1) ~= schema.split(3);
end
