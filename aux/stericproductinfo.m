%% STERICPRODUCTINFO - Inspect explicitly supplied native steric depth fields.
%
% Syntax
%   info = STERICPRODUCTINFO(product)
%
% Input arguments
%   product - Prepared steric product identifier.
%
% Output arguments
%   info - Native file metadata: product, path, fields, hasShallowField,
%       hasDeepField, hasDeep and hasHalosteric. A depth split requires
%       explicit full, shallow and deep grids. Empty grids are unsupported.
%       Supplied total and halosteric splits must agree. No correction
%       application policy is inferred from the product name.
%
% See also
%   steric2lonlatt, steric2gmsl, stericfilepath, stericfieldschema
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/05, OpenAI Codex

function info = stericproductinfo(product)
    arguments (Input)
        product (1, 1) StericProduct
    end
    info.product = char(product);
    info.path = stericfilepath(product);
    if ~isfile(info.path)
        error('ULMO:LoadData:FileNotFound', 'Missing steric input: %s', info.path);
    end
    % Inspect variable names and sizes without loading the large grid arrays.
    % Numerical geometry, support and layer closure are checked by steric2gmsl.
    metadata = whos('-file', info.path);
    info.fields = {metadata.name}';
    nonempty = {metadata(arrayfun(@(x) all(x.size > 0), metadata)).name};
    required = {'lon', 'lat', 'dates', 'stericSl'};
    if ~all(ismember(required, nonempty))
        error('ULMO:LoadData:VariableNotFound', 'Native coordinates, dates and stericSl are required.');
    end
    schema = stericfieldschema(metadata);
    if schema.incompleteSplit || schema.incompatibleSplits
        error('ULMO:stericproductinfo:InvalidSchema', ...
            'Supplied components require complete, compatible full/shallow/deep fields.');
    end
    info.hasShallowField = schema.available(1, 2);
    info.hasDeepField = schema.available(1, 3);
    info.hasDeep = schema.split(1);
    info.hasHalosteric = schema.available(3, 1);
end
