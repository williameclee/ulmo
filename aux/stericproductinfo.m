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
%       hasDeepField, hasDeep, hasFullDepth, nativeDepth and hasHalosteric.
%       nativeDepth is full-depth or shallow. hasFullDepth
%       describes the prepared product's coverage; the other flags describe
%       supported fields. For upper-ocean products, the unprefixed native fields
%       supply shallow coverage. A depth split requires
%       explicit full-depth, shallow and deep grids. Empty grids are unsupported.
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
%   2026/11/20, En-Chi Lee (williameclee@gmail.com)

function info = stericproductinfo(product)

    arguments (Input)
        product (1, 1) StericProduct
    end

    info.product = string(product);
    info.hasFullDepth = product.supportsFullDepth();
    info.nativeDepth = product.nativeDepth();
    info.path = stericfilepath(product);

    if ~isfile(info.path)
        error('ULMO:LoadData:FileNotFound', 'Missing steric input file: %s', info.path);
    end

    % Inspect variable names and sizes without loading the large grid arrays.
    % Numerical geometry, support and layer closure are checked by steric2gmsl.
    metadata = whos('-file', info.path);
    info.fields = {metadata.name}';
    nonempty = {metadata(arrayfun(@(x) all(x.size > 0), metadata)).name};
    required = {'lon', 'lat', 'dates', 'stericSl'};

    if ~all(ismember(required, nonempty))
        missingVars = required(~ismember(required, nonempty));
        error('ULMO:LoadData:VariableNotFound', ['Native coordinates, dates and stericSl are required. ', ...
                  'But file %s does not contain or has empty variables(s) for: %s'], ...
            info.path, strjoin(missingVars));
    end

    schema = stericfieldschema(metadata, info.nativeDepth);

    if schema.incompleteSplit || schema.incompatibleSplits
        error('ULMO:stericproductinfo:InvalidSchema', ...
        'Supplied components require complete, compatible full-depth/shallow/deep fields.');
    end

    info.hasShallowField = schema.available(1, 2);
    info.hasDeepField = schema.available(1, 3);
    info.hasDeep = schema.split(1);
    info.hasHalosteric = schema.available(3, schema.depths == info.nativeDepth);
end
