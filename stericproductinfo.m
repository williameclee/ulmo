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
%   steric2lonlatt, steric2gmsl, stericfilepath
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

function info = stericproductinfo(product)
    arguments (Input)
        product (1, 1) StericProduct
    end
    info.product = char(product);
    info.path = stericfilepath(product);
    if ~isfile(info.path)
        error('ULMO:LoadData:FileNotFound', 'Missing steric input: %s', info.path);
    end
    metadata = whos('-file', info.path);
    info.fields = {metadata.name}';
    nonempty = {metadata(arrayfun(@(x) all(x.size > 0), metadata)).name};
    required = {'lon', 'lat', 'dates', 'stericSl'};
    if ~all(ismember(required, nonempty))
        error('ULMO:LoadData:VariableNotFound', 'Native coordinates, dates and stericSl are required.');
    end
    names = ["stericSl", "shallowStericSl", "deepStericSl"; ...
             "thermostericSl", "shallowThermostericSl", "deepThermostericSl"; ...
             "halostericSl", "shallowHalostericSl", "deepHalostericSl"];
    available = ismember(names, nonempty);
    split = any(available(:, 2:3), 2);
    if any(split & ~all(available, 2)) || ...
            (available(3, 1) && split(1) ~= split(3))
        error('ULMO:stericproductinfo:InvalidSchema', ...
            'Supplied components require complete, compatible full/shallow/deep fields.');
    end
    info.hasShallowField = available(1, 2);
    info.hasDeepField = available(1, 3);
    info.hasDeep = split(1);
    info.hasHalosteric = available(3, 1);
end
