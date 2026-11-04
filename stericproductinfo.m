function info = stericproductinfo(product)
%STERICPRODUCTINFO Inspect depth support in the standard native steric file.
% An empty deep field explicitly marks unsupported depth coverage.
% A shallow/full pair supports deriving deep steric. An explicit deep field
% also supports deriving shallow steric from full minus deep. Missing halo
% fields are malformed inputs, not evidence that a product is shallow-only.
    arguments
        product {mustBeTextScalar}
    end
    product = char(product);
    info.product = product;
    info.path = fullfile(getenv('IFILES'), 'STERIC', product, ...
        [product, '-StericSeaLevel.mat']);
    if ~isfile(info.path)
        error('ULMO:LoadData:FileNotFound', 'Missing steric input: %s', info.path);
    end
    fields = whos('-file', info.path);
    info.fields = {fields.name}';
    nonempty = {fields(arrayfun(@(x) all(x.size > 0), fields)).name};
    required = {'lon', 'lat', 'dates', 'stericSl'};
    info.hasShallowField = ismember('shallowStericSl', nonempty);
    info.hasDeepField = ismember('deepStericSl', nonempty);
    if ismember('deepStericSl', info.fields)
        info.hasDeep = info.hasDeepField;
    else
        info.hasDeep = info.hasShallowField;
    end
    info.applyHaloCorrection = ~strcmpi(product, 'SIO');
    if info.applyHaloCorrection
        required{end+1} = 'halostericSl';
        if info.hasDeep && ~any(ismember( ...
                {'shallowHalostericSl', 'deepHalostericSl'}, nonempty))
            error('ULMO:LoadData:VariableNotFound', ...
                'Depth-resolved halosteric data required in %s.', info.path);
        end
    end
    missing = required(~ismember(required, nonempty));
    if ~isempty(missing)
        error('ULMO:LoadData:VariableNotFound', ...
            'Missing variables in %s: %s', info.path, strjoin(missing, ', '));
    end
end
