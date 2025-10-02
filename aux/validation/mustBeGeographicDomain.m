function mustBeGeographicDomain(domain)

    if isa(domain, "GeoDomain") || isa(domain, 'polyshape')
        return
    end

    if ismatrix(domain) && isnumeric(domain) && size(domain, 2) == 2
        return
    end

    if ischar(domain) || isstring(domain)

        if exist(domain, 'file') ~= 2
            error ('ULMO:notValidGeographicDomain', ...
            'String input for a geographic domain should be a valid function name.');
        end

        return
    end

    if iscell(domain)

        if numel(domain) > 3
            error ('ULMO:notValidGeographicDomain', ...
            'Cell array input for a geographic domain should have at most 3 elements.');
        end

        if ischar(domain{1}) || isstring(domain{1})

            if exist(domain{1}, 'file') ~= 2
                error ('ULMO:notValidGeographicDomain', ...
                'First element of cell array input for a geographic domain should be a valid function name.');
            end

        else
            error ('ULMO:notValidGeographicDomain', ...
            'First element of cell array input for a geographic domain should be a string.');
        end

        if any(~isnumeric(domain{2:end}))
            error ('ULMO:notValidGeographicDomain', ...
            'Second (and third) elements of cell array input for a geographic domain should be numeric.');
        end

    end

    error ('ULMO:notValidGeographicDomain', ...
    'Unrecognised format for a geographic domain.');

end
