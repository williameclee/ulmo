%% MASS2WEQ
% Converts mass to water equivalent
%
% Syntax
%   weq = mass2weq(mass, domain)
%   weq = mass2weq(__, density, unit)
%   weq = mass2weq(__, density, unit, inputUnit)
%   [weq, convFac] = mass2weq(__)
%
% Input arguments
%   mass - Mass
%       The unit is specified by the 'inputUnit' parameter.
%   domain - Domain over which the mass is defined
%       - GeoDomain object
%       - Numeric scalar specifying the area in m^2
%       - Nx2 numeric array of (lon, lat) vertices defining a polygon
%   density - Density of the material (default: 'water')
%       - Numeric scalar specifying the density in kg/m^3
%       - 'water' (1000 kg/m^3)
%       - 'seawater' (1028 kg/m^3)
%   unit - Unit of the output water equivalent (default: 'mm')
%       - 'mm' (millimeters)
%       - 'cm' (centimeters)
%       - 'm' (meters)
%   inputUnit - Unit of the input mass (default: 'Gton')
%       - 'Gton' (gigatons, 1 Gton = 1e12 kg)
%       - 'kg' (kilograms)
%
% Output arguments
%   weq - Mass in water equivalent
%   convFac - Conversion factor from mass to water equivalent
%
% Last modified by
%   2025/11/03, williameclee@arizona.edu (@williameclee)

function varargout = mass2weq(varargin)
    %% Initialisation
    % Parse inputs
    [massGton, domain, density, unit, inputUnit] = ...
        parseinputs(varargin{:});

    if isnumeric(density)
        waterDensity = density;
    else

        switch density
            case 'water'
                waterDensity = 1;
            case 'seawater'
                waterDensity = 1.028;
        end

    end

    switch unit
        case 'mm'
            unitConvFac = 1;
        case 'cm'
            unitConvFac = 10;
        case 'm'
            unitConvFac = 1000;
    end

    switch inputUnit
        case 'Gton'
            inputUnitConvFac = 1e12;
        case 'kg'
            inputUnitConvFac = 1;
    end

    %% Computing the mass in water equivalent
    if isa(domain, 'GeoDomain')
        area = domain.Area;
    elseif isnumeric(domain) && isscalar(domain)
        area = domain;
    elseif isnumeric(domain) && ismatrix(domain)
        sphArea = spharea(domain);
        area = sphArea * (4 * pi * 6370e3 ^ 2);
    end

    % Convert unit from Gton to kg/m^2 to water equivalent
    convFac = inputUnitConvFac / area / waterDensity * unitConvFac;

    if iscell(massGton)
        massWeq = cellfun(@(x) x * convFac, massGton, "UniformOutput", false);
    else
        massWeq = massGton * convFac;
    end

    %% Collecting outputs
    if iscell(massGton)
        varargout = [massWeq(:)', {convFac}];
    else
        varargout = {massWeq, convFac};
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addRequired(p, 'Mass', ...
        @(x) isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x))));
    addRequired(p, 'Domain', ...
        @(x) ischar(x) || isstring(x) || iscell(x) ...
        || isa(x, 'GeoDomain') || isnumeric(x));
    addOptional(p, 'Density', 'water', ...
        @(x) isnumeric(x) || ((ischar(x) || isstring(x)) && ...
        ismember(x, {'water', 'seawater'})));
    addOptional(p, 'Unit', 'mm', @(x) (ischar(x) || isstring(x)) ...
        && ismember(x, {'mm', 'cm', 'm'}));
    addOptional(p, 'InputUnit', 'Gton', @(x) (ischar(x) || isstring(x)) ...
        && ismember(x, {'Gton', 'kg'}));

    parse(p, varargin{:});

    massGton = p.Results.Mass;
    domain = p.Results.Domain;
    density = p.Results.Density;
    unit = p.Results.Unit;
    inputUnit = p.Results.InputUnit;

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    varargout = {massGton, domain, density, unit, inputUnit};
end
