%% GIA2SLEPT
% Computes the vertical displacement provided by a GIA model projected onto
% the Slepian basis.
%
% Syntax
%   slep = giaz2slept(model, [r, phi, theta, omega], L)
%   slep = giaz2slept(model, domain, L)
%   slept = giaz2slept(years, __)
%   slept = giaz2slept(time, __)
%   slept = giaz2slept(__, 'Name', Value)
%   [slept, total] = giaz2slept(__)
%
% Input arguments
%   model - Name of the GIA model
%       - Name of a model computed by H. Steffen. It can be specified as,
%           e.g. 'Steffen_ice6g_vm5a' or {'Steffen', 'ice6g', 'vm5a'}.
%       Other models are specified in the same way.
%       - The input can also be the path to the model file.
%		The default model is 'Steffen_ice6g_vm5a'.
%       Format: string or 1 x 3 cell array.
%   L - The bandwith of the spherical harmonic coefficients
%       - A scalar L, denoting the maximum degree of the SH coefficients.
%       - A 1-by-2 vector [Lmin, Lmax], denoting the minimum and maximum
%           degree of the SH coefficients.
%       The default value is 18.
%       Format: numeric scalar or 1 x 2 vector.
%   years - Number of years to calculate the GIA change for
%		Format: numeric scalar.
%   time - Vector of dates to calculate the GIA change for
%		Format: numeric (datenum), datetime, or duration vector.
%   r - The angular extent of the spherical cap radius in degrees
%   domain - The domain of interest
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A cell array of the domain of interest (string) and the buffer
%           value (scalar): {domain, buf}.
%       - An N-by-2 matrix of the domain of interest in [lon, lat] format.
%       - A GeoDomain object.
%       The default domain is the global ocean with a buffer of 1 degree.
%       Format: string, cell array, numeric matrix, or GeoDomain object.
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the centre of the tapers in degrees
%       The default values are 0.
%       Format: scalar.
%   BeQuiet - Whether to surpress output messages
%		The default option is true.
%       Format: logical scalar.
%   MakePlot - Whether to plot the GIA change
%		The default option is true if no output is requested.
%       Format: logical scalar.
%
% Output arguments
%   slep - The GIA displacement projected onto the Slepian basis
%       Format: (L+1)^2 x 1 double array.
%   slept - The GIA displacement time series projected onto the Slepian
%       basis
%       Format: (L+1)^2 x ndates double array.
%   total - Total GIA displacment over the domain
%       The unit is m^3.
%       Format: numeric 1 x ndates vector.
%
% See also
%   GIA2PLMT, GIAZ2PLMT
%
% Data source
%   Steffen, H. (2021). Surface Deformations from Glacial Isostatic
%       Adjustment Models with Laterally Homogeneous, Compressible Earth
%       Structure (1.0) [Dataset]. Zenodo. doi: 10.5281/zenodo.5560862.
%
% Last modified by
%   2025/06/03, williameclee@arizona.edu (@williameclee)

function varargout = giaz2slept(varargin)
    %% Initialisation
    % Parse inputs
    [model, L, dYear, domain, beQuiet, makePlot] = ...
        parseinputs(varargin{:});

    %% Loading the model
    % Get the yearly trend
    [wPlm, wUPlm, wLPlm] = giaz2plmt(model, L);
    [wSlep, ~, N, ~, G] = ...
        plm2slep_new(wPlm, domain, L, "BeQuiet", beQuiet);
    wSlept = wSlep * dYear;

    if ~isempty(wUPlm) && ~isempty(wLPlm)
        wUSlep = plm2slep_new(wUPlm, domain, L, "BeQuiet", beQuiet);
        wLSlep = plm2slep_new(wLPlm, domain, L, "BeQuiet", beQuiet);
        wUSlept = wUSlep * dYear;
        wLSlept = wLSlep * dYear;
    else
        wUSlept = [];
        wLSlept = [];
    end

    varargout = {wSlept, wUSlept, wLSlept};

    %% Getting the total
    if nargout >= 4
        truncation = round(N);

        eigfunInt = integratebasis_new( ...
            G, domain, truncation, "BeQuiet", beQuiet);
        eigfunInt = eigfunInt(:);
        eigfunInt = eigfunInt * (4 * pi * 6371e3 ^ 2) * 1e-12; % Convert to m^3
        total = eigfunInt' * wSlept(1:truncation, :);
        total = mass2weq(total(:), domain);

        varargout = {wSlept, wUSlept, wLSlept, total};
    end

    if nargout > 0 || ~makePlot
        return
    end

    wPlmt = wPlm;
    wPlmt(:, 3:4) = wPlm(:, 3:4) * dYear(end);
    plotdispmap(model, wPlmt, wSlept, dYear, domain, L);

end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Fallback values
    modelD = 'Steffen_ice6g_vm5a';
    domainD = GeoDomain('oceans', "Buffer", 1, "DefaultParams", true);
    LD = 18;
    phiD = 0;
    thetaD = 0;
    omegaD = 0;

    % Allow skipping the time argument
    if nargin > 0 && ...
            (ischar(varargin{1}) || isstring(varargin{1}) || ...
            iscell(varargin{1}))
        varargin(2:end + 1) = varargin;
        varargin{1} = [];
    end

    % Parsing inputs
    p = inputParser;
    addOptional(p, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isduration(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'Domain', domainD, ...
        @(x) ischar(x) || isstring(x) || iscell(x) ...
        || isa(x, 'GeoDomain') || isnumeric(x) || isempty(x));
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'phi', phiD, @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'theta', thetaD, @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'omega', omegaD, @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'BeQuiet', true, @istruefalse);
    addParameter(p, 'MakePlot', true, @istruefalse);

    parse(p, varargin{:});
    time = p.Results.Time(:);
    model = conddefval(p.Results.Model, modelD);
    domain = conddefval(p.Results.Domain, domainD);
    L = conddefval(p.Results.L, LD);
    phi = conddefval(p.Results.phi, phiD);
    theta = conddefval(p.Results.theta, thetaD);
    omega = conddefval(p.Results.omega, omegaD);
    beQuiet = logical(p.Results.BeQuiet);
    makePlot = logical(p.Results.MakePlot);

    % Format the time into change in years
    if isempty(time)
        dYear = 1;
    elseif isscalar(time) && isnumeric(time)
        dYear = time;
    elseif isvector(time) && isnumeric(time)
        dYear = (time - time(1)) / days(years(1));
    elseif isvector(time) && isdatetime(time)
        dYear = years(time - time(1));
    elseif isduration(time)
        dYear = years(time);
    end

    dYear = dYear(:)';

    % Change the domain to a GeoDomain object if appropriate
    if isnumeric(domain) && isscalar(domain)
        domain = [domain, phi, theta, omega];
    elseif ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    varargout = {model, L, dYear, domain, beQuiet, makePlot};

end

function plotdispmap(model, plm, slept, dYear, domain, L)
    %% Preparing the data
    % Get the change rate
    if ~isscalar(dYear)
        slept = squeeze(slept(:, end));
        dYear = dYear(end);
    end

    meshSize = 1;
    mesh = plm2xyz(plm, meshSize, "BeQuiet", true);
    plmLcl = slep2plm_new(slept, domain, L, "BeQuiet", true);
    [meshLcl, lon, lat] = plm2xyz(plmLcl, meshSize, "BeQuiet", true);

    dYear = dYear(end);
    mesh = mesh * dYear;
    meshLcl = meshLcl * dYear;

    coastLonlat = gshhscoastline('c', 'LonOrigin', 180, "BeQuiet", true);
    domainLonlat = domain.Lonlat('LonOrigin', 180);

    [cLim, cStep] = optimalclim(meshLcl, 'Percentile', 1);
    mesh = max(min(mesh, cLim(2)), cLim(1));
    meshLcl = max(min(meshLcl, cLim(2)), cLim(1));

    %% Plotting
    % Protect underscore in model name
    if iscell(model)
        model = strjoin(model, '_');
    end

    model = strrep(model, '_', '\_');

    figure(999)
    set(gcf, "NumberTitle", 'off', "Name", ...
        sprintf('Localised GIA change in %.1f year(s) (%s)', ...
        dYear, upper(mfilename)))
    clf

    subplot(2, 1, 1)
    title(sprintf('Model: %s (global)', model))
    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", 'Vertical displacement [mm]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, mesh, cLevels, 'LineStyle', 'none')
    plotqdm(coastLonlat, 'k');
    plotqdm(domainLonlat, 'k', 'LineWidth', 1);
    hold off

    subplot(2, 1, 2)
    title(sprintf('Model: %s (localised)', model))
    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", 'Vertical displacement [mm]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, meshLcl, cLevels, 'LineStyle', 'none')
    plotqdm(coastLonlat, 'k');
    plotqdm(domainLonlat, 'k', 'LineWidth', 1);
    hold off
end
