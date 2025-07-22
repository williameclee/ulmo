%% GIA2SLEPT
% Computes the GIA correction time series projected onto the Slepian basis.
%
% Syntax
%   [date, slept] = gia2slept(model, domain)
%   [date, slept] = gia2slept(model, domain, L)
%   [date, slept] = gia2slept(model, r, L, phi, theta, omega)
%   [date, slept] = gia2slept(days, __)
%   [date, slept] = gia2slept(time, __)
%   [date, slept] = gia2slept(__, 'Name', Value)
%   [date, slept, sleptU, sleptL, total] = gia2slept(__)
%
% Input arguments
%   model - Name of the GIA model
%       - 'Steffen_ice6g_vm5a': A model computed by H. Steffen using the
%           ice6g ice model and vm5a viscosity profile. Other models from
%           this dataset are also available and use the original naming
%           scheme. For example, 'Steffen_anu-ice_i72'.
%           This family of models can also be specified as a cell array,
%           e.g. {'Steffen', 'ice6g', 'vm5a'}.
%       - 'Paulson07': A model based on the ICE-5G ice load model of
%           Peltier (2004). Suitable for both Antarctica and Greenland. As
%           corrected by Geruo A and J. Wahr.
%           Please avoid using this model for oceans.
%       - 'Wangetal08': A model based on the older ICE-4G ice model, and
%           viscosity which varies laterally. Suitable for Greenland.
%       - 'IJ05_R2': A model based on the Ivins et al. (2013) ice model.
%           Suitable for Antarctica.
%       - 'IJ05': A model based on the Ivins and James (2005) ice model.
%           Suitable for Antarctica.
%       - 'W12a_v1': A 'best' model from Whitehouse et al. (2012). Suitable
%           only for Antarctica.
%       - 'LM17.3'
%       - 'Caron18' and 'Caron19'
%       The input can also be the path to the model file.
%		The default model is 'Steffen_ice6g_vm5a'.
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
%   L - The bandwith of the spherical harmonic coefficients
%       - A scalar L, denoting the maximum degree of the SH coefficients.
%       - A 1-by-2 vector [Lmin, Lmax], denoting the minimum and maximum
%           degree of the SH coefficients.
%   days - Number of days to calculate the GIA change for
%   time - Vector of dates to calculate the GIA change for
%		The input can be in DATENUM or DATETIME format.
%   L - Maximum degree of the GIA model
%		If empty, the model is not truncated.
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the centre of the tapers in degrees
%       The default values are 0.
%   BeQuiet - Whether to surpress output messages
%		The default option is false.
%
% Output arguments
%   time - Time vector
%   slept - GIA correction time series projected onto the Slepian basis
%   sleptU, sleptL - Upper and lower bounds of the GIA correction time
%       series, if available
%   total - Total GIA correction
%
% Notes
%   By default, if the output unit is geoid (POT), it is shifted vertically to represent the sea surface height change to conserve the ocean mass. To remove such effect, simply set the degree-0 coefficient to zero.
%
% Last modified by
%   2025/07/21, williameclee@arizona.edu (@williameclee)

function varargout = gia2slept(varargin)
    %% Initialisation
    % Parse inputs
    [date, model, domain, L, phi, theta, omega, units, beQuiet] = ...
        parseinputs(varargin{:});

    %% Loading the model
    warning('off', 'SLEPIAN:gia2plmt:noBoundsToReturn');
    [plm, plmU, plmL] = gia2plmt( ...
        [], model, L, "Unit", units, "BeQuiet", beQuiet);

    switch units
        case {'massdensity', 'SD'}
        case {'geoid', 'POT'}

            try
                oceanDomain = GeoDomain('alloceans', "Buffer", 0.5);
                zplm = giaz2plmt(model, L, "BeQuiet", beQuiet);
                plmLcl = localise(plm, oceanDomain, L);
                zplm = localise(zplm, oceanDomain, L);
                zOffset = plmLcl(1, 3) - zplm(1, 3) / 1e3;
                plm(1, 3) = plm(1, 3) - zOffset / oceanDomain.SphArea;
            catch
                warning('Failed to apply the degree-0 correction for geoid unit');
            end

    end

    hasBounds = ~isempty(plmU) && ~isempty(plmL);

    L = conddefval(L, max(plm(:, 1)));

    %% Computing the basis
    [falpha, falphaU, falphaL, N, G] = findslepianbasis( ...
        plm, plmU, plmL, domain, phi, theta, omega, L, hasBounds, beQuiet);

    %% Getting the trend
    if isempty(date) || isscalar(date)
        % If time is scalar, interpret it as the day change
        if isempty(date)
            date = datetime;
            deltaYear = 1;
        else
            deltaYear = date / days(years(1));
        end

    else
        deltaYear = years(date - date(1));
    end

    slept = deltaYear(:) * falpha(:)';

    if hasBounds
        sleptU = deltaYear(:) * falphaU(:)';
        sleptL = deltaYear(:) * falphaL(:)';
    end

    %% Getting the total
    truncation = ceil(N);

    if isa(domain, 'GeoDomain') || ismatrix(domain)
        eigfunINT = integratebasis_new( ...
            G, domain, truncation, "BeQuiet", beQuiet);
    else
        eigfunINT = integratebasis_new( ...
            G, domain, truncation, phi, theta, "BeQuiet", beQuiet);
    end

    eigfunINT = eigfunINT(1:truncation);

    switch units
        case {'massdensity', 'SD'}
            eigfunINT = eigfunINT * (4 * pi * 6370e3 ^ 2) / 1e12;
        case {'geoid', 'POT'}
            eigfunINT = eigfunINT * 1e3 / domain.SphArea;
    end

    total = slept(:, 1:truncation) * eigfunINT';

    %% Collecting outputs and plotting
    if ~hasBounds
        sleptU = nan;
        sleptL = nan;
    end

    if isscalar(date)
        varargout = {slept, sleptU, sleptL, total};
    else
        varargout = {date, slept, sleptU, sleptL, total};
    end

    if nargout > 0
        return
    end

    try
        plotgiamap(model, plm, falpha, deltaYear, domain, L);
    catch
        warning('Failed to plot the GIA map');
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Fallback values
    modelD = 'Steffen_ice6g_vm5a';
    domainD = {'greenland', 0.5};
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
        @(x) isnumeric(x) || isdatetime(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'Domain', domainD, ...
        @(x) ischar(x) || isstring(x) || iscell(x) ...
        || isa(x, 'GeoDomain') || isnumeric(x) || isempty(x));
    addOptional(p, 'L', [], ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'phi', phiD, @(x) isnumeric);
    addOptional(p, 'theta', thetaD, @(x) isnumeric);
    addOptional(p, 'omega', omegaD, @(x) isnumeric);
    addParameter(p, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'massdensity', 'geoid', 'SD', 'POT'})));
    addParameter(p, 'BeQuiet', 0.5, @(x) islogical(x) || isnumeric(x));

    parse(p, varargin{:});
    time = p.Results.Time(:);
    model = conddefval(p.Results.Model, modelD);
    domain = conddefval(p.Results.Domain, domainD);
    L = p.Results.L;
    phi = conddefval(p.Results.phi, phiD);
    theta = conddefval(p.Results.theta, thetaD);
    omega = conddefval(p.Results.omega, omegaD);
    units = p.Results.Unit;
    beQuiet = uint8(double(p.Results.BeQuiet) * 2);

    if isnumeric(time)
        time = datetime(time, 'ConvertFrom', 'datenum');
    end

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    varargout = {time, model, domain, L, phi, theta, omega, units, beQuiet};

end

function varargout = findslepianbasis(plm, plmU, plmL, domain, phi, theta, omega, ...
        L, hasBounds, beQuiet)
    falphaU = nan;
    falphaL = nan;

    if isa(domain, 'GeoDomain') || ismatrix(domain)
        [falpha, ~, N, ~, G] = plm2slep_new( ...
            plm, domain, L, "BeQuiet", beQuiet);

        if hasBounds
            falphaU = plm2slep_new( ...
                plmU, domain, L, "BeQuiet", beQuiet);
            falphaL = plm2slep_new( ...
                plmL, domain, L, "BeQuiet", beQuiet);
        end

    else
        [falpha, ~, N, ~, G] = plm2slep_new( ...
            plm, domain, L, phi, theta, omega, "BeQuiet", beQuiet);

        if hasBounds
            falphaU = plm2slep_new(plmU, domain, L, ...
                phi, theta, omega, "BeQuiet", beQuiet);
            falphaL = plm2slep_new(plmL, domain, L, ...
                phi, theta, omega, "BeQuiet", beQuiet);
        end

    end

    varargout = {falpha, falphaU, falphaL, N, G};

end

function plotgiamap(model, plm, slept, deltaYear, domain, L)
    %% Preparing the data
    meshSize = 1;
    mesh = plm2xyz(plm, meshSize, "BeQuiet", true);
    plmLcl = slep2plm_new(slept, domain, L, "BeQuiet", true);
    [meshLcl, lon, lat] = plm2xyz(plmLcl, meshSize, "BeQuiet", true);

    deltaYear = deltaYear(end);
    mesh = mesh * deltaYear;
    meshLcl = meshLcl * deltaYear;

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
        deltaYear, upper(mfilename)))
    clf

    subplot(1, 2, 1)
    title(sprintf('Model: %s (global)', model))
    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", 'GIA change [kg/m^2]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, mesh, cLevels, 'LineStyle', 'none')
    plotqdm(coastLonlat, 'k');
    plotqdm(domainLonlat, 'k', 'LineWidth', 1);
    hold off

    subplot(1, 2, 2)
    title(sprintf('Model: %s (localised)', model))
    loadcbar(cLim, cStep, ...
        "Title", 'GIA change [kg/m^2]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, meshLcl, cLevels, 'LineStyle', 'none')
    plotqdm(coastLonlat, 'k');
    plotqdm(domainLonlat, 'k', 'LineWidth', 1);
    hold off
end
