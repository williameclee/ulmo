%% INTEGRATEBASIS
% Accepts a Slepian basis and integrates the functions within a region.
% The integrals reported are in fractional sphere area. You may want to
% convert this to real sphere area via x(4*pi*R^2).
%
% Syntax
%   eigfunINT = INTEGRATEBASIS(eigfun, domain, J)
%   eigfunINT = INTEGRATEBASIS(eigfun, r, J, phi, theta)
%
% Input arguments
%   eigfun - The eigfun
%       This can be either the G you get from GLMALPHA or a cell array of
%       the functions as in LOCALIZATION.
%   domain - The region
%       Can be formatted multiple ways as in GLMALPHA
%   J - How many functions you want to do this for
%       The default is all of them
%   phi, theta - Longitude and colatitude of the centre in degrees
%       There is no omega here because a rotation of the polar cap
%       functions does not affect their integration
%
% Output arguments
%   eigfunINT - The integrals of the eigfun over the region
%       These are in real sphere area, depending on your radius.
%
% See also
%   PLM2AVG, PLM2AVGP
%
% Note
%   The function is not properly tested yet.
%
% Last modified by
%   2025/10/16, williameclee@arizona.edu (@williameclee)
%   2024/08/15, williameclee@arizona.edu (@williameclee)
%   2016/11/01, charig@email.arizona.edu (@harig00)

function varargout = integratebasis_new(varargin)
    [eigfun, domain, J, phi, theta, forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin{:});

    %% Initialisation
    % Sort out what CC is
    if ischar(eigfun) || isstring(eigfun)
        % Evaluate the specified expression
        eval(eigfun);
    elseif isnumeric(eigfun)
        % We must have a G matrix from glmalpha (already sorted)
        % Reorder them into a cell array
        [n, m] = size(eigfun);
        J = conddefval(J, m);
        L = sqrt(n) - 1;
        % This should be an integer
        if (floor(L) ~= L)
            error('Something fishy about your L');
        end

        [~, ~, ~, lmcosi, ~, ~, ~, ~, ~, ronm] = addmon(L);
        % Collect the eigenvector output into a format that PLM2XYZ knows
        % how to interpret

        for j = 1:J
            % Create the blanks
            cosi = lmcosi(:, 3:4);
            % Stick in the coefficients of the 1st eigentaper
            cosi(ronm) = eigfun(:, j);
            % Construct the full matrix
            CC2{j} = [lmcosi(:, 1:2) cosi];
        end

        eigfun = CC2;
    elseif iscell(eigfun)
        % We are all good, just get the size.
        m = size(eigfun, 2);
        J = conddefval(J, m);
        L = addmup(size(eigfun{1}, 1), 'r');
    else
        error('What format is CC?');
    end

    %% Find file
    domainType = domaintype(domain);

    if strcmp(domainType, 'geographic')
        dataFile = getoutputfile(domain, L);

        if ~forceNew && exist(dataFile, 'file')
            load(dataFile, 'eigfunINT');

            if ~beQuiet
                fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">integrated basis functions</a>.\n', ...
                    callchaintext(callChain), dataFile, dataFile);
            end

            if ~isempty(J)

                try
                    eigfunINT = eigfunINT(1:ceil(J));
                    varargout = {eigfunINT};
                    return
                catch
                end

            end

        end

    end

    %% Main
    switch domaintype(domain)
        case 'geographic'
            % Default parameters
            defval('buf', 0);
            defval('pars', 0);
            rotb = 0;

            % Rotate if needed
            try
                rotb = feval(domain.Domain, 'rotated');
            catch
            end

            if isscalar(rotb) && rotb
                % Return the coordinates and do the rotation
                try
                    [XY, lonc, latc] = domain.Lonlat;
                catch
                    [XY, lonc, latc] = feval(domain, pars, buf);
                end

                [thetap, phip, ~] = ...
                    rottp(deg2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), ...
                    deg2rad(-lonc), deg2rad(latc), 0);
                lonp = rad2deg(phip);
                latp = 90 - rad2deg(thetap);
                [latf, lonf] = flatearthpoly(latp, lonp);
                XY = [lonf, latf];
            else

                try
                    XY = domain.Lonlat;
                catch
                    XY = feval(domain.Domain, ...
                        domain.Upscale, domain.Buffer);
                end

            end

        case 'lonlat'
            % We have straight coordinates
            XY = domain;

        case 'polar'
            XY = [(1:360)', repmat(90 - domain, 360, 1)];
            [thetap, phip, ~] = ...
                rottp(deh2rad(90 - XY(:, 2)), deg2rad(XY(:, 1)), ...
                0, deg2rad(-theta), deg2rad(-phi));
            lonp = rad2deg(phip);
            latp = 90 - rad2deg(thetap);
            XY = [lonp, latp];
    end

    % Run parallel if possible
    if exist('eigfunINT', 'var')

        parfor h = length(eigfunINT) + 1:J
            eigfunINT(h) = plm2avg(eigfun{h}, XY);
        end

    else

        parfor h = 1:J
            eigfunINT(h) = plm2avg(eigfun{h}, XY);
        end

    end

    %% Returning requested outpput
    if strcmp(domainType, 'geographic') && saveData
        save(dataFile, 'eigfunINT');

        if ~beQuiet
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">integrated basis functions</a>.\n', ...
                callchaintext(callChain), dataFile, dataFile);
        end

    end

    varargout = {eigfunINT};
end

%% Subfunctions
function varargout = parseinputs(varargin)
    dfOpt.eigfun = '[~, eigfun] = localization(15, domain);';
    dfOpt.domain = 'africa';
    dfOpt.truncation = [];
    dfOpt.phi = 0;
    dfOpt.theta = 0;

    ip = inputParser;
    addOptional(ip, 'eigfun', dfOpt.eigfun);
    addOptional(ip, 'Domain', dfOpt.domain, ...
        @(x) ischar(x) || isstring(x) || iscell(x) ...
        || isa(x, 'GeoDomain') || isnumeric(x) || isempty(x));
    addOptional(ip, 'J', dfOpt.truncation);
    addOptional(ip, 'phi', dfOpt.phi);
    addOptional(ip, 'theta', dfOpt.theta);
    addOptional(ip, 'MoreRegionSpecs', {});
    addParameter(ip, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'SaveData', true, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'CallChain', {}, @iscell);
    parse(ip, varargin{:});

    eigfun = conddefval(ip.Results.eigfun, dfOpt.eigfun);
    domain = conddefval(ip.Results.Domain, dfOpt.domain);
    J = conddefval(ip.Results.J, dfOpt.truncation);
    phi = conddefval(ip.Results.phi, dfOpt.phi);
    theta = conddefval(ip.Results.theta, dfOpt.theta);
    moreRegionSpecs = ip.Results.MoreRegionSpecs;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = logical(ip.Results.BeQuiet);
    callChain = [ip.Results.CallChain, {mfilename}];

    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, moreRegionSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = GeoDomain( ...
            domain{1}, "Buffer", domain{2}, moreRegionSpecs{:});
    elseif iscell(domain)
        domain = GeoDomain(domain{:}, moreRegionSpecs{:});
    end

    if isnumeric(domain) && length(domain) == 3
        phi = domain(2);
        theta = domain(3);
        domain = domain(1);
    end

    varargout = ...
        {eigfun, domain, J, phi, theta, forceNew, saveData, beQuiet, callChain};
end

function domainType = domaintype(domain)

    if isnumeric(domain)

        if isscalar(domain)
            % We have a polar cap
            domainType = 'polar';
        else
            % We have straight coordinates
            domainType = 'lonlat';
        end

    elseif isa(domain, 'GeoDomain')
        % We have a geographic region
        domainType = 'geographic';
    else
        error('What format is your domain? %s', class(domain));
    end

end

function filePath = getoutputfile(domain, L)
    dataFolder = fullfile(getenv('IFILES'), 'EIGFUNINT');

    if ~exist(dataFolder, 'dir')
        mkdir(dataFolder);
    end

    fileName = sprintf('EIGFUNINT-%s-%d.mat', domain.Id, L);

    filePath = fullfile(dataFolder, fileName);
end
