%% GLMALPHA
% Returns an lm-by-alpha matrix with unit-normalized spherical harmonic
% coefficients of the BANDLIMITED or PASSBAND Slepian functions of the
% SINGLE or DOUBLE polar cap, or of a GEOGRAPHICAL region of interest.
% Only in the geographical case are the eigenvalues automatically sorted;
% if not, the column dimension is always block-ordered by virtue of the
% construction. The matrix G is orthogonal, G'*G is the identity. In column
% truncated form, G(:,1:J)*G(:,1:J)' is not the identity also, but rather
% a projection with eigenvalues 0 and 1.
%
% Synatx
%   [G, V, EL, EM, N] = glmalpha(r, L)
%   [G, V, EL, EM, N] = glmalpha(domain, L)
%   [G, V, EL, EM, N] = glmalpha(__, upscale)
%   [G, V, EL, EM, N] = glmalpha(__, 'Name', value)
%   [G, V, EL, EM, N, GM2AL, MTAP, IMTAP] = ...
%       glmalpha(r, L, numpoles, blox, upco, rescale, truncation, anti)
%   glmalpha(demoId)
%
% Input arguments
%   r - The angular extent of the spherical cap radius in degrees
%   domain - The domain of interest
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
%   sord - The type of the region
%       It takes on different values depending on the type of the domain.
%       - For the polar cap (i.e. the domain is a scalar), it takes on the
%           following values:
%           1: Single polar cap of radius TH.
%           2: Double polar cap, each of radius TH.
%           The default value is 1.
%       - For a geographical domain, it denotes the amount of upscaling
%           applied to the vertices of the domain.
%           The default value is 10.
%   The following options are only for axisymmetric polar caps:
%   blox - The ordering of the spherical harmonic coefficients
%       - 0: Standard (lm) row ordering, l = 0:L, m = -l:l as ADDMOUT.
%       - 1: Block-diagonal row ordering, m = [0 -1 1 -2 2 ... -L L],
%           l = abs(m):L.
%       The default value is 0.
%   upco - The fraction of unit radius for continuation
%       - +ve fraction of unit radius for upward continuation.
%       - -ve fraction of unit radius for downward continuation.
%       The default value is 0.
%   rescale - Whether to rescale the eigenfunctions
%       - 0: Not rescaled.
%       - 1: Rescaled to have a unit integral over the unit sphere (only
%           relevant for the down/upward continued functions).
%       The default value is 0.
%   There are three more options that hold only for the geographic
%   cases
%   truncation - The number of eigenfunctions saved and returned
%   anti - Whether to compute the domain or its complement.
%       - 0: Get exactly the region that you specify.
%       - 1: Get the opposite of the region you specify.
%       The default value is 0.
%   rotateBack - Whether to rotate the eigenfunctions back to the pole
%   There are a few more cosmetic options
%   ForceNew - Force the function to recomputed the data
%   SaveData - Save the data to a file
%   BeQuiet - Suppress the output messages
%
% Output arguments
%   G - The unitary matrix of localisation coefficients
%       Each column is the SH expansion of the corresponding eigenfunction
%   V - The eigenvalues of the Slepian functions
%       If the domain is geographic, V and G are sorted in descending order
%   EL - The vector with spherical harmonic degrees as first index of G
%   EM - The vector with spherical harmonic orders as first index of G
%   N - The Shannon number
%   GM2AL - The sum over all orders of the squared coefficients
%       The TOTAL power, NOT the power spectral density.
%   MTAP - The order of the eigentapers, if the region is axisymmetric
%   IMTAP - The rank within that particular order of the eigentapers
%
% Examples
%   The following examples illustrate the use of glmalpha:
%   >>  glmalpha('demo1') % Illustrates the block sorting and checks
%       unitarity
%   >>  glmalpha('demo2') % Makes a coupling kernel a la BCOUPLING
%   >>  glmalpha('demo3') % Calculates something and uses PLOTSLEP to plot
%   The following examples create the same localisation matrix a
%       geographic domain:
%   >> domain = GeoDomain('greenland, 'Upscale', 10, 'Buffer', 1);
%   >>  G = glmalpha(domain, 20);
%   >>  domain = {'greenland', 1};
%   >>  G = glmalpha(domain, 20, 10);
%
% See also
%   GLMALPHAPTO, ADDMOUT, ADDMON, KERNELC, GALPHA, DLMLMP, GLM2LMCOSI,
%   LOCALIZATION
%
% Notes
%   Region functions such as ANTARCTICA have a default behavior to indicate
%   if their eigenfunctions should be rotated (e.g. back to a pole). If you
%   want eigenfunctions for the region at the equator then rotate them back
%   after the fact using ROTATEGP.
%   Should be able to update this to retain the rank order per m as well as
%   the global ordering. Does this work for the whole-sphere? In that case,
%   should really want G to be the identity - all though of course,
%   anything works, too. You don't get necessarily the spherical harmonics
%   back...
%
% Last modified by
%   2025/06/02, williameclee@arizona.edu (@williameclee)
%   2017/12/01, fjsimons@alum.mit.edu (@fjsimons)
%   2016/06/27, charig@princeton.edu (@harig00)
%   2016/10/11, plattner@alumni.ethz.ch (@AlainPlattner)

function varargout = glmalpha_new(varargin)
    %% Initialisation & demos
    % demos
    if (ischar(varargin{1}) || isstring(varargin{1})) && ...
            contains(varargin{1}, 'demo')

        switch varargin{1}
            case 'demo1'
                glmalpha_demo1;
            case 'demo2'
                glmalpha_demo2;
            case 'demo3'
                glmalpha_demo3;
        end

        return

    end

    % Parse inputs
    [domain, L, sord, blox, upco, resc, truncation, anti, rotb, ...
         forceNew, saveData, beQuiet] = ...
        parseinputs(varargin);

    defval('mesg', 'GLMALPHA Check passed')
    % Hold all messages
    mesg = NaN;

    [lp, bp, maxL, ldim] = ldimension(L);
    truncation = conddefval(truncation, ldim);

    % Output file
    [dataPath, dataExists, GM2AL, MTAP, IMTAP, xver] = ...
        getoutputfile(domain, L, sord, blox, upco, resc, ...
        truncation, anti, bp, ldim, [], [], [], []);

    if dataExists && ~forceNew
        warning('off', 'MATLAB:load:variableNotFound')
        load(dataPath, 'G', 'V', 'EL', 'EM', 'N', 'GM2AL', 'MTAP', 'IMTAP')

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), dataPath)
        end

        varargout = {G, V, EL, EM, N, GM2AL, MTAP, IMTAP};

        if nargout > 0
            return
        end

        %% Plot eignvalue-weighted maps
        if ~beQuiet
            fprintf('%s generating eigenvalue-weighted map, this make take a while...\n', upper(mfilename))
        end

        plotvweightmap(G, V, domain)

        return

    end

    %% Computing the localisation matrix
    % Find row indices into G belonging to the orders
    [EM, EL, ~, blkm] = addmout(maxL);

    % For geographical or lonlat regions
    if ismatrix(domain) || isa(domain, 'GeoDomain')
        upscale = sord;
        [G, V, N] = glmalpha_geographic( ...
            maxL, domain, upscale, anti, rotb, ldim, bp, EL, EM, xver, ...
            beQuiet, forceNew, mesg);

        G = G(:, 1:truncation);
        V = V(1:truncation);

        try
            save(dataPath, '-v7.3', 'G', 'V', 'EL', 'EM', 'N')
        catch
            save(dataPath, 'G', 'V', 'EL', 'EM', 'N')
        end

    else
        % For AXISYMMETRIC REGIONS
        [G, V, EL, EM, N, GM2AL, MTAP, IMTAP] = glmalpha_axisymmetric( ...
            domain, sord, L, lp, bp, EM, EL, blkm, blox, upco, xver);

        if ~strcmp(dataPath, 'neveravailable') && saveData

            try
                save(dataPath, '-v7.3', ...
                    'G', 'V', 'EL', 'EM', 'N', 'GM2AL', 'MTAP', 'IMTAP')
            catch
                save(dataPath, ...
                    'G', 'V', 'EL', 'EM', 'N', 'GM2AL', 'MTAP', 'IMTAP')
            end

        end

    end

    %% Returning requested variables
    varargout = {G, V, EL, EM, N, GM2AL, MTAP, IMTAP};

    if nargout > 0
        return
    end

    %% Plot eignvalue-weighted maps
    if ~beQuiet
        fprintf('%s generating eigenvalue-weighted map, this make take a while...\n', ...
            upper(mfilename))
    end

    plotvweightmap(G, V, domain)
end

%% Subfunctions
function varargout = parseinputs(Inputs)
    domainD = 30;
    LD = 18;
    sordD = [];
    bloxD = 0;
    upcoD = 0;
    rescD = 0;
    JD = [];
    antiD = false;
    rotateBackD = true;

    p = inputParser;
    addOptional(p, 'Domain', domainD, ...
        @(x) isa(x, 'GeoDomain') || ... % GeoDomain object
        ischar(x) || isstring(x) || iscell(x) || ... % geogrpahic domain
        isnumeric(x) || isempty(x)); % polar cap or lonlat
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'sord', sordD, ...
        @(x) isnumeric(x) || iscell(x) || isempty(x));
    addOptional(p, 'blox', bloxD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'upco', upcoD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'resc', rescD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'J', JD, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'anti', antiD, @(x) islogical(x) || isnumeric(x));
    addOptional(p, 'RotateBack', rotateBackD, ...
        @(x) isnumeric(x) || islogical(x) || isempty(x));
    addParameter(p, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveData', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    parse(p, Inputs{:});

    domain = conddefval(p.Results.Domain, domainD);
    L = conddefval(p.Results.L, LD);
    sord = conddefval(p.Results.sord, sordD);
    blox = conddefval(p.Results.blox, bloxD);
    upco = conddefval(p.Results.upco, upcoD);
    resc = conddefval(p.Results.resc, rescD);
    J = conddefval(p.Results.J, JD);
    anti = logical(p.Results.anti);
    rotateBack = logical(conddefval(p.Results.RotateBack, rotateBackD));
    forceNew = p.Results.ForceNew;
    saveData = p.Results.SaveData;
    beQuiet = p.Results.BeQuiet;

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, "Upscale", sord);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2}, "Upscale", sord);
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{1}, "Upscale", sord, domain{2:end});
    end

    varargout = ...
        {domain, L, sord, blox, upco, resc, J, anti, rotateBack, ...
         forceNew, saveData, beQuiet};
end

function [dataPath, dataExists, GM2AL, MTAP, IMTAP, xver] = ...
        getoutputfile(domain, L, sord, blox, upco, resc, truncation, ...
        anti, bp, ldim, GM2AL, MTAP, IMTAP, xver)

    if upco == 0 && resc == 0
        xver = conddefval(xver, 0);

        if isnumeric(domain) && isscalar(domain)
            domainType = 'cap';
        elseif isa(domain, 'GeoDomain')
            domainType = 'GeoDomain';
        elseif ismatrix(domain)
            domainType = 'lonlat';
        else
            error('Unrecognised domain type.')
        end

        dataFolder = fullfile(getenv('IFILES'), 'GLMALPHA');
    else
        xver = conddefval(xver, 0);

        domainType = 'pto';
        dataFolder = fullfile(getenv('IFILES'), 'GLMALPHAPTO');
    end

    switch domainType
        case 'cap'
            sord = conddefval(sord, 1);

            if ~bp
                outputFile = sprintf('glmalpha-%i-%i-%i-%i.mat', ...
                    domain, L, sord, blox);
            else
                outputFile = sprintf('glmalphabl-%i-%i-%i-%i-%i.mat', ...
                    domain, L(1), L(2), sord, blox);
            end

            % Initialise ordering matrices
            MTAP = zeros([1, ldim]);
            IMTAP = zeros([1, ldim]);

        case 'GeoDomain'

            if ~bp
                outputFile = sprintf('glmalpha-%s-%i-%i.mat', ...
                    domain.Id, L, truncation);
            else
                outputFile = sprintf('glmalphabl-%s-%i-%i-%i.mat', ...
                    domain.Id, L(1), L(2), truncation);
            end

        case 'lonlat'

            if ~bp
                outputFile = sprintf('glmalpha-%s-%i-%i.mat', ...
                    hash(domain, 'sha1'), L, truncation);
            else
                outputFile = sprintf('glmalphabl-%s-%i-%i-%i.mat', ...
                    hash(domain, 'sha1'), L(1), L(2), truncation);
            end

        case 'pto'
            % Make a hash, who cares if it's human-readable?
            outputFile = sprintf('%s.mat', ...
                hash([domain, L, phi, theta, omega], 'sha1'));
            % For excessive verification of the upco'd case
    end

    if anti
        outputFile = [outputFile(1:end - 4), '-anti.mat'];
    end

    dataPath = fullfile(dataFolder, outputFile);
    dataExists = isfile(dataPath);
end

function [lp, bp, maxL, ldim] = ldimension(L)
    % Figure out if it's lowpass or bandpass
    lp = isscalar(L);
    bp = length(L) == 2;

    if ~(lp || bp)
        error('The degree range is either one or two numbers')
    end

    maxL = max(L);

    % The spherical harmonic dimension
    if lp
        ldim = (L + 1) ^ 2;
    elseif bp
        ldim = (L(2) + 1) ^ 2 - L(1) ^ 2;
    end

end

function plotvweightmap(G, V, domain)
    [mesh, lon, lat] = eigwmesh(G, V, 2);
    mesh = mesh / max(mesh(:));
    figName = 'Eigenvalue-weighted map of the Slepian functions';

    if isa(domain, 'GeoDomain')
        lonlatd = domain.Lonlat('RotateBack', true);
    elseif ismatrix(domain)
        lonlatd = domain;
    end

    figure(999)
    set(gcf, 'Name', sprintf('%s (%s)', figName, upper(mfilename)), ...
        'NumberTitle', 'off')
    clf
    contourf(lon, lat, mesh, 0:0.1:1, 'LineStyle', 'none')
    hold on

    try
        [latd, lond] = flatearthpoly(lonlatd(:, 2), lonlatd(:, 1), 180);
        lonlatd = [lond, latd];
        plot(lonlatd(:, 1), lonlatd(:, 2), 'k')
    catch
    end

    hold off

    axis equal
    axis tight
    grid on

    formatlonticks
    formatlatticks
end
