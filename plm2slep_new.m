%% PLM2SLEP
% Finds the spherical harmonic expansion coefficients into a SINGLE-CAP,
% potentially rotated, Slepian basis of a function whose real spherical
% harmonic expansion coefficients are given.
%
% Syntax
%   PLM2SLEP(demoId)
%       Runs a demo with the specified name.
%   [falpha, V, N] = PLM2SLEP(Plm, r, L, phi, theta, omega)
%       Finds the expansion coefficients of the function into the Slepian
%       basis of a polar cap.
%   [falpha, V, N] = PLM2SLEP(Plm, domain, L)
%       Finds the expansion coefficients of the function into the Slepian
%       basis of a geographic domain.
%   [falpha, V, N] = PLM2SLEP(__, nosort, truncation)
%       Finds the expansion with the specified sorting and the number of
%       eigenfunctions to use.
%   [__, MTAP, G] = PLM2SLEP(__)
%
% Input arguments
%   Plm - Standard-type real spherical harmonic expansion coefficients
%   r - Radius of the concentration region in degrees
%       Alternatively, phi, theta and omega can also be specified following
%       r.
%       The default value is 30 degrees.
%       Format: scalar or 1 x 4 vector.
%   domain - Geographic domain or a latitude-longitude pair
%       - A geographic domain (GeoDomain object).
%       - A string of the domain name.
%       - A cell array of the form {'domain name', buf}.
%       - A N-by-2 matrix of longitude-latitude vertices.
%   L - Bandwidth of the window
%       The default value is 18.
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the centre of the tapers in degrees
%       The default values are 0.
%   NoSort - Logical flag to sort the eigenvalues and eigenfunctions
%       - false: Will sort the output according to the global eigenvalue.
%       - true: Will not sort thus the "block" sorting is preserved,
%       The default value is false (i.e. sort).
%   Truncation - Number of largest eigenfunctions in which to expand
%       The default value is all the (L+1)^2 eigenfunctions (no truncation).
%
% Output arguments
%   falpha - Expansion coefficients of the function into the Slepian
%       basis
%   V - Eigenvalues of the Slepian functions in question
%   N - Shannon number
%   MTAP - Orders of the Slepian functions in question, if preserved
%   G - Matrix with the spherical harmonic expansion
%            coefficients of the Slepian functions used
%
% Examples
%   Check that single-order functions correctly transform back
%   >>  plm2slep('demo1')
%   Slight variation on the above with multiple same-order ones
%   >>  plm2slep('demo2')
%   Moon without the South Pole Aitken Basin
%   >>  plm2slep('demo3')
%
% See also
%   PTOSLEP, GLMALPHA, GLMALPHAPTO, SLEP2PLM
%
% Last modified by
%   2025/05/23, williameclee@arizona.edu (@williameclee)
%   2023/09/26, fjsimons@alum.mit.edu (@fjsimons)
%   2013/04/24, charig@princeton.edu (@harig00)

function varargout = plm2slep_new(varargin)
    %% Initialisation
    % Demos
    if ischar(varargin{1}) || isstring(varargin{1})
        demoId = varargin{1};

        switch demoId
            case 'demo1'
                plm2slep_demo1
            case 'demo2'
                plm2slep_demo2
            case 'demo3'
                plm2slep_demo3
            otherwise
                error('Unknown demo ''%s''', upper(demoId))
        end

        return
    end

    % Parse inputs
    [lmcosi, domain, L, phi, theta, omega, nosort, J, beQuiet, GVN, fmt, isError, callChain] = ...
        parseinputs(varargin{:});

    maxL = max(L);
    % The spherical harmonic dimension
    if isscalar(L)
        ldim = (L + 1) ^ 2;
    else
        ldim = (L(2) + 1) ^ 2 - L(1) ^ 2;
    end

    J = conddefval(J, ldim);

    %% Computing the projection
    % If it is the standard North-Polar cap or a geographic region, it's easy
    if ~isempty(GVN)
        G = GVN{1};
        V = GVN{2};
        N = GVN{3};
        MTAP = nan;
    elseif phi == 0 && theta == 0 && omega == 0
        % Get the Slepian basis; definitely not block-sorted as for the rotated
        % versions this will make no sense at all anymore
        % Glmalpha can handle a string, cell, or coordinates as TH, so this is ok
        [G, V, ~, ~, N, ~, MTAP, ~] = ...
            glmalpha_new(domain, L, "J", J, "BeQuiet", beQuiet, "CallChain", callChain);

    else
        % Need to get a complete GLMALPHA but for the rotated basis
        % Definitely, "single-order" has lost its meaning here, but the MTAP
        % will still identify what the order of the unrotated original was
        [G, V, ~, ~, N, ~, MTAP, ~] = ...
            glmalphapto(domain, L, phi, theta, omega);
    end

    if ~nosort
        % Sort by decreasing eigenvalue
        [V, vi] = sort(V, 'descend');
        G = G(:, vi);

        if ~isnan(MTAP)
            MTAP = MTAP(vi);
        end

        % If you don't do this, the eigenfunctions are ordered in the way
        % that they correspond to single-orders back when, unrotated, they
        % belonged to a polar cap, and the eigenvalues are sorted within
        % these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.
    end

    % Get the mapping from LMCOSI into not-block-sorted GLMALPHA
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ronm] = addmon(maxL);

    % Make sure that the requested L acts as truncation on lmcosi
    % or if we don't have enough, pad with zeros
    if ismatrix(lmcosi)

        if size(lmcosi, 2) == 4
            lmcosi = lmcosi(:, 3:4);
        end

        if size(lmcosi, 1) < addmup(maxL)
            % The lm part is not used anyway
            lmcosi(addmup(maxL), 1, :) = 0;
        elseif size(lmcosi, 1) > addmup(maxL)
            lmcosi = lmcosi(1:addmup(maxL), :);
        end

        % Perform the expansion of the signal into the Slepian basis
        falpha = G' * lmcosi(ronm(1:(maxL + 1) ^ 2));

    elseif ndims(lmcosi) == 3

        if strcmp(fmt, 'timefirst')
            lmcosi = permute(lmcosi, [2, 3, 1]);
        end

        if size(lmcosi, 2) == 4
            lmcosi = lmcosi(:, 3:4, :);
        end

        if size(lmcosi, 1) < addmup(maxL)
            % The lm part is not used anyway
            lmcosi(addmup(maxL), 1, :) = 0;
        else
            lmcosi = lmcosi(1:addmup(maxL), :, :);
        end

        lmcosi = reshape(lmcosi, ...
            [size(lmcosi, 1) * size(lmcosi, 2), size(lmcosi, 3)]);

        % Perform the expansion of the signal into the Slepian basis
        if ~isError
            falpha = G' * lmcosi(ronm(1:(maxL + 1) ^ 2), :);
        else
            falpha = sqrt((G .^ 2)' * lmcosi(ronm(1:(maxL + 1) ^ 2), :) .^ 2);
        end

        if strcmp(fmt, 'timefirst')
            falpha = permute(falpha, [2, 1]);
        end

    end

    %% Returning requested outputs
    varargout = {falpha, V, N, MTAP, G};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    dfOpt.Domain = 30;
    dfOpt.L = 18;
    dfOpt.Phi = 0;
    dfOpt.Theta = 0;
    dfOpt.Omega = 0;
    dfOpt.NoSort = false;
    dfOpt.Truncation = [];
    dfOpt.Upscale = 0;
    dfOpt.Format = 'timefirst';
    dfOpt.IsError = false;
    dfOpt.BeQuiet = 0.5;

    ip = inputParser;
    addRequired(ip, 'lmcosi', ...
        @(x) isnumeric(x) || ischar(x));
    addOptional(ip, 'Domain', dfOpt.Domain, ...
        @(x) ischar(x) || iscell(x) || isa(x, "GeoDomain") || ...
        isnumeric(x) || isempty(x));
    addOptional(ip, 'L', dfOpt.L, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'phi', dfOpt.Phi, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'theta', dfOpt.Theta, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'omega', dfOpt.Omega, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'NoSort', dfOpt.NoSort, ...
        @(x) isnumeric(x) || islogical(x) || isempty(x));
    addOptional(ip, 'J', dfOpt.Truncation, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'Upscale', dfOpt.Upscale, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'MoreDomainSpecs', {}, @iscell);
    addParameter(ip, 'GVN', {}, @(x) iscell(x) && length(x) == 3);
    addParameter(ip, 'Format', dfOpt.Format, ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'IsError', dfOpt.IsError, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', dfOpt.BeQuiet, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, "CallChain", {}, @iscell);
    parse(ip, varargin{:});

    lmcosi = ip.Results.lmcosi;
    domain = conddefval(ip.Results.Domain, dfOpt.Domain);
    L = conddefval(ip.Results.L, dfOpt.L);
    phi = conddefval(ip.Results.phi, dfOpt.Phi);
    theta = conddefval(ip.Results.theta, dfOpt.Theta);
    omega = conddefval(ip.Results.omega, dfOpt.Omega);
    nosort = conddefval(ip.Results.NoSort, dfOpt.NoSort);
    J = conddefval(ip.Results.J, dfOpt.Truncation);
    upscale = conddefval(ip.Results.Upscale, dfOpt.Upscale);
    moreRegionSpecs = ip.Results.MoreDomainSpecs;
    beQuiet = double(conddefval(ip.Results.BeQuiet, dfOpt.BeQuiet)) * 2;
    GVN = ip.Results.GVN;
    fmt = ip.Results.Format;
    isError = logical(ip.Results.IsError);
    callChain = [ip.Results.CallChain, {mfilename}];

    if isnumeric(domain) && isvector(domain)

        if length(domain) ~= 4
            error('If domain is a vector, it must have 4 elements: [R, phi, theta, omega]')
        end

        phi = domain(2);
        theta = domain(3);
        omega = domain(4);
        domain = domain(1);
    elseif ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, "Upscale", upscale, moreRegionSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = GeoDomain(domain{1}, ...
            "Buffer", domain{2}, "Upscale", upscale, moreRegionSpecs{:});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{1}, "Upscale", upscale, domain{2:end});
    end

    varargout = {lmcosi, domain, L, phi, theta, omega, nosort, J, beQuiet, GVN, fmt, isError, callChain};
end
