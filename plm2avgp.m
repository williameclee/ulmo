%% PLM2AVGP
% is a parallel version of PLM2AVG, for when you are running it many many times.
% Computes the integral and average value of a spherical harmonic
% function (lmcosi) within a specific area (dom).  This is done by
% creating an integration vector and multiplying by the unwrapped
% coefficients of the field (from lmcosi).
%
% Syntax
% 	[Int, A, miniK, XY] = PLM2AVGP(lmcosi, dom)
%	plm2avgp_new(demoCase)
%
% Inputs
% 	lmcosi - Standard-type real spherical harmonic expansion coefficients
% 	dom - A string name with an approved region such as 'africa', OR
% 	XY - its coordinates (such as supplied from 'africa' etc)
%	demoCase - 'demo1' or 'demo2' to run one of the two demo cases
% 		* demo1: Integrate the output from GEOBOXCAP over a region and
% 			see that you get something close to the region's area
%		* demo2: Compare this integration to summing the field which has
%			been expanded onto an equal area Fibonacci grid
%
% Outputs
% 	Int - Integral of the function within the region
% 	A - Average value of the function within the region
% 	miniK - Integration vector (also the first row/column from Klmlmp in KERNELC)
% 	XY - Coordinates of the region we integrated over
%
% See also
%   KERNELC, SPHAREA, FIBONACCI_GRID
%
% Last modified by
% 	2026/01/29, williameclee@arizona.edu
% 	2011/01/25, charig@princeton.edu

function varargout = plm2avgp(varargin)

    if nargin >= 1 && (ischar(varargin{1}) || isStringScalar(varargin{1})) && ...
            (strcmp(varargin{1}, 'demo1') || strcmp(varargin{1}, 'demo2'))
        demoCase = varargin{1};

        switch demoCase
            case 'demo1'
                varargout = plm2avgp_demo1();
            case 'demo2'
                plm2avgp_demo2();
                varargout = {};
        end

        return
    end

    [lmcosi, dom, beQuiet, callChain] = parseinputs(varargin);

    % Highest degree (bandwidth of the expansion)
    Lmax = lmcosi(end, 1);
    XY = dom.Lonlat;

    % Calculate northernmost and southernmost colatitude
    thN = 90 - max(XY(:, 2)); thN = thN * pi / 180;
    thS = 90 - min(XY(:, 2)); thS = thS * pi / 180;

    % Introduce and dimensionalize variables and arrays
    [dems, ~, mz, ~, ~, mzo, bigm, ~, ~, ~] = addmon(Lmax);
    dimK = (Lmax + 1) ^ 2;
    lenm = length(dems);

    % Calculate some Gauss-Legendre points
    intv = cos([thS thN]);
    % Degree of Gauss-Legendre integration
    ngl = 200;
    nGL = min(ngl, size(XY, 1) / 2);
    % These are going to be the required colatitudes - forget XY
    [w, x, N] = gausslegendrecof(nGL, [], intv);

    if beQuiet == 0
        fprintf('[ULMO>%s] Calculated %i Gauss-Legendre points and weights.\n', ...
            callchaintext([callChain, {mfilename}]), N);
    end

    % Some arrays needed in a bit
    % Note that dimK==sum(dubs)
    dubs = repmat(2, lenm, 1);
    dubs(mz) = 1;

    % Calculate all the Legendre functions themselves
    Xlm = nan(length(x), lenm);
    ind = 0;

    for l = 0:Lmax
        Xlm(:, ind + 1:ind + l + 1) = (legendre(l, x(:)', 'sch') * sqrt(2 * l + 1))';
        ind = ind + l + 1;
    end

    %disp('Legendre functions calculated')

    % No need for this to be longer, or to loop over anything
    % Normally it would help take out the first dimension of redundancy between
    % XlmXlmp and Klmlmp
    comdi = dubs';

    % In our ordering, the -1 precedes 1 and stands for the cosine term
    % This creates an indexing vector into Xlm that enlarges it from [0 01 012]
    % to [0 0-11 0-11-22]
    % Since we only have Xlm and not XlmXlm, we stop with coss.  In Kernelc,
    % this was again expanded into bigo to account for this redundancy in two
    % dimensions
    comdex = [1:lenm];
    coss = gamini(comdex, comdi);
    bigo = coss;

    % Get the longitudinal integration info for the domain
    defval('Nk', 10);
    % Now we may have multiple pairs
    phint = dphregion(acos(x) * 180 / pi, Nk, XY);
    phint = phint * pi / 180;

    % No need to initialize miniK since we can make it all at once
    % Calculate the longitudinal integrals for each l,m combination
    parfor lm1dex = 1:dimK
        m1 = bigm(lm1dex);

        if m1 <= 0
            I(:, lm1dex) = coscos(acos(x), m1, 0, phint);
        elseif m1 > 0
            I(:, lm1dex) = sincos(acos(x), m1, 0, phint);
        end

    end

    % COSCOS was reused here even though one of the m values is 0 because it is
    % useful for doing odd geographics in matrix form (i.e. if the continent
    % looks like a circle with hole in it)

    % Make miniK, really the first row/column of Klmlmp from KERNELC
    miniK = w(:)' * (Xlm(:, bigo) .* I);

    % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
    % did above here, taking the output of YLM and multiplying
    miniK = miniK / 4 / pi;

    % Compare with KERNELC
    defval('xver', 0)

    if beQuiet == 0
        % This then makes miniK(1) the fractional area on the sphere
        % Check by comparing to output from spharea, if you want
        % Note: SPHAREA defaults to 17 abcissas and weights, while this code uses 101.
        % So differences to be expected when continents are squiggly.
        % Check the first term which should equal the area on the unit sphere
        A1 = spharea(XY);
        A2 = areaint(XY(:, 2), XY(:, 1));
        fprintf('[ULMO>%s] Area check...  PLM2AVG A: %6.7f ; SPHAREA A: %6.7f ; AREAINT A: %6.7f', ...
            callchaintext([callChain, {mfilename}]), miniK(1), A1, A2)
    end

    % Now take the vector miniK and multiply with the vector for the function
    % you want (lmcosi) and you will get the integral of that field within the
    % region of interest
    thecofs = lmcosi(:, 3:4);
    theINT = miniK * thecofs(mzo);

    % To get the average value of the function in the region, divide by the
    % area, which is likely most accurate in miniK (i.e. more accurate than SPHAREA)
    A = theINT / miniK(1);

    varargout = {theINT, A, miniK, XY};

end

function varargout = parseinputs(Inputs)
    % Define fallback values
    dfOpt.lmcosi = [0, 0, 0, 0];
    dfOpt.domain = 'greenland';
    dfOpt.pars = [];
    dfOpt.BeQuiet = 0.5;

    ip = inputParser;
    addOptional(ip, 'lmcosi', dfOpt.lmcosi, ...
        @(x) isnumeric(x) && ismatrix(x) && size(x, 2) == 4);
    addOptional(ip, 'Domain', dfOpt.domain, ...
        @(x) ischar(x) || iscell(x) || isa(x, "GeoDomain") || ...
        isnumeric(x) || isempty(x));
    addOptional(ip, 'pars', dfOpt.pars, ...
        @(x) isnumeric(x) || ischar(x) || iscell(x) || isempty(x));
    addParameter(ip, 'BeQuiet', dfOpt.BeQuiet, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'CallChain', {}, @iscell);
    parse(ip, Inputs{:});

    lmcosi = conddefval(ip.Results.lmcosi, dfOpt.lmcosi);
    domain = conddefval(ip.Results.Domain, dfOpt.domain);
    pars = conddefval(ip.Results.pars, dfOpt.pars);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);
    callChain = [ip.Results.CallChain, {mfilename}];

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, "Upscale", pars);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2}, "Upscale", pars);
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{1}, "Upscale", pars, domain{2:end});
    end

    varargout = {lmcosi, domain, beQuiet, callChain};
end
