%% KERNELCP
% Calculates the localisation matrix for some domain on the sphere.
% Not for polar patches! And not good for near-polar patches! (See
% GRUNBAUM)
% NOT WITHOUT MODIFICATIONS FOR REGIONS CONTAINING THE NORTH POLE OR THE
% SOUTH POLE! (For that, see GLMALPHA). Unit normalisation as in YLM.
% This is the parallel version of KERNELC.
%
% Syntax
%   kernelcp(demoId)
%       Runs a demonstration of the KERNELCP function.
%   [Kernel, XY, K1, K] = kernelcp(Lmax, dom, pars, ngl, rotb)
%       Calculates the localisation kernel for a given domain.
%
% Input arguments
%   Lmax - Maximum angular degree (bandwidth)
%   dom - The domain on the sphere where the kernel is to be calculated
%       - 'patch': spherical patch, with parameters in the argument PARS.
%           Note: better use GRUNBAUM/PLM2ROT in this case.
%       - 'sqpatch': square patch, with parameters in argument PARS.
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A cell array of the domain of interest (string) and the buffer
%           value (scalar): {domain, buf}.
%       - An N-by-2 matrix of the domain of interest in [lon, lat] format.
%       - A GeoDomain object.
%       The default value is 'patch'.
%   pars - Specification of the domain
%       - For patch, it should be [th0, ph0, thR], where...
%           - th0: Colatitude of the cap center, in radians.
%           - ph0: Longitude of the cap center, in radians.
%           - thR: Radius of the cap, in radians.
%       - For sqpatch, it should be [thN, thS, phW, phE]
%       - For a geographic domain specified as a string or cell, it is the
%           amount of upscaling to apply to the vertices of the domain.
%       - A string with the file name where the data will be saved.
%   ngl - The degree of the Gauss-Legendre integration, or 'alternative'
%           for an alternative calculation method (not parallel)
%       The default value is 200.
%   RotateBack - Whether to rotate the kernel back to its original
%       orientation
%       - 0: No rotation
%       - 1: For, e.g. 'antarctica', 'contshelves', if you rotated
%           coordinates to make the integration procedure work, this option
%           makes sure that the kernel matrix reflects this. If not, you
%           have to apply counterrotation after diagonalizing in
%           LOCALIZATION.
%       The default value is 0.
%
% Output arguments
%   Kernel - The localization kernel whose eigenfunctions we want,
%       indexed as: degree  [0  1  1  1  2  2  2  2  2]
%                   order   [0  0 -1  1  0 -1  1 -2  2]
%       The function LOCALIZATION later reindexes this in LMCOSI fashion.
%       Note: you can use ADDMOUT and ADDMON to modify, and see, e.g.
%       PLOTSLEP and KLMLMP2ROT for some implementations
%   XY - The outlines of the region into which you are localising
%   K1 - An intermediate result useful when rotb = 1, see KLMLMP2ROT
%   K - An verification result useful when rotb=1, see KLMLMP2ROT
%
% Examples
%   >>  L = 19;
%   >>  [Kernel, XY] = kernelcp(L, 'australia');
%   >> Kernel2 = kernelcp(L, 'australia', [], 'alternative');
%   and then DIFER, EIG, PLOTSLEP, etc, to evaluate the difference
%   >>  kernelcp('demo1') % For an illustration of the Antarctica matrix
%   >>  kernelcp('demo2') % For an illustration of the Antarctica functions
%   >>  kernelcp('demo3') % For a show of Australia
%   >>  kernelcp('demo4') % For a demonstration of the rotation of the
%       kernel
%   >>  kernelcp('demo5') % For a demonstration of kernelcp
%
% See also
%   LOCALIZATION, SDWREGIONS, GLMALPHA, DLMLMP, KERNELC2D,
%   PLOTSLEP, PLM2AVG, KERNELC, LEGENDREPRODINT, DLMLMP
%
% Last modified by
%   2025/10/16, williameclee@arizona.edu (@williameclee)
%   2024/08/15, williameclee@arizona.edu (@williameclee)
%   2023/11/11, fjsimons@alum.mit.edu (@fjsimons)
%   2017/05/26, plattner@alumni.ethz.ch (@AlainPlattner)
%   2016/09/23, charig@princeton.edu (@harig00)

function varargout = kernelcp_new(varargin)
    % Demos
    if ischar(varargin{1}) || isstring(varargin{1}) && ...
            contains(varargin{1}, 'demo')

        switch varargin{1}
            case 'demo1'
                kernelcp_demo1
            case 'demo2'
                kernelcp_demo2
            case 'demo3'
                kernelcp_demo3
            case 'demo4'
                kernelcp_demo4
            case 'demo5'
                kernelcp_demo5
            otherwise
                error(['Unknown demo ', demoname])
        end

        return

    end

    % Parse inputs
    [Lmax, domain, pars, ngl, rotb, forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin);
    K1 = nan; %#ok<NASGU>
    K = nan; %#ok<NASGU>

    tic;

    %% Initialisation
    % Generic path name that I like
    [outputPath1, outputPath2, outputFolder1, ~] = ...
        getoutputfile(Lmax, domain, pars, rotb);

    % Check if the kernel has already been calculated
    vars = {'Klmlmp', 'XY', 'K1', 'K'};

    if ~forceNew && ~(isstring(ngl) || ischar(ngl)) && ...
            exist(outputPath1, 'file') && ...
            all(ismember(vars, who('-file', outputPath1)))
        % Check the KERNELC directory
        data = load(outputPath1, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">localisation kernel</a>.\n', ...
                callchaintext(callChain), outputPath1, outputPath1);
        end

        varargout = {data.Klmlmp, data.XY, data.K1, data.K};
        return
    elseif ~forceNew && ~(isstring(ngl) || ischar(ngl)) && ...
            exist(outputPath2, 'file') && ...
            all(ismember(vars, who('-file', outputPath2)))
        % Check if you have a file in the old KERNELCP directory
        data = load(outputPath2, vars{:});

        if beQuiet <= 1
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">localisation kernel</a>.\nConsider moving your kernel files back to the KERNELC directory\n', ...
                callchaintext(callChain), outputPath1, outputPath1);
        end

        varargout = {data.Klmlmp, data.XY, data.K1, data.K};
        return
    end

    %% Computing the kernel matrix
    if beQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO><a href="matlab: open(''%s'')">%s</a>] Computing localisation kernel, %s\n', ...
            mfilename("fullpath"), mfilename, templine);
    end

    if strcmp(domain, 'patch')
        % For future reference
        th0 = pars(1);
        ph0 = pars(2);
        thR = pars(3);

        if th0 == 0
            disp('Really, should be putting in the GRUNBAUM call here')
            error('Not for polar caps! Use GRUNBAUM or SDWCAP instead')
            % BUT IN COMPARING, NOTE THAT THE SIGN MAY BE OFF
        end

        if thR > th0
            disp('Really, should be putting in the GRUNBAUM call here')
            error('Not for near-polar caps! Use GRUNBAUM, SDWCAP, then rotate')
        end

        % Convert all angles to degrees for CAPLOC only
        [lon, lat] = caploc(rad2deg([ph0 pi / 2 - th0]), rad2deg(thR), 100, 1);
        % Northern and Southern points, in radians
        thN = (th0 - thR);
        thS = (th0 + thR);
        XY = [lon, lat];
    elseif strcmp(domain, 'sqpatch')
        thN = pars(1);
        thS = pars(2);
        phW = pars(3);
        phE = pars(4);
        XY = rad2deg( ...
            [phW, pi / 2 - thN; phW, pi / 2 - thS; phE, pi / 2 - thS; ...
             phE, pi / 2 - thN; phW, pi / 2 - thN]);
    else

        if isnumeric(domain)
            % This case when coordinates in degrees are input as matrix
            XY = domain;

            if isstring(pars) || ischar(pars)
                % Use the input to define the file name that will be created
                outputPath1 = sprintf('%s/WREG-%s-%i.mat', outputFolder1, pars, Lmax);
            end

        else

            try
                XY = domain.Lonlat;
            catch
                XY = feval(domain.Domain, domain.Upscale, domain.Buffer);
            end

        end

        thN = deg2rad(90 - max(XY(:, 2)));
        thS = deg2rad(90 - min(XY(:, 2)));
    end

    % Introduce and dimensionalize variables and arrays
    [dems, ~, ~, ~, ~, mzo, bigm, bigl] = addmon(Lmax);
    dimK = (Lmax + 1) ^ 2;
    lenm = length(dems);
    Klmlmp = NaN(dimK, dimK);

    if strcmp(ngl, 'alternative')
        % Hold on and see if we can go more quickly here using GEOBOXCAP
        fax = 2 ^ 5;
        % Try oversampling the Nyquist degree by a certain factor
        degres = 180 / sqrt(Lmax * (Lmax + 1)) / fax;
        % Calculate the mask function
        [~, ~, r, lon, lat] = geoboxcap(Lmax, domain, [], degres);
        % Prepare the reindexing arrays
        % We're not fully using the recursion here, so there is wastage
        % Perform the masked spherical harmonic transform
        wbar = waitbar(0, sprintf('%s: Loop over all degrees and orders', upper(mfilename)));
        % With the recursions as they are, we are not yet taking full
        % advantage of this method. See Mark Wieczorek's Fortran code which
        % presumably works better for this case.
        for l = 0:Lmax
            % Remember the normalization conventions
            theYplus = ylm(l, 0:l, deg2rad((90 - lat)), deg2rad(lon)) ...
                * 2 * sqrt(pi);
            theYmins = ylm(l, -1:-1:-l, deg2rad((90 - lat)), deg2rad(lon)) ...
                * 2 * sqrt(pi);

            for m = 0:l
                waitbar((addmup(l) + m) / addmup(Lmax), wbar);
                % Return the expansion coefficients in "standard" real
                % harmonics order
                lmcosiplus = xyz2plm((-1) ^ m * theYplus(:, :, m + 1) .* r, Lmax);
                % Reindex the coefficients to "standard" localization kernel order
                % and put them into where the positive m sits
                posm = addmoff(l - 1) + 2 * m + 1;
                % Redundant check that we are at the right order and degree
                % difer(bigl(posm)-l)
                % difer(bigm(posm)-m)
                Klmlmp(posm, :) = lmcosiplus(2 * size(lmcosiplus, 1) + mzo)';

                if m > 0
                    % Return the expansion coefficients in "standard" real
                    % harmonics order
                    lmcosimins = xyz2plm((-1) ^ m * theYmins(:, :, m) .* r, Lmax);
                    % Also do the negatives which come right before the positives
                    Klmlmp(posm - 1, :) = lmcosimins(2 * size(lmcosimins, 1) + mzo)';
                end

            end

        end

        delete(wbar)

        % NOTE : THIS PIECE OF THE CODE IS REPEATED VERBATIM BELOW
        % By whichever way you obtained the kernel, now check if you might want
        % to rotate it back so its eigenvectors are "right", right away,
        % e.g. for Antarctica or ContShelves without needing to rotate as
        % part of LOCALIZATION
        if rotb
            disp(' ')
            disp('The input coordinates were rotated. Kernel will be unrotated,')
            disp('so its eigenfunctions will show up in the right place')
            disp(' ')
            % Get the rotation parameters for this particular region
            [XY, lonc, latc] = feval(domain, pars);

            if nargout < 4
                % Rotate the kernel, properly
                [Klmlmp, ~] = klmlmp2rot(Klmlmp, lonc, latc);
            else
                % Some extra verification in here
                [Klmlmp, ~, ~] = klmlmp2rot(Klmlmp, lonc, latc);
            end

            % else
            %     [lonc, latc, K1, K] = deal(0);
        end

        % Do not save this way of calculating the kernels
        outputPath1 = 'neveravailable';
    else
        % Regular Gauss-Legendre method here (not alternative)

        % Calculating different Gauss-Legendre points for all possible product
        % degrees is not a good idea since they get multiplied by more
        % functions of theta
        intv = cos([thS thN]);
        nGL = max(ngl, 2 * Lmax);
        % These are going to be the required colatitudes - forget XY
        [w, x, N] = gausslegendrecof(nGL, [], intv);
        save(fullfile(getenv('IFILES'), 'gausslegendrecof.mat'), 'w', 'x', 'N')
        % fprintf('%s calculated %i Gauss-Legendre points and weights\n', upper(mfilename), N)

        % First calculate the Legendre functions themselves
        % Note that dimK==sum(dubs)
        % dubs = repmat(2, lenm, 1);
        % dubs(mz) = 1;
        % comdi = [];
        % First, calculate all the Legendre functions themselves
        Xlm = NaN(length(x), lenm);

        % Calculate the Legendre polynomials
        ind = 0;

        for l = 0:Lmax
            Xlm(:, ind + 1:ind + l + 1) = ...
                (legendre(l, x(:)', 'sch') * sqrt(2 * l + 1))';
            ind = ind + l + 1;
        end

        if strcmp(domain, 'patch')
            % Get the parameters of the dom
            phint = dphpatch(acos(x), thR, th0, ph0);
        elseif strcmp(domain, 'sqpatch')
            % Always the same longitudinal integration interval
            phint = repmat([phW, phE], length(x), 1);
        else
            % Now we may have multiple pairs
            % Changed "dom" to "XY" here CTH
            phint = dphregion(acosd(x), [], XY);
            phint = deg2rad(phint);
        end

        % The number of elements that will be calculated is
        % nel = (dimK ^ 2 + dimK) / 2;

        parfor lm1dex = 1:dimK
            l1 = bigl(lm1dex);
            m1 = bigm(lm1dex);
            % Can only use the loop variable once per index.  So for Klmlmp,
            % also use index=lm1dex
            index = lm1dex;
            ondex = 0;
            I = NaN([length(x), dimK - index + 1]);
            % Instead of counting up andex and undex, write expressions for
            % them analytically for the row of Klmlmp we are calculating
            % countdown = [dimK:-1:1];
            % andex = 1 + sum(countdown(1:(index - 1)));
            % undex = sum(countdown(1:index));
            % We know bigo, andex, and undex, so just calculated exactly which
            % parts of XlmXlmp you need for this specific iteration
            smalll1 = abs(l1);
            smallm1 = abs(m1);
            pos1 = 1/2 * (smalll1) ^ 2 +1/2 * smalll1 + smallm1 + 1;
            XlmXlmp = 0;

            for lm2dex = lm1dex:dimK
                l2 = bigl(lm2dex);
                m2 = bigm(lm2dex);
                smalll2 = abs(l2);
                smallm2 = abs(m2);
                pos2 = 1/2 * (smalll2) ^ 2 +1/2 * smalll2 + smallm2 + 1;
                ondex = ondex + 1;
                XlmXlmp(1:length(x), ondex) = Xlm(:, pos1) .* Xlm(:, pos2);
                % Now evaluate the longitudinal integrals at the GL points
                if m1 > 0 && m2 > 0
                    I(:, ondex) = sinsin(acos(x), m1, m2, phint);
                elseif m1 <= 0 && m2 <= 0
                    I(:, ondex) = coscos(acos(x), m1, m2, phint);
                elseif m1 > 0 && m2 <= 0 % Got rid of redundant, pars below here
                    I(:, ondex) = sincos(acos(x), m1, m2, phint);
                elseif m1 <= 0 && m2 > 0
                    I(:, ondex) = sincos(acos(x), m2, m1, phint);
                end

            end

            % Do the calculation and set as a temp variable
            temprow = (w(:)' * (XlmXlmp .* I));
            % Pad the temp variable with the appropriate zeros out front
            temprow = [zeros(1, (index - 1)) temprow];
            % Now we can distribute over the kernel.  We need to do it this way
            % because if you slice Klmlmp with the loop variable (lm1dex) then
            % all other indicies need to be constant, or ':', or 'end.'
            Klmlmp(lm1dex, :) = temprow;
        end %parfor

        % Symmetrize the Kernel
        Klmlmp = Klmlmp + Klmlmp' - eye(size(Klmlmp)) .* Klmlmp;

        % Verify that the first value correctly gives the area of
        % the patch
        if strcmp(domain, 'patch')
            parea = 2 * pi * (1 - cos(thR));
            apo = abs(parea - Klmlmp(1)) / parea;
            % fprintf( ...
            %     'Area of the patch approximated to within %5.2f %s\n', ...
            %     apo * 100, '%')

            if apo * 100 > 1
                error('Something wrong with the area element: radians/degrees ?')
            end

        elseif strcmp(domain, 'sqpatch')
            parea = (cos(thN) - cos(thS)) * (phE - phW);
            apo = abs(parea - Klmlmp(1)) / parea;
            % fprintf( ...
            %     'Area of the patch approximated to within %5.2f %s\n', ...
            %     apo * 100, '%')

            if apo * 100 > 1
                error('Something wrong with the area element: radians/degrees ?')
            end

        else
            % fprintf('Area of the domain approximated as %8.3e\n', ...
            %     Klmlmp(1))
        end

        % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
        % did above here, taking the output of YLM and multiplying
        Klmlmp = Klmlmp / (4 * pi);

        % This then makes Klmlmp(1) the fractional area on the sphere

        % Save this now
        if isstring(pars) || ischar(pars)
            domain = pars;
        end

        % This is where the save statement used to be
    end

    if beQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

    % NOTE : THIS PIECE OF THE CODE IS REPEATED VERBATIM ABOVE
    % By whichever way you obtained the kernel, now check if you might want
    % to rotate it back so its eigenvectors are "right", right away,
    % e.g. for Antarctica or ContShelves without needing to rotate as
    % part of LOCALIZATION
    if rotb && ...
            (((ischar(domain) || isstring(domain)) && ismember(domain, {'antarctica', 'contshelves'})) || ...
            isa(domain, 'GeoDomain') && ismember(domain.Domain, {'antarctica', 'contshelves'}))
        % Get the rotation parameters for this particular region
        if ischar(domain) || isstring(domain)
            [XY, lonc, latc] = feval(domain, pars);
        elseif isa(domain, 'GeoDomain')
            [XY, lonc, latc] = feval(domain.Domain, "Upscale", domain.Upscale, "Buffer", domain.Buffer);
        else
            error(sprintf('ULMO:%s:CannotRotate', upper(mfilename)), ...
                'Domain type %s not recognised for rotation', upper(class(domain)))
        end

        if nargout < 4
            % Rotate the kernel, properly
            [Klmlmp, K1] = klmlmp2rot(Klmlmp, lonc, latc);
            K = 0;
        else
            % Some extra verification in here
            [Klmlmp, K1, K] = klmlmp2rot(Klmlmp, lonc, latc);
        end

    else
        [lonc, latc, K1, K] = deal(0);
    end

    %% Saving and returning variables
    varargout = {Klmlmp, XY, K1, K};

    if strcmp(outputPath1, 'neveravailable') || ~saveData
        return
    end

    save(outputPath1, 'Lmax', 'Klmlmp', 'domain', 'ngl', 'XY', ...
        'lonc', 'latc', 'K1', 'K', '-v7.3')

    if beQuiet <= 1
        fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">localisation kernel</a>.\n', ...
            callchaintext(callChain), outputPath1, outputPath1);
    end

end

%% Subfunctions
function varargout = parseinputs(Inputs)
    % Define fallback values
    dfOpt.Lmax = 12;
    dfOpt.domain = 'greenland';
    dfOpt.pars = [];
    dfOpt.ngl = 200;
    dfOpt.rotb = false;
    dfOpt.ForceNew = false;
    dfOpt.SaveData = true;
    dfOpt.BeQuiet = 0.5;

    ip = inputParser;
    addOptional(ip, 'Lmax', dfOpt.Lmax, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(ip, 'Domain', dfOpt.domain, ...
        @(x) ischar(x) || iscell(x) || isa(x, "GeoDomain") || ...
        isnumeric(x) || isempty(x));
    addOptional(ip, 'pars', dfOpt.pars, ...
        @(x) isnumeric(x) || ischar(x) || iscell(x) || isempty(x));
    addOptional(ip, 'ngl', dfOpt.ngl, ...
        @(x) isnumeric(x) || ischar(x) || isempty(x));
    addOptional(ip, 'rotb', dfOpt.rotb, ...
        @(x) isnumeric(x) || islogical(x) || isempty(x));
    addParameter(ip, 'ForceNew', dfOpt.ForceNew, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'SaveData', dfOpt.SaveData, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'BeQuiet', dfOpt.BeQuiet, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'CallChain', {}, @iscell);
    parse(ip, Inputs{:});

    Lmax = conddefval(ip.Results.Lmax, dfOpt.Lmax);
    domain = conddefval(ip.Results.Domain, dfOpt.domain);
    pars = conddefval(ip.Results.pars, dfOpt.pars);
    ngl = conddefval(ip.Results.ngl, dfOpt.ngl);
    rotb = logical(conddefval(ip.Results.rotb, dfOpt.rotb));
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
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

    varargout = ...
        {Lmax, domain, pars, ngl, rotb, forceNew, saveData, beQuiet, callChain};
end

function varargout = getoutputfile(Lmax, domain, pars, rotb)
    outputFolder1 = fullfile(getenv('IFILES'), 'KERNELC');
    outputFolder2 = fullfile(getenv('IFILES'), 'KERNELCP');

    if isstring(domain) || ischar(domain)

        switch domain
            case 'sqpatch' % If the domain is a square patch
                domainType = 'sqpatch';
            case 'patch' % If the domain is a spherical patch
                domainType = 'patch';
            otherwise
                error('Unknown domain type')
        end

    elseif iscell(domain) || isa(domain, "GeoDomain")
        domainType = 'geographic';
    else
        domainType = 'lonlat';
    end

    switch domainType
        case 'sqpatch'
            defval('pars', deh2rad([30, 90, 10, 90]));
            pars = round(rad2deg(pars));
            outputName = sprintf('%s-%i-%i-%i-%i-%i.mat', ...
                domain, Lmax, pars(1), pars(2), pars(3), pars(4));
        case 'patch'
            defval('pars', deg2rad([90, 75, 30]));
            pars = round(rad2deg(pars));
            outputName = sprintf('%s-%i-%i-%i-%i.mat', domain, Lmax, ...
                pars(1), pars(2), pars(3));
        case 'geographic'

            if ismember(domain.Domain, {'antarctica', 'contshelves'}) && rotb
                outputName = sprintf('WREG-%s-%i-1.mat', domain.Id, Lmax);
            else
                outputName = sprintf('WREG-%s-%i.mat', domain.Id, Lmax);
            end

        case 'lonlat'
            outputName = sprintf('%s-%i.mat', hash(domain, 'sha1'), Lmax);
    end

    outputPath1 = fullfile(outputFolder1, outputName);
    outputPath2 = fullfile(outputFolder2, outputName);

    varargout = ...
        {outputPath1, outputPath2, ...
         outputFolder1, outputFolder2};
end
