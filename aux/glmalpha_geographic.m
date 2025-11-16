%% GLMALPHA_GEOGRAPHIC
% Finds G, V, and N for the given grographical domain.
% This is an auxiliary function separated from GLMALPHA.
%
% Notes
%   So far we only have wanted to remove small portions of the kernel.
%   Therefore at the moment here, we make the whole thing, and the apply
%   the bp after the fact. However, in the future if you want a bp at large
%   L, we should modify kernelcp to only make a partial kernel.
%
% Input arguments
%   L - The maximum degree and order of the spherical harmonic expansion
%   domain - The domain of the spherical harmonic expansion
%       - A GeoDomain object.
%       - A lonlat matrix.
%   upscale - How many times to upscale the data
%       This corresponds to the argument SORD in GLMALPHA.
%   anti - Whether to use the anti-domain
%   rotb - Whether to rotate back the domain
%       This is only available for an GeoDomain object, where the domain
%       function supports the argument 'rotated'.
%   ldim - The dimension of the kernel
%   bp - is a band-pass L
%
% See also
%   GLMALPHA, KERNELC, KERNELCP, LOCALIZATION
%
% Last modified by
% 	2025/11/03, williameclee@arizona.edu

function [G, V, N] = glmalpha_geographic( ...
        L, domain, upscale, isAntiDomain, rotateBack, ldim, isBp, ...
        EL, EM, xver, beQuiet, forceNew, mesg, callChain)
    %% Computing the kernel
    % See if we can run this calculation in parallel
    canRunParallel = license('test', 'distrib_computing_toolbox') && ...
        (verLessThan('matlab', '8.2') && ...
        parpool('processes').NumWorkers > 0 || ...
        ~verLessThan('matlab', '8.2'));
    % Actually computing the kernel
    if canRunParallel
        Kernel = kernelcp_new(L, domain, upscale, "ForceNew", forceNew, ...
            "BeQuiet", beQuiet, "CallChain", callChain);
    else
        Kernel = kernelc(L, domain, upscale);
    end

    if isAntiDomain
        % Get the complimentary region
        Kernel = eye(size(Kernel)) - Kernel;
    end

    if isBp
        % Remove the beginning section of the kernel
        rem = isBp * L(1) ^ 2;
        Kernel(:, 1:rem) = [];
        Kernel(1:rem, :) = [];
    end

    %% Calculating the eigenfunctions/values
    if ~isBp % low pass
        [G, V] = eig(Kernel);
    else % band pass
        [Gbp, V] = eig(Kernel);
        G(L(1) ^ 2 + 1:end, :) = Gbp;
    end

    % Sort the eigenfunctions/eigenvalues
    [V, isrt] = sort(sum(real(V), 1));
    V = fliplr(V);
    V = V(:);
    G = G(:, fliplr(isrt));

    [~, ~, ~, ~, ~, ~, ems, els, R1, R2] = addmon(L);
    % This indexes the orders of G back as 0 -101 -2-1012 etc
    G = G(R1, :);
    % Check indexing
    difer(els(R1) - EL, [], [], mesg)
    difer(ems(R1) - EM, [], [], mesg)

    % Calculate Shannon number and compare this with the theory
    N = sum(V);

    % Is the Shannon number right? Need the area of the region
    if ~isBp % low pass
        difer(ldim * Kernel(1) - N, [], [], mesg)
    else % band pass
        difer(ldim * spharea(domain) - N, [], [], mesg)
    end

    % Check if the expansion of a basis function is indeed either 1 or 0
    if ~isBp && xver == 1
        disp('Excessive verification')
        % Is the area right? Don't be too demanding
        difer(Kernel(1) - abs(isAntiDomain - spharea(domain)), 4, [], mesg)
        % This is a bit double up... but it's only for excessive verification
        [V1, C, ~, ~, ~, ~, GG] = localization(L, domain, upscale);
        difer(V(:) - V1(:), [], [], mesg)
        % A test by expansion and orthogonality
        for index = 1:length(C)
            salpha = G' * C{index}(R2);
            % Only one of these functions should get "hit"
            difer(sum(abs(salpha) > 1e-8) - 1, [], [], mesg)
        end

        % Yet another way, see LOCALIZATION
        [~, ~, ~, ~, ~, ~, ~, ~, R1, ~] = addmon(L);
        % It's only a matter of indexing and ordering
        difer(GG(R1, :) - G)
    end

    %% Rotation checking
    % This is only available for an GeoDomain object and not a lonlat
    % matrix.
    if isnumeric(domain)

        if rotateBack
            warning('Rotation is not supported for lonlat matrices')
        end

        return
    end

    try
        rotateBack = feval(domain.Domain, 'rotated');
    catch
        rotateBack = false;
    end

    % Now do the rotation
    if isscalar(rotateBack) && islogical(rotateBack) && rotateBack
        % Get the rotation parameters to rotate G. Note, the region
        % rotation angles that we return from the functions (lonc, latc)
        % are the same regardless of if we did a buffer, as they pertain
        % to the original region
        try
            [~, lonc, latc] = feval(domain.Domain);
            G = rotateGp(G, lonc, latc);
        catch
        end

    end

end
