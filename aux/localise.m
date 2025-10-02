%% LOCALISE
% Localise spherical harmonic coefficients to a domain.
% This goal is very similar to that of Slepian functions', but this function don't truncate anything.
%
% Syntax
%   Plm = LOCALISE(Plm, domain, L)
%   Plm = LOCALISE(Plm, domain, L, "Inverse", true)
%   Plm = LOCALISE(Plm, "K", K)
%   [Plm, K] = LOCALISE(Plm, __)
%
% Input arguments
%   Plm - Spherical harmonic coefficients
%       The coefficients are in the form lmcosi, but the first two columns (degree and order) can be omitted.
%       Can be three-dimensional, where the third dimension is the data.
%   domain - Geographic domain
%       A geographic domain (GeoDomain object).
%   L - Bandwidth of the window
%       The default value is the degree of the data.
%   "Inverse" - Invert the kernel
%       The default value is false.
%   "K" - Kernel matrix
%       The default value is computed using KERNELCP_NEW.
%       For large L, loading the kernel matrix can be time-consuming.
%
% Output arguments
%   Plm - Localised spherical harmonic coefficients
%   K - Kernel matrix
%
% See also
%   KERNELCP_NEW
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/08/03, williameclee@arizona.edu (@williameclee)

function [plm, K] = localise(plm, varargin)
    ip = inputParser;
    addRequired(ip, 'Plm', ...
        @(x) isnumeric(x) && (size(x, 2) == 4 || size(x, 2) == 2));
    addOptional(ip, 'domain', @(x) isa(x, 'GeoDomain'));
    addOptional(ip, 'L', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'Inverse', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'K', [], @ismatrix);
    addParameter(ip, 'KernelOrder', [], @isvector);
    addParameter(ip, 'IsError', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', true, @(x) islogical(x) || isnumeric(x));
    parse(ip, plm, varargin{:});
    plm = ip.Results.Plm;
    domain = ip.Results.domain;
    L = ip.Results.L;
    isInverted = ip.Results.Inverse;
    K = ip.Results.K;
    j2 = ip.Results.KernelOrder;
    isError = ip.Results.IsError;
    beQuiet = ip.Results.BeQuiet;

    is3d = ndims(plm) == 3;

    if is3d
        nData = size(plm, 3);
    end

    includesDO = size(plm, 2) == 4;

    if includesDO

        if isempty(L)
            L = max(plm(:, 1));
        end

        if is3d
            plm = plm(:, 3:4, :);
        else
            plm = plm(:, 3:4);
        end

    end

    if isempty(j2)
        j2 = kernelorder(L);
    end

    nInput = size(plm, 1);
    nTarget = addmup(L);

    if nInput < nTarget

        if is3d
            plm(nTarget, 2, 1) = 0;
        else
            plm(nTarget, 2) = 0;
        end

    elseif nInput > nTarget

        if is3d
            plm = plm(1:nTarget, :, :);
        else
            plm = plm(1:nTarget, :);
        end

    end

    if is3d
        plm = reshape(plm, [nInput * size(plm, 2), size(plm, 3)]);
        plm = plm(j2, :);
    else
        plm = plm(j2);
    end

    if isempty(K)
        K = kernelcp_new(L, domain, "BeQuiet", beQuiet);
    end

    if isInverted
        K = eye(size(K)) - K;
    end

    if ~isError
        ofun = K * plm;
    else
        ofun = sqrt((K .^ 2) * (plm .^ 2));
    end

    if is3d
        plm = zeros([nTarget * 2, nData]);
        plm(j2, :) = ofun;
        plm = reshape(plm, [nTarget, 2, nData]);
    else
        plm = zeros([nTarget, 2]);
        plm(j2) = ofun;
    end

    if includesDO

        [order, degree] = addmon(L);

        if is3d
            plm = cat(2, repmat(degree, [1, 1, nData]), repmat(order, [1, 1, nData]), plm);
        else
            plm = [degree, order, plm];
        end

    end

end
