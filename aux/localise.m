%% LOCALISE - Localise spherical harmonic coefficients to a domain.
% This goal is very similar to that of Slepian functions', but this 
% function don't truncate anything.
%
% Syntax
%   Plm = LOCALISE(Plm, domain, L)
%   Plm = LOCALISE(Plm, domain, L, "Inverse", true)
%   Plm = LOCALISE(Plm, "K", K)
%   coeffs = LOCALISE(Plm, "L", L, "K", K, "OutputIndices", indices)
%   [Plm, K] = LOCALISE(Plm, __)
%
% Input arguments
%   Plm - Spherical harmonic coefficients
%       The coefficients are in the form lmcosi, but the first two columns 
%       (degree and order) can be omitted.
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
%   OutputIndices (name-value) - Selected rows in the kernel's harmonic ordering.
%       Positive integer vector, with entries no greater than (L+1)^2.
%       Preserves the requested order and repeated indices. Nonempty indices
%       compute only the selected coefficients and return a matrix without
%       degree/order columns. Empty indices retain the full lmcosi output.
%       Default selection: [] (all coefficients in the original format).
%   IsError (name-value) - Whether to propagate input standard deviations.
%       Uses squared kernel entries and assumes uncorrelated input errors.
%       Default option: false.
%
% Output arguments
%   Plm - Localised spherical harmonic coefficients
%       With OutputIndices, size is (number of indices, number of epochs),
%       including a column vector for one epoch, regardless of input metadata.
%   K - Kernel matrix
%       Remains the full kernel even when OutputIndices is supplied.
%
% See also
%   KERNELCP_NEW
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified
%   2026/10/26, OpenAI Codex

function [plm, K] = localise(plm, varargin)
    ip = inputParser;
    addRequired(ip, 'Plm', ...
        @(x) isnumeric(x) && (size(x, 2) == 4 || size(x, 2) == 2));
    addOptional(ip, 'domain', @(x) isa(x, 'GeoDomain'));
    addOptional(ip, 'L', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'Inverse', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'K', [], @ismatrix);
    addParameter(ip, 'KernelOrder', [], @isvector);
    addParameter(ip, 'OutputIndices', [], @isnumeric);
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
    outputIndices = ip.Results.OutputIndices;
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

    if ~isempty(outputIndices)
        assert(isreal(outputIndices) && isvector(outputIndices) ...
            && all(isfinite(outputIndices)) ...
            && all(outputIndices >= 1 & outputIndices == fix(outputIndices) ...
            & outputIndices <= (L + 1)^2), ...
            'ULMO:localise:InvalidOutputIndices', ...
            'OutputIndices must be a vector of integer indices from 1 to (L+1)^2.');
    end

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
        plm = reshape(plm, [], size(plm, 3));
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

    if isempty(outputIndices)
        selectedKernel = K;
    else
        selectedKernel = K(outputIndices, :);
    end

    if ~isError
        ofun = selectedKernel * plm;
    else
        ofun = sqrt((selectedKernel .^ 2) * (plm .^ 2));
    end

    if ~isempty(outputIndices)
        plm = ofun;
        return
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
