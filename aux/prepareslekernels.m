%% PREPARESLEKERNELS - Prepares matrices for the fixed-ocean elastic SLE.
%
% Syntax
%   kernels = prepareslekernels(L, O)
%   kernels = prepareslekernels(__, Frame="CF", RotationFeedback=true)
%
% Input arguments
%   L - Maximum harmonic degree. Integer >= 1 (>= 2 with rotation).
%   O - Ocean localization matrix in kernel order.
%       Size: ((L+1)^2, (L+1)^2). Unit: dimensionless.
%   Frame (name-value) - Reference frame, CF or CM. Default frame: CF.
%   LoveNumSource (name-value) - Love-number source. Default source: ISSM.
%   RotationFeedback (name-value) - Whether to include rotation.
%       Default option: true.
%   KernelOrder (name-value) - Packing indices. Default order: kernelorder(L).
%
% Output arguments
%   kernels - Struct of ordinary numeric matrices, with no function handles.
%       O: ocean matrix; response: sparse G+R; KF: sparse forcing matrix.
%       oceanFunction: O(:,1); oceanArea: O(1,1); kernelOrder: packing indices.
%       waterDensity: conversion from metres to kg/m^2; config: model settings.
%
% Notes
%   S and F use kg/m^2, so G+R and KF are dimensionless. KF=P*(G+R)
%   -e00*e00.'/O(1,1), where P=I-e00*O(1,:)/O(1,1).
%   The system A=I-P*(G+R)*O is applied from these matrices without storing
%   another full dense matrix. No idempotence of the truncated O is assumed.
%
% See also
%   SOLVESLE, SOLVESLESYSTEM, PREPARELOWDEGREEMAP, LOVENUMBER
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function kernels = prepareslekernels(L, O, options)
    arguments (Input)
        L (1, 1) double {mustBeInteger, mustBePositive}
        O (:, :) double {mustBeReal}
        options.Frame (1, 1) string {mustBeMember(options.Frame, ["CF", "CM"])} = "CF"
        options.LoveNumSource (1, 1) string = "ISSM"
        options.RotationFeedback (1, 1) logical = true
        options.KernelOrder (:, 1) double = kernelorder(L)
    end
    M = (L + 1)^2;
    if ~isequal(size(O), [M, M]) || ~isfinite(O(1, 1)) || O(1, 1) <= 0
        error('ULMO:prepareslekernels:InvalidOceanKernel', ...
            'O must have size (L+1)^2 by (L+1)^2 and positive finite area.');
    end
    % Check in blocks to avoid a full logical temporary for a large kernel.
    for first = 1:256:M
        if any(~isfinite(O(:, first:min(first+255, M))), 'all')
            error('ULMO:prepareslekernels:InvalidOceanKernel', 'O must be finite.');
        end
    end
    if options.RotationFeedback && L < 2
        error('ULMO:prepareslekernels:RotationDegree', 'Rotation requires L >= 2.');
    end
    frame = char(options.Frame);
    loveNumSrc = char(options.LoveNumSource);
    doRotation = options.RotationFeedback;
    kernelOrder = options.KernelOrder;
    expectedOrder = kernelorder(L);
    if ~isequal(kernelOrder(:), expectedOrder(:))
        error('ULMO:prepareslekernels:InvalidOrder', 'Invalid packing indices.');
    end
    [~, degree] = addmon(L);
    degrees = [degree, degree];
    degrees = degrees(kernelOrder);
    llnGeoids = lovenumber(degrees, 'LLN geoid', frame, "Source", loveNumSrc);
    llnVlms = lovenumber(degrees, 'LLN VLM', frame, "Source", loveNumSrc);
    llnFactors = 1 + llnGeoids - llnVlms;
    tlnGeoids = lovenumber(degrees, 'TLN geoid', frame, "Source", loveNumSrc);
    tlnGeoid2 = lovenumber(2, 'TLN geoid', frame, "Source", loveNumSrc);
    tlnVlms = lovenumber(degrees, 'TLN VLM', frame, "Source", loveNumSrc);
    tlnFactors = 1 + tlnGeoids - tlnVlms;
    llnGeoid2 = lovenumber(2, 'LLN geoid', frame, "Source", loveNumSrc);

    % Many values from Adhikari et al. (2016)
    WATER_DENSITY = 1000; % to be consistent with other functions(mainly MASS2WEQ) and Adhikari et al. (2016)
    EARTH_DENSITY = 5517; % to be consistent with PLM2POT
    EARTH_RADIUS = 6371e3;
    EARTH_ANGULAR_VELOCITY = 7.2921e-5;
    EQUATORIAL_INERTIA = 8.0077e37;
    POLAR_INERTIA = 8.0345e37;
    SECULAR_LOVE_NUMBER = 0.942; % 'Secular (fluid) Love number '
    CHANDLER_WOBBLE_FREQUENCY = EARTH_ANGULAR_VELOCITY ...
        * (POLAR_INERTIA - EQUATORIAL_INERTIA) / EQUATORIAL_INERTIA ...
        * (1 - tlnGeoid2 / SECULAR_LOVE_NUMBER); % Should use the elastic Earth value
    GRAVITY = 9.81;

    % Construct the SLE kernel
    % Gravitational potential
    response = spdiags(3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1), 0, ...
        numel(degrees), numel(degrees));

    if doRotation
        cmnFac = (4 * pi / 3) * (1 + llnGeoid2) * EARTH_RADIUS ^ 6 * EARTH_ANGULAR_VELOCITY ^ 2 / GRAVITY;
        c00c20Fac = 4/3 / POLAR_INERTIA;
        c21s21Fac = EARTH_ANGULAR_VELOCITY / 5 / EQUATORIAL_INERTIA / CHANDLER_WOBBLE_FREQUENCY;

        % C00
        response(1, 1) = response(1, 1) - cmnFac * c00c20Fac;
        response(1, 5) = response(1, 5) + cmnFac * c00c20Fac / sqrt(5);
        % C20
        response(5, 1) = response(5, 1) + cmnFac * c00c20Fac / sqrt(5) * tlnFactors(5);
        response(5, 5) = response(5, 5) - cmnFac * c00c20Fac / 5 * tlnFactors(5);
        % C21 and S21
        response(6, 6) = response(6, 6) + cmnFac * c21s21Fac * tlnFactors(6);
        response(7, 7) = response(7, 7) + cmnFac * c21s21Fac * tlnFactors(7);
    end

    response = WATER_DENSITY * response;
    if any(~isfinite(nonzeros(response)))
        error('ULMO:prepareslekernels:InvalidResponse', 'Earth response must be finite.');
    end
    KF = response;
    KF(1, :) = KF(1, :) - sparse((O(1, :) * response) / O(1, 1));
    KF(1, 1) = KF(1, 1) - 1 / O(1, 1);
    config = rmfield(options, 'KernelOrder');
    config.L = L;
    config.Version = 1;
    kernels = struct('O', O, 'response', response, 'KF', KF, ...
        'oceanFunction', O(:, 1), 'oceanArea', O(1, 1), ...
        'kernelOrder', kernelOrder, 'waterDensity', WATER_DENSITY, 'config', config);
end
