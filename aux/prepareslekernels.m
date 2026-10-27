%% PREPARESLEKERNELS - Prepares matrices for the fixed-ocean elastic SLE.
%
% Syntax
%   [earthResponseKernel, forcingKernel, waterDensity] = prepareslekernels(L, oceanKernel)
%   [earthResponseKernel, forcingKernel] = prepareslekernels(__, ...
%       Frame="CF", RotationFeedback=true)
%
% Input arguments
%   L - Maximum harmonic degree. Integer >= 1 (>= 2 with rotation).
%   oceanKernel - Ocean localization matrix in kernel order.
%       Size: ((L+1)^2, (L+1)^2). Unit: dimensionless.
%   Frame (name-value) - Reference frame, CF or CM. Default frame: CF.
%   LoveNumSource (name-value) - Love-number source. Default source: ISSM.
%   RotationFeedback (name-value) - Whether to include rotation.
%       Default option: true.
%   KernelOrder (name-value) - Packing indices. Default order: kernelorder(L).
%
% Output arguments
%   earthResponseKernel - Sparse elastic Earth response G+R, including rotation.
%       Maps total surface mass density to relative sea-level load before the
%       ocean mass adjustment. Size: ((L+1)^2, (L+1)^2). Unit: dimensionless.
%   forcingKernel - Sparse mass-constrained forcing matrix KF.
%       Maps forcing loads to the SLE right-hand side.
%       Size: ((L+1)^2, (L+1)^2). Unit: dimensionless.
%   waterDensity - Density used to convert sea-level metres to surface load.
%       Size: scalar. Unit: kg/m^3.
%
% Notes
%   S and F use kg/m^2, so G+R and KF are dimensionless. KF=P*(G+R)
%   -e00*e00.'/oceanKernel(1,1), where P=I-e00*oceanKernel(1,:)/oceanKernel(1,1).
%   The system A=I-P*(G+R)*oceanKernel is applied from these matrices without storing
%   another full dense matrix. No idempotence of the truncated oceanKernel is assumed.
%
% See also
%   SOLVESLE, SOLVESLESYSTEM, PREPARELOWDEGREEMAP, LOVENUMBER
%
% Created by
%   2026/10/27, En-Chi Lee (williameclee@arizona.edu)

function [earthResponseKernel, forcingKernel, waterDensity] = prepareslekernels(L, oceanKernel, options)
    arguments (Input)
        L (1, 1) double {mustBeInteger, mustBePositive}
        oceanKernel (:, :) double {mustBeReal}
        options.Frame (1, 1) string {mustBeMember(options.Frame, ["CF", "CM"])} = "CF"
        options.LoveNumSource (1, 1) string = "ISSM"
        options.RotationFeedback (1, 1) logical = true
        options.KernelOrder (:, 1) double = kernelorder(L)
    end
    numCoeffs = (L + 1)^2;
    if ~isequal(size(oceanKernel), [numCoeffs, numCoeffs]) || ~isfinite(oceanKernel(1, 1)) || oceanKernel(1, 1) <= 0
        error('ULMO:prepareslekernels:InvalidOceanKernel', ...
            'oceanKernel must have size (L+1)^2 by (L+1)^2 and positive finite area.');
    end
    % Check in blocks to avoid a full logical temporary for a large kernel.
    for first = 1:256:numCoeffs
        if any(~isfinite(oceanKernel(:, first:min(first+255, numCoeffs))), 'all')
            error('ULMO:prepareslekernels:InvalidOceanKernel', 'oceanKernel must be finite.');
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
    waterDensity = 1000; % to be consistent with other functions(mainly MASS2WEQ) and Adhikari et al. (2016)
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
    earthResponseKernel = spdiags(3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1), 0, ...
        numel(degrees), numel(degrees));

    if doRotation
        cmnFac = (4 * pi / 3) * (1 + llnGeoid2) * EARTH_RADIUS ^ 6 * EARTH_ANGULAR_VELOCITY ^ 2 / GRAVITY;
        c00c20Fac = 4/3 / POLAR_INERTIA;
        c21s21Fac = EARTH_ANGULAR_VELOCITY / 5 / EQUATORIAL_INERTIA / CHANDLER_WOBBLE_FREQUENCY;

        % C00
        earthResponseKernel(1, 1) = earthResponseKernel(1, 1) - cmnFac * c00c20Fac;
        earthResponseKernel(1, 5) = earthResponseKernel(1, 5) + cmnFac * c00c20Fac / sqrt(5);
        % C20
        earthResponseKernel(5, 1) = earthResponseKernel(5, 1) + cmnFac * c00c20Fac / sqrt(5) * tlnFactors(5);
        earthResponseKernel(5, 5) = earthResponseKernel(5, 5) - cmnFac * c00c20Fac / 5 * tlnFactors(5);
        % C21 and S21
        earthResponseKernel(6, 6) = earthResponseKernel(6, 6) + cmnFac * c21s21Fac * tlnFactors(6);
        earthResponseKernel(7, 7) = earthResponseKernel(7, 7) + cmnFac * c21s21Fac * tlnFactors(7);
    end

    earthResponseKernel = waterDensity * earthResponseKernel;
    if any(~isfinite(nonzeros(earthResponseKernel)))
        error('ULMO:prepareslekernels:InvalidResponse', 'Earth response must be finite.');
    end
    forcingKernel = earthResponseKernel;
    forcingKernel(1, :) = forcingKernel(1, :) - sparse((oceanKernel(1, :) * earthResponseKernel) / oceanKernel(1, 1));
    forcingKernel(1, 1) = forcingKernel(1, 1) - 1 / oceanKernel(1, 1);
end
