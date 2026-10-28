%% PREPARESLEKERNELS - Prepares matrices for the fixed-ocean elastic SLE.
%
% Syntax
%   [earthRespKer, forcKer, waterDen] = prepareslekernels(L, oceanKernel)
%   [earthRespKer, forcKer] = prepareslekernels(__, Frame="CF", RotationFeedback=true)
%
% Input arguments
%   L - Maximum harmonic degree.
%       Integer >= 1 (>= 2 with rotation).
%   oceanKernel - Ocean localization matrix in kernel order.
%       Size: ((L+1)^2, (L+1)^2).
%       Unit: dimensionless.
%   Frame (name-value) - Reference frame, CF or CM.
%       Default frame: CF.
%   LoveNumSource (name-value) - Love-number source.
%       Default source: ISSM.
%   RotationFeedback (name-value) - Whether to include rotation.
%       Default option: true.
%   KernelOrder (name-value) - Packing indices.
%       Default order: kernelorder(L).
%
% Output arguments
%   earthRespKer - Sparse elastic Earth response G+R, including rotation.
%       Maps total surface mass density to RSL load before the ocean mass 
%       adjustment.
%       Size: ((L+1)^2, (L+1)^2).
%       Unit: dimensionless.
%   forcKer - Sparse mass-constrained forcing matrix KF.
%       Maps forcing loads to the SLE right-hand side.
%       Size: ((L+1)^2, (L+1)^2).
%       Unit: dimensionless.
%   waterDen - Density used to convert sea-level metres to surface load.
%       Size: scalar. Unit: kg/m^3.
%
% Notes
%   The kernels defining the mass-density SLE are
%       S = KS * S + KF * F,
%       KS = (I - e00 * e00^T * O / O00) * (G + R) * O,
%       KF = (I - e00 * e00^T * O / O00) * (G + R) - e00 * e00^T / O00.
%   where:
%     - O: Ocean localisation kernel (ocnKer)
%     - G+R: Earth's response kernel to load (earthRespKer)
%     - KS: RSL load's response kernel to itself
%     - KF: RSL load's response kernel to forcing (forcKer)
%     - O00: ocean area fraction
%     - e00: unit vector with its first entry equal to 1
%
% See also
%   SOLVESLE, SOLVESLESYSTEM, PREPARELOWDEGREEMAP, LOVENUMBER
%
% Created by
%   2026/10/28, En-Chi Lee (williameclee@arizona.edu)

function [earthRespKer, forcKer, waterDen] = prepareslekernels(L, ocnKer, options)

    arguments (Input)
        L (1, 1) double {mustBeInteger, mustBePositive}
        ocnKer (:, :) double {mustBeReal}
        options.Frame {mustBeTextScalar, mustBeMember(options.Frame, ["CF", "CM"])} = "CF"
        options.LoveNumSource {mustBeTextScalar} = "ISSM"
        options.RotationFeedback (1, 1) logical = true
        options.KernelOrder (:, 1) double = kernelorder(L)
    end

    arguments (Output)
        earthRespKer (:, :) double {mustBeReal}
        forcKer (:, :) double {mustBeReal}
        waterDen (1, 1) double {mustBeReal}
    end

    %% Initialisation
    nCoeffs = (L + 1) ^ 2;

    if ~isequal(size(ocnKer), [nCoeffs, nCoeffs]) || ~isfinite(ocnKer(1, 1)) || ocnKer(1, 1) <= 0
        error('ULMO:prepareslekernels:InvalidOceanKernel', ...
        'oceanKernel must have size (L+1)^2 by (L+1)^2 and positive finite area.');
    end

    % Check in blocks to avoid a full logical temporary for a large kernel.
    blockSize = 256;

    for first = 1:blockSize:nCoeffs

        if any(~isfinite(ocnKer(:, first:min(first + 255, nCoeffs))), 'all')
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

    %% Data preparation
    % Load Love numbers
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

    % Physical constants (many from Adhikari et al., 2016)
    waterDen = 1000;
    % to be consistent with other functions (mainly MASS2WEQ) and Adhikari et al. (2016)
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

    %% Kernal construction
    % Gravitational response
    earthRespKer = spdiags(3 / EARTH_DENSITY * llnFactors ./ (2 * degrees + 1), 0, ...
        numel(degrees), numel(degrees));

    % Rotational response
    if doRotation
        cmnFac = (4 * pi / 3) * (1 + llnGeoid2) * EARTH_RADIUS ^ 6 * EARTH_ANGULAR_VELOCITY ^ 2 / GRAVITY;
        c00c20Fac = 4/3 / POLAR_INERTIA;
        c21s21Fac = EARTH_ANGULAR_VELOCITY / 5 / EQUATORIAL_INERTIA / CHANDLER_WOBBLE_FREQUENCY;

        % C00
        earthRespKer(1, 1) = earthRespKer(1, 1) - cmnFac * c00c20Fac;
        earthRespKer(1, 5) = earthRespKer(1, 5) + cmnFac * c00c20Fac / sqrt(5);
        % C20
        earthRespKer(5, 1) = earthRespKer(5, 1) + cmnFac * c00c20Fac / sqrt(5) * tlnFactors(5);
        earthRespKer(5, 5) = earthRespKer(5, 5) - cmnFac * c00c20Fac / 5 * tlnFactors(5);
        % C21/S21
        earthRespKer(6, 6) = earthRespKer(6, 6) + cmnFac * c21s21Fac * tlnFactors(6);
        earthRespKer(7, 7) = earthRespKer(7, 7) + cmnFac * c21s21Fac * tlnFactors(7);
    end

    earthRespKer = earthRespKer * waterDen; % RSL -> load

    if any(~isfinite(nonzeros(earthRespKer)))
        error('ULMO:prepareslekernels:InvalidResponse', 'Earth response must be finite.');
    end

    % Forcing kernel
    forcKer = earthRespKer;
    forcKer(1, :) = forcKer(1, :) - sparse((ocnKer(1, :) * earthRespKer) / ocnKer(1, 1));
    forcKer(1, 1) = forcKer(1, 1) - 1 / ocnKer(1, 1);
end
