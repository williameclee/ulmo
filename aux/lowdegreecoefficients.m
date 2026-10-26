%% LOWDEGREECOEFFICIENTS - Finds locations and degrees of jointly estimated harmonics.
%
% Syntax
%   [locs, kernelIds, degrees] = lowdegreecoefficients(c20, c30, c21s21)
%
% Output arguments
%   locs - Rows in lmcosi and cosine/sine selectors (1/2).
%   kernelIds - Indices in the real spherical harmonic kernel ordering.
%   degrees - Harmonic degree used for unit conversion.
%
% Notes
%   Order: C10, C11, S11, optional C20, optional C30, optional C21, S21.
%
% Created by
%   2026/10/26, En-Chi Lee (williameclee@gmail.com)

function [locs, kernelIds, degrees] = lowdegreecoefficients(c20, c30, c21s21)

    arguments (Input)
        c20 (1, 1) logical = false
        c30 (1, 1) logical = false
        c21s21 (1, 1) logical = false
    end

    locs = [2, 1; 3, 1; 3, 2];
    kernelIds = [2, 3, 4];
    degrees = [1, 1, 1];

    if c20
        locs(end + 1, :) = [4, 1];
        kernelIds(end + 1) = 5;
        degrees(end + 1) = 2;
    end

    if c30
        locs(end + 1, :) = [7, 1];
        kernelIds(end + 1) = 10;
        degrees(end + 1) = 3;
    end

    if c21s21
        locs = [locs; 5, 1; 5, 2];
        kernelIds = [kernelIds, 6, 7];
        degrees = [degrees, 2, 2];
    end

end
