%% LOWDEGREECOEFFICIENTS - Maps jointly estimated coefficients to harmonic indices.
%
% Syntax
%   [locs, kernelIds, degrees] = lowdegreecoefficients()
%   [locs, kernelIds, degrees] = lowdegreecoefficients(c20, c30, c21s21)
%
% Input arguments
%   c20 (optional) - Include C20.
%       Scalar logical.
%       Default option: false.
%   c30 (optional) - Include C30.
%       Scalar logical.
%       Default option: false.
%   c21s21 (optional) - Include C21 and S21.
%       Scalar logical.
%       Default option: false.
%
% Output arguments
%   locs - Harmonic location of each coefficient.
%       Each row contains an lmcosi row and cosine/sine selector (1/2).
%       Add 2 to the selector to obtain the column in an lmcosi array.
%       K = 3 + c20 + c30 + 2*c21s21.
%       Size: (K, 2).
%   kernelIds - Indices in the real harmonic kernel ordering.
%       Size: (1, K).
%   degrees - Harmonic degrees for unit conversion.
%       Size: (1, K).
%
% Notes
%   Output order: C10, C11, S11, optional C20, optional C30, then C21/S21.
%   Degree-1 coefficients are always included.
%
% See also
%   SOLVEDEGREE1, KERNELORDER, CONVERTGRAVITY
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
