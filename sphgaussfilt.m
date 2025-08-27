%% SPHGAUSSFILT
% Applys a Gaussian filter to a spherical harmonic field.
%
% Syntax
%   Plm = sphgaussfilt(Plm, sigma)
%   Plm = sphgaussfilt(Plm, sigma, unit)
%   Plm = sphgaussfilt(Plm, sigma, unit, lmin)
%
% Input arguments
%   Plm - Spherical harmonic coefficients
%   sigma - Standard deviation of the Gaussian filter
%   unit - Unit of the standard deviation
%       'degree' - The standard deviation is in degrees
%       'radian' - The standard deviation is in radians
%       'km' - The standard deviation is in kilometres
%       The default unit is 'degree'
%   lmin - Minimum spherical harmonic degree to apply the filter
%       If given, the filter is applied only to spherical harmonic degrees 
%       at or above this threshold.
%       The default value is none, i.e. the filter is applied to all 
%       spherical harmonic degrees.
%
% Output arguments
%   Plm - Filtered spherical harmonic coefficients
%
% Last modified by
%   2025/07/30, williameclee@arizona.edu (@williameclee)

function Plm = sphgaussfilt(Plm, sigma, unit, options)

    arguments
        Plm (:, 4, :) {isnumeric}
        sigma (1, 1) {mustBeNumeric, mustBePositive} = 5
        unit (1, :) {mustBeMember(unit, {'degree', 'radian', 'km'})} = 'degree'
        options.Lmin (1, 1) {mustBeNumeric, mustBeNonnegative} = Inf
    end

    lmin = options.Lmin;

    if ismatrix(Plm)
        l = Plm(:, 1);
    elseif ndims(Plm) == 3
        l = Plm(:, 1, :);
    end

    switch unit
        case 'degree'
            sigma = 180 / sqrt(2) / sigma;
            gfilter = exp(-l .* (l + 1) / sigma ^ 2/2);
        case 'radian'
            gfilter = exp(-l .* (l + 1) / sigma ^ 2/2);
        case 'km' % modified from Adhikari et al. (2019)
            % sigma = sigma/6371;
            gfilter = exp(- l .* (l + 1) * (sigma / 6371) ^ 2 / (4 * log(2)));
    end

    if lmin < Inf
        gfilter(l < lmin) = 1;
    end

    if ismatrix(Plm)
        Plm(:, end - 1:end) = Plm(:, end - 1:end) .* gfilter;
    elseif ndims(Plm) == 3
        Plm(:, end - 1:end, :) = Plm(:, end - 1:end, :) .* gfilter;
    else
        error('Invalid input dimensions');
    end

end
