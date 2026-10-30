%% WRITESHNC - Saves compact harmonic coefficients and uncertainties to NetCDF.
% Writes a temporary file and replaces the destination only on success.
%
% Syntax
%   writeshnc(path, coefficients, uncertainty, labels, dates, unit)
%
% Input arguments
%   path - Output NetCDF file.
%   coefficients - Double array, size (M, 2, N).
%   uncertainty - Same shape as coefficients, or [] to omit the variable.
%   labels - Degree/order pairs, size (M, 2).
%   dates - Unzoned UTC datetimes, size (N, 1).
%   unit - GRAV, POT or SD.
%
% See also
%   CREATESHNC, READSHNC
%
% Created by
%   2026/10/30

function writeshnc(path, coefficients, uncertainty, labels, dates, unit)
    arguments (Input)
        path {mustBeTextScalar}
        coefficients (:, 2, :) double
        uncertainty double
        labels (:, 2) double
        dates (:, 1) datetime
        unit {mustBeTextScalar}
    end
    assert(size(coefficients, 1) == size(labels, 1) && size(coefficients, 3) == numel(dates) && ...
        (isempty(uncertainty) || isequal(size(coefficients), size(uncertainty))), ...
        'ULMO:writeshnc:InvalidDimensions', 'Coefficient, label, date and uncertainty sizes must match.');
    folder = fileparts(path);
    if isempty(folder), folder = pwd; end
    temporary = string(tempname(folder)) + ".nc";
    cleanup = onCleanup(@() removetemporary(temporary));
    createshnc(temporary, labels, dates, unit, 'Uncertainty', ~isempty(uncertainty));
    ncwrite(temporary, 'coefficients', coefficients);
    if ~isempty(uncertainty), ncwrite(temporary, 'uncertainty', uncertainty); end
    movefile(temporary, path, 'f');
end

function removetemporary(path)
    if isfile(path), delete(path); end
end
