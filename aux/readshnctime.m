%% READSHNCTIME - Reads the time coordinate of a spherical-harmonic NetCDF file.
% Returns sorted, unique, unzoned UTC datetimes without reading coefficients.
%
% Syntax
%   dates = readshnctime(path)
%
% Input arguments
%   path - NetCDF file path with a 'time' argument.
%
% Output arguments
%   dates - UTC timestamps. Size: (N, 1), datetime.
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/10, En-Chi Lee (williameclee@gmail.com)

function dates = readshnctime(path)

    arguments (Input)
        path {mustBeTextScalar, mustBeFile}
    end

    arguments (Output)
        dates (:, 1) datetime
    end

    % Validate the shared time encoding and ordering without loading coefficients.
    times = ncread(path, 'time');
    origin = "seconds since 2000-01-01 00:00:00";
    units = ncreadatt(path, 'time', 'units');
    assert(any(strcmpi(units, [origin, origin + " UTC"])), ...
        'ULMO:readshnc:InvalidTime', 'Unsupported time coordinate units.');
    assert(~isempty(times) && all(isfinite(times)) && all(diff(times) > 0), ...
        'ULMO:readshnc:InvalidTime', 'Time must be finite, sorted and unique.');
    % Return a column of unzoned datetimes representing UTC.
    dates = datetime(2000, 1, 1) + seconds(times(:));
end
