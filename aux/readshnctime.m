%% READSHNCTIME - Reads the time coordinate of a spherical-harmonic NetCDF file.
% Returns sorted, unique, unzoned UTC datetimes without reading coefficients.
%
% Syntax
%   dates = readshnctime(path)
%
% Input arguments
%   path - NetCDF file written using CREATESHNC.
%
% Output arguments
%   dates - UTC timestamps. Size: (N, 1), datetime.
%
% Created by
%   2026/10/30, En-Chi Lee (williameclee@gmail.com)

function dates = readshnctime(path)

    arguments (Input)
        path {mustBeTextScalar, mustBeFile}
    end

    arguments (Output)
        dates (:, 1) datetime
    end

    times = ncread(path, 'time');
    assert(strcmpi(ncreadatt(path, 'time', 'units'), 'seconds since 2000-01-01 00:00:00'), ...
        'ULMO:readshnc:InvalidTime', 'Unsupported time coordinate units.');
    assert(~isempty(times) && all(isfinite(times)) && all(diff(times) > 0), ...
        'ULMO:readshnc:InvalidTime', 'Time must be finite, sorted and unique.');
    dates = datetime(2000, 1, 1) + seconds(times(:));
end
