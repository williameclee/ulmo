%% CREATESHNC - Creates a compact spherical-harmonic NetCDF file.
% Coefficients and optional uncertainties have MATLAB dimension order
% (lm, cs, time); degree/order labels are stored once as uint16 vectors.
%
% Syntax
%   createshnc(path, labels, dates, unit)
%   createshnc(__, "Uncertainty", true, ...
%       "Attributes", {"equatorRadius", radius, "source", "ESAESM2"})
%
% Input arguments
%   path - New NetCDF file; existing files are rejected.
%   labels - Ascending degree/order pairs. Size: (M, 2).
%   dates - Sorted unique unzoned timestamps. Size: (N, 1).
%   unit - Coefficient unit: GRAV, POT or SD.
%   Uncertainty (name-value) - Whether to create an uncertainty variable.
%       Default option: false.
%   Attributes (name-value) - Additional global attributes as a cell vector
%       {name1, value1, name2, value2, ...}. Names accept char/string scalars.
%       coefficientUnit is reserved for unit. Default attributes: {}.
%
% Output file content
%   degree, order - uint16 labels on dimension lm.
%   cs - uint8 coordinate: 0 cosine, 1 sine.
%   time - Double seconds since 2000-01-01 00:00:00 UTC.
%   coefficients, uncertainty (optional) - Double arrays, size (M, 2, N).
%       Missing values: NaN. Use NCWRITE to fill time slabs after creation.
%   coefficientUnit (global attribute) - Coefficient unit supplied by unit.
%   Additional global attributes are taken from Attributes.
%
% See also
%   READSHNC, READSHNCTIME, PREPAREESM
%
% Created by
%   2026/10/30

function createshnc(path, labels, dates, unit, options)
    arguments (Input)
        path {mustBeTextScalar}
        labels (:, 2) double {mustBeInteger, mustBeNonnegative, mustBeFinite}
        dates (:, 1) datetime
        unit {mustBeTextScalar, mustBeMember(unit, ["GRAV", "POT", "SD"])}
        options.Uncertainty (1, 1) logical = false
        options.Attributes (1, :) cell {mustBeAttributes} = {}
    end
    assert(~isfile(path), 'ULMO:createshnc:FileExists', 'File already exists: %s', path);
    assert(~isempty(labels) && all(labels <= double(intmax('uint16')), 'all'), ...
        'ULMO:createshnc:InvalidLabels', 'Labels must fit uint16.');
    assert(~isempty(dates) && all(isfinite(dates)) && isempty(dates.TimeZone) && ...
        all(diff(dates) > seconds(0)), 'ULMO:createshnc:InvalidDates', ...
        'Dates must be nonempty, finite, unzoned, sorted and unique.');
    nRows = size(labels, 1);
    nTimes = numel(dates);
    names = {'degree', 'order'};
    for column = 1:2
        nccreate(path, names{column}, 'Dimensions', {'lm', nRows}, ...
            'Datatype', 'uint16', 'Format', 'netcdf4');
        ncwrite(path, names{column}, uint16(labels(:, column)));
    end
    nccreate(path, 'cs', 'Dimensions', {'cs', 2}, 'Datatype', 'uint8');
    ncwrite(path, 'cs', uint8([0; 1]));
    ncwriteatt(path, 'cs', 'description', '0=cosine, 1=sine');
    nccreate(path, 'time', 'Dimensions', {'time', nTimes}, 'Datatype', 'double');
    ncwrite(path, 'time', seconds(dates - datetime(2000, 1, 1)));
    ncwriteatt(path, 'time', 'units', 'seconds since 2000-01-01 00:00:00');
    ncwriteatt(path, 'time', 'calendar', 'proleptic_gregorian');
    names = "coefficients";
    if options.Uncertainty
        names(end + 1) = "uncertainty";
    end
    for name = names
        nccreate(path, name, 'Dimensions', {'lm', nRows, 'cs', 2, 'time', nTimes}, ...
            'Datatype', 'double', 'FillValue', NaN, 'DeflateLevel', 4, 'Shuffle', true, ...
            'ChunkSize', [min(nRows, 1024), 2, min(nTimes, 8)]);
    end
    ncwriteatt(path, '/', 'coefficientUnit', char(unit));
    for k = 1:2:numel(options.Attributes)
        ncwriteatt(path, '/', options.Attributes{k}, options.Attributes{k + 1});
    end
end

function mustBeAttributes(attributes)
    assert(mod(numel(attributes), 2) == 0, 'ULMO:createshnc:InvalidAttributes', ...
        'Attributes must contain name/value pairs.');
    names = strings(1, numel(attributes) / 2);
    for k = 1:numel(names)
        name = attributes{2 * k - 1};
        mustBeTextScalar(name);
        names(k) = string(name);
    end
    assert(all(~ismissing(names) & strlength(names) > 0) && ...
        numel(unique(names)) == numel(names) && ~any(names == "coefficientUnit"), ...
        'ULMO:createshnc:InvalidAttributes', ...
        'Attribute names must be nonempty and unique; coefficientUnit is reserved.');
end
