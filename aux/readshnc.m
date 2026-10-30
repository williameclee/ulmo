%% READSHNC - Reads selected times and harmonic rows from compact NetCDF data.
% Reads time first, then only the coefficient/uncertainty slab in the requested
% range. No coefficient reads occur for an empty time selection.
%
% Syntax
%   [coefficients, uncertainty, dates, labels] = readshnc(path, timelim, degree)
%
% Input arguments
%   path - NetCDF file written using CREATESHNC.
%   timelim (optional) - Inclusive two-element datetime/datenum range.
%       Default range: [] (all times).
%   degree (optional) - Maximum degree to read; never pads stored data.
%       Default degree: [] (native degree).
%   EndExclusive (name-value) - Whether to exclude the range's upper boundary.
%       Default option: false. Used for epoch averaging windows.
%
% Output arguments
%   coefficients, uncertainty - Double arrays, size (M, 2, N).
%       Missing uncertainty variables produce NaNs for the selected slab.
%   dates - Selected UTC timestamps. Size: (N, 1), datetime.
%   labels - Selected degree/order pairs converted to double. Size: (M, 2).
%
% See also
%   CREATESHNC, READSHNCTIME
%
% Created by
%   2026/10/30

function [coefficients, uncertainty, dates, labels] = readshnc(path, timelim, degree, options)
    arguments (Input)
        path {mustBeTextScalar, mustBeFile}
        timelim {mustBeTimeRange} = []
        degree double = []
        options.EndExclusive (1, 1) logical = false
    end
    arguments (Output)
        coefficients (:, 2, :) double
        uncertainty (:, 2, :) double
        dates (:, 1) datetime
        labels (:, 2) double
    end
    dates = readshnctime(path);
    keep = true(size(dates));
    if ~isempty(timelim)
        if isnumeric(timelim), timelim = datetime(timelim, 'ConvertFrom', 'datenum'); end
        keep = dates >= timelim(1) & dates <= timelim(2);
        if options.EndExclusive, keep = keep & dates < timelim(2); end
    end
    ids = find(keep);
    info = ncinfo(path);
    variable = info.Variables(strcmp({info.Variables.Name}, 'coefficients'));
    assert(isscalar(variable) && isequal({variable.Dimensions.Name}, {'lm', 'cs', 'time'}) && ...
        variable.Size(2) == 2 && variable.Size(3) == numel(dates), ...
        'ULMO:readshnc:InvalidDimensions', 'Expected coefficients(lm, cs, time).');
    errorVariable = info.Variables(strcmp({info.Variables.Name}, 'uncertainty'));
    if ~isempty(errorVariable)
        assert(isequal({errorVariable.Dimensions.Name}, {'lm', 'cs', 'time'}) && ...
            isequal(errorVariable.Size, variable.Size), 'ULMO:readshnc:InvalidDimensions', ...
            'Uncertainty dimensions must match coefficients.');
    end
    nRows = variable.Size(1);
    if ~isempty(degree)
        validateattributes(degree, {'double'}, {'scalar', 'integer', 'nonnegative', 'finite'});
        nRows = min(nRows, (degree + 1) * (degree + 2) / 2);
    end
    labels = [double(ncread(path, 'degree', 1, nRows)), double(ncread(path, 'order', 1, nRows))];
    coefficients = nan(nRows, 2, numel(ids));
    uncertainty = coefficients;
    dates = dates(keep);
    if isempty(ids), return; end
    start = [1, 1, ids(1)];
    count = [nRows, 2, numel(ids)];
    coefficients = reshape(ncread(path, 'coefficients', start, count), count);
    if ~isempty(errorVariable)
        uncertainty = reshape(ncread(path, 'uncertainty', start, count), count);
    end
end
