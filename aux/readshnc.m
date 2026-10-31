%% READSHNC - Reads selected times and harmonic rows from compact NetCDF data.
% Reads time first, then only the coefficient/uncertainty slab in the requested
% range. Returns lmcosi arrays by default. No coefficient reads occur for an
% empty time selection.
%
% Syntax
%   [plmt, uncs, dates] = readshnc(path, timelim, Lmax)
%   [coeffs, uncs, dates, labels] = readshnc(__, "OutputFormat", "compact")
%
% Input arguments
%   path - NetCDF file written using CREATESHNC.
%   timelim (optional) - Inclusive two-element datetime/datenum range.
%       Default range: [] (all times).
%   Lmax (optional) - Maximum degree to read; never pads stored data.
%       Default degree: [] (native degree).
%   OutputFormat (name-value) - Returned coefficient/uncertainty layout.
%       'lmcosi': degree, order, cosine, sine columns; 'compact': cosine/sine.
%       Default format: 'lmcosi'.
%   EndExclusive (name-value) - Whether to exclude the range's upper boundary.
%       Default option: false. Used for epoch averaging windows.
%
% Output arguments
%   coeffs, uncs - Double arrays, size (M, 4, N) for lmcosi,
%       or (M, 2, N) for compact. Lmcosi uncertainties carry valid degree/order
%       columns; unavailable errors are NaN only in the cosine/sine columns.
%   dates - Selected UTC timestamps. Size: (N, 1), datetime.
%   labels - Selected degree/order pairs converted to double. Size: (M, 2).
%
% See also
%   CREATESHNC, READSHNCTIME
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)

function [coeffs, uncs, dates, labels] = readshnc(path, timelim, Lmax, options)

    arguments (Input)
        path {mustBeTextScalar, mustBeFile}
        timelim {mustBeTimeRange} = []
        Lmax double = []
        options.OutputFormat ...
            {mustBeTextScalar, mustBeMember(options.OutputFormat, ["lmcosi", "compact"])} = 'lmcosi'
        options.EndExclusive (1, 1) logical = false
    end

    arguments (Output)
        coeffs (:, :, :) double
        uncs (:, :, :) double
        dates (:, 1) datetime
        labels (:, 2) double
    end

    % Read the small time coordinate first to determine the required data slab
    dates = readshnctime(path);
    keep = true(size(dates));

    if ~isempty(timelim)

        if isnumeric(timelim)
            timelim = datetime(timelim, 'ConvertFrom', 'datenum');
        end

        keep = dates >= timelim(1) & dates <= timelim(2);

        if options.EndExclusive
            keep = keep & dates < timelim(2);
        end

    end

    % Check named dimensions before interpreting coefficient or uncertainty axes
    ids = find(keep);
    fileInfo = ncinfo(path);
    varInfo = fileInfo.Variables(strcmp({fileInfo.Variables.Name}, 'coefficients'));
    assert(isscalar(varInfo) && isequal({varInfo.Dimensions.Name}, {'lm', 'cs', 'time'}) && ...
        varInfo.Size(2) == 2 && varInfo.Size(3) == numel(dates), ...
        'ULMO:readshnc:InvalidDimensions', 'Expected coefficients(lm, cs, time).');
    uncInfo = fileInfo.Variables(strcmp({fileInfo.Variables.Name}, 'uncertainty'));

    if ~isempty(uncInfo)
        assert(isequal({uncInfo.Dimensions.Name}, {'lm', 'cs', 'time'}) && ...
            isequal(uncInfo.Size, varInfo.Size), 'ULMO:readshnc:InvalidDimensions', ...
        'Uncertainty dimensions must match coefficients.');
    end

    % For degree-major triangular ordering, degrees through L occupy the first
    % (L+1)*(L+2)/2 rows. Clamp the request to the stored degree.
    nRows = varInfo.Size(1);

    if ~isempty(Lmax)
        validateattributes(Lmax, {'double'}, {'scalar', 'integer', 'nonnegative', 'finite'});
        nRows = min(nRows, (Lmax + 1) * (Lmax + 2) / 2);
    end

    % Read labels once and leave unavailable uncertainties as NaN
    labels = [double(ncread(path, 'degree', 1, nRows)), double(ncread(path, 'order', 1, nRows))];
    coeffs = nan(nRows, 2, numel(ids));
    uncs = coeffs;
    dates = dates(keep);

    % Sorted times make the selected interval contiguous. Empty selections skip
    % coefficient reads; reshape restores dimensions dropped for singleton slabs.
    if ~isempty(ids)
        start = [1, 1, ids(1)];
        count = [nRows, 2, numel(ids)];
        coeffs = reshape(ncread(path, 'coefficients', start, count), count);

        if ~isempty(uncInfo)
            uncs = reshape(ncread(path, 'uncertainty', start, count), count);
        end

    end

    % Expand labels only for the public lmcosi layout, including uncertainty arrays
    if ~strcmp(options.OutputFormat, 'lmcosi')
        return
    end

    repeatedLabels = repmat(labels, 1, 1, numel(ids));
    coeffs = cat(2, repeatedLabels, coeffs);
    uncs = cat(2, repeatedLabels, uncs);

end
