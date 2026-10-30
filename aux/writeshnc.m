%% WRITESHNC - Saves lmcosi or compact harmonic arrays to NetCDF.
% Degree/order are stored once; both input layouts produce the same schema.
%
% Syntax
%   writeshnc(path, plmt, uncs, dates, unit)
%   writeshnc(path, coeffs, uncs, dates, unit, "Labels", labels)
%   writeshnc(__, "Attributes", {"equatorRadius", radius, "source", "ESAESM2"})
%
% Input arguments
%   path - Output NetCDF file.
%   coeffs - Double array, size (M, 4, N) for lmcosi (degree, order,
%       cosine, sine), or (M, 2, N) for compact cosine/sine data.
%   uncs - Same layout/shape as coefficients, or [] to omit the variable.
%       Lmcosi uncertainty labels must match the coefficient labels.
%   dates - Unzoned UTC datetimes, size (N, 1).
%   unit - GRAV, POT or SD.
%   Labels (name-value) - Degree/order pairs, size (M, 2).
%       Required for compact input; inferred from lmcosi input. If supplied
%       with lmcosi, must agree. Labels must be constant across time.
%       Default labels: [] (infer from lmcosi).
%   Attributes (name-value) - Additional global attributes as alternating
%       name/value cells; passed to CREATESHNC. Default attributes: {}.
%
% Notes
%   The function writes a temporary file and replaces the destination only
%   on success.
%
% See also
%   CREATESHNC, READSHNC
%
% Created by
%   2026/10/30, En-Chi Lee (williameclee@gmail.com)

function writeshnc(path, coeffs, uncs, dates, unit, options)

    arguments (Input)
        path {mustBeTextScalar}
        coeffs (:, :, :) double
        uncs double
        dates (:, 1) datetime
        unit {mustBeTextScalar}
        options.Labels (:, 2) double = zeros(0, 2)
        options.Attributes (1, :) cell = {}
    end

    assert(~isempty(dates) && ismember(size(coeffs, 2), [2, 4]) && ...
        size(coeffs, 3) == numel(dates) && ...
        (isempty(uncs) || isequal(size(coeffs), size(uncs))), ...
        'ULMO:writeshnc:InvalidDimensions', ...
    'Coefficient, date and uncertainty sizes must match.');
    labels = options.Labels;

    if size(coeffs, 2) == 4
        embeddedLabels = coeffs(:, 1:2, 1);
        assert(all(coeffs(:, 1:2, :) == embeddedLabels, 'all') && ...
            (isempty(labels) || isequal(labels, embeddedLabels)) && ...
            (isempty(uncs) || all(uncs(:, 1:2, :) == embeddedLabels, 'all')), ...
            'ULMO:writeshnc:InvalidLabels', ...
        'Lmcosi labels must agree and remain constant across time.');
        labels = embeddedLabels;
        coeffs = coeffs(:, 3:4, :);

        if ~isempty(uncs)
            uncs = uncs(:, 3:4, :);
        end

    end

    assert(~isempty(labels) && size(coeffs, 1) == size(labels, 1), ...
        'ULMO:writeshnc:InvalidLabels', 'Compact input requires matching Labels.');
    folder = fileparts(path);

    if isempty(folder)
        folder = pwd;
    end

    tmpPath = string(tempname(folder)) + ".nc";
    cleanup = onCleanup(@() removetemporary(tmpPath));
    createshnc(tmpPath, labels, dates, unit, 'Uncertainty', ~isempty(uncs), ...
        'Attributes', options.Attributes);
    ncwrite(tmpPath, 'coefficients', coeffs);

    if ~isempty(uncs)
        ncwrite(tmpPath, 'uncertainty', uncs);
    end

    movefile(tmpPath, path, 'f');
end

function removetemporary(path)

    if isfile(path)
        delete(path);
    end

end
