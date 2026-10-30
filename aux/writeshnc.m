%% WRITESHNC - Saves lmcosi or compact harmonic arrays to NetCDF.
% Writes a temporary file and replaces the destination only on success.
% Degree/order are stored once; both input layouts produce the same schema.
%
% Syntax
%   writeshnc(path, plmt, uncertainty, dates, unit)
%   writeshnc(path, coefficients, uncertainty, dates, unit, "Labels", labels)
%   writeshnc(__, "Attributes", {"equatorRadius", radius, "source", "ESAESM2"})
%
% Input arguments
%   path - Output NetCDF file.
%   coefficients - Double array, size (M, 4, N) for lmcosi (degree, order,
%       cosine, sine), or (M, 2, N) for compact cosine/sine data.
%   uncertainty - Same layout/shape as coefficients, or [] to omit the variable.
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
% See also
%   CREATESHNC, READSHNC
%
% Created by
%   2026/10/30

function writeshnc(path, coefficients, uncertainty, dates, unit, options)
    arguments (Input)
        path {mustBeTextScalar}
        coefficients (:, :, :) double
        uncertainty double
        dates (:, 1) datetime
        unit {mustBeTextScalar}
        options.Labels (:, 2) double = zeros(0, 2)
        options.Attributes (1, :) cell = {}
    end
    assert(~isempty(dates) && ismember(size(coefficients, 2), [2, 4]) && ...
        size(coefficients, 3) == numel(dates) && ...
        (isempty(uncertainty) || isequal(size(coefficients), size(uncertainty))), ...
        'ULMO:writeshnc:InvalidDimensions', 'Coefficient, date and uncertainty sizes must match.');
    labels = options.Labels;
    if size(coefficients, 2) == 4
        embeddedLabels = coefficients(:, 1:2, 1);
        assert(all(coefficients(:, 1:2, :) == embeddedLabels, 'all') && ...
            (isempty(labels) || isequal(labels, embeddedLabels)) && ...
            (isempty(uncertainty) || all(uncertainty(:, 1:2, :) == embeddedLabels, 'all')), ...
            'ULMO:writeshnc:InvalidLabels', 'Lmcosi labels must agree and remain constant across time.');
        labels = embeddedLabels;
        coefficients = coefficients(:, 3:4, :);
        if ~isempty(uncertainty), uncertainty = uncertainty(:, 3:4, :); end
    end
    assert(~isempty(labels) && size(coefficients, 1) == size(labels, 1), ...
        'ULMO:writeshnc:InvalidLabels', 'Compact input requires matching Labels.');
    folder = fileparts(path);
    if isempty(folder), folder = pwd; end
    temporary = string(tempname(folder)) + ".nc";
    cleanup = onCleanup(@() removetemporary(temporary));
    createshnc(temporary, labels, dates, unit, 'Uncertainty', ~isempty(uncertainty), ...
        'Attributes', options.Attributes);
    ncwrite(temporary, 'coefficients', coefficients);
    if ~isempty(uncertainty), ncwrite(temporary, 'uncertainty', uncertainty); end
    movefile(temporary, path, 'f');
end

function removetemporary(path)
    if isfile(path), delete(path); end
end
