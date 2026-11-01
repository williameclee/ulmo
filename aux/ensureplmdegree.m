%% ENSUREPLMDEGREE - Pads or truncates a Plm matrix to ensure it has the correct degree.
%
% Syntax
%   plm = ENSUREPLMDEGREE(plm, L, fmt)
%
% Input arguments
%   plm - Four-column [degree, order, cosine, sine] coefficients or
%       two-column [cosine, sine] coefficients, with optional epochs.
%   L - Nonnegative integer target spherical-harmonic degree.
%   fmt - 'traditional': [harmonic, column, epoch];
%       'timefirst': [epoch, harmonic, column] (default).
%
% Output arguments
%   plm - Resized coefficients, preserving all epochs and the column layout.
%       Added cosine/sine coefficients are zero; four-column inputs also
%       receive degree/order metadata for the padded rows.
%
% Created by
%	2025/08/14, williameclee@arizona.edu
%
% Last modified
%	2026/10/24, williameclee@gmail.com

function plm = ensureplmdegree(plm, L, fmt)

    arguments
        plm (:, :, :) double
        L (1, 1) {mustBeFinite, mustBeInteger, mustBeNonnegative}
        fmt (1, :) char {mustBeMember(fmt, {'timefirst', 'traditional'})} = 'timefirst'
    end

    isTraditional = ismember(size(plm, 2), [2, 4]);
    isTimefirst = ismember(size(plm, 3), [2, 4]) && ~isTraditional;

    if strcmpi(fmt, 'traditional') && ~isTraditional

        if isTimefirst
            fmt = 'timefirst';
            warning('ULMO:InvalidInputSize', ...
                'PLM appears to be in TIMEFIRST format based on its size ([%s]). Proceeding as such.', ...
                strjoin(string(size(plm)), ' x '))
        else
            error('ULMO:InvalidInputSize', ...
                'When using TRADITIONAL format, PLM must have 2 or 4 columns ([addmup(L) x (2 or 4) x n]), got %s instead.', ...
                strjoin(string(size(plm)), ' x '))
        end

    elseif strcmpi(fmt, 'timefirst') && ~isTimefirst

        if isTraditional
            fmt = 'traditional';
            warning('ULMO:InvalidInputSize', ...
                'PLM appears to be in TRADITIONAL format based on its size ([%s]). Proceeding as such.', ...
                strjoin(string(size(plm)), ' x '))
        else
            error('ULMO:InvalidInputSize', ...
                'When using TIMEFIRST format, PLM must have 2 or 4 columns ([n x addmup(L) x (2 or 4)]), got %s instead.', ...
                strjoin(string(size(plm)), ' x '))
        end

    end

    % Preserve the original array without permuting when no resizing is needed.
    harmonicDim = 1 + strcmpi(fmt, 'timefirst');

    if size(plm, harmonicDim) == addmup(L)
        return
    end

    if strcmpi(fmt, 'timefirst')
        plm = permute(plm, [2, 3, 1]); % timefirst -> traditional
    end

    if size(plm, 1) > addmup(L)
        % Truncate
        plm = plm(1:addmup(L), :, :);
    elseif size(plm, 2) == 2
        % Coefficient-only arrays have no degree/order columns to populate.
        plm(addmup(L), 2, size(plm, 3)) = 0;
    else
        [order, degree] = addmon(L);
        plm(1:addmup(L), 1, :) = repmat(degree(:), [1, 1, size(plm, 3)]);
        plm(1:addmup(L), 2, :) = repmat(order(:), [1, 1, size(plm, 3)]);
    end

    if strcmpi(fmt, 'timefirst')
        plm = permute(plm, [3, 1, 2]); % traditional -> timefirst
    end

end
