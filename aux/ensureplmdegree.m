%% ENSUREPLMDEGREE
% Pad or truncate a plm matrix to ensure it has the correct degree.
% The input format is assumed to be strictly timefirst (t, lmcosi).
%
% See also
%	SOLVESLE, SOLVEDEGREE1
%
% Notes
%	This is a helper function with limited documentation.%
%
% Authored by
%	2025/08/14, williameclee@arizona.edu (@williameclee)

function plm = ensureplmdegree(plm, L, fmt)

    arguments
        plm (:, :, :) double
        L (1, 1) {mustBeFinite, mustBeInteger, mustBeNonnegative}
        fmt (1, :) char {mustBeMember(fmt, {'timefirst', 'traditional'})} = 'timefirst'
    end

    if strcmpi(fmt, 'traditional') && size(plm, 2) ~= 4

        if size(plm, 3) == 4
            fmt = 'timefirst';
            warning('ULMO:InvalidInputSize', ...
                'PLM appears to be in TIMEFIRST format based on its size ([%s]). Proceeding as such.', ...
                strjoin(string(size(plm)), ' x '))
        else
            error('ULMO:InvalidInputSize', ...
                'When using TRADITIONAL format, PLM must have 4 columns ([addmup(L) x 4 x (n)]), got %s instead.', ...
                strjoin(string(size(plm)), ' x '))
        end

    elseif strcmpi(fmt, 'timefirst') && size(plm, 3) ~= 4

        if size(plm, 2) == 4
            fmt = 'traditional';
            warning('ULMO:InvalidInputSize', ...
                'PLM appears to be in TRADITIONAL format based on its size ([%s]). Proceeding as such.', ...
                strjoin(string(size(plm)), ' x '))
        else
            error('ULMO:InvalidInputSize', ...
                'When using TIMEFIRST format, PLM must have 4 columns ([n x addmup(L) x 4]), got %s instead.', ...
                strjoin(string(size(plm)), ' x '))
        end

    end

    if ndims(plm) == 3 & strcmpi(fmt, 'timefirst')
        plm = permute(plm, [2, 3, 1]); % timefirst -> traditional
    end

    if size(plm, 1) == addmup(L)
        % Do nothing, correct size
    elseif size(plm, 1) > addmup(L)
        % Truncate
        plm = plm(:, 1:addmup(L), :);
    else
        [order, degree] = addmon(L);
        plm(1:addmup(L), 1, :) = repmat(degree(:), [1, 1, size(plm, 3)]);
        plm(1:addmup(L), 2, :) = repmat(order(:), [1, 1, size(plm, 3)]);
    end

    if ndims(plm) == 3 & strcmpi(fmt, 'timefirst')
        plm = permute(plm, [3, 1, 2]); % traditional -> timefirst
    end

end
