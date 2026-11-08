%% MUSTBETIMERANGE - Validates the bounds of a time range.
%
% Syntax
%   mustBeTimeRange(timelim)
%
% Input arguments
%   timelim - Empty for the full record, or two datetime values or serial
%       date numbers in nondecreasing order. Infinite bounds are permitted.
%
% See also
%   ssh2slept, ssh2lonlatt, readshnc
%
% Created by
%   2025/10/01, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function mustBeTimeRange(timelim)
    eid = 'ULMO:notValidTimeRange';

    if isempty(timelim)
        return
    end

    if ~(isnumeric(timelim) || isdatetime(timelim))
        error(eid, 'Time range must be a numeric or datetime array.');
    elseif numel(timelim) ~= 2
        error(eid, 'Time range must be a two-element array.');
    elseif ~(timelim(2) >= timelim(1))
        error(eid, 'The end date must be on or after the start date.');
    end

end
