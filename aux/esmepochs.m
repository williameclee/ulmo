%% ESMEPOCHS - Defines averaging windows for an ESM time step.
% Keeps epoch selection separate from coefficient averaging and caching.
%
% Syntax
%   epochs = esmepochs(timestep, nativeDates)
%
% Input arguments
%   timestep - Averaging time step. Currently supports 'midmonth'.
%   nativeDates - Available six-hour timestamps, in unzoned UTC.
%       Size: (N, 1), datetime.
%
% Output arguments
%   epochs - Start and exclusive end of each averaging window.
%       Size: (M, 2), datetime. Midpoints label the averaged output.
%
% Notes
%   Calendar windows span every month containing or between native dates.
%   A future GRACE option can supply real/fallback windows here, using the
%   same N-by-2 contract as GRACEEPOCHS, without changing the averaging code.
%   GRACE epoch discovery is not implemented by this helper yet.
%
% See also
%   ESM2PLMT, AVERAGEESMEPOCHS
%
% Created by
%   2026/10/30

function epochs = esmepochs(timestep, nativeDates)
    arguments (Input)
        timestep (1, 1) string
        nativeDates (:, 1) datetime
    end
    assert(~isempty(nativeDates) && ~any(isnat(nativeDates)) && ...
        isempty(nativeDates.TimeZone), 'ULMO:esmepochs:InvalidDates', ...
        'Provide nonempty, valid unzoned UTC timestamps.');
    switch timestep
        case "midmonth"
            first = dateshift(min(nativeDates), 'start', 'month');
            last = dateshift(max(nativeDates), 'start', 'month');
            starts = (first:calmonths(1):last)';
            epochs = [starts, starts + calmonths(1)];
        otherwise
            error('ULMO:esm2plmt:InvalidTimeStep', ...
                'Unsupported ESM averaging time step: %s.', timestep);
    end
end
