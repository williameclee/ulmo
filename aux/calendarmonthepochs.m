%% CALENDARMONTHEPOCHS - Returns calendar-month windows spanning an inclusive date range.
%
% Syntax
%   epochs = calendarmonthepochs(trange)
%
% Input arguments
%   trange - Inclusive [first, last] unzoned timestamps or MATLAB datenums.
%       Equal endpoints select one month; an endpoint on a month boundary includes that month.
%       Size: two-element vector, with first <= last.
%
% Output arguments
%   epochs - First day and exclusive end of each calendar month from the
%       first through the last supplied date, including intervening months.
%       Size: (N, 2).
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)

function epochs = calendarmonthepochs(trange)

    arguments (Input)
        trange {mustBeA(trange, {'datetime', 'numeric'}), mustBeVector}
    end

    arguments (Output)
        epochs (:, 2) datetime
    end

    % Validate the inclusive endpoints and normalise legacy datenum input
    assert(numel(trange) == 2, 'ULMO:calendarmonthepochs:InvalidRange', ...
    'Provide the inclusive first and last dates as a two-element vector.');

    if isnumeric(trange)
        assert(all(isfinite(trange)), 'ULMO:calendarmonthepochs:InvalidDates', ...
        'Dates must be finite.');
        trange = datetime(trange, 'ConvertFrom', 'datenum');
    end

    assert(~any(isnat(trange) | isinf(trange)) && isempty(trange.TimeZone), ...
        'ULMO:calendarmonthepochs:InvalidDates', 'Provide finite unzoned timestamps.');
    assert(trange(1) <= trange(2), 'ULMO:calendarmonthepochs:InvalidRange', ...
    'The first date must not exceed the last date.');
    % Include each endpoint's month, even when the last date is a month boundary
    first = dateshift(trange(1), 'start', 'month');
    last = dateshift(trange(2), 'start', 'month');
    starts = (first:calmonths(1):last)';
    % Use calendar increments so February and unequal month lengths are handled
    epochs = [starts, starts + calmonths(1)];
end
