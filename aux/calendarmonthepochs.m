%% CALENDARMONTHEPOCHS - Returns calendar-month windows spanning an inclusive date range.
%
% Syntax
%   epochs = calendarmonthepochs(timerange)
%
% Input arguments
%   timerange - Inclusive [first, last] unzoned timestamps or MATLAB datenums.
%       Size: two-element vector, with first <= last. Equal endpoints select
%       one month; an endpoint on a month boundary includes that month.
%
% Output arguments
%   epochs - First day and exclusive end of each calendar month from the
%       first through the last supplied date, including intervening months.
%       Size: (N, 2), datetime.
%
% See also
%   AVERAGEEPOCHS
%
% Created by
%   2026/10/30

function epochs = calendarmonthepochs(timerange)
    arguments (Input)
        timerange {mustBeA(timerange, {'datetime', 'numeric'}), mustBeVector}
    end
    arguments (Output)
        epochs (:, 2) datetime
    end
    assert(numel(timerange) == 2, 'ULMO:calendarmonthepochs:InvalidRange', ...
        'Provide the inclusive first and last dates as a two-element vector.');
    if isnumeric(timerange)
        assert(all(isfinite(timerange)), 'ULMO:calendarmonthepochs:InvalidDates', ...
            'Dates must be finite.');
        timerange = datetime(timerange, 'ConvertFrom', 'datenum');
    end
    assert(~any(isnat(timerange) | isinf(timerange)) && isempty(timerange.TimeZone), ...
        'ULMO:calendarmonthepochs:InvalidDates', 'Provide finite unzoned timestamps.');
    assert(timerange(1) <= timerange(2), 'ULMO:calendarmonthepochs:InvalidRange', ...
        'The first date must not exceed the last date.');
    first = dateshift(timerange(1), 'start', 'month');
    last = dateshift(timerange(2), 'start', 'month');
    starts = (first:calmonths(1):last)';
    epochs = [starts, starts + calmonths(1)];
end
