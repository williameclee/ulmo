%% CALENDARMONTHEPOCHS - Returns calendar-month windows spanning supplied dates.
%
% Syntax
%   epochs = calendarmonthepochs(dates)
%
% Input arguments
%   dates - Nonempty vector of unzoned timestamps or MATLAB datenums.
%       Dates need not be sorted or regularly sampled.
%
% Output arguments
%   epochs - First day and exclusive end of each calendar month from the
%       earliest through the latest supplied date, including intervening months.
%       Size: (N, 2), datetime.
%
% See also
%   AVERAGEEPOCHS
%
% Created by
%   2026/10/30

function epochs = calendarmonthepochs(dates)
    arguments (Input)
        dates {mustBeA(dates, {'datetime', 'numeric'}), mustBeVector}
    end
    arguments (Output)
        epochs (:, 2) datetime
    end
    if isnumeric(dates)
        assert(all(isfinite(dates)), 'ULMO:calendarmonthepochs:InvalidDates', ...
            'Dates must be finite.');
        dates = datetime(dates, 'ConvertFrom', 'datenum');
    end
    assert(~isempty(dates) && ~any(isnat(dates) | isinf(dates)) && ...
        isempty(dates.TimeZone), 'ULMO:calendarmonthepochs:InvalidDates', ...
        'Provide nonempty, finite unzoned timestamps.');
    first = dateshift(min(dates), 'start', 'month');
    last = dateshift(max(dates), 'start', 'month');
    starts = (first:calmonths(1):last)';
    epochs = [starts, starts + calmonths(1)];
end
