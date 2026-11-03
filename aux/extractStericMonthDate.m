%% EXTRACTSTERICMONTHDATE - Extracts a calendar-month midpoint from a source name.
% Recognises year_YYYY_month_MM or _YYYYMMDD in the filename stem and
% returns the midpoint of the identified calendar month.
%
% Syntax
%   date = extractStericMonthDate(name)
%
% Input arguments
%   name - Source filename or path.
%       The file need not exist.
%
% Output arguments
%   date - Datetime at the midpoint between the first day
%       of the identified month and the first day of the following month.
%       Odd-length months have a midpoint at noon.
%
% See also
%   saveSourceStericMonth
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function date = extractStericMonthDate(name)

    arguments (Input)
        name {mustBeTextScalar}
    end

    arguments (Output)
        date (1, 1) datetime
    end

    [~, name] = fileparts(name);
    parts = regexp(name, 'year_(\d{4})_month_(\d{2})', 'tokens', 'once');

    if isempty(parts)
        parts = regexp(name, '_(\d{4})(\d{2})\d{2}', 'tokens', 'once');
    end

    assert(~isempty(parts), 'ULMO:extractStericMonthDate:InvalidFileName', ...
        'Cannot identify month from file name %s.', name);
    year = str2double(parts{1});
    month = str2double(parts{2});
    assert (month >= 1 && month <= 12, 'ULMO:extractStericMonthDate:InvalidFileName', ...
        ['Expected month to be between 1 and 12. ', ...
     'But extracted month %.0f from file name %s.'], ...
        month, name);
    start = datetime(year, month, 1);
    date = start + (start + calmonths(1) - start) / 2;
end
