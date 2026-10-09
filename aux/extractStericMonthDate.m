%% extractStericMonthDate
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
        parts = regexp(name, '_(20\d{2})(\d{2})\d{2}', 'tokens', 'once');
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
