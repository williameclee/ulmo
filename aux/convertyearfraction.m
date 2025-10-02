%% CONVERTYEARFRACTION
% Convert a year fraction to a datetime object.
%
% Authored by:
%   2025/07/21, williameclee@arizona.edu (@williameclee)

function dt = convertyearfraction(fyear)
    calYear = floor(fyear);
    calFDday = fyear - calYear;

    isLeap = eomday(calYear, 2) == 29; % Check if the year is a leap year

    dt = datetime(calYear, 1, 1, "Format", "yyyy/MM/dd hh:mm:ss") ...
        + calFDday .* (365 + isLeap); % Convert to datetime

end
