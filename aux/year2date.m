function time = year2date(timeInYears)
    theYear = floor(timeInYears);
    timeInDays = (timeInYears - theYear) .* ...
        (365 + double(isleap(theYear)));
    time = datetime(theYear, 1, 1) + days(timeInDays);
end
