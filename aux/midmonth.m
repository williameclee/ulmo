function midMonth = midmonth(range)

    arguments
        range (1, 2) datetime
    end

    start = datetime(year(range(1)), month(range(1)), 1);
    last = datetime(year(range(2)), month(range(2)), 1) + calmonths(1);
    dmonths = ceil((last - start) / days(28));
    startMonth = start + calmonths(0:dmonths - 1);
    endMonth = start + calmonths(1:dmonths);
    midMonth = startMonth + (endMonth - startMonth) / 2;
    midMonth = midMonth(midMonth >= range(1) & midMonth <= range(2));
end
