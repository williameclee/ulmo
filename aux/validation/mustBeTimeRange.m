function mustBeTimeRange(timelim)
    eid = 'ULMO:notValidTimeRange';

    if isempty(timelim)
        return
    end

    if ~(isnumeric(timelim) || isdatetime(timelim))
        error(eid, 'Time range must be a numeric or datetime array.');
    end

    if numel(timelim) ~= 2
        error(eid, 'Time range must be a two-element array.');
    end

    if ~(timelim(2) > timelim(1))
        error(eid, 'Time range must be increasing.');
    end

end
