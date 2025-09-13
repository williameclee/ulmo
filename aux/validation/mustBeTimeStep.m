function mustBeTimeStep(timeStep)
    eid = 'ULMO:notValidTimeStep';

    if isempty(timeStep)
        return
    end

    if (ischar(timeStep) || isstring(timeStep))

        if ~ismember(timeStep, {'midmonth'})
            error(eid, 'Time step string must be ''midmonth''. Got "%s".', timeStep);
        end

    elseif ~isnumeric(timeStep) || isduration(timeStep)

        if ~(isscalar(timeStep) && timeStep > 0)
            error(eid, 'Time step must be a positive scalar.');
        end

    else
        error(eid, 'Time step must be a string, numeric, or duration scalar. Got %s.', class(timeStep));
    end

end
