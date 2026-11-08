%% MUSTBETIMESTEP - Validates a time step for gridded sea-level data.
%
% Syntax
%   MUSTBETIMESTEP(timeStep)
%
% Input arguments
%   timeStep - Native sampling ([]), a positive finite scalar number of
%       days or duration, or a text scalar naming 'midmonth' or 'GRACE'.
%       Named steps are case-insensitive.
%
% See also
%   normalisetimeconfig, processfieldtime
%
% Created by
%   2025/10/01, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/07, En-Chi Lee (williameclee@gmail.com)

function mustBeTimeStep(timeStep)
    eid = 'ULMO:notValidTimeStep';

    if isempty(timeStep)
        return
    end

    if (ischar(timeStep) || isstring(timeStep))
        mustBeTextScalar(timeStep)

        if ~strcmpi(timeStep, {'midmonth', 'GRACE'})
            error(eid, 'Time step string must be ''midmonth'' or ''GRACE''. Got "%s".', timeStep);
        end

    elseif isnumeric(timeStep) || isduration(timeStep)

        if ~(isscalar(timeStep) && isfinite(timeStep) && timeStep > 0)
            error(eid, 'Time step must be a positive scalar.');
        end

    else
        error(eid, 'Time step must be a string, numeric, or duration scalar. Got %s.', class(timeStep));
    end

end
