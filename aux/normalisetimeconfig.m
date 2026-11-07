%% NORMALISETIMECONFIG - Normalises temporal processing and projection identity.
%
% Syntax
%   [tstep, method, useMean, config] = NORMALISETIMECONFIG(tstep, mthd)
%
% Input arguments
%   tstep - Native ([]), positive numeric days or duration, 'midmonth',
%       or 'GRACE'. Numeric days become durations; named steps become lower
%       case scalar strings. Duration display formats are preserved.
%   mthd - Temporal interpolation method or 'mean'; converted to lower
%       case. Methods are validated by the processing function when used.
%
% Output arguments
%   tstep, mthd - Normalised arguments in the conventions above.
%   useMean - Whether named epochs or fixed-step means are requested.
%   config - Projection-cache token: 'native', a named step, or
%       'T<seconds>seconds-<method>' with 17 significant digits for seconds.
%
% See also
%   processfieldtime
%
% Created by
%   2026/11/06, En-Chi Lee

function [tstep, mthd, useMean, config] = normalisetimeconfig(tstep, mthd)

    arguments (Input)
        tstep {mustBeTimeStep}
        mthd {mustBeTextScalar}
    end

    mthd = lower(mthd);

    if ~isempty(tstep) && isnumeric(tstep)
        tstep = days(tstep);
    elseif ischar(tstep) || isstring(tstep)
        tstep = lower(string(tstep));
    end

    % Both 'midmonth' and 'GRACE' timesteps default to using epoch averages
    useMean = ~isempty(tstep) && (isstring(tstep) || strcmpi(mthd, 'mean'));

    if isempty(tstep)
        config = 'native';
    elseif isduration(tstep)
        config = sprintf('T%sseconds-%s', sprintf('%.17g', seconds(tstep)), mthd);
    else
        config = char(tstep);
    end

end
