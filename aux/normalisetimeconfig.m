%% NORMALISETIMECONFIG - Normalises a time step and generates its projection cache token.
%
% Syntax
%   [tstep, mthd, useMean, config] = NORMALISETIMECONFIG(tstep, mthd)
%
% Input arguments
%   tstep - Native sampling ([]), a positive scalar number of days or
%       duration, or one of the names 'midmonth' and 'GRACE'.
%   mthd - Interpolation method or 'mean' for fixed-step averaging.
%
% Output arguments
%   tstep - Normalised time step.
%       Numeric days become a duration and named steps become lower-case scalar strings.
%   mthd - Input method in lower case.
%   useMean - Whether to average over named epochs or fixed-step windows.
%   config - Projection cache name.
%       Native sampling gives 'native'. Named epochs use their lower-case
%       name. Fixed steps use 'T<seconds>seconds-<method>', with seconds
%       written to 17 significant digits.
%
% See also
%   processfieldtime, mustBeTimeStep
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

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

    % Both named time steps use epoch averages.
    useMean = ~isempty(tstep) && (isstring(tstep) || strcmpi(mthd, 'mean'));

    if isempty(tstep)
        config = 'native';
    elseif isduration(tstep)
        config = sprintf('T%sseconds-%s', sprintf('%.17g', seconds(tstep)), mthd);
    else
        config = char(tstep);
    end

end
