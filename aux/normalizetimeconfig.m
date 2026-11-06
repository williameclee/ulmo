%% NORMALIZETIMECONFIG - Normalize temporal processing and projection identity.
%
% Syntax
%   [timestep, method, useMean, configuration] = NORMALIZETIMECONFIG(timestep, method)
%
% Input arguments
%   timestep - Native ([]), positive numeric days or duration, 'midmonth',
%       or 'GRACE'. Numeric days become durations; named steps become lower
%       case scalar strings. Duration display formats are preserved.
%   method - Temporal interpolation method or 'mean'; converted to lower
%       case. Methods are validated by the processing function when used.
%
% Output arguments
%   timestep, method - Normalized arguments in the conventions above.
%   useMean - Whether named epochs or fixed-step means are requested.
%   configuration - Projection-cache token: 'native', a named step, or
%       'T<seconds>seconds-<method>' with 17 significant digits for seconds.
%
% Notes
%   This does not construct grid-cache filenames or identify source products.
%   Legacy grid path builders retain their own suffixes. Monthly-source
%   equivalence must be resolved by the caller using its source convention.
%
% See also
%   processfieldtime, temporalmeanepochs, interptemporal
%
% Created by
%   2026/11/06, En-Chi Lee

function [timestep, method, useMean, configuration] = normalizetimeconfig(timestep, method)
    arguments (Input)
        timestep {mustBeTimeStep}
        method {mustBeTextScalar}
    end
    method = lower(method);
    if ~isempty(timestep) && isnumeric(timestep), timestep = days(timestep); end
    if ischar(timestep) || isstring(timestep)
        assert(isscalar(string(timestep)), 'ULMO:normalizetimeconfig:TimeStep', 'Supply one named time step.');
        timestep = lower(string(timestep));
    end
    useMean = ~isempty(timestep) && (isstring(timestep) || strcmpi(method, 'mean'));
    if isempty(timestep)
        configuration = 'native';
    elseif isduration(timestep)
        configuration = sprintf('T%sseconds-%s', sprintf('%.17g', seconds(timestep)), method);
    else
        configuration = char(timestep);
    end
end
