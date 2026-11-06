%% PROCESSFIELDTIME - Apply the existing temporal convention to named fields.
%
% Syntax
%   [data, dates] = PROCESSFIELDTIME(data, dates, timestep, method, limits)
%
% Input arguments
%   data - Nonempty scalar structure of fields with time on axis three.
%       Each field has size (nX, nY, nNativeTime); scalar series use
%       (1, 1, nNativeTime). Spatial sizes may differ between fields.
%   dates - Strictly increasing native datetime timestamps.
%       Size: (nNativeTime, 1).
%   timestep - Native ([]), positive duration or numeric days, 'midmonth',
%       or 'GRACE'. Named epochs always average.
%   method - Temporal interpolation method, or 'mean' for fixed-step means.
%       Native selection ignores this argument.
%   limits - Inclusive datetime limits selecting processed timestamps, or [].
%       Averaging selects window centres without clipping native samples.
%       Native selection ignores limits; callers select returned dates later.
%       Size: two-element vector or empty.
%
% Output arguments
%   data - Processed fields, retaining names, spatial dimensions and units.
%       Size of each field: (nX, nY, nTime).
%   dates - Processed datetime column, (nTime, 1).
%
% Notes
%   Each field uses its own finite samples. This function applies no residual,
%   depth, uncertainty, source completeness or global-mean policy.
%
% See also
%   normalizetimeconfig, interptemporal, temporalmeanepochs, averageepochs
%
% Created by
%   2026/11/06, En-Chi Lee

function [data, dates] = processfieldtime(data, dates, timestep, method, limits)
    arguments (Input)
        data (1, 1) struct
        dates datetime
        timestep {mustBeTimeStep}
        method {mustBeTextScalar}
        limits = []
    end
    [timestep, method, useMean] = normalizetimeconfig(timestep, method);
    if isempty(timestep), dates = dates(:); return; end
    names = fieldnames(data);
    assert(~isempty(names), 'ULMO:processfieldtime:EmptyFields', 'Supply at least one field.');
    if useMean, epochs = temporalmeanepochs(dates, timestep, limits); end
    for iName = 1:numel(names)
        if useMean
            [data.(names{iName}), outputDates] = averageepochs(dates, data.(names{iName}), epochs);
        else
            [data.(names{iName}), outputDates] = interptemporal(dates, data.(names{iName}), ...
                timestep, method, TimeRange=limits, BeQuiet=true);
        end
    end
    dates = outputDates(:);
end
