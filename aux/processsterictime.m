%% PROCESSSTERICTIME - Apply the steric loader's temporal convention to fields.
%
% Syntax
%   [data, dates] = PROCESSSTERICTIME(data, dates, timestep, method, limits)
%
% Input arguments
%   data - Nonempty scalar structure of fields with time on axis three.
%       Each field has size (nX, nY, nNativeTime); scalar series use (1, 1,
%       nNativeTime). Spatial sizes can differ between fields.
%   dates - Strictly increasing native datetime timestamps.
%       Size: (nNativeTime, 1).
%   timestep - Native ([]), positive duration, 'midmonth' or 'GRACE' selection.
%       Callers convert numeric day counts to durations before calling.
%   method - Temporal interpolation method, or mean for fixed-step averaging.
%       Named epochs always average; native selection ignores this argument.
%   limits - Inclusive datetime limits selecting processed timestamps, or [].
%       Averaging selects window centres without clipping native samples.
%       Native selection ignores these limits; callers select outputs later.
%       Size: two-element vector, or empty for the complete record.
%
% Output arguments
%   data - Processed fields with time on axis three.
%       Field names and spatial dimensions are unchanged.
%       Each field has size (nX, nY, nTime); input units are preserved.
%   dates - Processed datetime column.
%       Size: (nTime, 1).
%
% Notes
%   Each field uses its own finite samples. Process grids and GMSL separately.
%
% See also
%   steric2lonlatt, steric2slept, interptemporal, temporalmeanepochs
%
% Created by
%   2026/11/05, En-Chi Lee

function [data, dates] = processsterictime(data, dates, timestep, method, limits)
    % Grids and means enter separately but share temporal conventions. Epoch
    % averaging uses each array's own finite observations; it deliberately does
    % not match a GMSL to each individual cell's temporal support.
    if isempty(timestep), dates = dates(:); return; end
    names = fieldnames(data);
    useMean = ischar(timestep) || isstring(timestep) || strcmpi(method, 'mean');
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
