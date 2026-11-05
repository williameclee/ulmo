%% PROCESSSTERICTIME - Apply the steric loader's temporal convention to fields.
%
% Syntax
%   [data, dates] = PROCESSSTERICTIME(data, dates, timestep, method, limits)
%
% Input arguments
%   data - Structure of arrays with time on axis three.
%   dates - Native datetime column, (nTime, 1).
%   timestep - Native ([]), duration, midmonth or GRACE selection.
%   method - Temporal interpolation method, or mean for fixed-step averaging.
%   limits - Datetime range for early temporal selection, or [].
%
% Output arguments
%   data - Processed fields with time on axis three.
%   dates - Processed datetime column.
%
% Notes
%   Each field uses its own finite samples. Process grids and GMSL separately.
%
% See also
%   steric2lonlatt, interptemporal, temporalmeanepochs
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
