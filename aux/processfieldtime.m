%% PROCESSFIELDTIME - Apply the existing temporal convention to named fields.
%
% Syntax
%   [data, dates] = PROCESSFIELDTIME(data, dates, timestep, mthd, limits)
%
% Input arguments
%   data - Nonempty scalar structure of fields with time on axis three.
%       Each field has size (nX, nY, nNativeTime); scalar series use
%       (1, 1, nNativeTime).
%   dates - Strictly increasing native datetime timestamps.
%       Size: (nNativeTime, 1).
%   timestep - Native ([]), positive duration or numeric days, 'midmonth',
%       or 'GRACE'. Named epochs always average.
%   mthd - Temporal interpolation method, or 'mean' for fixed-step means.
%   limits - Inclusive datetime limits selecting processed timestamps, or [].
%       Averaging selects window centres without clipping native samples.
%       Size: two-element vector or empty.
%
% Output arguments
%   data - Processed fields, retaining names, spatial dimensions and units.
%       Size of each field: (nX, nY, nTime).
%   dates - Processed datetime column, (nTime, 1).
%
% See also
%   normalisetimeconfig, interptemporal, temporalmeanepochs, averageepochs
%
% Created by
%   2026/11/06, En-Chi Lee

function [data, dates] = processfieldtime(data, dates, timestep, mthd, limits)

    arguments (Input)
        data (1, 1) struct
        dates datetime
        timestep {mustBeTimeStep}
        mthd {mustBeTextScalar}
        limits = []
    end

    [timestep, mthd, useMean] = normalisetimeconfig(timestep, mthd);

    if isempty(timestep) % No interpolation needed
        dates = dates(:);
        return
    end

    names = fieldnames(data);
    assert(~isempty(names), 'ULMO:processfieldtime:EmptyFields', 'Supply at least one field.');

    if useMean
        epochs = temporalmeanepochs(dates, timestep, limits);
    end

    for iName = 1:numel(names)

        if useMean
            [data.(names{iName}), outputDates] = averageepochs(dates, data.(names{iName}), epochs);
        else
            [data.(names{iName}), outputDates] = interptemporal(dates, data.(names{iName}), ...
                timestep, mthd, TimeRange = limits, BeQuiet = true);
        end

    end

    dates = outputDates(:);
end
