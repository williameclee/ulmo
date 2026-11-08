%% PROCESSFIELDTIME - Averages or interpolates gridded fields in time.
%
% Syntax
%   [data, dates] = PROCESSFIELDTIME(data, dates, timestep, mthd)
%   [data, dates] = PROCESSFIELDTIME(data, dates, timestep, mthd, limits)
%
% Input arguments
%   data - Scalar structure of fields with time on the third dimension.
%       Size of each field: (nX, nY, nNativeTime). Spatial sizes may differ.
%   dates - Strictly increasing datetime vector with one date per native sample.
%   timestep - Native sampling ([]), a positive scalar number of days or
%       duration, or one of the names 'midmonth' and 'GRACE'.
%   mthd - Interpolation method or 'mean' for fixed-step averaging.
%       Named steps always use epoch means. Native sampling returns all
%       fields and dates directly, regardless of the method or limits.
%   limits (optional) - Inclusive datetime limits for output dates.
%       For means, limits select window centres. Each selected window
%       uses all its available native samples.
%       Size: two-element vector or empty.
%       Default range: [] (full native date range).
%
% Output arguments
%   data - Processed fields with their original names, spatial dimensions and units.
%       Size of each field: (nX, nY, nTime).
%   dates - Processed datetime column vector.
%       Size: (nTime, 1).
%
% See also
%   normalisetimeconfig, interptemporal, temporalmeanepochs, averageepochs
%
% Created by
%   2026/11/06, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/07, En-Chi Lee (williameclee@gmail.com)

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
