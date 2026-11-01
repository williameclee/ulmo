%% TEMPORALMEANEPOCHS - Build averaging windows for gridded sea level data.
% GRACE uses the reference solution windows. Midmonth uses calendar
% months. Fixed steps are centered on dates starting at the first sample.
% The time range selects window midpoints without clipping the samples.
%
% Created by
%	2026/11/01, En-Chi Lee (williameclee@arizona.edu)

function epochs = temporalmeanepochs(nativeDates, timestep, timelim)

    if isnumeric(nativeDates)
        nativeDates = datetime(nativeDates, 'ConvertFrom', 'datenum');
    end

    nativeRange = [min(nativeDates), max(nativeDates)];

    if isempty(timelim)
        timelim = nativeRange;
    end

    if strcmpi(timestep, 'GRACE')
        timelim = [max(nativeRange(1), timelim(1)), min(nativeRange(2), timelim(2))];

        if timelim(1) > timelim(2)
            epochs = NaT(0, 2);
            return
        end

        % Use the same gap boundaries as the full product before selecting midpoints
        if nativeRange(1) == nativeRange(2)
            monthStart = dateshift(nativeRange(1), 'start', 'month');
            nativeRange = [monthStart, monthStart + calmonths(1)];
        end

        epochs = graceepochs({'CSR', 'RL06', 60}, nativeRange);
        epochCentres = mean(epochs, 2);
    elseif strcmpi(timestep, 'midmonth')
        epochCentres = midmonth(nativeRange);
        epochStarts = dateshift(epochCentres(:), 'start', 'month');
        epochs = [epochStarts, epochStarts + calmonths(1)];
    else
        epochCentres = (nativeRange(1):timestep:nativeRange(2))';
        epochs = [epochCentres - timestep / 2, epochCentres + timestep / 2];
    end

    keep = epochCentres >= timelim(1) & epochCentres <= timelim(2);
    epochs = epochs(keep, :);
end
