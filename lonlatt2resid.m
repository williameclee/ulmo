%% LONLATT2RESID
% Calculates the trend of a time series in grid format.
%
% Authored by
%	2025/07/28, williameclee@arizona.edu (@williameclee)

function trend = lonlatt2resid(t, data)
    % Initialisation
    if isdatetime(t)
        t = years(t - t(1));
    elseif isduration(t)
        t = years(t);
    end

    t = t(:);

    nTime = length(t);
    [nLon, nLat, ~] = size(data);

    flatData = reshape(data, [nLon * nLat, nTime])';

    % Solve linear regression
    X = [ones(size(t)), t];
    B = X \ flatData;
    flatTrend = B(2, :);

    % Deal with NaNs
    hasMissing = find(any(isnan(flatData), 1) & ~all(isnan(flatData), 1));

    for iGrid = hasMissing
        isValidData = ~isnan(flatData(:, iGrid));
        squeezedX = X(isValidData, :);
        squeezedData = flatData(isValidData, iGrid);

        B = squeezedX \ squeezedData;
        flatTrend(iGrid) = B(2);
    end

    trend = reshape(flatTrend, [nLon, nLat]);
end
