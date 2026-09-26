%% FILLGRACEEPOCHS
% Fills runs of missing solution months with aligned synthetic calendar
% windows.
%
% Real solution windows remain unchanged. Each missing-month run begins at
% the latest preceding solution end and ends at the earliest following start.
% Calendar month boundaries inside those bounds partition the run.
%
% Syntax
%   epochs = FILLGRACEEPOCHS(epochs)
%   epochs = FILLGRACEEPOCHS(epochs, timelim)
%   [epochs, synthetic] = FILLGRACEEPOCHS(__)
%
% Input arguments
%   epochs - Start and end timestamps of actual solution windows
%       Data type: DATETIME | DOUBLE (datenum)
%       Dimension: [M x 2]
%   timelim - Time range of interest
%       A two-element vector [startTime, endTime] specifying the time range
%       to restrict output epochs. If empty, defaults to the range spanned
%       by the midpoints of the input epochs: [min(labels), max(labels)].
%       The default is [] (full span of input epochs).
%       Data type: DATETIME | DOUBLE (datenum) | []
%       Dimension: [1 x 2] | [2 x 1] | []
%
% Output arguments
%   epochs - Chronologically sorted array of actual and synthetic windows
%       Start time is in the first column and end time is in the second
%       column.
%       Data type: DATETIME
%       Dimension: [N x 2]
%   synthetic - Logical flag indicating synthetic (gap-filled) windows
%       true for synthetic calendar windows inserted during missing months;
%       false for actual input solution windows.
%       Data type: LOGICAL
%       Dimension: [N x 1]
%
% See also
%   GRACEEPOCHS, AVERAGEEPOCHS, SSH2LONLATT, STERIC2LONLATT
%
% Created by
%   2026/09/25, En-Chi Lee (williameclee@arizona.edu)

function [epochs, synthetic] = fillgraceepochs(epochs, timelim)

    arguments (Input)
        epochs (:, 2) {mustBeA(epochs, {'datetime', 'numeric'})}
        timelim {mustBeTimeRange} = []
    end

    arguments (Output)
        epochs (:, 2) datetime
        synthetic (:, 1) logical
    end

    if isnumeric(epochs)
        epochs = datetime(epochs, 'ConvertFrom', 'datenum');
    end

    labels = mean(epochs, 2);

    if isempty(timelim)
        timelim = [min(labels), max(labels)];
    elseif isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    if ~isdatetime(timelim) || numel(timelim) ~= 2 || any(isnat(timelim))
        error('ULMO:Epochs:InvalidRange', 'timelim must contain two valid dates.');
    end

    % Build whole runs before filtering, even when the request is inside a gap.
    first = dateshift(min([labels; timelim(:)]), 'start', 'month');
    last = dateshift(max([labels; timelim(:)]), 'start', 'month');
    months = (first:calmonths(1):last)';
    realMonths = dateshift(labels, 'start', 'month');
    missing = ~ismember(months, realMonths);
    runStarts = find(diff([false; missing]) == 1);
    runEnds = find(diff([missing; false]) == -1);
    extra = NaT(0, 2);

    for k = 1:numel(runStarts)
        a = months(runStarts(k));
        b = months(runEnds(k)) + calmonths(1);
        before = realMonths < a;
        after = realMonths >= b;
        if any(before), a = max(epochs(before, 2)); end
        if any(after), b = min(epochs(after, 1)); end
        if b <= a, continue; end
        boundaries = months(runStarts(k) + 1:runEnds(k));
        boundaries = boundaries(boundaries > a & boundaries < b);
        edges = [a; boundaries; b];
        extra = [extra; edges(1:end - 1), edges(2:end)]; %#ok<AGROW>
    end

    synthetic = [false(size(epochs, 1), 1); true(size(extra, 1), 1)];
    epochs = [epochs; extra];
    labels = mean(epochs, 2);
    keep = labels >= min(timelim) & labels <= max(timelim);
    epochs = epochs(keep, :);
    synthetic = synthetic(keep);
    [~, order] = sort(labels(keep));
    epochs = epochs(order, :);
    synthetic = synthetic(order);
end
