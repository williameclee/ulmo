%% AVERAGEEPOCHS
% Averages 3D gridded temporal samples within specified solution or epoch
% windows [start, end).
%
% Overlapping solution windows deliberately share samples without temporal
% interpolation. If no native samples fall within a window, the output grid
% for that epoch will contain NaNs.
%
% Syntax
%   avgs = AVERAGEEPOCHS(nativeDates, data, epochs)
%   [avgs, dates] = AVERAGEEPOCHS(nativeDates, data, epochs)
%
% Input arguments
%   nativeDates - Timestamps corresponding to the time dimension of data
%       Data type: DATETIME | DOUBLE (datenum)
%       Dimension: [M x 1] | [1 x M]
%   data - 3D numeric array of gridded spatial data over time
%       The first two dimensions are spatial dimensions and the third
%       dimension corresponds to nativeDates.
%       Data type: NUMERIC
%       Dimension: [nX x nY x M]
%   epochs - Start and end timestamps defining each averaging window [start, end)
%       Data type: DATETIME | DOUBLE (datenum)
%       Dimension: [N x 2]
%
% Output arguments
%   avgs - 3D array of time-averaged grids for each epoch window
%       Integer data produces double averages. Floating-point data keeps its precision.
%       Dimension: [nX x nY x N]
%   dates - Midpoint timestamps of each epoch window
%       Data type: DATETIME
%       Dimension: [N x 1]
%
% See also
%   GRACEEPOCHS, FILLGRACEEPOCHS, SSH2LONLATT, STERIC2LONLATT,
%   INTERPTEMPORAL
%
% Created by
%   2026/09/25, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/01, En-Chi Lee (williameclee@gmail.com)

function [avgs, dates] = averageepochs(nativeDates, data, epochs)

    arguments (Input)
        nativeDates {mustBeA(nativeDates, {'datetime', 'numeric'}), mustBeVector}
        data {mustBeNumeric}
        epochs (:, 2) {mustBeA(epochs, {'datetime', 'numeric'})}
    end

    arguments (Output)
        avgs {mustBeNumeric}
        dates (:, 1) datetime
    end

    if isnumeric(nativeDates)
        nativeDates = datetime(nativeDates, 'ConvertFrom', 'datenum');
    end

    if isnumeric(epochs)
        epochs = datetime(epochs, 'ConvertFrom', 'datenum');
    end

    if any(isnat(epochs), 'all') || any(epochs(:, 2) <= epochs(:, 1))
        error('ULMO:Epochs:InvalidWindows', ...
        'Epochs must be N-by-2 finite start/end times with end > start.');
    end

    if size(data, 3) ~= numel(nativeDates)
        error('ULMO:Epochs:DimensionMismatch', ...
        'Mesh time dimension must match native dates.');
    end

    dates = epochs(:, 1) + (epochs(:, 2) - epochs(:, 1)) / 2;
    outputSize = [size(data, 1), size(data, 2), size(epochs, 1)];
    if isinteger(data)
        avgs = nan(outputSize);
    else
        avgs = nan(outputSize, 'like', data);
    end

    for k = 1:size(epochs, 1)
        selected = nativeDates(:) >= epochs(k, 1) & nativeDates(:) < epochs(k, 2);

        if any(selected)
            avgs(:, :, k) = mean(data(:, :, selected), 3, 'omitmissing');
        end

    end

end
