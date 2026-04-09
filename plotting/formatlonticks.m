%% FORMATLONTICKS - Format longitude tick labels for a plot
%
% Syntax
%   formatlonticks
%       Formats the x-axis tick labels of the current plot as longitude
%       ticks.
%   formatlonticks(lon)
%       Formats the x-axis tick labels of the current plot as longitude
%       ticks with the specified values.
%   formatlonticks(ax)
%       Formats the x-axis tick labels of the specified axes as longitude
%       ticks.
%   lonL = formatlonticks(__)
%       Returns the formatted longitude tick labels as a cell array of
%       character vectors.
%
% Input arguments
%   lon - Longitude ticks in degrees
%   ax - Axes object to format
%
% Output arguments
%   lonL - Formatted longitude tick labels
%
% Last modified by
%   2024/10/25, williameclee@arizona.edu (@williameclee)

function varargout = formatlonticks(varargin)

    if nargin > 0 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
        ax = varargin{1};
        varargin(1) = [];
    else
        ax = nan;
    end

    ip = inputParser;
    addOptional(ip, 'lon', [], @(x) isnumeric(x) || iscell(x));
    parse(ip, varargin{:});
    lon = ip.Results.lon;

    if iscell(lon)
        lon = cell2mat(lon);
    elseif isempty(lon)

        if isnumeric(ax) && isnan(ax)
            lon = xticks;
        else
            lon = xticks(ax);
        end

    end

    if isnumeric(ax) && isnan(ax)
        xticks(lon)
    else
        xticks(ax, lon)
    end

    if isscalar(lon)
        lonL = formatlontick(lon, nan);
    elseif isnumeric(lon) || iscell(lon)

        if iscell(lon)
            lon = cell2mat(lon);
        end

        minDiff = min(diff(lon));

        if minDiff >= 1
            sigNum = 0;
        else
            sigNum = ceil(-log10(minDiff));
        end

        lonL = arrayfun(@(x) formatlontick(x, sigNum), lon, 'UniformOutput', false);
    else
        error('Invalid input type %s', class(lon));
    end

    if nargout > 0
        varargout = {lonL};
        return
    end

    clear varargout

    if isnumeric(ax) && isnan(ax)
        xticklabels(lonL)
    else
        xticklabels(ax, lonL)
    end

end

%% Subfunctions
function lonL = formatlontick(lon, sigNum)

    lon = mod(lon + 180, 360) - 180;

    if ~isnan(sigNum)
        lonL = sprintf(['%0.', num2str(sigNum), 'f'], abs(lon));
    end

    lonL = num2str(lonL);

    if lon > 0
        lonL = [lonL, '°E'];
    elseif lon < 0
        lonL = [lonL, '°W'];
    else
        lonL = '0°';
    end

end
