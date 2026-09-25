%% FORMATLATTICKS - Format latitude tick labels for a plot. 
% The tick labels will always have the same number of significant digits.
%
% Syntax
%   formatlatticks
%       Formats the y-axis tick labels of the current plot as latitude
%       ticks.
%   formatlatticks(lat)
%       Formats the y-axis tick labels of the current plot as latitude
%       ticks with the specified values.
%   formatlatticks(ax)
%       Formats the y-axis tick labels of the specified axes as latitude
%       ticks.
%   latL = formatlatticks(__)
%       Returns the formatted latitude tick labels as a cell array of
%       character vectors.
%
% Input arguments
%   lat - Latitude ticks in degrees
%   ax - Axes object to format
%
% Output arguments
%   latL - Formatted latitude tick labels
%
% Last modified by
%   2024/10/11, williameclee@arizona.edu (@williameclee)

function varargout = formatlatticks(varargin)

    if nargin > 0 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
        ax = varargin{1};
        varargin(1) = [];
    else
        ax = nan;
    end

    ip = inputParser;
    addOptional(ip, 'lat', [], @(x) isnumeric(x) || iscell(x));
    parse(ip, varargin{:});
    lat = ip.Results.lat;

    if iscell(lat)
        lat = cell2mat(lat);
    elseif isempty(lat)

        if isnumeric(ax) && isnan(ax)
            lat = yticks;
        else
            lat = yticks(ax);
        end

    end

    if isnumeric(ax) && isnan(ax)
        yticks(lat)
    else
        yticks(ax, lat)
    end

    if isscalar(lat)
        latL = formatlattick(lat, nan);
    elseif isnumeric(lat) || iscell(lat)

        if iscell(lat)
            lat = cell2mat(lat);
        end

        minDiff = min(diff(lat));

        if minDiff >= 1
            sigNum = 0;
        else
            sigNum = ceil(-log10(minDiff));
        end

        latL = arrayfun(@(x) formatlattick(x, sigNum), lat, 'UniformOutput', false);
    else
        error('Invalid input type %s', class(lat));
    end

    if nargout > 0
        varargout = {latL};
        return
    end

    clear varargout

    if isnumeric(ax) && isnan(ax)
        yticklabels(latL)
    else
        yticklabels(ax, latL)
    end

end

%% Subfunctions
function latL = formatlattick(lat, sigNum)
    sgn = sign(lat);

    if ~isnan(sigNum)
        lat = sprintf(['%0.', num2str(sigNum), 'f'], abs(lat));
    end

    latL = num2str(lat);

    if sgn > 0
        latL = [latL, '°N'];
    elseif sgn < 0
        latL = [latL, '°S'];
    else
        latL = '0°';
    end

end
