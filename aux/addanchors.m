%% ADDANCHORS
% Add anchors to horizontal or vertical segments of a polygon to ensure that the segments are not too long.
% This is useful when trying to plot the polygon on a map, as the segments are not displayed as great circle arcs.
%
% Syntax
%   p = addanchors(p)
%   lonlat = addanchors(lonlat)
%   [lon, lat] = addanchors(lon, lat)
%
% Input arguments
%   p - A polyshape object
%   lonlat - A 2-column matrix of longitudes and latitudes
%   lon, lat - Vectors of longitudes and latitudes
%
% Output arguments
%   p - A polyshape object
%   lonlat - A 2-column matrix of longitudes and latitudes
%   lon, lat - Vectors of longitudes and latitudes
%
% Notes
%   This function was originally part of the slepian_ulmo package.
%
% Last modified by
%   2024/10/15, En-Chi Lee (williameclee@arizona.edu)

function varargout = addanchors(varargin)
    isPoly = false;

    if nargin == 1

        if isa(varargin{1}, 'polyshape')
            isPoly = true;
            varargin{1} = varargin{1}.Vertices;
            varargin{1} = closecoastline(varargin{1});
        elseif size(varargin{1}, 1) == 2
            varargin{1} = varargin{1}';
        elseif size(varargin{1}, 2) ~= 2
            error('Invalid input argument')
        end

        lon = varargin{1}(:, 1);
        lat = varargin{1}(:, 2);
    elseif nargin == 2
        lon = varargin{1};
        lat = varargin{2};
    else
        error('Invalid number of input arguments')
    end

    lon = lon(:);
    lat = lat(:);

    maxSep = 5;
    tol = 1;

    diffLat = abs(diff(lat));
    diffLon = abs(diff(lon));

    needsAnchor = find( ...
        (diffLat > maxSep & diffLon <= tol) | ...
        (diffLon > maxSep & diffLat <= tol));
    difLength = max(diffLat(needsAnchor), diffLon(needsAnchor));

    lata = lat;
    lona = lon;

    if isempty(needsAnchor)

        if nargout == 1

            if isPoly
                varargout = {p};
            else
                varargout = {[lon, lat]};
            end

        else
            varargout = {lon, lat};
        end

        return
    end

    for i = length(needsAnchor):-1:1
        latInterp = interp1([0, 1], ...
            [lata(needsAnchor(i)), lata(needsAnchor(i) + 1)], ...
            linspace(0, 1, round(difLength(i) / tol)));
        lata = [lata(1:needsAnchor(i) - 1); latInterp'; lata(needsAnchor(i) + 1:end)];
        lonInterp = interp1([0, 1], ...
            [lona(needsAnchor(i)), lona(needsAnchor(i) + 1)], ...
            linspace(0, 1, round(difLength(i) / tol)));
        lona = [lona(1:needsAnchor(i) - 1); lonInterp'; lona(needsAnchor(i) + 1:end)];
    end

    if nargout == 0
        figName = sprintf('Polygon with more anchors (%s)', upper(mfilename));
        figure(998)
        set(gcf, 'Name', figName, 'NumberTitle', 'off')
        clf
        plotqdm(lona, lata)
        hold on
        plot(lona, lata, 'k.-')
        hold off
    elseif nargout == 1

        if isPoly
            varargout = {polyshape(lona, lata, "Simplify", false)};
        else
            varargout = {[lona, lata]};
        end

    else
        varargout = {lona, lata};
    end

end
