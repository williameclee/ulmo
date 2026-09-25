%% EARTHQUAKES
% Finds the longitude and latitude coordinates of the mask around megathrust earthquakes that may leave signals in GRACE data.
%
% Syntax
%   XY = earthquakes(upscale, buf)
%   XY = earthquakes(upscale, buf, epicenter)
%   earthquakes(__)
%
% Input arguments
%   upscale - How many times to upscale the boundary of the mask
%       The minimum value is 3, so that the mask looks circular.
%       The default value is 0 (but actually 3).
%   buf - The radius of the mask in degrees
%       The default value is 10.
%   epicenter - The epicentres of the earthquakes
%       The epicenter should be specified as a 2-by-2 matrix, where the 
%       first column is the longitude and the second column is the 
%       latitude.
%       The default value is the epicentres of 2004 Indian Ocean earthquake 
%       and the 2011 Tōhoku earthquake.
%
% Output arguments
%   XY - Closed-curved coordinates of the mask
%       The mask is defined as the buffer around the epicentres of the
%       earthquakes.
%       The coordinates are in the form of a closed-curved polygon.
%   map - Map of the mask boundary
%       When there is no output argument, a quick and dirty map will be
%       displayed.
%
% Examples
%   >>  XY = earthquakes(3, 10);
%   >>  epicentres = [95.854, 3.316; 142.269, 38.322];
%   >>  XY = earthquakes(3, 10, epicentres);
%
% Last modified by
%   2024/08/09, williameclee@arizona.edu (@williameclee)

function XY = earthquakes(varargin)
    %% Initilisation
    % Parse inputs
    p = inputParser;
    addOptional(p, 'Upscale', 0, ...
        @(x) isnumeric(x) || isempty(x));
    addOptional(p, 'Buffer', 10, ...
        @(x) isnumeric(x));
    addOptional(p, 'Epicenter', [], ...
        @(x) (isnumeric(x) && length(x) == 2) || isempty(x));
    parse(p, varargin{:});

    upscale = p.Results.Upscale;
    buf = p.Results.Buffer;
    epicentre = p.Results.Epicenter;

    if isempty(upscale) || upscale == 1
        upscale = 0;
    end

    if isempty(epicentre)
        epicentre = ...
            [95.854, 3.316; ... % 2004 Indian Ocean earthquake
             142.269, 38.322]; % 2011 Tōhoku earthquake
    end

    clear p % to avoid confusion with the mask polygon

    %% Computing the mask
    p = cell(1, size(epicentre, 1));

    for i = 1:size(epicentre, 1)
        [maskLat, maskLon] = bufferm(epicentre(i, 2), epicentre(i, 1), ...
            buf, "outPlusInterior");
        p{i} = polyshape(maskLon, maskLat);
    end

    p = union(p{:});
    XY = p.Vertices;

    XY = closecoastline(XY);

    if upscale > 3
        XY = bezier(XY, upscale);
    else
        XY = bezier(XY, 3);
    end

    %% Returning or plotting data
    if nargout > 0
        return
    end

    load coastlines %#ok<LOAD>

    figName = sprintf('Earthquakes (%s)', upper(mfilename));

    figure(999)
    set(gcf, 'Name', figName, 'NumberTitle', 'off')
    clf

    plotqdm(XY)
    hold on
    plot(coastlon, coastlat, 'k')
    plot(XY(:, 1), XY(:, 2), 'b.-')
    hold off
end
