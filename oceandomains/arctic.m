%% ARCTIC
% Finds the longitude and latitude coordinates of the Arctic Ocean.
%
% Syntax
%   XY = arctic(upscale, buf)
%   XY = arctic(upscale, buf, latlim, morebuffers)
%   XY = arctic(__, 'Name', value)
%   [XY, p] = arctic(__)
%   [XY, lonc, latc] = arctic(__)
%	flag = arctic('rotated')
%
% Input arguments
%   upscale - How many times to upscale the data
%       This option should be used carefully, as spline interpolation do
%       not work well at sharp corners, and the original sampling rate
%       should be sufficient for most applications.
%       The default value is 0 (no upscaling).
%   buffer - The (negative) buffer from the coastlines in degrees
%       The default value is 0 (no buffer).
%   latlim - The latitudes of the polar caps in degrees
%       The inclination angle can be specified as a scalar or a 1-by-2
%       vector.
%       - If a scalar, the inclination angle is applied to both the north
%           and south polar caps.
%       - If a 1-by-2 vector, the first element is the inclination angle of
%           the north polar cap, and the second element is the inclination
%           angle of the south polar cap.
%       The default value is 90 (no polar caps).
%   morebuffers - Additional buffers to apply to the coastlines
%       The additional buffers must be specified as a cell array of
%       domain names and buffer values.
%       The default value is an empty cell array.
%   LonOrigin - The longitude origin of the data
%       The domain will be contained within the range
%       [LonOrigin - 180, LonOrigin + 180].
%       The default value is 180 (i.e. the range of longitude is [0, 360]).
%   RotateBack - A flag to indicate whether the coordinates should be
%       rotated back to their original position
%       The default value is true (coordinates rotated to equator).
%   ForceNew - Force the function to reload the data
%       The default value is false.
%   BeQuiet - Suppress the output messages
%       The default value is false.
%
% Output arguments
%   XY - Closed-curved coordinates of the domain
%		By default, the coordinates are rotated 90 degrees about the y-axis to centre at the equator.
%   p - Polygon of the domain boundary
%       Note that this is an ordinary polyshape object, and not a geoshape
%       or geopolyshape.
%   lonc, latc - Longitude and latitude by which the coordinates are rotated
%       The values are 0 and 90, respectively (i.e. rotated to the equator).
%   map - Map of the domain boundary
%       When there is no output argument, a quick and dirty map will be
%       displayed.
%	flag - A flag to indicate that the coordinates are rotated
%       The only possible value is 1.
%
% Examples
%   The 'default' domain used for most applications is given by
%   XY = arctic('Buffer', 1);
%
% Notes
%   The function is written intentionally so that it is compatible with
%   other region functions from slepian_alpha, i.e.
%   XY = arctic(upscale, buf)
%   works as intended.
%
% Data source
%   The coastline data is based on the GSHHG coastline data:
%       Wessel, P., & Smith, W. H. F. (1996).
%       doi: 10.1029/96JB00104
%   The ocean boundaries are based on IHO's 'Limits of oceans and seas':
%       International Hydrographic Organization & Sieger, R. (2012).
%       doi: 10.1594/PANGAEA.777975
%
% See also
%   OCEANPOLY, GSHHSCOASTLINE, BUFFER4OCEANS
%
% Last modified
%   2026/02/12, williameclee@arizona.edu (@williameclee)
%     - Added varible check before loading
%   2024/08/15, williameclee@arizona.edu (@williameclee)

function varargout = arctic(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');

    % First check: rotation
    if nargin == 1 && strcmpi(varargin{1}, 'rotated')
        varargout = {true};
        return
    end

    % Parse the inputs
    lonOriginD = 180;
    [upscale, latlim, buf, moreBufs, lonOrigin, rotateBack, ...
         forceNew, saveData, beQuiet] = ...
        parseoceaninputs(varargin, "DefaultLonOrigin", lonOriginD);
    oceanParts = ...
        {'Arctic Ocean, eastern part', ...
     'Arctic Ocean, western part'};

    %% Check if the data file exists
    [dataFile, ~, dataExists] = oceanfilename(mfilename, ...
        'Upscale', upscale, 'Latlim', latlim, ...
        'Buffer', buf, 'MoreBuffers', moreBufs, 'RotateBack', rotateBack);

    if dataExists && ~forceNew && all(ismember({'XY', 'p'}, who('-file', dataFile)))
        load(dataFile, 'XY', 'p')

        % Make sure the requested data exists
        if beQuiet < 2
            fprintf('%s loaded %s\n', upper(mfilename), dataFile)
        end

        if lonOrigin ~= lonOriginD
            [Y, X] = flatearthpoly(XY(:, 2), XY(:, 1), lonOrigin);
            p = polyshape(X, Y);
            XY = poly2xy(p);
        end

        varargout = returncoastoutputs(nargout, XY, p);

        if nargout == 3
            varargout = {XY, 0, 90};
        end

        return
    end

    %% Compute the ocean boundary
    % Find the ocean boundary (not accounting for the coastlines)
    [oceanPoly, oceanLatlim, oceanLonlim] = oceanpoly( ...
        oceanParts, latlim, lonOrigin, 'BeQuiet', beQuiet);
    % Find the coastline
    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LonOrigin', lonOrigin, 'BeQuiet', beQuiet);
    % Crop out more buffers
    [~, coastPoly] = buffer4oceans(coastPoly, ...
        'MoreBuffers', moreBufs, 'LonOrigin', lonOrigin);
    % Manually remove small holes
    coastPoly = manualadjustments(coastPoly, buf, lonOrigin);
    % Subtract the land regions from the ocean boundary
    p = subtract(oceanPoly, coastPoly);

    % Turn the polygon into a well-defined curve
    XY = poly2xy(p, upscale);

    % Rotate the coordinates if not asked to rotate back
    if ~rotateBack
        XY = rotateXY(XY);
    end

    %% Save and return data
    varargout = returncoastoutputs(nargout, XY, p);

    if nargout == 3
        varargout = {XY, 0, 90};
    end

    if ~saveData
        return
    end

    save(dataFile, 'XY', 'p', '-v7.3')

    if beQuiet < 2
        fprintf('%s saved %s\n', upper(mfilename), dataFile)
    end

end

%% Subfunctions
% Manually remove small holes in the coastline
function coastPoly = manualadjustments(coastPoly, buf, lonOrigin)
    coastPoly = addlandregion(coastPoly, ...
        'Latlim', [65, 66; 62, 63], ...
        'Lonlim', [23, 26; 316, 318], ...
        'LongitudeOrigin', lonOrigin);

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [65, 67; 63, 64; 59, 62; 66, 68; 74, 75; 77, 78], ...
            'Lonlim', [278, 282; 284 286; 292, 294; 185, 195; 255, 256; 255, 256], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 3
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [59, 63; 58, 61], ...
            'Lonlim', [300, 310; 270, 280], ...
            'LongitudeOrigin', lonOrigin);
    end

end

function lonlatR = rotateXY(lonlat)
    col = deg2rad(90 - lonlat(:, 2));
    lon = deg2rad(lonlat(:, 1));

    [colR, lonR] = rottp(col, lon, deg2rad(0), deg2rad(-90), deg2rad(0));

    lonlatR = [rad2deg(lonR), 90 - rad2deg(colR)];
    lonlatP = simplify(polyshape(lonlatR));
    lonlatR = poly2xy(lonlatP);
end
