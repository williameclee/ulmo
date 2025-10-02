%% BUFFER4OCEANS
% Buffers the coastlines
%
% Last modified by
%   2025/06/03, williameclee@arizona.edu (@williameclee)

function [coastXY, coastPoly] = buffer4oceans(coast, varargin)
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    % Parse inputs
    p = inputParser;
    addRequired(p, 'Coast');
    addOptional(p, 'MoreBuffers', []);
    addOptional(p, 'LonOrigin', 0);
    parse(p, coast, varargin{:});
    coast = p.Results.Coast;
    moreBufs = p.Results.MoreBuffers;
    lonOrigin = p.Results.LonOrigin;

    %% Buffer continents with different buffers
    if isa(coast, 'polyshape')
        coastPoly = coast;
    else
        coastPoly = polyshape(coast);
    end

    if isempty(moreBufs)
        coastXY = boundary(coastPoly);
        return
    end

    for i = 1:length(moreBufs) / 2

        try
            moreCoastXY = feval(moreBufs{i * 2 - 1}, [], moreBufs{i * 2}, "RotateBack", true);
        catch
            moreCoastXY = feval(moreBufs{i * 2 - 1}, [], moreBufs{i * 2});
        end

        [moreCoastY, moreCoastX] = flatearthpoly( ...
            moreCoastXY(:, 2), moreCoastXY(:, 1), lonOrigin);
        moreCoastX = moreCoastX - 360 * floor(min(moreCoastX) / 360);
        moreCoastPoly = polyshape(moreCoastX, moreCoastY);
        coastPoly = union(coastPoly, moreCoastPoly);
    end

    [coastX, coastY] = boundary(coastPoly);
    coastXY = [coastX, coastY];

end
