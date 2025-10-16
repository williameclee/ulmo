function [p, lon, lat] = domain2geopoly(varargin)

    parseStart = 2;

    switch class(varargin{1})
        case 'GeoDomain'
            domain = varargin{1};
            lonlat = domain.Lonlat;
        case 'cell'
            domain = varargin{1};
            lonlat = feval(domain{:});
        case {'char', 'string'}

            if ismember(varargin{1}, {'coast', 'coasts', 'coastline', 'coastlines'})
                lonlat = gshhscoastline('c');
            else
                domain = varargin{1};
                lonlat = feval(domain);
            end

        case 'double'

            if size(varargin{1}, 2) == 2
                lonlat = varargin{1};
            elseif size(varargin{1}, 2) == 1

                if ismatrix(varargin{2}) && size(varargin{2}, 2) == 1
                    lonlat = [varargin{1}, varargin{2}];
                end

                parseStart = 3;

            end

    end

    if nargin >= parseStart
        p = inputParser;
        addOptional(p, 'Reducem', 0.5, ...
            @(x) isnumeric(x) && isscalar(x) && x > 0);
        addParameter(p, 'LonOrigin', 180, @isnumeric);
        addParameter(p, 'Inverse', false, @istruefalse);
        parse(p, varargin{parseStart:end});
        reducemFactor = p.Results.Reducem;
        lonOrigin = p.Results.LonOrigin;
        isInversed = p.Results.Inverse;
    else
        reducemFactor = 0.5;
        lonOrigin = 180;
        isInversed = true;
    end

    lat = lonlat(:, 2);
    lon = lonlat(:, 1);
    [lat, lon] = reducem(lat, lon, reducemFactor);

    if isInversed
        [lat, lon] = flatearthpoly(lat, lon, lonOrigin);
        worldBoundaryPoly = polyshape([lonOrigin - 180, -90; lonOrigin - 180, 90; lonOrigin + 180, 90; lonOrigin + 180, -90]);
        p = subtract(worldBoundaryPoly, polyshape([lon, lat]));
        p = simplify(p);
        lat = p.Vertices(:, 2);
        lon = p.Vertices(:, 1);
    end

    [lat, lon] = flatearthpoly(lat, lon, lonOrigin);
    [lat, lon] = addanchors(lat, lon);
    p = geoshape(lat, lon, 'Geometry', 'polygon');
end
