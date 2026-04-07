%% ICESHEETPOLY - Gets approximate
%
% Syntax
%   [icePoly, iceLonlat] = ICESHEETPOLY(icesheet)
%
% Input arguments
%   icesheet - Name of the ice sheet
%       The following ice sheets are recognised (case-insensitive):
%       - "Fennoscandian", "Fennoscandian Ice Sheet", "FIS"
%       - "Laurentide", "Laurentide Ice Sheet", "LIS"
%       - "Greenland", "Greenland Ice Sheet", "GIS"
%       - "Patagonian", "Patagonian Ice Sheet", "PIS"
%       - "Antarctic", "Antarctic Ice Sheet", "AIS"
%       - "West Antarctic", "West Antarctic Ice Sheet", "WAIS"
%       - "East Antarctic", "East Antarctic Ice Sheet", "EAIS"
%
% Output arguments
%   icePoly - Polyshape object of the ice sheet
%   iceLonlat - Longitude-latitude vertices of the ice sheet polygon
%
% Author
%   2026/04/07, En-Chi Lee (williameclee@arizona.edu)

function [icePoly, iceLonlat] = icesheetPoly(icesheet)

    arguments (Input)
        icesheet {mustBeTextScalar}
    end

    arguments (Output)
        icePoly (1, 1) polyshape
        iceLonlat (:, 2) double
    end

    switch lower(icesheet)
        case {"fennoscandian", "fennoscandian ice sheet", "fis"}
            iceLonlat = [0, 45; 30, 45; 95, 70; 95, 90; 0, 90; NaN, NaN; 360, 68; 344, 58; 344, 47; 360, 47];
        case {"laurentide", "laurentide ice sheet", "lis"}
            iceLonlat = [240, 90; 185, 70; 185, 50; 215, 50; 230, 30; 290, 30; 315, 45; 315, 55; 300, 60; 300, 70; 286, 77; 300, 90];
        case {"greenland", "greenland ice sheet", "gis"}
            iceLonlat = [300, 60; 300, 70; 286, 77; 300, 90; 355, 90; 355, 75; 340, 70; 325, 60; 315, 55];
        case {"patagonian", "patagonian ice sheet", "pis"}
            iceLonlat = [275, -35; 275, -60; 300, -60; 300, -35];
        case {"antarctic", "antarctic ice sheet", "ais"}
            iceLonlat = [0, -90; 360, -90; 360, -60; 0, -60];
        case {"west antarctic", "west antarctic ice sheet", "wais"}
            iceLonlat = [185, -60; 166, -78; 153, -83; 162, -85; 200, -87; 240, -85; 290, -85; 306, -83; 318, -83; 340, -74; 340, -60];
        case {"east antarctic", "east antarctic ice sheet", "eais"}
            iceLonlat = [0, -90; 0, -60; 185, -60; 166, -78; 153, -83; 162, -85; 200, -87; 240, -85; 290, -85; 306, -83; 318, -83; 340, -74; 340, -60; 360, -60; 360, -90];
        otherwise
            error("Unrecognised ice sheet name: %s", icesheet)
    end

    icePoly = polyshape(iceLonlat(:, 1), iceLonlat(:, 2));
end
