% Integrate the output from GEOBOXCAP which puts 1 in the region and 0
% elsewhere.  This should integrate to the area of the region.  In
% practice, since the SH expansion of the GEOBOXCAP mask has rings, this
% comparison has a good amount of error.
% This is FRACTIONAL area on the unit sphere.
function vout = plm2avgp_demo1(varargin)
    dom = 'australia';
    L = 60;
    [~, ~, ~, ~, ~, lmcosi] = geoboxcap(L, dom, act = 1);
    [Int, A, miniK, XY] = plm2avgp(lmcosi, dom);
    plotplm(lmcosi, [], [], 4);
    colorbar
    fprintf("\nIntegration check... This should equal the area of the region with error mostly from GEOBOXCAP\n");
    fprintf("PLM2AVG Int: %6.7f ; SPHAREA A: %6.7f ; ERROR: %6.2f%%\n", ...
        Int, spharea(XY), (spharea(XY) - Int) / spharea(XY) * 100)

    % Provide output where requested
    vout = {Int, A, miniK, XY};
end
