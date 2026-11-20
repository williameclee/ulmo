%% STERICPRODUCT - Identifies prepared steric products and their depth coverage.
% Full-depth products represent the water column from the surface to the
% seafloor. The prepared SIO, NOAA, JAMSTEC, ISAS and IPRC products cover only
% the upper ocean. In particular, the prepared ISAS integral stops at 2,000 m
% even though its source profiles include deeper levels.
% Product coverage does not imply observations at every depth in every cell.
% Component and layer availability must also be checked in the supplied file.
%
% Last modified
%   2026/11/20, En-Chi Lee (williameclee@gmail.com)

classdef StericProduct < uint8

    enumeration
        SIO (0)
        EN4c14 (1)
        CMEMS (2)
        NOAA (3)
        ORAS5 (4)
        JAMSTEC (5)
        ISAS (6)
        IPRC (7)
        IAPv4 (8)
    end

    methods

        %% NATIVEDEPTH - Returns the depth range of the prepared native integral.
        function depth = nativeDepth(product)

            if product.supportsFullDepth()
                depth = "full-depth";
            else
                depth = "shallow";
            end

        end

        %% SUPPORTSFULLDEPTH - Tests the depth coverage of a prepared product.
        function supported = supportsFullDepth(product)
            supported = ismember(product, ...
                [StericProduct.EN4c14, StericProduct.CMEMS, StericProduct.ORAS5, StericProduct.IAPv4]);
        end

    end

end
