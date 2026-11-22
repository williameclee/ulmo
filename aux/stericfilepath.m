%% STERICFILEPATH - Builds the filename for a native or derived steric product.
%
% Syntax
%   outPath = STERICFILEPATH(product)
%   outPath = STERICFILEPATH(product, timeStep, meshSize, lonOrigin, ...
%       intpMthd)
%
% Input arguments
%   product - Prepared product identifier accepted by StericProduct.
%   timeStep (optional) - Temporal cache selector.
%       Positive numeric scalar or duration, 'midmonth', or 'GRACE'.
%       Uses the same temporal cache key as STERIC2SLEPT.
%       Unit: days for numeric values.
%       Default selector: [] (native timestamps).
%   meshSize (optional) - Spatial grid spacing.
%       Unit: degrees.
%       Default spacing: [] (no spacing suffix).
%   lonOrigin (optional) - Longitude centre.
%       Unit: degrees east.
%       Default centre: [] (no origin suffix).
%   intpMthd (optional) - Interpolation method.
%       Used for fixed-step temporal processing and spatial interpolation.
%       Default method: linear.
%
% Output arguments
%   outPath - MAT file path under the product folder in IFILES/STERIC.
%
% See also
%   steric2lonlatt, steric2slept, steric2gmsl, normalisetimeconfig
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/22, En-Chi Lee (williameclee@gmail.com)

function outPath = stericfilepath( ...
        product, timeStep, meshSize, lonOrigin, intpMthd)

    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty}
        timeStep {mustBeTimeStep} = []
        meshSize {mustBePositive} = []
        lonOrigin {mustBeFinite, mustBeReal} = []
        intpMthd (1, 1) string = "linear"
    end

    outFolder = fullfile(getenv("IFILES"), 'STERIC', char(product));

    [timeStep, intpMthd, ~, timeConfig] = normalisetimeconfig(timeStep, intpMthd);
    timeStepStr = '';

    if ~isempty(timeStep)
        timeStepStr = ['-', timeConfig];
    end

    spaceIntpStr = '';

    if ~isempty(meshSize) && ~isempty(lonOrigin)
        spaceIntpStr = ['-H', num2str(meshSize), 'O', num2str(lonOrigin)];
    elseif ~isempty(meshSize)
        spaceIntpStr = ['-H', num2str(meshSize)];
    elseif ~isempty(lonOrigin)
        spaceIntpStr = ['-O', num2str(lonOrigin)];
    end

    intpMethodStr = '';

    if ~isempty(spaceIntpStr) && ~isduration(timeStep)
        % A fixed-step key already includes the interpolation method.
        intpMethodStr = sprintf('-%s', intpMthd);
    end

    outFile = sprintf('%s-StericSeaLevel%s%s%s.mat', ...
        product, timeStepStr, spaceIntpStr, intpMethodStr);

    outPath = fullfile(outFolder, outFile);
end
