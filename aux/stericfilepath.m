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
%       Numeric steps use the MATLAB %d format. Duration steps use their
%       displayed value and unit with spaces removed.
%       Unit: days for numeric values.
%       Default selector: [] (native timestamps).
%   meshSize (optional) - Spatial grid spacing.
%       Unit: degrees.
%       Default spacing: [] (no spacing suffix).
%   lonOrigin (optional) - Longitude centre.
%       Unit: degrees east.
%       Default centre: [] (no origin suffix).
%   intpMthd (optional) - Interpolation method.
%       Used only when a temporal or spatial selector is supplied.
%       Default method: "".
%
% Output arguments
%   outPath - MAT file path under the product folder in IFILES/STERIC.
%
% See also
%   steric2lonlatt, steric2gmsl
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/20, En-Chi Lee (williameclee@gmail.com)

function outPath = stericfilepath( ...
        product, timeStep, meshSize, lonOrigin, intpMthd)

    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty}
        timeStep {mustBeTimeStep} = []
        meshSize {mustBePositive} = []
        lonOrigin {mustBeFinite, mustBeReal} = []
        intpMthd (1, 1) string = ""
    end

    outFolder = fullfile(getenv("IFILES"), 'STERIC', char(product));

    timeStepStr = '';

    if ~isempty(timeStep)

        if ischar(timeStep) || isstring(timeStep)
            timeStepStr = sprintf('-T%s', lower(timeStep));
        elseif isduration(timeStep)
            timeStepStr = sprintf('-T%s', erase(sprintf('%s', timeStep), ' '));
        elseif isnumeric(timeStep)
            timeStepStr = sprintf('-T%d', timeStep);
        else
            error(sprintf('ULMO:%s:InvalidTimeStep', upper(mfilename)), ...
                'Unrecognised time step input for saving: %s', class(timeStep));
        end

        % Earlier linear fixed-step caches could contain bin averages.
        if isduration(timeStep) && strcmpi(intpMthd, 'linear')
            timeStepStr = [timeStepStr, '-interp'];
        end

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

    if ~(isempty(timeStepStr) && isempty(spaceIntpStr))
        intpMethodStr = sprintf('-%s', lower(intpMthd));
    end

    outFile = sprintf('%s-StericSeaLevel%s%s%s.mat', ...
        product, timeStepStr, spaceIntpStr, intpMethodStr);

    outPath = fullfile(outFolder, outFile);
end
