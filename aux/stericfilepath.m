%% STERICFILEPATH - Generate native and interpolated steric product paths.
%
% Syntax
%   nativePath = STERICFILEPATH(product)
%   outputPath = STERICFILEPATH(product, timeStep, meshSize, lonOrigin, ...
%       intpMthd, truncBeforeIntp, timelim, component)
%
% Input arguments
%   product - Prepared StericProduct identifier.
%   timeStep, meshSize, lonOrigin, intpMthd (optional) - Temporal and spatial
%       cache selectors, as in STERIC2LONLATT. Default selectors: empty.
%   truncBeforeIntp (optional) - Whether to include timelim in the filename.
%       Default flag: false.
%   timelim (optional) - Datetime or datenum limits. Default limits: empty.
%   component (optional) - Component suffix, including its leading hyphen.
%       Default suffix: empty.
%
% Output arguments
%   outPath - MAT filename under the IFILES/STERIC/product directory.
%       The function does not create directories or read/write data.
%
% See also
%   steric2lonlatt, steric2gmsl
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)

function outPath = stericfilepath( ...
        product, timeStep, meshSize, lonOrigin, intpMthd, truncBeforeIntp, timelim, component)

    arguments (Input)
        product StericProduct {mustBeScalarOrEmpty}
        timeStep {mustBeTimeStep} = []
        meshSize {mustBePositive} = []
        lonOrigin {mustBeFinite, mustBeReal} = []
        intpMthd {mustBeTextScalar} = ""
        truncBeforeIntp (1, 1) {mustBeNumericOrLogical} = false
        timelim {mustBeTimeRange} = []
        component (1, :) char = ''
    end

    outputFolder = fullfile(getenv("IFILES"), 'STERIC', char(product));

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

    timelimStr = '';

    if truncBeforeIntp && ~isempty(timelim)

        if isnumeric(timelim)
            timelim = datetime(timelim, 'ConvertFrom', 'datenum');
        end

        timeFmt = 'yyyyMMdd';

        if any(timeofday(timelim) ~= 0)
            timeFmt = 'yyyyMMdd''T''HHmmssSSSSSSSSS';
        end

        timelimStr = sprintf('-%s_%s', string(timelim(1), timeFmt), string(timelim(2), timeFmt));
    end

    outputFile = sprintf('%s-StericSeaLevel%s%s%s%s%s.mat', ...
        product, timeStepStr, spaceIntpStr, intpMethodStr, timelimStr, component);

    outPath = fullfile(outputFolder, outputFile);
end
