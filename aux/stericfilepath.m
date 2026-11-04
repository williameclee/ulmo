%% STERICFILEPATH - Generate native and interpolated steric product paths.
%
% Syntax
%   outPath = STERICFILEPATH(product)
%   outPath = STERICFILEPATH(product, timeStep, meshSize, lonOrigin, ...
%       intpMthd, truncBeforeIntp, timelim, component)
%
% Input arguments
%   product - Prepared product: SIO, EN4c14, CMEMS, NOAA, or ORAS5.
%   timeStep (optional) - Temporal cache selector.
%       Positive numeric scalar or duration, 'midmonth', or 'GRACE'.
%       Numeric values are encoded as integers; durations retain their units.
%       Unit: days for numeric values.
%       Default selector: [] (native timestamps).
%   meshSize (optional) - Spatial grid spacing encoded in the filename.
%       Unit: degrees.
%       Default spacing: [] (no spacing suffix).
%   lonOrigin (optional) - Longitude centre encoded in the filename.
%       Unit: degrees east.
%       Default centre: [] (no origin suffix).
%   intpMthd (optional) - Interpolation method encoded in the filename.
%       Used only when a temporal or spatial selector is supplied.
%       Default method: "" (empty text).
%   truncBeforeIntp (optional) - Whether to include timelim in the filename.
%       Default flag: false.
%   timelim (optional) - Datetime or datenum limits.
%       Size: (1, 2) or (2, 1), or empty.
%       Unit: days for datenum values.
%       Default limits: [] (no time-range suffix).
%   component (optional) - Component suffix, including its leading hyphen.
%       Default suffix: '' (empty text).
%
% Output arguments
%   outPath - MAT filename under the IFILES/STERIC/product directory.
%
% Notes
%   This function only constructs a path; it does not create directories or
%   read or write data. It preserves STERIC2LONLATT cache naming, including
%   the extra '-interp' suffix for duration steps with linear interpolation.
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
