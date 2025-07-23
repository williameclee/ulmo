function varargout = steric2lonlatt(varargin)
    [product, timelim, timeStep, meshSize, lonOrigin, intpMthd, timeFmt, unit, ...
         forceNew, beQuiet, saveData] = ...
        parseinputs(varargin);

    outputPath = outputpath(product, timeStep, meshSize, lonOrigin, intpMthd);

    if exist(outputPath, 'file') && ~forceNew
        load(outputPath, 'lon', 'lat', 'dates', 'stericSl');

        if beQuiet <= 1
            fprintf('%s loaded %s\n', upper(mfilename), outputPath);
        end

        [stericSl, dates, lon, lat] = ...
            formatoutput(squeeze(stericSl), dates, lon, lat, timelim, timeFmt, unit);
        varargout = {stericSl, dates, lon, lat};

        return
    end

    %% Compute data
    inputPath = outputpath(product, [], [], [], []);

    if ~exist(inputPath, 'file')
        error(sprintf('%s:LoadData:MethodNotImplemented', upper(mfilename)), ...
            'Input file %s not found', inputPath);
    end

    load(inputPath, 'lon', 'lat', 'dates', 'stericSl');
    stericSl = squeeze(stericSl);

    if beQuiet <= 1
        fprintf('%s loaded %s\n', upper(mfilename), outputPath);
    end

    if ~isempty(timeStep)
        [stericSl, dates] = ...
            interptemporal(dates, stericSl, timeStep, intpMthd, beQuiet);
    end

    if ~(isempty(meshSize) && isempty(lonOrigin))
        [stericSl, lon, lat] = ...
            interpspatial(lon, lat, stericSl, meshSize, lonOrigin, intpMthd, beQuiet);
    end

    if saveData
        save(outputPath, 'lon', 'lat', 'dates', 'stericSl', '-v7.3');

        if beQuiet <= 1
            fprintf('%s saved %s\n', upper(mfilename), outputPath);
        end

    end

    [stericSl, dates, lon, lat] = ...
        formatoutput(stericSl, dates, lon, lat, timelim, timeFmt, unit);

    varargout = {stericSl, dates, lon, lat};

end

%% Subfunctions
% Interpolation
function [meshIntp, datesIntp] = ...
        interptemporal(dates, mesh, timeStep, intpMthd, beQuiet)

    if beQuiet <= 1
        fprintf('%s: Interpolating temporally, this may take a while...\n', ...
            upper(mfilename));
    end

    if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
        startDate = datetime(year(dates(1)), month(dates(1)), 1, 0, 0, 0);
        endDate = datetime(year(dates(end)), month(dates(end)), 1, 0, 0, 0) + calmonths(1);
        dmonths = ceil((endDate - startDate) / days(28));
        startDates = startDate + calmonths(0:dmonths - 1);
        endDates = startDate + calmonths(1:dmonths);
        datesIntp = startDates + (endDates - startDates) / 2;
        datesIntp = datesIntp(datesIntp >= dates(1) & datesIntp <= dates(end));
    else

        if mean(diff(dates)) > timeStep
            warning(sprintf('%s:InterpolationStepTooSmall', upper(mfilename)), ...
                'The interpolation time step (%s) is smaller than the mean data resolution (%s)', timeStep, mean(diff(dates)));
        end

        datesIntp = dates(1):timeStep:dates(end);
    end

    meshFlat = reshape(mesh, [prod(size(mesh, 1:2)), size(mesh, 3)])';
    meshIntp = interp1(dates, meshFlat, datesIntp, intpMthd)';
    meshIntp = reshape(meshIntp, [size(mesh, 1:2), length(datesIntp)]);
end

function [meshIntp, lonIntp, latIntp] = ...
        interpspatial(lon, lat, mesh, meshSize, lonOrigin, intpMthd, beQuiet)

    if beQuiet <= 1
        fprintf('%s: Interpolating spatially, this may take a while...\n', ...
            upper(mfilename));
    end

    if isempty(meshSize)
        meshSize = 1/2;
    end

    if isempty(lonOrigin)
        lonOrigin = 180;
    end

    ogLonOrigin = 200;
    lonIntp = (-180:meshSize:180) + lonOrigin;
    latIntp = -90:meshSize:90;
    [lonnIntp, lattIntp] = meshgrid(lonIntp, latIntp);
    lonnIntp(lonnIntp < ogLonOrigin - 180) = lonnIntp(lonnIntp < ogLonOrigin - 180) + 360;
    lonnIntp(lonnIntp > ogLonOrigin + 180) = lonnIntp(lonnIntp > ogLonOrigin + 180) - 360;
    lonPad = [lon(end) - 360; lon; lon(1) + 360];
    meshPad = cat(1, mesh(end, :, :), mesh, mesh(1, :, :));
    [lonn, latt] = meshgrid(lonPad, lat);

    meshPad = permute(meshPad, [2, 1, 3]);

    meshIntp = nan([length(latIntp), length(lonIntp), size(mesh, 3)], 'single');

    for iDate = 1:size(meshIntp, 3)
        meshIntp(:, :, iDate) = ...
            interp2(lonn, latt, squeeze(meshPad(:, :, iDate)), ...
            lonnIntp, lattIntp, intpMthd);
        % meshIntp(:, :, iDate) = ...
        %     interp2(lonn, latt, squeeze(meshPad(:, :, iDate)), ...
        %     mod(lonnIntp, 360), lattIntp, intpMthd);
    end

    meshIntp = permute(meshIntp, [2, 1, 3]);
end

% Format output
function varargout = formatoutput(stericSl, dates, lon, lat, timelim, timeFmt, unit)

    if ~isempty(timelim)
        isValidTime = (dates >= timelim(1)) & (dates <= timelim(2));
        stericSl = stericSl(:, :, isValidTime);
        dates = dates(isValidTime);
    end

    if strcmpi(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    if strcmpi(unit, 'mm')
        stericSl = stericSl * 1e3;
    end

    varargout = {stericSl, dates, lon, lat};
end

% Parse input arguments
function varargout = parseinputs(inputs)
    ip = inputParser;
    addOptional(ip, 'Product', 'SIO', ...
        @(x) ischar(validatestring(x, {'SIO', 'EN4'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || ((isnumeric(x) || isdatetime(x)) && length(x) == 2));
    addOptional(ip, 'TimeStep', [], ...
        @(x) isempty(x) || ((isduration(x) || isnumeric(x)) && isscalar(x)) || ischar(validatestring(x, {'midmonth'})));
    addOptional(ip, 'MeshSize', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addOptional(ip, 'LonOrigin', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'InterpolationMethod', 'linear', @ischar);
    addParameter(ip, 'TimeFormat', 'datetime', ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'Unit', 'm', ...
        @(x) ischar(validatestring(x, {'mm', 'm'})));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
    parse(ip, inputs{:});

    product = upper(ip.Results.Product);
    timelim = ip.Results.TimeRange;

    timeStep = ip.Results.TimeStep;
    meshSize = ip.Results.MeshSize;
    lonOrigin = ip.Results.LonOrigin;
    intpMthd = ip.Results.InterpolationMethod;

    timeFmt = lower(ip.Results.TimeFormat);
    unit = lower(ip.Results.Unit);

    forceNew = logical(ip.Results.ForceNew);
    beQuiet = uint8(double(ip.Results.BeQuiet) * 2);
    saveData = logical(ip.Results.SaveData);

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    if ~isempty(timeStep) && isnumeric(timeStep)
        timeStep = days(timeStep);
    end

    if ~strcmpi(product, 'SIO')
        error(sprintf('%s:LoadData:ProductNotImplemented', upper(mfilename)), ...
            'Unsupported product: %s. Currently only SIO is supported', upper(product));
    end

    varargout = ...
        {product, timelim, timeStep, meshSize, lonOrigin, intpMthd, timeFmt, unit, ...
         forceNew, beQuiet, saveData};

end

% Find output path
function outputPath = outputpath(product, timeStep, meshSize, lonOrigin, intpMthd)
    outputFolder = fullfile(getenv("IFILES"), 'HoMAGE', product);

    timeStepStr = '';

    if ~isempty(timeStep)

        if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
            timeStepStr = '-Tmidmonth';
        elseif isduration(timeStep)
            timeStepStr = sprintf('-T%s', erase(sprintf('%s', timeStep), ' '));
        elseif isnumeric(timeStep)
            timeStepStr = sprintf('-T%d', timeStep);
        else
            error(sprintf('%s:InvalidTimeStep', upper(mfilename)), ...
                'Unrecognised time step input for saving: %s', class(timeStep));
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
        intpMethodStr = ['-', lower(intpMthd)];
    end

    outputFile = sprintf('%s_StericSeaLevel%s%s%s.mat', ...
        product, timeStepStr, spaceIntpStr, intpMethodStr);

    outputPath = fullfile(outputFolder, outputFile);
end
