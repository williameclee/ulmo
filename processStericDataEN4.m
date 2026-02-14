function processStericDataEN4(inputFolder, outputFolder, options)

    arguments (Input)
        inputFolder (1, :) char
        outputFolder (1, :) char
        options.ForceNew (1, 1) logical = false
    end

    forceNew = options.ForceNew;
    tlim = [datetime(1990, 1, 1), datetime(2010, 12, 31)];

    zipPattern = 'EN.4.2.2.analyses.c14.*.zip';
    zipFiles = dir(fullfile(inputFolder, zipPattern));
    zipFiles = {zipFiles.name};

    if ~isempty(zipFiles)
        recycle('on');

        for iFile = 1:length(zipFiles)
            zipFile = zipFiles{iFile};
            zipPath = fullfile(inputFolder, zipFile);
            unzip(zipPath, inputFolder);
            delete(zipPath);
        end

    end

    inputPattern = 'EN.4.2.2.f.analysis.c14.*.nc';
    inputFiles = dir(fullfile(inputFolder, inputPattern));
    inputFiles = {inputFiles.name};

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    parfor iFile = 1:length(inputFiles)
        inputFile = inputFiles{iFile};
        inputPath = fullfile(inputFolder, inputFile);
        computedensity(inputPath, outputFolder, ForceNew = forceNew);
    end

    outputPattern = 'EN4c14-M*.mat';
    outputFiles = dir(fullfile(outputFolder, outputPattern));
    outputFiles = {outputFiles.name};

    climPath = fullfile(outputFolder, sprintf('EN4c14-C%s_%s.mat', datetime(tlim, "Format", 'yyyyMM')));

    computeclimatology(tlim, outputFolder, outputFiles, climPath, ...
        ForceNew = forceNew);

    parfor iFile = 1:length(outputFiles)
        outputFile = outputFiles{iFile};
        outputPath = fullfile(outputFolder, outputFile);
        computesteric(outputPath, climPath, ForceNew = true);
    end

    aggregatedPath = fullfile(outputFolder, 'EN4c14-steric.nc');

    aggregatesteric(outputFolder, outputFiles, aggregatedPath, ...
        ForceNew = true);
end

%% Subfunctions
function computedensity(inputPath, outputFolder, options)

    arguments
        inputPath (1, :) char
        outputFolder (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    date = datetime(1800, 1, 1) + days(ncread(inputPath, 'time'));

    outputFile = sprintf('EN4c14-M%s.mat', datetime(date, "Format", 'yyyyMM'));
    outputPath = fullfile(outputFolder, outputFile);

    vars = {'salinity', 'consTemp', 'density', 'lon', 'lat', 'depth', 'date'};

    if ~options.ForceNew && exist(outputPath, 'file') && ...
            all(ismember(vars, who('-file', outputPath)))

        if ~options.BeQuiet
            fprintf('[ULMO>%s] Checked %s of %s, data already exist.\n', ...
                callchaintext([options.CallChain, {mfilename}]), filehref(outputPath, 'steric data'), datetime(date, "Format", 'yyyy/MM'));
        end

        return
    end

    lon = single(ncread(inputPath, 'lon'));
    lat = single(ncread(inputPath, 'lat'));
    depth = single(ncread(inputPath, 'depth'));
    potTemp = single(ncread(inputPath, 'temperature')) - 273.15; % K -> °C
    salinityPsu = single(ncread(inputPath, 'salinity'));

    lon = lon(:)';
    lat = lat(:);
    potTemp = permute(potTemp, [2, 1, 3]);
    salinityPsu = permute(salinityPsu, [2, 1, 3]);

    pres = gsw_p_from_z(repmat(-depth(:)', [length(lat), 1]), lat);

    salinity = nan(size(salinityPsu), 'single');

    for iDepth = 1:length(depth)
        layerSalinityPsu = salinityPsu(:, :, iDepth);
        salinity(:, :, iDepth) = gsw_SA_from_SP(layerSalinityPsu, pres(iDepth), mod(lon, 360), lat);
    end

    consTemp = gsw_CT_from_pt(salinity, potTemp);

    density = nan(size(potTemp), 'single');

    for iDepth = 1:length(depth)
        density(:, :, iDepth) = gsw_rho( ...
            squeeze(salinity(:, :, iDepth)), squeeze(consTemp(:, :, iDepth)), pres(iDepth));
    end

    try
        save(outputPath, vars{:}, '-append');
    catch
        save(outputPath, vars{:}, '-v7.3');
    end

    if ~options.BeQuiet
        fprintf('[ULMO>%s] Computed %s of %s.\n', callchaintext([options.CallChain, {mfilename}]), filehref(outputPath, 'density data'), datetime(date, "Format", 'yyyy/MM'));
    end

end

function computeclimatology(tlim, inputFolder, inputFiles, outputPath, options)

    arguments
        tlim (1, 2) datetime
        inputFolder (1, :) char
        inputFiles (1, :) cell
        outputPath (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
    end

    vars = {'consTempClim', 'salinityClim', 'densityClim'};

    if ~options.ForceNew && exist(outputPath, 'file') && ...
            all(ismember(vars, who('-file', outputPath)))

        if ~options.BeQuiet
            fprintf('Skipped climatology %s\n', outputPath);
        end

        return
    end

    for iFile = 1:length(inputFiles)
        inputFile = inputFiles{iFile};
        inputPath = fullfile(inputFolder, inputFile);

        if ~exist(inputPath, 'file')
            error('Input file %s does not exist.', inputFile);
        end

        if any(~ismember({'salinity', 'consTemp', 'density', 'date'}, who('-file', inputPath)))
            error('Input file %s is missing some variables.', inputFile);
        end

        load(inputPath, 'salinity', 'consTemp', 'density', 'date');

        if date < tlim(1) || date > tlim(2)

            if ~options.BeQuiet
                fprintf('Skipped %s for climatology (out of time range)\n', inputFile);
            end

            continue
        end

        if ~exist('salinityClim', 'var')
            salinityClim = zeros(size(salinity), 'single');
            salinityCnt = zeros(size(salinity), 'uint16');
        end

        if ~exist('consTempClim', 'var')
            consTempClim = zeros(size(consTemp), 'single');
            consTempCnt = zeros(size(consTemp), 'uint16');
        end

        if ~exist('densityClim', 'var')
            densityClim = zeros(size(density), 'single');
            densityCnt = zeros(size(density), 'uint16');
        end

        isValidSalinity = ~isnan(salinity);
        salinityClim(isValidSalinity) = ...
            salinityClim(isValidSalinity) + salinity(isValidSalinity);
        salinityCnt(isValidSalinity) = salinityCnt(isValidSalinity) + 1;

        isValidConsTemp = ~isnan(consTemp);
        consTempClim(isValidConsTemp) = ...
            consTempClim(isValidConsTemp) + consTemp(isValidConsTemp);
        consTempCnt(isValidConsTemp) = consTempCnt(isValidConsTemp) + 1;

        isValidDensity = ~isnan(density);
        densityClim(isValidDensity) = ...
            densityClim(isValidDensity) + density(isValidDensity);
        densityCnt(isValidDensity) = densityCnt(isValidDensity) + 1;
    end

    salinityClim = salinityClim ./ single(salinityCnt);
    salinityClim(salinityCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    consTempClim = consTempClim ./ single(consTempCnt);
    consTempClim(consTempCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    densityClim = densityClim ./ single(densityCnt);
    densityClim(densityCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    save(outputPath, vars{:}, '-v7.3');

    if ~options.BeQuiet
        fprintf('Computed climatology %s\n', outputPath);
    end

end

function computesteric(dataPath, climatologyPath, options)

    arguments
        dataPath (1, :) char
        climatologyPath (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
    end

    if ~exist(dataPath, 'file')
        error('Data file %s does not exist.', dataPath);
    end

    if ~exist(climatologyPath, 'file')
        error('Climatology file %s does not exist.', climatologyPath);
    end

    if any(~ismember({'density', 'depth'}, who('-file', dataPath)))
        error('Data file %s is missing some variables.', dataPath);
    end

    if ~ismember({'densityClim'}, who('-file', climatologyPath))
        error('Climatology file %s is missing some variables.', climatologyPath);
    end

    vars = ...
        {'stericSl', 'shallowStericSl', 'deepStericSl', ...
         'thermostericSl', 'shallowThermostericSl', 'deepThermostericSl', ...
         'halostericSl', 'shallowHalostericSl', 'deepHalostericSl'};

    % Check if stericSl already exists and the file is younger than the climatology
    if ~options.ForceNew && all(ismember(vars, who('-file', dataPath))) && ...
            (dir(dataPath).datenum > dir(climatologyPath).datenum)

        if ~options.BeQuiet
            fprintf('Skipped %s, already computed\n', dataPath);
        end

        return
    end

    load(dataPath, 'density', 'salinity', 'consTemp', 'lat', 'depth');
    load(climatologyPath, 'densityClim', 'consTempClim', 'salinityClim');

    % Compute halosteric density
    if ismember('haloDensity', who('-file', dataPath))
        load(dataPath, 'haloDensity');
    else
        pres = gsw_p_from_z(repmat(-depth(:)', [length(lat), 1]), lat);

        thermoDensity = nan(size(density), 'single');
        haloDensity = nan(size(density), 'single');

        for iDepth = 1:length(depth)
            thermoDensity(:, :, iDepth) = gsw_rho( ...
                squeeze(salinityClim(:, :, iDepth)), squeeze(consTemp(:, :, iDepth)), pres(iDepth));
            haloDensity(:, :, iDepth) = gsw_rho( ...
                squeeze(salinity(:, :, iDepth)), squeeze(consTempClim(:, :, iDepth)), pres(iDepth));
        end

    end

    % Integrate steric sea level
    layerTop = [0; (depth(1:end - 1) + depth(2:end)) / 2];
    layerBottom = [layerTop(2:end); 5500]; % TODO: Use actual depths from https://www.metoffice.gov.uk/hadobs/en4/depths.txt
    layerThickness = abs(layerTop - layerBottom);
    layerThickness(end) = layerThickness(end - 1) ...
        + (layerThickness(end - 1) - layerThickness(end - 2)); % Extrapolate bottom layer thickness

    thermostericSls = (densityClim ./ thermoDensity - 1) .* ...
        reshape(layerThickness, 1, 1, []);
    halostericSls = (densityClim ./ haloDensity - 1) .* ...
        reshape(layerThickness, 1, 1, []);
    thermostericSl = sum(thermostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
    halostericSl = sum(halostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

    stericSls = (densityClim ./ density - 1) .* ...
        reshape(layerThickness, 1, 1, []);
    stericSl = sum(stericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

    isShallow = layerTop < 2000;
    shallowStericSls = (densityClim(:, :, isShallow) ./ density(:, :, isShallow) - 1) .* ...
        reshape(layerThickness(isShallow), 1, 1, []);
    shallowStericSl = sum(shallowStericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
    shallowThermostericSls = (densityClim(:, :, isShallow) ./ thermoDensity(:, :, isShallow) - 1) .* ...
        reshape(layerThickness(isShallow), 1, 1, []);
    shallowThermostericSl = sum(shallowThermostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
    shallowHalostericSls = (densityClim(:, :, isShallow) ./ haloDensity(:, :, isShallow) - 1) .* ...
        reshape(layerThickness(isShallow), 1, 1, []);
    shallowHalostericSl = sum(shallowHalostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

    isDeep = layerTop >= 2000;
    deepStericSls = (densityClim(:, :, isDeep) ./ density(:, :, isDeep) - 1) .* ...
        reshape(layerThickness(isDeep), 1, 1, []);
    deepStericSl = sum(deepStericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
    deepThermostericSls = (densityClim(:, :, isDeep) ./ thermoDensity(:, :, isDeep) - 1) .* ...
        reshape(layerThickness(isDeep), 1, 1, []);
    deepThermostericSl = sum(deepThermostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
    deepHalostericSls = (densityClim(:, :, isDeep) ./ haloDensity(:, :, isDeep) - 1) .* ...
        reshape(layerThickness(isDeep), 1, 1, []);
    deepHalostericSl = sum(deepHalostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

    try
        save(dataPath, vars{:}, '-append');
    catch
        save(dataPath, vars{:}, '-v7.3', '-append');
    end

    if ~options.BeQuiet
        fprintf('Computed steric sea level for %s\n', dataPath);
    end

end

function aggregatesteric(inputFolder, inputFiles, outputPath, options)

    arguments
        inputFolder (1, :) char
        inputFiles (1, :) cell
        outputPath (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
    end

    [~, ~, outputExt] = fileparts(outputPath);

    matStericVars = ...
        {'stericSl', 'shallowStericSl', 'deepStericSl', ...
         'thermostericSl', 'shallowThermostericSl', 'deepThermostericSl', ...
         'halostericSl', 'shallowHalostericSl', 'deepHalostericSl'};

    switch outputExt
        case '.mat'
            vars = ...
                [{'dates', 'lat', 'lon'}, matStericVars];

            if ~options.ForceNew && exist(outputPath, 'file') && ...
                    all(ismember(vars, who('-file', outputPath)))

                if ~options.BeQuiet
                    fprintf('Steric sea level already exists in %s\n', outputPath);
                end

                return
            end

        case '.nc'

            if ~options.ForceNew && exist(outputPath, 'file')
                info = ncinfo(outputPath);
                varNamesInFile = {info.Variables.Name};

                if all(ismember( ...
                        {'time', 'lat', 'lon', 'steric', 'steric_shallow', 'steric_deep', ...
                         'thermosteric', 'thermosteric_shallow', 'thermosteric_deep', ...
                         'halosteric', 'halosteric_shallow', 'halosteric_deep'}, ...
                        varNamesInFile))

                    if ~options.BeQuiet
                        fprintf('Steric sea level already exists in %s\n', outputPath);
                    end

                    return
                end

            end

        otherwise
            error('Unsupported output file extension: %s', outputExt);
    end

    for iFile = 1:length(inputFiles)
        inputFile = inputFiles{iFile};
        inputPath = fullfile(inputFolder, inputFile);

        if ~exist(inputPath, 'file')
            error('Input file %s does not exist.', inputFile);
        end

        if any(~ismember([{'date', 'lat', 'lon'}, matStericVars], who('-file', inputPath)))
            missingVars = setdiff([{'dates', 'lat', 'lon'}, matStericVars], who('-file', inputPath));
            error('Input file %s is missing variables: %s\n', inputFile, strjoin(missingVars, ', '));
        end

        load(inputPath, 'date', matStericVars{:});

        if ~exist('lat', 'var') || ~exist('lon', 'var')
            load(inputPath, 'lat', 'lon');
        end

        if ~exist('stericSls', 'var')
            stericSls = nan([size(stericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            shallowStericSls = nan([size(shallowStericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            deepStericSls = nan([size(deepStericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            thermostericSls = nan([size(thermostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            shallowThermostericSls = nan([size(shallowThermostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            deepThermostericSls = nan([size(deepThermostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            halostericSls = nan([size(halostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            shallowHalostericSls = nan([size(shallowHalostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            deepHalostericSls = nan([size(deepHalostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            dates = NaT([1, length(inputFiles)]);
        end

        stericSls(:, :, iFile) = stericSl;
        shallowStericSls(:, :, iFile) = shallowStericSl;
        deepStericSls(:, :, iFile) = deepStericSl;
        thermostericSls(:, :, iFile) = thermostericSl;
        shallowThermostericSls(:, :, iFile) = shallowThermostericSl;
        deepThermostericSls(:, :, iFile) = deepThermostericSl;
        halostericSls(:, :, iFile) = halostericSl;
        shallowHalostericSls(:, :, iFile) = shallowHalostericSl;
        deepHalostericSls(:, :, iFile) = deepHalostericSl;
        dates(iFile) = date;
    end

    stericSl = stericSls;
    shallowStericSl = shallowStericSls;
    deepStericSl = deepStericSls;
    thermostericSl = thermostericSls;
    shallowThermostericSl = shallowThermostericSls;
    deepThermostericSl = deepThermostericSls;
    halostericSl = halostericSls;
    shallowHalostericSl = shallowHalostericSls;
    deepHalostericSl = deepHalostericSls;

    switch outputExt
        case '.mat'
            save(outputPath, vars{:}, '-v7.3');
        case '.nc'

            if exist(outputPath, 'file')
                delete(outputPath);
            end

            % Steric
            write2nc(outputPath, 'steric', permute(stericSl, [2, 1, 3]), ...
                Dimensions = {'lon', 'lat', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Total steric sea level'});
            write2nc(outputPath, 'steric_shallow', shallowStericSl, ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow total steric sea level above 2000 m'});
            write2nc(outputPath, 'steric_deep', deepStericSl, ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep total steric sea level below 2000 m'});
            write2nc(outputPath, 'thermosteric', permute(thermostericSl, [2, 1, 3]), ...
                Dimensions = {'lon', 'lat', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Thermosteric sea level'});
            write2nc(outputPath, 'thermosteric_shallow', shallowThermostericSl, ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow thermosteric sea level above 2000 m'});
            write2nc(outputPath, 'thermosteric_deep', deepThermostericSl, ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep thermosteric sea level below 2000 m'});
            write2nc(outputPath, 'halosteric', permute(halostericSl, [2, 1, 3]), ...
                Dimensions = {'lon', 'lat', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Halosteric sea level'});
            write2nc(outputPath, 'halosteric_shallow', shallowHalostericSl, ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow halosteric sea level above 2000 m'});
            write2nc(outputPath, 'halosteric_deep', deepHalostericSl, ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep halosteric sea level below 2000 m'});
            nccreate(outputPath, 'lat', ...
                Dimensions = {'lat', length(lat)}, ...
                Datatype = 'single', Format = 'netcdf4');
            ncwrite(outputPath, 'lat', lat);
            ncwriteatt(outputPath, 'lat', "Units", 'degrees north');
            nccreate(outputPath, 'lon', ...
                Dimensions = {'lon', length(lon)}, ...
                Datatype = 'single', Format = 'netcdf4');
            ncwrite(outputPath, 'lon', lon);
            ncwriteatt(outputPath, 'lon', "Units", 'degrees east');
            nccreate(outputPath, 'time', ...
                Dimensions = {'time', length(dates)}, ...
                Datatype = 'single', Format = 'netcdf4');
            timeDays = single(days(dates - datetime(1800, 1, 1)));
            ncwrite(outputPath, 'time', timeDays);
            ncwriteatt(outputPath, 'time', "Units", 'days since 1800-01-01');
        otherwise
            error('Unsupported output file extension: %s', outputExt);
    end

    if ~options.BeQuiet
        fprintf('Aggregated steric sea level to %s\n', outputPath);
    end

end

function write2nc(outputPath, varname, var, options)

    arguments (Input)
        outputPath (1, :) char
        varname (1, :) char
        var (:, :, :) single
        options.Dimensions (1, :) cell = NaN
        options.Attributes (1, :) cell = NaN
    end

    if iscell(options.Dimensions) && length(options.Dimensions) ~= ndims(var)
        error ('Name of dimensions provided does not match variable dimensions (length = %d vs %d).', ...
            length(options.Dimensions), ndims(var));
    elseif iscell(options.Dimensions)
        % Pre-allocate dimensions cell array: name and size for each dimension
        dimensions = cell(1, 2 * length(options.Dimensions));

        for iDim = 1:length(options.Dimensions)
            idxName  = 2 * iDim - 1;
            idxSize  = 2 * iDim;
            dimensions{idxName} = options.Dimensions{iDim};
            dimensions{idxSize} = size(var, iDim);
        end

    else
        % Pre-allocate dimensions cell array when names are autogenerated
        dimensions = cell(1, 2 * ndims(var));

        for iDim = 1:ndims(var)
            idxName  = 2 * iDim - 1;
            idxSize  = 2 * iDim;
            dimensions{idxName} = sprintf('dim%d', iDim);
            dimensions{idxSize} = size(var, iDim);
        end

    end

    nccreate(outputPath, varname, ...
        "Dimensions", dimensions, ...
        "Datatype", class(var), "Format", 'netcdf4');
    ncwrite(outputPath, varname, var);

    if iscell(options.Attributes)

        if mod(length(options.Attributes), 2) ~= 0
            error('Attributes must be provided as name-value pairs (got length = %s).', ...
                length(options.Attributes));
        end

        for iAttr = 1:2:length(options.Attributes)
            attrName = options.Attributes{iAttr};
            attrValue = options.Attributes{iAttr + 1};
            ncwriteatt(outputPath, varname, attrName, attrValue);
        end

    end

end
