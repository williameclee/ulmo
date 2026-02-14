function processStericDataEN4(inputFolder, outputFolder, aggregatePath, climatologyTimeRange, options)

    arguments (Input)
        inputFolder (1, :) char
        outputFolder (1, :) char
        aggregatePath (1, :) char = fullfile(outputFolder, 'EN4c14-steric.nc')
        climatologyTimeRange (1, 2) datetime = [datetime(1990, 1, 1), datetime(2010, 12, 31)]
        options.ForceNew (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    tlim = climatologyTimeRange;

    forceNew = options.ForceNew;
    callChain = [options.CallChain, {mfilename}];

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
        computedensity(inputPath, outputFolder, ForceNew = forceNew, CallChain = callChain);
    end

    outputPattern = 'EN4c14-M*.mat';
    outputFiles = dir(fullfile(outputFolder, outputPattern));
    outputFiles = {outputFiles.name};

    climPath = fullfile(outputFolder, sprintf('EN4c14-C%s_%s.mat', datetime(tlim, "Format", 'yyyyMM')));

    computeclimatology(tlim, outputFolder, outputFiles, climPath, ...
        ForceNew = forceNew, CallChain = callChain);

    parfor iFile = 1:length(outputFiles)
        outputFile = outputFiles{iFile};
        outputPath = fullfile(outputFolder, outputFile);
        computesteric(outputPath, climPath, ForceNew = forceNew, CallChain = callChain);
    end

    aggregatesteric(outputFolder, outputFiles, aggregatePath, ...
        ForceNew = forceNew, CallChain = callChain);
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

    callChain = [options.CallChain, {mfilename}];

    date = datetime(1800, 1, 1) + days(ncread(inputPath, 'time'));

    outputFile = sprintf('EN4c14-M%s.mat', datetime(date, "Format", 'yyyyMM'));
    outputPath = fullfile(outputFolder, outputFile);

    vars = {'salinity', 'consTemp', 'density', 'lon', 'lat', 'depth', 'date'};

    if ~options.ForceNew && exist(outputPath, 'file') && ...
            all(ismember(vars, who('-file', outputPath)))

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exist.\n', ...
                callchaintext(callChain), datetime(date, "Format", 'yyyy/MM'), filehref(outputPath, 'density data'));
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
        cprintf('[ULMO>%s] Computed %s %s.\n', callchaintext(callChain), ...
            datetime(date, "Format", 'yyyy/MM'), filehref(outputPath, 'density data'));
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
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    vars = {'consTempClim', 'salinityClim', 'densityClim'};

    if ~options.ForceNew && exist(outputPath, 'file') && ...
            all(ismember(vars, who('-file', outputPath)))

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s, already exist.\n', ...
                callchaintext(callChain), filehref(outputPath, 'climatology data'));
        end

        return
    end

    numClimFiles = 0;

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
                cprintf('[ULMO>%s] Skipped %s for climatology (out of time range).\n', ...
                    callchaintext(callChain), filehref(inputPath, 'data'));
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

        numClimFiles = numClimFiles + 1;
    end

    assert(numClimFiles > 0, 'No files found for climatology in the specified time range.');

    salinityClim = salinityClim ./ single(salinityCnt);
    salinityClim(salinityCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    consTempClim = consTempClim ./ single(consTempCnt);
    consTempClim(consTempCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    densityClim = densityClim ./ single(densityCnt);
    densityClim(densityCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    save(outputPath, vars{:}, '-v7.3');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s from %d files.\n', ...
            callchaintext(callChain), filehref(outputPath, 'climatology data'), numClimFiles);
    end

end

function computesteric(dataPath, climatologyPath, options)

    arguments
        dataPath (1, :) char
        climatologyPath (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

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
            cprintf('[ULMO>%s] Skipped computing %s, already exist and is newer than climatology.\n', ...
                callchaintext(callChain), filehref(dataPath, 'steric sea level data'));
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

        haloDensity = nan(size(density), 'single');

        for iDepth = 1:length(depth)
            haloDensity(:, :, iDepth) = gsw_rho( ...
                squeeze(salinity(:, :, iDepth)), squeeze(consTempClim(:, :, iDepth)), pres(iDepth));
        end

    end

    % Compute thermosteric density
    if ismember('thermoDensity', who('-file', dataPath))
        load(dataPath, 'thermoDensity');
    else
        pres = gsw_p_from_z(repmat(-depth(:)', [length(lat), 1]), lat);

        thermoDensity = nan(size(density), 'single');

        for iDepth = 1:length(depth)
            thermoDensity(:, :, iDepth) = gsw_rho( ...
                squeeze(salinityClim(:, :, iDepth)), squeeze(consTemp(:, :, iDepth)), pres(iDepth));
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
        cprintf('[ULMO>%s] Computed %s.\n', ...
            callchaintext(callChain), filehref(dataPath, 'steric sea level data'));
    end

end

function aggregatesteric(inputFolder, inputFiles, outputPath, options)

    arguments
        inputFolder (1, :) char
        inputFiles (1, :) cell
        outputPath (1, :) char
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

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
                    cprintf('[ULMO>%s] Skipped computing %s, already exist.\n', ...
                        callchaintext(callChain), filehref(outputPath, 'steric sea level data'));
                end

                return
            end

        case '.nc'

            if ~options.ForceNew && exist(outputPath, 'file')

                try
                    info = ncinfo(outputPath);
                    varNamesInFile = {info.Variables.Name};

                    if all(ismember( ...
                            {'time', 'lat', 'lon', 'steric', 'steric_shallow', 'steric_deep', ...
                             'thermosteric', 'thermosteric_shallow', 'thermosteric_deep', ...
                             'halosteric', 'halosteric_shallow', 'halosteric_deep'}, ...
                            varNamesInFile))

                        if ~options.BeQuiet
                            cprintf('[ULMO>%s] Skipped computing %s, already exist.\n', ...
                                callchaintext(callChain), filehref(outputPath, 'steric sea level data'));
                        end

                        return
                    end

                catch ME
                    warning('Error checking existing NetCDF file: %s\nProceeding to create new file. Error details:\n%s', ...
                        outputPath, getReport(ME));
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
            missingVars = setdiff([{'date', 'lat', 'lon'}, matStericVars], who('-file', inputPath));
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

            % Coordinates
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

            % Variables
            write2nc(outputPath, 'steric', single(stericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Total steric sea level'});
            write2nc(outputPath, 'steric_shallow', single(shallowStericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow total steric sea level above 2000 m'});
            write2nc(outputPath, 'steric_deep', single(deepStericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep total steric sea level below 2000 m'});
            write2nc(outputPath, 'thermosteric', single(thermostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Thermosteric sea level'});
            write2nc(outputPath, 'thermosteric_shallow', single(shallowThermostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow thermosteric sea level above 2000 m'});
            write2nc(outputPath, 'thermosteric_deep', single(deepThermostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep thermosteric sea level below 2000 m'});
            write2nc(outputPath, 'halosteric', single(halostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Halosteric sea level'});
            write2nc(outputPath, 'halosteric_shallow', single(shallowHalostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow halosteric sea level above 2000 m'});
            write2nc(outputPath, 'halosteric_deep', single(deepHalostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep halosteric sea level below 2000 m'});

            try
                ncinfo(outputPath);
            catch ME
                error('Error creating NetCDF file: %s\nDetails:\n%s', outputPath, getReport(ME));
            end

        otherwise
            error('Unsupported output file extension: %s', outputExt);
    end

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Aggregated %s.\n', ...
            callchaintext(callChain), filehref(outputPath, 'steric sea level data'));
    end

end
