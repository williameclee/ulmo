%% PROCESSSTERICDATACMEMS - Process steric sea level from CMEMS data
%
% Author
%	2026/02/15, williameclee@arizona.edu (@williameclee)

function processStericDataCmems(inputFolder, outputFolder, aggregatePath, climatologyTimeRange, options)

    arguments (Input)
        inputFolder (1, :) char
        outputFolder (1, :) char
        aggregatePath (1, :) char = fullfile(outputFolder, 'CMEMS-steric.nc');
        climatologyTimeRange (1, 2) datetime = [datetime(1990, 1, 1), datetime(2010, 12, 31)];
        options.ForceNew (1, 1) logical = false;
        options.BeQuiet (1, 1) logical = false;
        options.CallChain (1, :) cell = {};
    end

    tlim = climatologyTimeRange;

    forceNew = options.ForceNew;
    beQuiet = options.BeQuiet;
    callChain = [options.CallChain, {mfilename}];

    %% Main processing steps
    % Break aggregated input into individual files
    outputFiles = breakAggregatedInput(inputFolder, outputFolder, ...
        ForceNew = forceNew, BeQuiet = beQuiet, CallChain = callChain);

    % Compute density for each file
    parfor iFile = 1:length(outputFiles)
        outputFile = outputFiles{iFile};
        outputPath = fullfile(outputFolder, outputFile);
        computeDensity(outputPath, ForceNew = forceNew, BeQuiet = beQuiet, CallChain = callChain);
    end

    % Compute climatology
    climPath = fullfile(outputFolder, sprintf('CMEMS-C%s_%s.mat', datetime(tlim, "Format", 'yyyyMM')));

    computeclimatology(tlim, outputFolder, outputFiles, climPath, ...
        ForceNew = forceNew, CallChain = callChain);

    % % Step 4: Compute steric sea level anomalies
    % parfor iFile = 1:length(outputFiles)
    %     outputFile = outputFiles{iFile};
    %     outputPath = fullfile(outputFolder, outputFile);
    %     computesteric(outputPath, climPath, ForceNew = forceNew, CallChain = callChain);
    % end

    % % Step 5: Aggregate steric data
    % aggregatesteric(outputFolder, outputFiles, aggregatePath, ...
    %     ForceNew = forceNew, CallChain = callChain);
end

%% Subfunctions
function outputFiles = breakAggregatedInput(inputFolder, outputFolder, options)

    arguments (Input)
        inputFolder {mustBeTextScalar, mustBeFolder}
        outputFolder {mustBeTextScalar}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    arguments (Output)
        outputFiles (1, :) cell
    end

    forceNew = options.ForceNew;
    callChain = [options.CallChain, {mfilename}];

    % Define input file paths for CMEMS data
    salinityFilePattern = "cmems_mod_glo_phy_my_0.083deg_P1M-m_so*.nc";
    tempFilePattern = "cmems_mod_glo_phy_my_0.083deg_P1M-m_thetao*.nc";
    salinityFile = dir(fullfile(inputFolder, salinityFilePattern));
    tempFile = dir(fullfile(inputFolder, tempFilePattern));
    salinityFile = salinityFile.name;
    tempFile = tempFile.name;

    assert(~isempty(salinityFile), 'No salinity file found matching pattern: %s', salinityFilePattern);
    assert(~isempty(tempFile), 'No temperature file found matching pattern: %s', tempFilePattern);

    salinityPath = fullfile(inputFolder, salinityFile);
    tempPath = fullfile(inputFolder, tempFile);

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Read time variable from both files and ensure they match
    timeSalinity = ncread(salinityPath, 'time');
    timeTemp = ncread(tempPath, 'time');

    if ~isequal(timeSalinity, timeTemp)
        error('Time steps in salinity and temperature files do not match.');
    end

    dates = datetime(1950, 1, 1) + hours(timeSalinity);

    lat = ncread(salinityPath, 'latitude');
    lon = ncread(salinityPath, 'longitude');
    depth = ncread(salinityPath, 'depth');

    outputFiles = cell(length(dates), 1);

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Found %d time steps in CMEMS data.\n', callchaintext(callChain), length(dates));
    end

    % Process each time step and save to individual .mat files
    for iDate = 1:length(dates)
        date = dates(iDate);

        % Save the data to a .mat file
        outputFile = sprintf('CMEMS-M%s.mat', datetime(date, "Format", 'yyyyMM'));
        outputPath = fullfile(outputFolder, outputFile);

        outputFiles{iDate} = outputFile;

        if exist(outputPath, 'file') && ~forceNew && ...
                all(ismember({'salinityPsu', 'potTemp', 'lat', 'lon', 'depth', 'date'}, who('-file', outputPath)))

            if ~options.BeQuiet
                cprintf('[ULMO>%s] Skipped saving %s %s, already exists.\n', ...
                    callchaintext(callChain), datetime(date, "Format", 'yyyyMM'), filehref(outputPath, 'T-S data'));
            end

            continue;
        end

        % Read salinity and temperature for the current time step
        salinityPsu = single(ncread(salinityPath, 'so', [1, 1, 1, iDate], [Inf, Inf, Inf, 1]));
        potTemp = single(ncread(tempPath, 'thetao', [1, 1, 1, iDate], [Inf, Inf, Inf, 1]));

        salinityPsu = permute(salinityPsu, [2, 1, 3]); % (lat, lon, z)
        potTemp = permute(potTemp, [2, 1, 3]); % (lat, lon, z)

        save(outputPath, 'salinityPsu', 'potTemp', 'lat', 'lon', 'depth', 'date', '-v7.3');

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Saved %s %s (%d/%d).\n', ...
                callchaintext(callChain), datetime(date, "Format", 'yyyyMM'), filehref(outputPath, 'T-S data'), ...
                iDate, length(dates));
        end

    end

end

function computeDensity(dataPath, options)

    arguments (Input)
        dataPath (1, :) char {mustBeFile}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    % Load data from .mat file
    inputVars = {'salinityPsu', 'potTemp', 'lat', 'lon', 'depth', 'date'};

    if any(~ismember(inputVars, who('-file', dataPath)))
        missinginputVars = setdiff(inputVars, who('-file', dataPath));
        error('Input file %s is missing required variables: %s', ...
            dataPath, strjoin(missinginputVars, ', '));
    elseif ismember('density', who('-file', dataPath)) && ~options.ForceNew

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exists.\n', ...
                callchaintext(callChain), datetime(date, "Format", 'yyyy/MM'), filehref(dataPath, 'density data'));
        end

        return
    end

    load(dataPath, inputVars{:});

    % Convert salinity and temperature to absolute salinity and in-situ density
    pres = gsw_p_from_z(repmat(-depth(:)', [length(lat), 1]), lat); %#ok<USENS> - actually loaded through the inputVars list
    salinity = nan(size(salinityPsu), 'single'); %#ok<USENS> - actually loaded through the inputVars list

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

    % Save density data back to the same .mat file
    save(dataPath, 'density', '-append');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', callchaintext(callChain), ...
            datetime(date, "Format", 'yyyy/MM'), filehref(dataPath, 'density data'));
    end

end
