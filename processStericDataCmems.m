%% PROCESSSTERICDATACMEMS - Process steric sea level from CMEMS data
%
% Author
%	2026/02/15, williameclee@arizona.edu (@williameclee)

function processStericDataCmems(inputFolder, outputFolder, aggregatePath, climatologyTimeRange, options)
    % Process steric sea level data from CMEMS data

    arguments (Input)
        inputFolder (1, :) char
        outputFolder (1, :) char
        aggregatePath (1, :) char = fullfile(outputFolder, 'CMEMS-steric.nc');
        climatologyTimeRange (1, 2) datetime = [datetime(1990, 1, 1), datetime(2010, 12, 31)];
        options.ForceNew (1, 1) logical = false;
        options.CallChain (1, :) cell = {};
    end

    tlim = climatologyTimeRange;

    forceNew = options.ForceNew;
    callChain = [options.CallChain, {mfilename}];

    % Step 1: Break aggregated input into individual files
    outputFiles = breakAggregatedInput(inputFolder, outputFolder, ...
        ForceNew = forceNew, CallChain = callChain);

    % % Step 2: Compute density for each file

    % % Step 3: Compute climatology
    % climPath = fullfile(outputFolder, sprintf('CMEMS-C%s_%s.mat', datetime(tlim, "Format", 'yyyyMM')));

    % computeclimatology(tlim, outputFolder, outputFiles, climPath, ...
    %     ForceNew = forceNew, CallChain = callChain);

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

    cprintf('[ULMO>%s] Found %d time steps in CMEMS data.\n', callchaintext(callChain), length(dates));

    % Process each time step and save to individual .mat files
    for iDate = 1:length(dates)
        date = dates(iDate);

        % Read salinity and temperature for the current time step
        salinityPsu = single(ncread(salinityPath, 'so', [1, 1, 1, iDate], [Inf, Inf, Inf, 1]));
        potTemp = single(ncread(tempPath, 'thetao', [1, 1, 1, iDate], [Inf, Inf, Inf, 1]));

        salinityPsu = permute(salinityPsu, [2, 1, 3]); % (lat, lon, z)
        potTemp = permute(potTemp, [2, 1, 3]); % (lat, lon, z)

        % Save the data to a .mat file
        outputFile = sprintf('CMEMS-M%s.mat', datetime(date, "Format", 'yyyyMM'));
        outputPath = fullfile(outputFolder, outputFile);

        outputFiles{iDate} = outputFile;

        if exist(outputPath, 'file') && ~forceNew && ...
                all(ismember({'salinityPsu', 'potTemp', 'lat', 'lon', 'depth', 'date'}, who('-file', outputPath)))
            cprintf('[ULMO>%s] Skipped saving %s %s, already exists.\n', ...
                callchaintext(callChain), datetime(date, "Format", 'yyyyMM'), filehref(outputPath, 'T-S data'));
            continue;
        end

        save(outputPath, 'salinityPsu', 'potTemp', 'lat', 'lon', 'depth', 'date', '-v7.3');
        cprintf('[ULMO>%s] Saved %s %s (%d/%d).\n', ...
            callchaintext(callChain), datetime(date, "Format", 'yyyyMM'), filehref(outputPath, 'T-S data'), ...
            iDate, length(dates));
    end

end

function computedensity(inputPath, outputFolder, options)
    % Compute density from CMEMS data

    arguments
        inputPath (1, :) char
        outputFolder (1, :) char
        options.ForceNew (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    % Read data from netCDF file
    date = datetime(1950, 1, 1) + seconds(ncread(inputPath, 'time'));
    lat = ncread(inputPath, 'latitude');
    lon = ncread(inputPath, 'longitude');
    depth = ncread(inputPath, 'depth');
    salinityPsu = ncread(inputPath, 'so'); % Practical Salinity Units
    potTemp = ncread(inputPath, 'thetao'); % Potential Temperature

    % Convert salinity and temperature to absolute salinity and in-situ density
    pres = gsw_p_from_z(repmat(-depth(:)', [length(lat), 1]), lat);
    salinity = nan(size(salinityPsu));
    density = nan(size(potTemp));

    for iDepth = 1:length(depth)
        salinity(:, :, iDepth) = ...
            gsw_SA_from_SP(squeeze(salinityPsu(:, :, iDepth)), pres(:, iDepth), mod(lon, 360), lat);

        density(:, :, iDepth) = gsw_rho( ...
            squeeze(salinity(:, :, iDepth)), squeeze(potTemp(:, :, iDepth)), pres(:, iDepth));
    end

    % Save density data
    outputFile = sprintf('CMEMS-M%s.mat', datetime(date, "Format", 'yyyyMM'));
    outputPath = fullfile(outputFolder, outputFile);

    if exist(outputPath, 'file')
        save(outputPath, 'density', '-v7.3', '-append');
    else
        save(outputPath, 'depth', 'date', 'lat', 'lon', 'density', '-v7.3');
    end

end

% The computeclimatology, computesteric, and aggregatesteric functions can be reused as-is from processStericDataEN4.m
