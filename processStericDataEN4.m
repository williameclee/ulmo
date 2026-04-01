%% PROCESSSTERICDATAEN4 - Process steric sea level from EN4 data
%
% Author
%	2026/02/14, williameclee@arizona.edu (@williameclee)
%
% Last modified
%   2026/02/17, williameclee@arizona.edu (@williameclee)
%     - Extracted auxiliary functions to their own files for reusability

function processStericDataEN4(inputFolder, outputFolder, aggregatePath, climatologyTimeRange, options)

    arguments (Input)
        inputFolder (1, :) char
        outputFolder (1, :) char
        aggregatePath (1, :) char = fullfile(outputFolder, 'EN4c14-steric.nc')
        climatologyTimeRange (1, 2) datetime = [datetime(1990, 1, 1), datetime(2010, 12, 31)]
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    tlim = climatologyTimeRange;

    forceNew = options.ForceNew;
    beQuiet = options.BeQuiet;
    callChain = [options.CallChain, {mfilename}];

    % Unzip if needed
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

    %% Main processing steps
    % Compute density for each file
    parfor iFile = 1:length(inputFiles)
        inputFile = inputFiles{iFile};
        inputPath = fullfile(inputFolder, inputFile);
        computeDensity(inputPath, outputFolder, ForceNew = forceNew, CallChain = callChain);
    end

    outputPattern = 'EN4c14-M*.mat';
    outputFiles = dir(fullfile(outputFolder, outputPattern));
    outputFiles = {outputFiles.name};

    % Compute climatology
    climPath = fullfile(outputFolder, sprintf('EN4c14-C%s_%s.mat', datetime(tlim, "Format", 'yyyyMM')));

    computeStericClimatology(tlim, outputFolder, outputFiles, climPath, ...
        ForceNew = forceNew, BeQuiet = beQuiet, CallChain = callChain);

    % Compute halosteric and thermosteric densities
    parfor iFile = 1:length(outputFiles)
        outputFile = outputFiles{iFile};
        outputPath = fullfile(outputFolder, outputFile);
        computeStericDensityVar(outputPath, climPath, ...
            ForceNew = forceNew, BeQuiet = beQuiet, CallChain = callChain);
    end

    % Compute steric sea level anomalies
    parfor iFile = 1:length(outputFiles)
        outputFile = outputFiles{iFile};
        outputPath = fullfile(outputFolder, outputFile);
        computeStericSeaLevel(outputPath, climPath, Bottom = 5500, ...
            ForceNew = forceNew, CallChain = callChain);
    end

    aggregateStericSeaLevel(outputFolder, outputFiles, aggregatePath, ...
        ForceNew = forceNew, CallChain = callChain);
end

%% Subfunctions
function computeDensity(inputPath, outputFolder, options)

    arguments (Input)
        inputPath {mustBeFile}
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
