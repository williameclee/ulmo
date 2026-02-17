%% COMPUTESTERICCLIMATOLOGY - Computes climatology of temperature, salinity, and density
%
% Last modified
%   2026/02/16, williameclee@arizona.edu (@williameclee)
%     - Extracted from PROCESSSTERICDATAEN4 for reusability

function computeStericClimatology(tlim, inputFolder, inputFiles, outputPath, options)

    arguments (Input)
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

    if ~options.BeQuiet
        msg = 'this may take a while...';
        cprintf('[ULMO>%s] Computing %s from %d files, %s\n', ...
            callchaintext(callChain), filehref(outputPath, 'climatology data'), length(inputFiles), msg);
        tic;
    end

    cnt = 0;

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
                fprintf(repmat('\b', 1, cnt));
                cnt = cprintf('[ULMO>%s] Skipped %s for climatology, out of time range (%d/%d).\n', ...
                    callchaintext(callChain), filehref(inputPath, 'data'), iFile, length(inputFiles));
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

        if ~options.BeQuiet
            fprintf(repmat('\b', 1, cnt));
            cnt = cprintf('[ULMO>%s] Processed %s for climatology (%d/%d).\n', ...
                callchaintext(callChain), filehref(inputPath, 'data'), iFile, length(inputFiles));
        end

    end

    fprintf(repmat('\b', 1, cnt));

    assert(numClimFiles > 0, 'No files found for climatology in the specified time range.');

    salinityClim = salinityClim ./ single(salinityCnt);
    salinityClim(salinityCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    consTempClim = consTempClim ./ single(consTempCnt);
    consTempClim(consTempCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    densityClim = densityClim ./ single(densityCnt);
    densityClim(densityCnt == 0) = nan; %#ok<NASGU> - actually saved through VARS variable

    if ~options.BeQuiet
        fprintf(repmat('\b', 1, length(msg) + 1));
        cprintf('took %.1f seconds.\n', toc);
    end

    save(outputPath, vars{:}, '-v7.3');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Saved %s.\n', ...
            callchaintext(callChain), filehref(outputPath, 'climatology data'));
    end

end
