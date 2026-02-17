%% COMPUTESTERICSEALEVEL - Computes sea level anomalies from density and climatology
%
% Last modified
%   2026/02/17, williameclee@arizona.edu (@williameclee)
%     - Extracted from PROCESSSTERICDATAEN4 for reusability

function computeStericSeaLevel(dataPath, climatologyPath, options)

    arguments (Input)
        dataPath {mustBeTextScalar, mustBeFile}
        climatologyPath {mustBeTextScalar, mustBeFile}
        options.Bottom (1, 1) double = 6000
        options.HasDeepLayer (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    bottom = options.Bottom;
    hasDeepLayer = options.HasDeepLayer;
    callChain = [options.CallChain, {mfilename}];

    if ~exist(dataPath, 'file')
        error('Data file %s does not exist.', dataPath);
    elseif ~exist(climatologyPath, 'file')
        error('Climatology file %s does not exist.', climatologyPath);
    elseif any(~ismember({'density', 'depth'}, who('-file', dataPath)))
        error('Data file %s is missing some variables.', dataPath);
    elseif ~ismember({'densityClim'}, who('-file', climatologyPath))
        error('Climatology file %s is missing some variables.', climatologyPath);
    end

    vars = ...
        {'stericSl', 'thermostericSl', 'halostericSl'};
    depthVars = ...
        {'shallowStericSl', 'deepStericSl', ...
         'shallowThermostericSl', 'deepThermostericSl', ...
         'shallowHalostericSl', 'deepHalostericSl'};

    % Check if stericSl already exists and the file is younger than the climatology
    if ~options.ForceNew && ...
            all(ismember(vars, who('-file', dataPath))) && ...
            (~hasDeepLayer || all(ismember(depthVars, who('-file', dataPath)))) && ...
            (dir(dataPath).datenum > dir(climatologyPath).datenum)

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s, already exist and is newer than climatology.\n', ...
                callchaintext(callChain), filehref(dataPath, 'steric sea level data'));
        end

        return
    end

    load(dataPath, 'date', 'density', 'salinity', 'consTemp', 'lat', 'depth');
    load(climatologyPath, 'densityClim', 'consTempClim', 'salinityClim');

    % Compute halosteric density
    if ismember('haloDensity', who('-file', dataPath))
        load(dataPath, 'haloDensity');
    else
        pres = gsw_p_from_z(repmat(-depth(:)', [length(lat), 1]), lat);

        haloDensity = nan(size(density), 'single');

        for iDepth = 1:length(depth)
            haloDensity(:, :, iDepth) = gsw_rho( ...
                squeeze(salinity(:, :, iDepth)), squeeze(consTempClim(:, :, iDepth)), pres(1, iDepth));
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
                squeeze(salinityClim(:, :, iDepth)), squeeze(consTemp(:, :, iDepth)), pres(1, iDepth));
        end

    end

    % Integrate steric sea level
    layerTop = [0; (depth(1:end - 1) + depth(2:end)) / 2];
    layerBottom = [layerTop(2:end); bottom];
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

    if hasDeepLayer
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
    end

    try
        save(dataPath, vars{:}, '-append');
    catch
        save(dataPath, vars{:}, '-v7.3', '-append');
    end

    if hasDeepLayer

        try
            save(dataPath, depthVars{:}, '-append');
        catch
            save(dataPath, depthVars{:}, '-v7.3', '-append');
        end

    end

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', ...
            callchaintext(callChain), datetime(date, "Format", 'yyyyMM'), filehref(dataPath, 'steric sea level data'));
    end

end
