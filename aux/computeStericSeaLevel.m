%% COMPUTESTERICSEALEVEL - Computes sea level anomalies from density and climatology
%
% Last modified
%   2026/03/02, williameclee@arizona.edu (@williameclee)
%     - Corrected ocean bottom thickness calculation
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
        missingDataVars = setdiff({'density', 'depth'}, who('-file', dataPath));
        error('Data file %s is missing required variables: %s', ...
            dataPath, strjoin(missingDataVars, ', '));
    elseif any(~ismember({'densityClim', 'consTempClim', 'salinityClim'}, who('-file', climatologyPath)))
        missingClimVars = setdiff({'densityClim', 'consTempClim', 'salinityClim'}, who('-file', climatologyPath));
        error('Climatology file %s is missing required variables: %s', ...
            climatologyPath, strjoin(missingClimVars, ', '));
    end

    vars = ...
        {'stericSl', 'thermostericSl', 'halostericSl'};
    depthVars = ...
        {'shallowStericSl', 'deepStericSl', ...
         'shallowThermostericSl', 'deepThermostericSl', ...
         'shallowHalostericSl', 'deepHalostericSl'};
    inputVars = {'density', 'haloDensity', 'thermoDensity', 'depth'};
    inputClimVars = {'densityClim'};

    ddata = load(dataPath, 'date');

    % Check if date variable exists and is in datetime format
    if ~isfield(ddata, 'date') || ~isa(ddata.date, 'datetime')
        error('Input file %s is missing required variable "date" in datetime format.', dataPath);
    end

    % Check if stericSl already exists and the file is younger than the climatology
    if ~options.ForceNew && ...
            all(ismember(vars, who('-file', dataPath))) && ...
            (~hasDeepLayer || all(ismember(depthVars, who('-file', dataPath)))) && ...
            (dir(dataPath).datenum > dir(climatologyPath).datenum)

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exist and is newer than climatology.\n', ...
                callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'steric sea level data'));
        end

        return
    elseif any(~ismember(inputVars, who('-file', dataPath)))
        missingDataVars = setdiff(inputVars, who('-file', dataPath));
        error('Data file %s is missing required variables: %s', ...
            dataPath, strjoin(missingDataVars, ', '));
    elseif any(~ismember(inputClimVars, who('-file', climatologyPath)))
        missingClimVars = setdiff(inputClimVars, who('-file', climatologyPath));
        error('Climatology file %s is missing required variables: %s', ...
            climatologyPath, strjoin(missingClimVars, ', '));
    end

    data = load(dataPath, inputVars{:});
    cdata = load(climatologyPath, inputClimVars{:});

    % Integrate steric sea level
    layerTop = [0; (data.depth(1:end - 1) + data.depth(2:end)) / 2];
    layerBottom = [layerTop(2:end); bottom];
    layerThickness = abs(layerTop - layerBottom);
    % layerThickness(end) = layerThickness(end - 1) ...
    %     + (layerThickness(end - 1) - layerThickness(end - 2)); % Extrapolate bottom layer thickness

    thermostericSls = (cdata.densityClim ./ data.thermoDensity - 1) .* ...
        reshape(layerThickness, 1, 1, []);
    halostericSls = (cdata.densityClim ./ data.haloDensity - 1) .* ...
        reshape(layerThickness, 1, 1, []);
    thermostericSl = sum(thermostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
    halostericSl = sum(halostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

    stericSls = (cdata.densityClim ./ data.density - 1) .* ...
        reshape(layerThickness, 1, 1, []);
    stericSl = sum(stericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

    if hasDeepLayer
        isShallow = layerTop < 2000;
        shallowStericSls = (cdata.densityClim(:, :, isShallow) ./ data.density(:, :, isShallow) - 1) .* ...
            reshape(layerThickness(isShallow), 1, 1, []);
        shallowStericSl = sum(shallowStericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
        shallowThermostericSls = (cdata.densityClim(:, :, isShallow) ./ data.thermoDensity(:, :, isShallow) - 1) .* ...
            reshape(layerThickness(isShallow), 1, 1, []);
        shallowThermostericSl = sum(shallowThermostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
        shallowHalostericSls = (cdata.densityClim(:, :, isShallow) ./ data.haloDensity(:, :, isShallow) - 1) .* ...
            reshape(layerThickness(isShallow), 1, 1, []);
        shallowHalostericSl = sum(shallowHalostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable

        isDeep = layerTop >= 2000;
        deepStericSls = (cdata.densityClim(:, :, isDeep) ./ data.density(:, :, isDeep) - 1) .* ...
            reshape(layerThickness(isDeep), 1, 1, []);
        deepStericSl = sum(deepStericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
        deepThermostericSls = (cdata.densityClim(:, :, isDeep) ./ data.thermoDensity(:, :, isDeep) - 1) .* ...
            reshape(layerThickness(isDeep), 1, 1, []);
        deepThermostericSl = sum(deepThermostericSls, 3, 'omitnan'); %#ok<NASGU> - actually saved through VARS variable
        deepHalostericSls = (cdata.densityClim(:, :, isDeep) ./ data.haloDensity(:, :, isDeep) - 1) .* ...
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
            callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'steric sea level data'));
    end

end
