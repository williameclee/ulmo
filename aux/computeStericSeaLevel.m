%% COMPUTESTERICSEALEVEL - Computes sea level anomalies from density and climatology
%
% Last modified
%   2026/09/25, williameclee@arizona.edu (@williameclee)

function computeStericSeaLevel(dataPath, climatologyPath, options)
    %% Validation and checks
    arguments (Input)
        dataPath {mustBeTextScalar, mustBeFile}
        climatologyPath {mustBeTextScalar, mustBeFile}
        options.Bottom double {mustBePositive, mustBeFinite} = 6000
        options.HasDeepLayer (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    bottom = options.Bottom;
    hasDeepLayer = options.HasDeepLayer;
    callChain = [options.CallChain, {mfilename}];

    ddata = load(dataPath, 'date');
    % Check if date variable exists and is in datetime format
    if ~isfield(ddata, 'date') || ~isa(ddata.date, 'datetime')
        error('Input file %s is missing required variable "date" in datetime format.', dataPath);
    end

    % The needed variables
    inputVars = {'density', 'haloDensity', 'thermoDensity', 'depth'};
    inputClimVars = {'densityClim'};
    outputVars = ...
        {'stericSl', 'thermostericSl', 'halostericSl'};
    outputDepthVars = ...
        {'shallowStericSl', 'deepStericSl', ...
         'shallowThermostericSl', 'deepThermostericSl', ...
         'shallowHalostericSl', 'deepHalostericSl'};

    if ~options.ForceNew && ...
            all(ismember(outputVars, who('-file', dataPath))) && ...
            (~hasDeepLayer || all(ismember(outputDepthVars, who('-file', dataPath)))) && ...
            (dir(dataPath).datenum > dir(climatologyPath).datenum)
        % Check if steric sea level variables already exist and the file is younger than the climatology
        % If so, no need to recompute
        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exist and is newer than climatology.\n', ...
                callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'steric sea level data'));
        end

        return
    elseif any(~ismember(inputVars, who('-file', dataPath)))
        % Make sure all required time-dependent input variables exist in the data file
        missingDataVars = setdiff(inputVars, who('-file', dataPath));
        error('Data file %s is missing required variables: %s', ...
            dataPath, strjoin(missingDataVars, ', '));
    elseif ~ismember('densityClim', who('-file', climatologyPath))
        % Make sure all required climatology variables exist in the climatology file
        error('Climatology file %s is missing density data', ...
            climatologyPath);
    end

    % Load variables from .mat files
    data = load(dataPath, inputVars{:});
    cdata = load(climatologyPath, inputClimVars{:});

    %% Main computation
    depth = double(data.depth);

    if isvector(depth)
        depth = repmat(depth(:)', size(data.density, 1), 1);
    end

    assert(size(depth, 1) == size(data.density, 1) && ...
        size(depth, 2) == size(data.density, 3), 'Depth grid does not match density.');
    assert(all(isfinite(depth), 'all') && all(depth >= 0, 'all') && ...
        all(diff(depth, 1, 2) > 0, 'all'), 'Depth must be finite and increasing.');

    if isscalar(bottom)
        bottom = repmat(bottom, size(depth, 1), 1);
    else
        bottom = bottom(:);
    end

    assert(numel(bottom) == size(depth, 1), 'Bottom must be scalar or one value per latitude.');

    % Calculate layer thicknesses
    layerTop = [zeros(size(depth, 1), 1), (depth(:, 1:end - 1) + depth(:, 2:end)) / 2];
    layerBottom = [layerTop(:, 2:end), bottom];
    assert(all(layerBottom >= layerTop, 'all'), 'Bottom lies above the final layer top.');
    layerThk = layerBottom - layerTop;
    shallowThk = max(0, min(layerBottom, 2000) - min(layerTop, 2000));
    deepThk = max(0, layerBottom - max(layerTop, 2000));

    % Integrate steric sea level
    denPtrbtn = cdata.densityClim ./ data.density - 1;
    thermoDenPtrbtn = cdata.densityClim ./ data.thermoDensity - 1;
    haloDenPtrbtn = cdata.densityClim ./ data.haloDensity - 1;

    stericSl = integrateLayers(denPtrbtn, layerThk); %#ok<NASGU>
    thermostericSl = integrateLayers(thermoDenPtrbtn, layerThk); %#ok<NASGU>
    halostericSl = integrateLayers(haloDenPtrbtn, layerThk); %#ok<NASGU>

    if hasDeepLayer
        shallowStericSl = integrateLayers(denPtrbtn, shallowThk); %#ok<NASGU>
        deepStericSl = integrateLayers(denPtrbtn, deepThk); %#ok<NASGU>
        shallowThermostericSl = integrateLayers(thermoDenPtrbtn, shallowThk); %#ok<NASGU>
        deepThermostericSl = integrateLayers(thermoDenPtrbtn, deepThk); %#ok<NASGU>
        shallowHalostericSl = integrateLayers(haloDenPtrbtn, shallowThk); %#ok<NASGU>
        deepHalostericSl = integrateLayers(haloDenPtrbtn, deepThk); %#ok<NASGU>
    end

    try
        save(dataPath, outputVars{:}, '-append');
    catch
        save(dataPath, outputVars{:}, '-v7.3', '-append');
    end

    if hasDeepLayer

        try
            save(dataPath, outputDepthVars{:}, '-append');
        catch
            save(dataPath, outputDepthVars{:}, '-v7.3', '-append');
        end

    end

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', ...
            callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), ...
            filehref(dataPath, 'steric sea level data'));
    end

end

%% Subfunctions
function intgVal = integrateLayers(val, layerThk)
    weights = reshape(layerThk, size(layerThk, 1), 1, []);
    valid = isfinite(val) & weights > 0;
    terms = val .* weights;
    terms(~valid) = 0;
    intgVal = sum(terms, 3);
    intgVal(~any(valid, 3)) = NaN;
end
