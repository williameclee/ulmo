%% COMPUTESTERICDENSITYVAR - Computes density from temperature and salinity, with one fixed as climatology
%
% Author
%   2026/03/02, williameclee@arizona.edu (@williameclee)

function computeStericDensityVar(dataPath, climatologyPath, options)

    arguments (Input)
        dataPath {mustBeTextScalar, mustBeFile}
        climatologyPath {mustBeTextScalar, mustBeFile}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    ddata = load(dataPath, 'date');
    % Check if date variable exists and is in datetime format
    if ~isfield(ddata, 'date') || ~isa(ddata.date, 'datetime')
        error('Input file %s is missing required variable "date" in datetime format.', dataPath);
    end

    inputVars = {'salinity', 'consTemp', 'lat', 'pres'};
    inputClimVars = {'consTempClim', 'salinityClim'};
    outputVars = {'haloDensity', 'thermoDensity'};

    if all(ismember(outputVars, who('-file', dataPath))) && ~options.ForceNew

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exists.\n', ...
                callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'density variable data'));
        end

        return
    elseif any(~ismember(inputVars, who('-file', dataPath)))
        missingInputVars = setdiff(inputVars, who('-file', dataPath));
        error('Input file %s is missing required variables: %s', ...
            dataPath, strjoin(missingInputVars, ', '));
    elseif any(~ismember(inputClimVars, who('-file', climatologyPath)))
        missingClimVars = setdiff(inputClimVars, who('-file', climatologyPath));
        error('Climatology file %s is missing required variables: %s', ...
            climatologyPath, strjoin(missingClimVars, ', '));
    end

    load(dataPath, inputVars{:});
    load(climatologyPath, inputClimVars{:});

    % Compute halosteric density
    if ~ismember('haloDensity', who('-file', dataPath))
        haloDensity = nan(size(density), 'single');

        for iPres = 1:length(pres)
            haloDensity(:, :, iPres) = gsw_rho( ...
                squeeze(salinity(:, :, iPres)), squeeze(consTempClim(:, :, iPres)), pres(1, iPres)); %#ok<USENS>
        end

    end

    % Compute thermosteric density
    if ~ismember('thermoDensity', who('-file', dataPath))
        thermoDensity = nan(size(density), 'single');

        for iPres = 1:length(pres)
            thermoDensity(:, :, iPres) = gsw_rho( ...
                squeeze(salinityClim(:, :, iPres)), squeeze(consTemp(:, :, iPres)), pres(1, iPres)); %#ok<USENS>
        end

    end

    % Save density data back to the same .mat file
    save(dataPath, outputVars{:}, '-append');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', callchaintext(callChain), ...
            datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'density variable data'));
    end

end
