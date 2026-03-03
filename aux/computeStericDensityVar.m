%% COMPUTESTERICDENSITYVAR - Computes density from temperature and salinity, with one fixed as climatology
%
% Author
%   2026/03/02, williameclee@arizona.edu (@williameclee)

function computeStericDensityVar(dataPath, climatologyPath, options)
	%% Validation and checks
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

	% The needed variables or variables to be computed
    inputVars = {'salinity', 'consTemp', 'lat', 'pres'};
    inputClimVars = {'consTempClim', 'salinityClim'};
    outputVars = {'haloDensity', 'thermoDensity'};

    if ~options.ForceNew && all(ismember(outputVars, who('-file', dataPath))) && ...
            (dir(dataPath).datenum > dir(climatologyPath).datenum)
        % Check if density variables already exist and the file is younger than the climatology
        % If so, no need to recompute
        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exists.\n', ...
                callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), ...
                filehref(dataPath, 'thermo- and halosteric density data'));
        end

        return
    elseif any(~ismember(inputVars, who('-file', dataPath)))
        % Make sure all required time-dependent input variables exist in the data file
        missingInputVars = setdiff(inputVars, who('-file', dataPath));
        error('Input file %s is missing required variables: %s', ...
            dataPath, strjoin(missingInputVars, ', '));
    elseif any(~ismember(inputClimVars, who('-file', climatologyPath)))
        % Make sure all required climatology variables exist in the climatology file
        missingClimVars = setdiff(inputClimVars, who('-file', climatologyPath));
        error('Climatology file %s is missing required variables: %s', ...
            climatologyPath, strjoin(missingClimVars, ', '));
    end

	% Load variables from .mat files
    data = load(dataPath, inputVars{:});
    cdata = load(climatologyPath, inputClimVars{:});

    % Compute halosteric density
    haloDensity = nan(size(data.salinity), 'single');

    for iPres = 1:size(data.pres, 2)
        haloDensity(:, :, iPres) = gsw_rho( ...
            squeeze(data.salinity(:, :, iPres)), squeeze(cdata.consTempClim(:, :, iPres)), ...
            data.pres(1, iPres));
    end

    % Compute thermosteric density
    thermoDensity = nan(size(data.salinity), 'single');

    for iPres = 1:size(data.pres, 2)
        thermoDensity(:, :, iPres) = gsw_rho( ...
            squeeze(cdata.salinityClim(:, :, iPres)), squeeze(data.consTemp(:, :, iPres)), ...
            data.pres(1, iPres));
    end

    % Save density data back to the same .mat file
    save(dataPath, outputVars{:}, '-append');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', callchaintext(callChain), ...
            datetime(ddata.date, "Format", 'yyyy/MM'), ...
            filehref(dataPath, 'thermo- and halosteric density data'));
    end

end
