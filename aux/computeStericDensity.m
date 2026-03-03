%% COMPUTESTERICDENSITY - Computes density from temperature and salinity
%
% Last modified
%   2026/02/15, williameclee@arizona.edu (@williameclee)
%     - Extracted from PROCESSSTERICDATAEN4 for reusability

function computeStericDensity(dataPath, options)
    %% Validation and checks
    arguments (Input)
        dataPath {mustBeTextScalar, mustBeFile}
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

    % The needed variables
    inputVars = {'salinity', 'consTemp', 'lat', 'lon', 'pres', 'date'};

    if ~options.ForceNew && ismember('density', who('-file', dataPath))
        % Check if density variable already exists
        % If so, no need to recompute
        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exists.\n', ...
                callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'density data'));
        end

        return
    elseif any(~ismember(inputVars, who('-file', dataPath)))
        % Make sure all required input variables exist in the data file
        missingInputVars = setdiff(inputVars, who('-file', dataPath));
        error('Input file %s is missing required variables: %s', ...
            dataPath, strjoin(missingInputVars, ', '));
    end

    % Load variables from .mat files
    data = load(dataPath, inputVars{:});

    %% Main computation
    % Compute density
    density = nan(size(data.consTemp), 'single');

    for iDepth = 1:size(data.pres, 2)
        density(:, :, iDepth) = gsw_rho( ...
            squeeze(data.salinity(:, :, iDepth)), squeeze(data.consTemp(:, :, iDepth)), ...
            data.pres(1, iDepth));
    end

    % Save density data back to the same .mat file
    save(dataPath, 'density', '-append');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', callchaintext(callChain), ...
            datetime(data.date, "Format", 'yyyy/MM'), filehref(dataPath, 'density data'));
    end

end
