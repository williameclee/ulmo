function computeStericDensity(dataPath, options)

    arguments (Input)
        dataPath (1, :) char {mustBeFile}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    % Load data from .mat file
    inputVars = {'salinity', 'consTemp', 'lat', 'lon', 'pres', 'date'};
    ddata = load(dataPath, 'date');

    if any(~ismember(inputVars, who('-file', dataPath)))
        missinginputVars = setdiff(inputVars, who('-file', dataPath));
        error('Input file %s is missing required variables: %s', ...
            dataPath, strjoin(missinginputVars, ', '));
    elseif ismember('density', who('-file', dataPath)) && ~options.ForceNew

        if ~options.BeQuiet
            cprintf('[ULMO>%s] Skipped computing %s %s, already exists.\n', ...
                callchaintext(callChain), datetime(ddata.date, "Format", 'yyyy/MM'), filehref(dataPath, 'density data'));
        end

        return
    end

    data = load(dataPath, inputVars{:});

    density = nan(size(data.consTemp), 'single');

    for iDepth = 1:size(data.pres, 2)
        density(:, :, iDepth) = gsw_rho( ...
            squeeze(data.salinity(:, :, iDepth)), squeeze(data.consTemp(:, :, iDepth)), data.pres(1, iDepth));
    end

    % Save density data back to the same .mat file
    save(dataPath, 'density', '-append');

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Computed %s %s.\n', callchaintext(callChain), ...
            datetime(data.date, "Format", 'yyyy/MM'), filehref(dataPath, 'density data'));
    end

end
