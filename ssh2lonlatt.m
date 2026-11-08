%% SSH2LONLATT
% Loads sea surface height (SSH) anomaly data and interpolate to regular
% latitude-longitude mesh.
%
% Syntax
%   [ssh, sigma, dates] = SSH2LONLATT(product)
%   [ssh, sigma, dates] = SSH2LONLATT(product, timestep, meshsize)
%   [ssh, sigma, dates] = SSH2LONLATT(product, timestep, meshsize, timelim)
%   [ssh, sigma, dates] = SSH2LONLATT(product, timestep, meshsize, timelim)
%   [ssh, sigma, dates] = SSH2LONLATT(__, "Name", Value)
%   [ssh, sigma, dates, lon, lat, partialMonths] = SSH2LONLATT(__)
%
% Input arguments
%	product - Name of the SSH product to load
%		- 'MEaSUREs': The MEaSUREs Gridded Sea Surface Height Anomalies
%           Version 2205 product.
%       - 'NASASSH' and 'NASASSHMonthly': Native NASA-SSH maps and monthly means.
%       Other products: 'ORAS5', 'DT2024', 'CMEMS'.
%       The default product is 'MEaSUREs'.
%   timestep - Temporal interpolation time step
%       - Numeric or duration scalar, in units of days.
%       - 'midmonth': mean of native samples in each calendar month,
%         dated at the midpoint.
%       - 'GRACE': mean within CSR RL06 degree-60 GRACE solution windows.
%         Missing months use calendar boundaries inside the gap, with outer
%         bounds aligned to adjacent real solutions. Dates are window midpoints.
%         Monthly inputs are averaged by native timestamps, not reconstructed
%         at daily resolution.
%       When not specified, no temporal interpolation is performed.
%   meshsize - Spatial interpolation grid size in degrees
%       - Real scalar, in units of degrees.
%       When not specified, no spatial interpolation is performed.
%   timelim - Time range of interest
%       - Datetime or numeric vector, in units of datenum
%       When not specified, the full time range of the data is returned.
%   LonOrigin - Longitudinal centre of the interpolated output mesh
%       - Real scalar, in units of degrees.
%       When not specified, the midpoint of the input longitude range is
%       used. This option is only relevant when spatial interpolation is
%       performed.
%   Interpolation - Interpolation method, for both temporal and spatial
%       - An INTERP1 method, or 'mean' to average native samples.
%         For 'midmonth' and 'GRACE', this option controls spatial
%         interpolation only. Fixed-step means use
%         half-open windows centered on the output dates. Spatial
%         interpolation uses 'linear' when 'mean' is requested.
%       The default interpolation method is 'linear'. This option is only
%       relevant when temporal or spatial interpolation is performed.
%   OutputFormat - Format of the output mesh
%       - 'meshgrid' or 'ndgrid' for (lat, lon, t) or (lon, lat, t)
%           ordering.
%       The default format is 'meshgrid'.
%   TimeFormat - Format of the output time vector
%       - 'datetime' or 'datenum'.
%       The default format is 'datetime'.
%   Unit - Unit of the output steric sea level data
%       - 'm' (metres) or 'mm' (millimetres, equivalent to kg/m^2).
%       The default unit is 'm'.
%	ForceNew - Logical flag to force reprocess of the data
%		The default option is false.
%	SaveData - Logical flag to cache the interpolated grid over the full data time range
%		The default option is true.
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is 'soft true', only printing important
%       messages.
%
% Output arguments
%	ssh - SSH mesh
%		The dimensions and units depend on the options specified.
%       The default dimensions are (lat, lon, t) in units of metres.
%	sigma - Uncertainty or error mesh of the SSH data
%		The units and dimensions are the same as for SSH.
%	dates - The time stamps of the data
%		The data type depends on the TimeFormat option.
%	lon - Longitude of the mesh
%		The dimensions (i.e. row or column vector) depend on the
%       OutputFormat option.
%	lat - Latitude of the mesh
%		The dimensions (i.e. row or column vector) depend on the
%       OutputFormat option.
%
%   partialMonths - Logical column aligned with dates for native NASASSHMonthly
%       maps. Empty for other products or interpolated grids and epochs.
%
% Data sources
%	MEaSUREs Gridded Sea Surface Height Anomalies Version 2205
% 		https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205
%
% Created by
%	2025/05/19, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%	2026/11/07, En-Chi Lee (williameclee@gmail.com)

function [ssh, sshSigma, dates, lon, lat, partialMonths] = ssh2lonlatt(product, timestep, meshsize, timelim, options)

    arguments (Input)
        product SSHProduct {mustBeScalarOrEmpty} = 'MEaSUREs'
        timestep {mustBeTimeStep} = []
        meshsize = []
        timelim {mustBeTimeRange} = []
        options.LonOrigin {mustBeFinite, mustBeReal} = []
        options.Interpolation {mustBeTextScalar} = 'linear'
        options.TimeFormat (1, :) DateFormat = 'datetime'
        options.OutputFormat (1, :) MeshFormat = 'meshgrid'
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, {'mm', 'm'})} = 'm'
        options.ForceNew (1, 1) {mustBeNumericOrLogical} = false
        options.BeQuiet (1, 1) {mustBeNumericOrLogical} = 0.5
        options.SaveData (1, 1) {mustBeNumericOrLogical} = true
        options.CallChain cell = {}
    end

    arguments (Output)
        ssh (:, :, :) {mustBeReal}
        sshSigma (:, :, :) {mustBeReal}
        dates (:, 1) {mustBeVector}
        lon (:, :) {mustBeReal, mustBeVector}
        lat (:, :) {mustBeReal, mustBeVector}
        partialMonths (:, 1) logical
    end

    lonOrigin = options.LonOrigin;
    intpMthd = options.Interpolation;
    timeFmt = char(options.TimeFormat);
    outputFmt = char(options.OutputFormat);
    unit = lower(options.Unit);
    forceNew = logical(options.ForceNew);
    beQuiet = uint8(double(options.BeQuiet) * 2);
    saveData = logical(options.SaveData);
    callChain = [options.CallChain, {mfilename}];

    if numel(meshsize) == 2 && isempty(lonOrigin)
        lonOrigin = meshsize(2);
        meshsize = meshsize(1);
    elseif numel(meshsize) > 2
        error('ULMO:InvalidInput', ...
            'meshsize must be a scalar or a 2-element vector of (h, lonOrigin), got %d elements.', ...
            numel(meshsize));
    end

    if meshsize < 0
        error('ULMO:InvalidMeshSize', ...
            'Mesh size must be a positive real scalar, got %f.', meshsize);
    elseif isnan(meshsize)
        warning('ULMO:InvalidMeshSize', ...
        'Mesh size is NaN, argument ignored.');
        meshsize = [];
    end

    % Preserve the legacy path builder's method spelling.
    [timestep, ~, useMean] = normalisetimeconfig(timestep, intpMthd);

    if ~isempty(timelim) && isnumeric(timelim)
        timelim = datetime(timelim, 'ConvertFrom', 'datenum');
    end

    if ~isempty(timelim) && timelim(1) > timelim(2)
        warning('The start time (%s) is later than the end time (%s), flipping the order.', timelim(1), timelim(2));
        timelim = fliplr(timelim);
    end

    %% Checking for existing data
    spatialIntpMthd = intpMthd;

    if strcmpi(spatialIntpMthd, 'mean')
        spatialIntpMthd = 'linear';
    end

    vars = {'sshs', 'dates', 'lon', 'lat'};
    nativeMonthly = product == SSHProduct.NASASSHMonthly && isempty(timestep) && ...
        isempty(meshsize) && isempty(lonOrigin);
    partialMonths = false(0, 1);
    inPath = outputpath(product);
    if nativeMonthly
        if nargout >= 6
            vars{end+1} = 'partialMonths';
        end
    else
        intpOutPath = outputpath(product, intpMthd, timestep, meshsize, lonOrigin);
        intpDataExists = ~forceNew && isfile(intpOutPath) && ...
            all(ismember(vars, who('-file', intpOutPath)));
        if intpDataExists
            inPath = intpOutPath;
        end
    end

    if ~exist(inPath, 'file')
        error('ULMO:InputDataNotFound', ...
            'Input data not found at %s.', inPath);
    elseif ~all(ismember(vars, who('-file', inPath)))
        notFoundVars = vars(~ismember(vars, who('-file', inPath)));
        error('ULMO:InputDataIncomplete', ...
            'Input data at %s is missing variables: %s.', ...
            inPath, strjoin(notFoundVars, ', '));
    end

    if ismember('sshErrors', who('-file', inPath))
        data = load(inPath, vars{:}, 'sshErrors');
    else
        data = load(inPath, vars{:});
        data.sshErrors = [];
    end

    if beQuiet <= 1
        fprintf('[ULMO>%s] Loaded %s.\n', callchaintext(callChain), filehref(inPath, 'SSH product'));
    end

    if ~nativeMonthly && ~intpDataExists
        nativeDates = data.dates;
        if isnumeric(nativeDates), nativeDates = datetime(nativeDates, 'ConvertFrom', 'datenum'); end

        if useMean
            data = preparegrid(data, meshsize, lonOrigin, spatialIntpMthd, beQuiet, callChain);

            % The cached epochs span the product. Apply timelim only to the output.
            fields = struct('height', data.sshs);
            if ~isempty(data.sshErrors), fields.errorSquared = data.sshErrors.^2; end
            [fields, data.dates] = processfieldtime(fields, nativeDates, timestep, intpMthd, []);
            data.sshs = fields.height;
            if ~isempty(data.sshErrors)
                % Preserve RMS uncertainty; these are not independent errors.
                data.sshErrors = sqrt(fields.errorSquared);
            end
        else
            if ~isempty(timestep)
                fields = struct('height', data.sshs);
                if ~isempty(data.sshErrors), fields.error = data.sshErrors; end
                [fields, data.dates] = processfieldtime(fields, nativeDates, timestep, intpMthd, []);
                data.sshs = fields.height;
                if ~isempty(data.sshErrors), data.sshErrors = fields.error; end
            end

            data = preparegrid(data, meshsize, lonOrigin, spatialIntpMthd, beQuiet, callChain);
        end

        if saveData
            save(intpOutPath, '-struct', 'data', 'sshs', 'sshErrors', 'dates', 'lon', 'lat', '-v7.3');

            if beQuiet <= 1
                fprintf('[ULMO>%s] Saved %s.\n', callchaintext(callChain), filehref(intpOutPath, 'SSH product'));
            end
        end
    end

    if nativeMonthly && nargout >= 6
        partialMonths = logical(data.partialMonths(:));
        assert(numel(partialMonths) == numel(data.dates), 'ULMO:ssh2lonlatt:PartialMonths', ...
            'Expected one partial-month flag per native date.');
        if ~isempty(timelim)
            partialMonths = partialMonths(data.dates >= timelim(1) & data.dates <= timelim(2));
        end
    end

    if nargout == 0
        timeFmt = 'datetime';
    end

    [ssh, sshSigma, dates, lon, lat] = formatoutput( ...
        data.sshs, data.sshErrors, data.dates, data.lon, data.lat, ...
        timelim, timeFmt, outputFmt, unit);

    if nargout > 0
        return
    end

    plotsealeveltseries(dates, ssh, lon, lat, product, unit, 'Global mean sea surface height');

end

%% Subfunctions
function data = preparegrid(data, meshsize, lonOrigin, intpMthd, beQuiet, callChain)

    if ~isempty(meshsize) || ~isempty(lonOrigin)

        if isempty(meshsize)
            meshsize = 1/6;
        end

        if isempty(lonOrigin)
            lonOrigin = 180;
        end

        nativeLon = data.lon;
        nativeLat = data.lat;
        [data.sshs, data.lon, data.lat] = interpspatial( ...
            data.lon, data.lat, data.sshs, meshsize, lonOrigin, intpMthd, beQuiet, callChain);

        if ~isempty(data.sshErrors)
            data.sshErrors = interpspatial( ...
                nativeLon, nativeLat, data.sshErrors, meshsize, lonOrigin, intpMthd, beQuiet, callChain);
        end

    end

end

% Interpolation
function [meshIntp, lonIntp, latIntp] = interpspatial( ...
        lon, lat, mesh, meshSize, lonOrigin, intpMthd, beQuiet, callChain)

    lonIntp = (-180:meshSize:180) + lonOrigin;
    latIntp = -90:meshSize:90;
    [lonnIntp, lattIntp] = meshgrid(lonIntp, latIntp);

    if isequal(lonIntp(:), lon(:)) && isequal(latIntp(:), lat(:))
        meshIntp = mesh;
        return
    end

    if beQuiet == 0
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO>%s] Interpolating spatially, %s\n', ...
            callchaintext(callChain), templine);
    end

    % lonnIntp(lonnIntp < ogLonOrigin - 180) = lonnIntp(lonnIntp < ogLonOrigin - 180) + 360;
    % lonnIntp(lonnIntp > ogLonOrigin + 180) = lonnIntp(lonnIntp > ogLonOrigin + 180) - 360;

    lon = lon - 360 * floor(lon(1) / 360);

    if lon(1) ~= 0

        if lon(end) - lon(1) == 360
            lon = lon(1:end - 1);
            mesh = mesh(1:end - 1, :, :);
        end

        [lon, idx] = sort(mod(lon, 360));
        mesh = mesh(idx, :, :);
    end

    if lon(end) - lon(1) ~= 360
        lon = [lon(end) - 360; lon; lon(1) + 360];
        mesh = cat(1, mesh(end, :, :), mesh, mesh(1, :, :));
    end

    [lonn, latt] = meshgrid(lon, lat);
    mesh = permute(mesh, [2, 1, 3]);

    meshIntp = nan([length(latIntp), length(lonIntp), size(mesh, 3)], 'single');

    for iDate = 1:size(meshIntp, 3)
        meshIntp(:, :, iDate) = ...
            interp2(lonn, latt, squeeze(mesh(:, :, iDate)), ...
            mod(lonnIntp, 360), lattIntp, intpMthd);
    end

    meshIntp = permute(meshIntp, [2, 1, 3]);

    if beQuiet == 0
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end

% Format the output
function [sshs, sshErrors, dates, lon, lat] = ...
        formatoutput(sshs, sshErrors, dates, lon, lat, timelim, timeFmt, outputFmt, unit)

    if ~isempty(timelim)
        isValidTime = ...
            (dates >= timelim(1)) & (dates <= timelim(2));
        sshs = sshs(:, :, isValidTime);

        if ~isempty(sshErrors)
            sshErrors = sshErrors(:, :, isValidTime);
        end

        dates = dates(isValidTime);
    end

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

    switch outputFmt
        case 'meshgrid'
            lon = lon(:)';
            lat = lat(:);
            sshs = permute(sshs, [2, 1, 3]);

            if ~isempty(sshErrors)
                sshErrors = permute(sshErrors, [2, 1, 3]);
            end

        case 'ndgrid'
            lon = lon(:);
            lat = lat(:)';
    end

    switch unit
        case 'm'
        case 'mm'
            sshs = sshs * 1e3;

            if ~isempty(sshErrors)
                sshErrors = sshErrors * 1e3;
            end

    end

end

% Get output path for interpolated data
function outPath = outputpath( ...
        product, intpMthd, timeStep, meshSize, lonOrigin)

    arguments (Input)
        product (1, :) char
        intpMthd (1, :) char = []
        timeStep {mustBeTimeStep} = []
        meshSize {mustBePositive} = []
        lonOrigin {mustBeFinite, mustBeReal} = []
    end

    arguments (Output)
        outPath (1, :) char
    end

    root = getenv('IFILES');
    assert(isfolder(root), 'ULMO:ssh2lonlatt:MissingRoot', 'Set IFILES to the data directory.');
    outFolder = fullfile(root, 'SSH');

    if isempty(intpMthd) && isempty(timeStep) && isempty(meshSize) && isempty(lonOrigin)
        outPath = fullfile(outFolder, sprintf('%s.mat', product));
        return
    end

    timeStepStr = '';

    if ~isempty(timeStep)

        if isstring(timeStep) || ischar(timeStep)
            timeStepStr = sprintf('-T%s', lower(timeStep));
        elseif isduration(timeStep)
            timeStepStr = sprintf('-T%s', erase(sprintf('%s', timeStep), ' '));
        elseif isnumeric(timeStep)
            timeStepStr = sprintf('-T%d', timeStep);
        else
            error('ULMO:ssh2lonlatt:InvalidTimeStep', ...
                'Unrecognised time step input for saving: %s', class(timeStep));
        end

        % Earlier linear fixed-step caches could contain bin averages.
        if isduration(timeStep) && strcmpi(intpMthd, 'linear')
            timeStepStr = [timeStepStr, '-interp'];
        elseif strcmpi(timeStep, 'midmonth') && ~strcmpi(intpMthd, 'mean')
            timeStepStr = [timeStepStr, '-mean'];
        end

    end

    spatialIntpStr = '';

    if ~isempty(meshSize) && ~isempty(lonOrigin)
        spatialIntpStr = ['-H', num2str(meshSize), 'O', num2str(lonOrigin)];
    elseif ~isempty(meshSize)
        spatialIntpStr = ['-H', num2str(meshSize)];
    elseif ~isempty(lonOrigin)
        spatialIntpStr = ['-O', num2str(lonOrigin)];
    end

    outFile = sprintf('%s%s%s-%s.mat', ...
        product, timeStepStr, spatialIntpStr, intpMthd);
    outPath = fullfile(outFolder, outFile);
end
