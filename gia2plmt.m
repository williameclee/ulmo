%% GIA2PLMT
% Reads a GIA model and converts it to the spheircal harmonic format.
%
% Syntax
%   Plmt = gia2plmt(model)
%		Returns the GIA geoid change in one year.
%   Plmt = gia2plmt(days, model)
%		Returns the GIA geoid change in the given number of days.
%   Plmt = gia2plmt(time, model)
%		Returns the GIA geoid changes at the given times.
%   Plmt = gia2plmt(time, model, L)
%		Returns the GIA geoid changes truncated to degree L.
%   Plmt = gia2plmt(__, 'Name', Value)
%   [Plmt, PlmtU, PlmtL] = gia2plmt(__)
%		Also returns the upper and lower bounds of the GIA model, if
%       available.
%
% Input arguments
%   model - Name of the GIA model
%       - 'Steffen_ice6g_vm5a': A model computed by H. Steffen using the
%           ice6g ice model and vm5a viscosity profile. Other models from
%           this dataset are also available and use the original naming
%           scheme. For example, 'Steffen_anu-ice_i72'.
%           This family of models can also be specified as a cell array,
%           e.g. {'Steffen', 'ice6g', 'vm5a'}.
%       - 'LM17.3': A model based on the data from the LM17.3 dataset.
%       - 'Paulson07': A model based on the ICE-5G ice load model of
%           Peltier (2004). Suitable for both Antarctica and Greenland. As
%           corrected by Geruo A and J. Wahr.
%           Please avoid using this model for oceans.
%       - 'Wangetal08': A model based on the older ICE-4G ice model, and
%           viscosity which varies laterally. Suitable for Greenland.
%       - 'IJ05_R2': A model based on the Ivins et al. (2013) ice model.
%           Suitable for Antarctica.
%       - 'IJ05': A model based on the Ivins and James (2005) ice model.
%           Suitable for Antarctica.
%       - 'W12a_v1': A 'best' model from Whitehouse et al. (2012). Suitable
%           only for Antarctica.
%       The input can also be the path to the model file.
%		The default model is 'Steffen_ice6g_vm5a'.
%   days - Number of days to calculate the GIA change for
%   time - Vector of dates to calculate the GIA change for
%		The input can be in DATENUM or DATETIME format.
%   L - Maximum degree of the GIA model
%		If empty, the model is not truncated.
%   OutputField - Field to output
%       - 'massdensity': Surface mass density
%       - 'geoid': Geoid change
%       The default field is 'massdensity'.
%   OutputFormat - Format of the output
%       - 'timefirst': The first dimension is time
%       - 'traditional': The first dimension is degree
%		The default format is 'timefirst'.
%   BeQuiet - Whether to surpress output messages
%		The default option is false.
%
% Output arguments
%   Plmt - GIA change in spherical harmonic format
%   PlmtU, PlmtL - Upper and lower bounds of the GIA change, if available
%
% See also
%   CORRECT4GIA
%
% Data source
%   Paulson, A., S. Zhong, and J. Wahr. (2007). Inference of mantle
%       viscosity from GRACE and relative sea level data, Geophys. J. Int.
%       171, 497–508. doi: 10.1111/j.1365-246X.2007.03556.x
%   Geruo, A., Wahr, J. & Zhong, S. (2013). Computations of the
%       viscoelastic response of a 3-D compressible Earth to surface
%       loading: An application to Glacial Isostatic Adjustment in
%       Antarctica and Canada. Geophys. J. Int. 192, 557-572.
%   Ivins, E. R., T. S. James, J. Wahr, E. J. O. Schrama, F. W. Landerer,
%       and K. M. Simon. (2013). Antarctic contribution to sea level rise
%       observed by GRACE with improved GIA correction, Journal of
%       Geophysical Research: Solid Earth, vol. 118, 3126-3141, doi:
%       10.1002/jgrb.50208
%   Ivins, E. R., T. S. James. (2005). Antarctic glacial isostatic
%       adjustment: a new assessment, Antarctic Science, vol. 17(4),
%       541-553,  doi: 10.1017/S0954102005002968
%   Wang, H., P. Wu, and W. van der Wal. (2018). Using postglacial sea
%       level, crustal velocities and gravity-rate-of-change to constrain
%       the influence of thermal effects on mantle lateral heterogeneities,
%       Journal of Geodynamics 46, 104-117. doi: 10.1016/j.jog.2008.03.003
%   Whitehouse, P. L., Bentley, M. J., Milne, G. A., King, M. A., Thomas,
%       I. D. (2012). A new glacial isostatic adjustment model for
%       Antarctica: calibrated and tested using observations of relative
%       sea-level change and present-day uplift rates. Geophysical Journal
%       International 190, 1464-1482. doi:10.1111/j.1365-246X.2012.05557.x
%   Steffen, H. (2021). Surface Deformations from Glacial Isostatic
%       Adjustment Models with Laterally Homogeneous, Compressible Earth
%       Structure (1.0) [dataset]. Zenodo. doi: 10.5281/zenodo.5560862.
%   Steffen, H., Li, T., Wu, P., Gowan, E. J., Ivins, E., Lecavalier, B.,
%       Tarasov, L., Whitehouse, P. L. (2021). LM17.3 - a global vertical
%       land motion model of glacial isostatic adjustment [dataset].
%       PANGAEA, doi: 10.1594/PANGAEA.932462
%
% Last modified by
%   2025/11/16, williameclee@arizona.edu (@williameclee)

function varargout = gia2plmt(varargin)
    %% Initialisation
    [time, model, L, unit, outputFmt, beQuiet, callChain] = ...
        parseinputs(varargin{:});

    %% Loading the model
    % Load this data (saved as lmcosiM)

    inputPath = finddatafile(model, unit);
    inputUnit = unit;

    if exist(inputPath, 'file')

        if beQuiet
            warning('off', 'MATLAB:load:variableNotFound');
        end

        data = load(inputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');

        if ~beQuiet
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s GIA model</a>\n', ...
                callchaintext(callChain), inputPath, inputPath, upper(model));
        end

    else

        switch unit
            case 'SD'
                inputUnit = 'POT';
            case 'POT'
                inputUnit = 'SD';
        end

        altInputPath = finddatafile(model, inputUnit);

        if ~exist(altInputPath, 'file')
            error('ULMO:LoadData:FileNotFound', ...
                'GIA model %s not found at expected location:\n%s or %s', upper(model), inputPath, altInputPath);
        end

        data = load(altInputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');

        if ~beQuiet
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s GIA model</a>\n', ...
                callchaintext(callChain), inputPath, inputPath, upper(model));
        end

    end

    if strcmp(unit, 'SD')
        data.lmcosiM(1, 3) = 0; % Ensure the mean is zero
    end

    hasBounds = isfield(data, 'lmcosiU') && isfield(data, 'lmcosiL');

    %% Some additional processing
    % Convert unit
    if ~strcmp(unit, inputUnit)
        data.lmcosiM = convertgravity(data.lmcosiM, inputUnit, unit);

        if hasBounds
            data.lmcosiU = convertgravity(data.lmcosiU, inputUnit, unit);
            data.lmcosiL = convertgravity(data.lmcosiL, inputUnit, unit);
        end

    end

    % Convert to SSH (i.e. find spatially constant term)
    if strcmp(unit, 'POT')

        if abs(data.lmcosiM(1, 3)) < 1e-7
            data.lmcosiM(1, 3) = 0;
        end

        try
            Linput = max(data.lmcosiM(:, 1));
            vlmLmcosi = giaz2plmt(model, Linput);
            vlmLmcosi(:, 3:4) = vlmLmcosi(:, 3:4) * 1e-3; % mm/yr -> m/yr
            rslLmcosi = data.lmcosiM;
            rslLmcosi(:, 3:4) = data.lmcosiM(:, 3:4) - vlmLmcosi(1:size(data.lmcosiM, 1), 3:4);
            oceanDomain = GeoDomain('alloceans', "Buffer", 0.5);
            rslLmcosi = localise(rslLmcosi, oceanDomain, 60, "BeQuiet", true);
            data.lmcosiM(1, 3) = data.lmcosiM(1, 3) - rslLmcosi(1, 3);
        catch ME

            switch model
                case 'Paulson07'
                    data.lmcosiM(1, 3) = -1.4e-4;
                    warning('ULMO:LoadData:FunctionExternalFailuse', ...
                        '%s failed to load GIAZ2PLMT for model %s, cannot convert POT to SSH, genetic value of %.5f is used\n(Error: %s)', ...
                        upper(mfilename), upper(model), data.lmcosiM(1, 3), ME.message);
                case 'A13'
                    data.lmcosiM(1, 3) = -1.3e-4;
                    warning('ULMO:LoadData:FunctionExternalFailuse', ...
                        '%s failed to load GIAZ2PLMT for model %s, cannot convert POT to SSH, genetic value of %.5f is used\n(Error: %s)', ...
                        upper(mfilename), upper(model), data.lmcosiM(1, 3), ME.message);
                otherwise
                    warning('ULMO:LoadData:FunctionExternalFailuse', ...
                        '%s failed to load GIAZ2PLMT for model %s, cannot convert POT to SSH\n(Error: %s)', ...
                        upper(mfilename), upper(model), ME.message);
            end

        end

    end

    % Save result to inputpath if previously not found
    if ~exist(inputPath, 'file')

        if ~exist(fileparts(inputPath), 'dir')
            mkdir(fileparts(inputPath));
        end

        save(inputPath, '-struct', 'data', 'lmcosiM');

        if hasBounds
            save(inputPath, '-struct', 'data', 'lmcosiU', 'lmcosiL', '-append');
        end

        if ~beQuiet
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s GIA model</a>\n', ...
                callchaintext(callChain), inputPath, inputPath, upper(model));
        end

    end

    % Truncate the model to the desired degree
    if ~isempty(L)

        if size(data.lmcosiM, 1) < addmup(L)

            if ~beQuiet
                warning('ULMO:LoadData:InsufficientResolution', ...
                    'Model %s resolution lower than the requested degree %d', model, L);
            end

            [data.lmcosiM(1:addmup(L), 2), data.lmcosiM(1:addmup(L), 1)] = addmon(L);
        else
            data.lmcosiM = data.lmcosiM(1:addmup(L), :);
        end

        if hasBounds

            if size(data.lmcosiU, 1) < addmup(L) || size(data.lmcosiL, 1) < addmup(L)
                [data.lmcosiU(1:addmup(L), 2), data.lmcosiU(1:addmup(L), 1)] = addmon(L);
                [data.lmcosiL(1:addmup(L), 2), data.lmcosiL(1:addmup(L), 1)] = addmon(L);
            else
                data.lmcosiU = data.lmcosiU(1:addmup(L), :);
                data.lmcosiL = data.lmcosiL(1:addmup(L), :);
            end

        end

    end

    % Time
    if isempty(time) || isscalar(time)
        % If time is scalar, interpret it as the day change
        if isempty(time)
            deltaYear = 1;
        else
            deltaYear = time / days(years(1));
        end

        GIAt = data.lmcosiM;
        GIAt(:, 3:4) = deltaYear * GIAt(:, 3:4);

        if hasBounds
            GIAtU = data.lmcosiU;
            GIAtU(:, 3:4) = deltaYear * GIAtU(:, 3:4);
            GIAtL = data.lmcosiL;
            GIAtL(:, 3:4) = deltaYear * GIAtL(:, 3:4);
        end

    else
        % Otherwise, interpret it as a vector of dates (in datenum format)
        % Reference the date string to the first date
        deltaYear = (time - time(1)) / days(years(1));

        GIAt = plm2plmt(data.lmcosiM, deltaYear);

        if hasBounds
            GIAtU = plm2plmt(data.lmcosiU, deltaYear);
            GIAtL = plm2plmt(data.lmcosiL, deltaYear);
        end

    end

    %% Collecting outputs
    switch outputFmt
        case 'timefirst'
            % Do nothing
        case 'traditional'
            GIAt = permute(GIAt, [2, 3, 1]);
    end

    if ~hasBounds
        GIAtU = [];
        GIAtL = [];

        if nargout > 1 && ~beQuiet
            warning(sprintf('ULMO:%s:NoBoundsToReturn', upper(mfilename)), ...
                'Upper and lower bounds are not available for model %s', upper(model));
        end

    else

        switch outputFmt
            case 'timefirst'
                % Do nothing
            case 'traditional'
                GIAtU = permute(GIAtU, [2, 3, 1]);
                GIAtL = permute(GIAtL, [2, 3, 1]);
        end

    end

    varargout = {GIAt, GIAtU, GIAtL};

    if nargout > 0
        return
    end

    %% Plotting
    plotgiamap(GIAt, time, deltaYear, model, unit)
end

%% Subfunctions
function varargout = parseinputs(varargin)
    dfOpts.model = 'Caron18';

    % Allow skipping the time argument
    if nargin > 0 && ...
            (ischar(varargin{1}) || isstring(varargin{1}) || iscell(varargin{1}))
        varargin(2:end + 1) = varargin;
        varargin{1} = [];
    end

    ip = inputParser;
    addOptional(ip, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isempty(x) || isduration(x));
    addOptional(ip, 'Model', dfOpts.model, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(ip, 'L', [], ...
        @(x) (isnumeric(x) && isscalar(x) && x > 0) || isempty(x));
    addParameter(ip, 'BeQuiet', false, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'Unit', 'SD', ...
        @(x) ischar(validatestring(x, {'massdensity', 'geoid', 'SD', 'POT', 'geo'})));
    addParameter(ip, 'OutputFormat', 'timefirst', ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'CallChain', {}, ...
        @(x) iscell(x));

    parse(ip, varargin{:});
    time = ip.Results.Time;
    model = conddefval(ip.Results.Model, dfOpts.model);
    L = ip.Results.L;
    beQuiet = ip.Results.BeQuiet;
    unit = ip.Results.Unit;
    outputFmt = ip.Results.OutputFormat;
    callChain = [ip.Results.CallChain, {mfilename}];

    if isdatetime(time)
        time = datenum(time); %#ok<DATNM>
    elseif isduration(time)
        time = days(time);
    end

    if iscell(model)

        if ~strcmpi(model{1}, 'Steffen')
            error('Unrecognised model name %s', upper(model{1}));
        elseif ~ismember(model{2}, {'anu-ice', 'ice6g', 'ice7g'})
            error('Unrecognised ice model %s', upper(model{2}));
        end

        model = sprintf('%s_%s_%s', model{1}, model{2}, model{3});

    end

    if ~isempty(regexp(lower(model), '^ice[-]?6g[_-]?d$', 'once'))
        model = 'ICE6GD';
    elseif startsWith(lower(model), {'paulson', 'caron'})
        model(1) = upper(model(1));
    end

    if L ~= round(L)
        warning( ...
            sprintf('ULMO:%s:InvalidInput:NonIntegerDegree', upper(mfilename)), ...
            'The degree L must be an integer. Rounding to %d', round(L));
        L = round(L);
    end

    switch lower(unit)
        case {'massdensity', 'sd'}
            unit = 'SD';
        case {'geoid', 'pot', 'geo'}
            unit = 'POT';
    end

    varargout = {time, model, L, unit, outputFmt, beQuiet, callChain};

end

function inputPath = finddatafile(model, unit)

    if isfile(model)
        inputPath = model;
        return
    end

    if ~isempty(getenv('GIA'))
        inputFolder = getenv('GIA');
    elseif ~isempty(getenv('IFILES'))
        inputFolder = fullfile(getenv('IFILES'), 'GIA');
    else
        error('GIA folder not found')
    end

    if strncmp(model, 'Morrow', 6)
        inputFolder = fullfile(inputFolder, model(1:6));
    elseif strncmpi(model, 'Steffen', 7)
        inputFolder = fullfile(inputFolder, 'Steffen21');
    elseif strcmp(model, 'LM17.3')
        inputFolder = fullfile(inputFolder, 'LM17.3');

        if exist(inputFolder, 'dir') ~= 7
            mkdir(inputFolder)
        end

    else
        inputFolder = fullfile(inputFolder, model);
    end

    % And the appropriate name
    inputPath = fullfile(inputFolder, sprintf('%s_%s.mat', model, unit));
end

function plmt = plm2plmt(plm, deltaYear)
    plmt = zeros([length(deltaYear), size(plm)]);

    plmt(:, :, 1:2) = repmat(reshape( ...
        plm(:, 1:2), [1, length(plm), 2]), [length(deltaYear), 1, 1]);
    plmt(:, :, 3:4) = ...
        deltaYear(:) .* reshape(plm(:, 3:4), [1, length(plm), 2]);
    plmt = squeeze(plmt);

end

function plotgiamap(GIAt, time, deltaYear, model, outputField)
    % Get the change rate
    if isempty(deltaYear)
        deltaYear = 1;
    else
        deltaYear = deltaYear(end);
    end

    if isempty(time) || isscalar(time)
        GIAchange = GIAt;
    else
        GIAchange = squeeze(GIAt(end, :, :) - GIAt(1, :, :));
        GIAchange(:, 1:2) = GIAt(1, :, 1:2);
    end

    [GIAmesh, lon, lat] = plm2xyz(GIAchange, "BeQuiet", true);

    % Get the coastlines
    coastLonlat = gshhscoastline('c', 'LonOrigin', 180, "BeQuiet", true);

    [cLim, cStep] = optimalclim(GIAmesh, 'Percentile', 1);
    GIAmesh = max(min(GIAmesh, cLim(2)), cLim(1));

    figure(999)
    % Protect underscore in model name
    model = strrep(model, '_', '\_');
    set(gcf, "NumberTitle", "off", "Name", ...
        sprintf('GIA change in %.1f year(s) (%s)', deltaYear, upper(mfilename)))
    clf

    title(sprintf('Model: %s', model))

    switch outputField
        case {'massdensity', 'SD'}
            cLabel = 'Surface mass density [kg/m^2]';
        case {'geoid', 'POT'}
            cLabel = 'Geoid rate [m/s]';
    end

    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", cLabel, ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, GIAmesh, cLevels, "LineStyle", 'none');
    plotqdm(coastLonlat, 'k');
    hold off
end
