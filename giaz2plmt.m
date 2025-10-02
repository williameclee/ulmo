%% GIAZ2PLMT
% Reads the vertical displacement provided by a GIA model and converts it
% to the spheircal harmonic format in mm.
%
% Syntax
%   plm = gia2plmt(model, L)
%		Returns the GIA vertical displacement in one year, truncated to
%       degree L.
%   plm = gia2plmt(years, __)
%		Returns the GIA vertical displacement in the given number of years.
%   plmt = gia2plmt(time, __)
%		Returns the GIA vertical displacements at the given dates.
%   plmt = gia2plmt(__, 'Name', Value)
%
% Input arguments
%   model - Name of the GIA model
%       - Name of a model computed by H. Steffen. It can be specified as,
%           e.g. 'Steffen_ice6g_vm5a' or {'Steffen', 'ice6g', 'vm5a'}.
%       - LM17.3 is also supported.
%       Other models are specified in the same way.
%       - The input can also be the path to the model file.
%		The default model is 'Steffen_ice6g_vm5a'.
%       Format: string or 1 x 3 cell array.
%   years - Number of years to calculate the GIA change for
%		Format: numeric scalar.
%   time - Vector of dates to calculate the GIA change for
%		Format: numeric (datenum), datetime, or duration.
%   L - Maximum degree of the soherical harmonic expansion
%		The default option is 18.
%       Format: numeric scalar.
%   BeQuiet - Whether to surpress output messages
%		The default option is true.
%       Format: logical scalar.
%   MakePlot - Whether to plot the GIA change
%		The default option is true if no output is requested.
%       Format: logical scalar.
%
% Output arguments
%   plm - GIA change in the given duration in spherical harmonic format
%       Format: addmup(L) x 4 double array.
%   plmt - GIA change at the given dates in spherical harmonic format
%       Format: addmup(L) x 4 x ndates double array.
%
% See also
%   GIA2PLMT, GIAZ2SLEPT
%
% Data source
%   Steffen, H. (2021). Surface Deformations from Glacial Isostatic
%       Adjustment Models with Laterally Homogeneous, Compressible Earth
%       Structure (1.0) [Dataset]. Zenodo. doi: 10.5281/zenodo.5560862.
%   Steffen, H., Li, T., Wu, P., Gowan, E. J., Ivins, E., Lecavalier, B.,
%       Tarasov, L., Whitehouse, P. L. (2021). LM17.3 - a global vertical
%       land motion model of glacial isostatic adjustment [dataset].
%       PANGAEA, doi: 10.1594/PANGAEA.932462
%
% Last modified by
%   2025/08/03, williameclee@arizona.edu (@williameclee)

function varargout = giaz2plmt(varargin)
    %% Initialisation
    % Parse the inputs
    [model, L, dYear, beQuiet, makePlot] = parseinputs(varargin{:});

    % Loading the model
    if contains(lower(model), 'steffen')
        wPlmt = findsteffendata(model);
        wUPlmt = [];
        wLPlmt = [];
    elseif strcmpi(model, 'lm17.3')
        wPlmt = lm17_vup;
        wPlmt = wPlmt(wPlmt(:, 1) <= L, :);
        wUPlmt = [];
        wLPlmt = [];
    elseif contains(lower(model), 'caron')
        % Load the Caron model
        inputFolder = fullfile(getenv('IFILES'), 'GIA', capitalise(model));
        inputPath = fullfile(inputFolder, [capitalise(model), '_VLM.mat']);
        load(inputPath, 'lmcosiM', 'lmcosiU', 'lmcosiL');
        wPlmt = lmcosiM(lmcosiM(:, 1) <= L, :);
        wUPlmt = lmcosiU(lmcosiU(:, 1) <= L, :);
        wLPlmt = lmcosiL(lmcosiL(:, 1) <= L, :);
    elseif strcmpi(model, 'ice6gd') || strcmpi(model, 'ice-6g_d')
        inputFolder = fullfile(getenv('IFILES'), 'GIA', 'ICE-6G_D');
        inputPath = fullfile(inputFolder, 'ICE-6G_D_VLM.mat');
        load(inputPath, 'wSph');
        wPlmt = wSph;
        wUPlmt = [];
        wLPlmt = [];
    else
        error('ULMO:LoadData:FileNotFound', 'Unrecognised model name %s', upper(model));
    end

    if (size(wPlmt, 1) < addmup(L) || ...
            (~isempty(wUPlmt) && size(wUPlmt, 1) < addmup(L)) || ...
            (~isempty(wLPlmt) && size(wLPlmt, 1) < addmup(L))) && ~beQuiet
        warning('Model %s resolution lower than the requested degree %d', model, L);
    end

    if size(wPlmt, 1) < addmup(L)
        [order, degree] = addmon(L);
        wPlmt(1:addmup(L), 1) = degree;
        wPlmt(1:addmup(L), 2) = order;
    end

    if ~isempty(wUPlmt) && size(wUPlmt, 1) < addmup(L)
        [order, degree] = addmon(L);
        wUPlmt(1:addmup(L), 1) = degree;
        wUPlmt(1:addmup(L), 2) = order;
    end

    if ~isempty(wLPlmt) && size(wLPlmt, 1) < addmup(L)
        [order, degree] = addmon(L);
        wLPlmt(1:addmup(L), 1) = degree;
        wLPlmt(1:addmup(L), 2) = order;
    end

    %% Converting the model
    if isscalar(dYear)
        wPlmt(:, 3:4) = dYear .* wPlmt(:, 3:4);
    else
        wPlmt = plm2plmt(wPlmt, dYear);
    end

    if ~isempty(wUPlmt) && ~isempty(wLPlmt)

        if isscalar(dYear)
            wUPlmt(:, 3:4) = dYear .* wUPlmt(:, 3:4);
            wLPlmt(:, 3:4) = dYear .* wLPlmt(:, 3:4);
        else
            wUPlmt = plm2plmt(wUPlmt, dYear);
            wLPlmt = plm2plmt(wLPlmt, dYear);
        end

    end

    %% Collecting and displaying outputs
    varargout = {wPlmt, wUPlmt, wLPlmt};

    if nargout > 0 || ~makePlot
        return
    end

    plotdispmap(wPlmt, dYear, model)
end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Allow skipping the time argument
    if nargin > 0 && ...
            (ischar(varargin{1}) || isstring(varargin{1}) || iscell(varargin{1}))
        varargin(2:end + 1) = varargin;
        varargin{1} = [];
    end

    % Parse the inputs
    modelD = 'Steffen_ice6g_vm5a';
    LD = 18;
    p = inputParser;
    addOptional(p, 'Time', [], ...
        @(x) isnumeric(x) || isdatetime(x) || isduration(x) || isempty(x));
    addOptional(p, 'Model', modelD, ...
        @(x) ischar(x) || (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(p, 'L', LD, ...
        @(x) isnumeric(x) || isempty(x));
    % Quiet by default because this function itself does not display anything
    addParameter(p, 'BeQuiet', true, @(x) istruefalse(x));
    addParameter(p, 'MakePlot', true, @(x) istruefalse(x));

    parse(p, varargin{:});
    time = p.Results.Time;
    model = conddefval(p.Results.Model, modelD);
    L = conddefval(p.Results.L, LD);
    beQuiet = logical(p.Results.BeQuiet);
    makePlot = logical(p.Results.MakePlot);

    % Format the time into change in years
    if isempty(time)
        dYear = 1;
    elseif isscalar(time) && isnumeric(time)
        dYear = time;
    elseif isvector(time) && isnumeric(time)
        dYear = (time - time(1)) / days(years(1));
    elseif isvector(time) && isdatetime(time)
        dYear = years(time - time(1));
    elseif isduration(time)
        dYear = years(time);
    end

    % Format the model name
    if iscell(model)

        if ~strcmpi(model{1}, 'Steffen')
            error('Unrecognised model name %s, only STEFFEN is supported at the moment', ...
                upper(model{1}));
        elseif ~ismember(model{2}, {'anu-ice', 'ice6g', 'ice7g'})
            error('Unrecognised ice model %s', upper(model{2}));
        end

        % Make sure the model name is capitalised
        model{1} = capitalise(model{1});
        % Assemble the model name
        model = sprintf('%s_%s_%s', model{1}, model{2}, model{3});

    end

    varargout = {model, L, dYear, beQuiet, makePlot};

end

function plm = findsteffendata(model)
    %% Finding the model data
    if isfile(model)
        % If the model is a path, load it directly
        inputPath = model;
    else
        % Otherwise, find the model in the GIA folder
        % Find if the GIA folder exists
        if ~isempty(getenv('GIA'))
            inputFolder = getenv('GIA');
        elseif ~isempty(getenv('IFILES'))
            inputFolder = fullfile(getenv('IFILES'), 'GIA');
        else
            error( ...
                ['GIA folder not found', newline, ...
                 'It should be kept at $IFILES/GIA/ or $GIA/', newline, ...
             'Use SETENV to specify the environment path']);
        end

        inputFolder = fullfile(inputFolder, 'Steffen21');

        % Get the appropriate file name
        model = replace(model, 'Steffen_', '');
        model = replace(model, 'steffen_', '');
        model = [model, '_vup.mat'];
        inputFile = model;

        inputPath = fullfile(inputFolder, inputFile);

        % Make sure the file exists
        if exist(inputPath, 'file') ~= 2
            error('ULMO:LoadData:FileNotFound', ...
                'Model %s not found at %s', upper(model), inputPath);
        end

    end

    %% Loading the model
    load(inputPath, 'z', 'plm');

    if exist('plm', 'var')
        return
    else
        plm = xyz2plm_new(flip(z'), 96, "BeQuiet", true);
        save(inputPath, 'plm', '-append');
    end

end

function plmt = plm2plmt(plm, multplicationFactor)
    % Put the time dimension in the third dimension
    tShape = [1, 1, length(multplicationFactor)];
    multplicationFactor = reshape(multplicationFactor, tShape);

    % Preallocate the output
    plmt = zeros([size(plm), length(multplicationFactor)]);
    % Copy the degree and order
    plmt(:, 1:2, :) = repmat(plm(:, 1:2), tShape);
    % Multiply the coefficients by the year change
    plmt(:, 3:4, :) = multplicationFactor .* plm(:, 3:4);
end

function plotdispmap(wPlmt, dYear, model)
    % Get the change rate
    if ~isscalar(dYear)
        wPlmt = squeeze(wPlmt(:, :, end));
        dYear = dYear(end);
    end

    meshSize = 1;
    [wXY, lon, lat] = plm2xyz(wPlmt, meshSize, "BeQuiet", true);

    % Get the coastlines
    coastLonlat = gshhscoastline('c', 'LonOrigin', 180, "BeQuiet", true);

    % Get the colour scheme and limits
    [cLim, cStep] = optimalclim(wXY, 'Percentile', 1);
    wXY = max(min(wXY, cLim(2)), cLim(1));

    % Protect underscore in model name
    model = strrep(model, '_', '\_');

    figTitle = sprintf('GIA vertical displacement in %s year(s)', num2str(dYear));

    figure(999)
    clf
    set(gcf, "NumberTitle", "off", "Name", ...
        sprintf('%s (%s)', figTitle, upper(mfilename)))
    title([figTitle, newline, sprintf('Model: %s', model)])

    [~, cLevels] = loadcbar(cLim, cStep, ...
        "Title", 'Vertical displacement [mm]', ...
        "Colormap", 'temperature anomaly');

    hold on
    contourf(lon, lat, wXY, cLevels, "LineStyle", 'none');
    plotqdm(coastLonlat, 'k');
    hold off
end
