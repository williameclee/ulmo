%% GSHHSCOASTLINE Retrieves the GSHHS coastline data and formats it
%
% Last modified
%   2026/02/12, williameclee@arizona.edu (@williameclee)
%     - Added variable check before loading

function varargout = gshhscoastline(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    % Parse inputs
    [dataQuality, latlim, lonlim, minLandArea, ...
         upscale, buf, ~, lonOrigin, ...
         gshhsFileName, gshhsFileExists, ...
         forcenew, saveData, beQuiet] = parseinputs(varargin);

    %% Reteiving the original data
    vars = {'GshhsCoasts', 'gshhsCoastXY', 'gshhsCoastPoly'};

    if gshhsFileExists && ~forcenew && all(ismember(vars, who('-file', gshhsFileName)))
        load(gshhsFileName, vars{:})

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), gshhsFileName);
        end

    else
        GshhsCoasts = gshhsstruct('DataQuality', dataQuality, ...
            'Upscale', upscale, 'Buffer', buf, ...
            'Quiet', true, 'ForceNew', forcenew);
    end

    if ~beQuiet
        fprintf('%s generating the coasline polygon, this may take a while...\n', ...
            upper(mfilename))
    end

    segmentsToInclude = 1:length(GshhsCoasts);
    segmentsToInclude = ...
        segmentsToInclude(all(segmentsToInclude ~= [57, 77, 78]'));
    gshhsCoastPoly = union([GshhsCoasts(segmentsToInclude).Polygon]);

    gshhsCoastXY = poly2xy(gshhsCoastPoly); %#ok<NASGU>

    % Save the data if requested
    if saveData

        if gshhsFileExists
            save(gshhsFileName, vars{:}, '-append')

            if ~beQuiet
                fprintf('%s updated %s\n', upper(mfilename), gshhsFileName)
            end

        else
            save(gshhsFileName, vars{:})

            if ~beQuiet
                fprintf('%s saved %s\n', upper(mfilename), gshhsFileName)
            end

        end

    end

    %% Cropping the data to the desired limits
    % Leave only the land with a minimum area
    GshhsCoasts = GshhsCoasts([GshhsCoasts.Area] >= minLandArea);

    % Crop the data to the desired limits
    gshhsCoastPoly = croptolims(gshhsCoastPoly, ...
        latlim, lonlim, lonOrigin);

    % Make sure the coastline is closed
    gshhsCoastXY = closecoastline(gshhsCoastPoly.Vertices);

    %% Returning requested data
    varargout = {gshhsCoastXY, gshhsCoastPoly, GshhsCoasts};

    if nargout > 0
        return
    end

    figure(10)
    set(gcf, "Name", sprintf('Coastline (%s)', upper(mfilename)), ...
        "NumberTitle", "off")
    clf

    plotqdm(gshhsCoastXY, 'k.-')

end

%% Subfunctions
function varargout = parseinputs(inputArguments)
    %% Defining the input parser
    ip = inputParser;
    addOptional(ip, 'dataQuality', 'c', ...
        @(x) ischar(x) && ismember(x, 'cfhil'));
    addOptional(ip, 'LatLim', [-90, 90], ...
        @(x) isnumeric(x) && length(x) == 2);
    addOptional(ip, 'LonLim', [0, 360], ...
        @(x) isnumeric(x) && length(x) == 2);
    addOptional(ip, 'MinLandArea', 30, @isnumeric);
    addOptional(ip, 'Upscale', 0, @isnumeric);
    addOptional(ip, 'Buffer', 0, @isnumeric);
    addOptional(ip, 'Tolerence', 0.2, @isnumeric);
    addParameter(ip, 'LonOrigin', 180, @isnumeric);
    addParameter(ip, 'ForceNew', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'SaveData', true, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    parse(ip, inputArguments{:});

    %% Assigning the inputs
    dataQuality = ip.Results.dataQuality;
    latlim = ip.Results.LatLim;
    lonlim = ip.Results.LonLim;
    minLandArea = ip.Results.MinLandArea;
    upscale = ip.Results.Upscale;
    buf = ip.Results.Buffer;
    tol = ip.Results.Tolerence;
    lonOrigin = ip.Results.LonOrigin;
    forcenew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = logical(ip.Results.BeQuiet);

    %% Additional parameters
    [gshhsFileName, ~, gshhsFileExists] = gshhsfilename( ...
        'DataQuality', dataQuality, ...
        'Upscale', upscale, 'Buffer', buf, ...
        'MinLandArea', minLandArea, 'Tolerence', tol);

    %% Returning the inputs
    varargout = ...
        {dataQuality, latlim, lonlim, minLandArea, ...
         upscale, buf, tol, lonOrigin, ...
         gshhsFileName, gshhsFileExists, ...
         forcenew, saveData, beQuiet};
end

function p = croptolims(p, latlim, lonlim, lonOrigin)
    %% Shifting
    XY = poly2xy(p);
    [Y, X] = ...
        flatearthpoly(XY(:, 2), XY(:, 1), lonOrigin);
    X = X - floor(min(X(:)) / 360) * 360;
    p = polyshape([X, Y]);

    %% Cropping
    bbox = polyshape( ...
        [lonlim(1), lonlim(2), lonlim(2), lonlim(1), lonlim(1)], ...
        [latlim(1), latlim(1), latlim(2), latlim(2), latlim(1)]);
    p = intersect(p, bbox);
end
