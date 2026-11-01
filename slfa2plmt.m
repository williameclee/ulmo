%% SLFA2PLMT
% Reads sea level fingerprints (SLF) SH coefficients
%
% Syntax
%   [time, plmLandload, plmRsl, plmGeoid, plmBedrock] = ...
%       slf2plmt(pcenter, rframe, rotation)
%   [__] = slf2plmt(__, L, timerange)
%   [__] = slf2plmt(__, 'Name', Value)
%
% Input arguments
%   pcenter - The GRACE data centre from which the SLF data is derived
%       Must be 'CSR', 'JPL', or 'GFZ'.
%       The default data centre is 'CSR'.
%   rframe - The reference frame of the SLF data
%       Must be 'CF' (centre-of-figure) or 'CM' (centre-of-mass).
%       The default reference frame is 'CF'.
%   rotation - Whether the rotation effect is included
%       Must be true or false.
%       The default is true.
%   L - The maximum degree of the SH coefficients to be returned
%       The default value is the maximum degree of the data files (60).
%   timerange - The time range of the SH coefficients to be returned
%       Must be a two-element DATETIME or DATENUM array.
%       The default time range is the entire time range of the data files.
%   ForceNew - Whether to force the recomputation of the SH coefficients
%       The default option is false.
%   BeQuiet - Whether to suppress the output
%       The default option is false.
%   SaveData - Whether to save the SH coefficients to a .mat file
%       The default option is true.
%
% Output arguments
%   time - The time stamps of the SH coefficients
%   plmLandload - The SH coefficients of the land load
%   plmRsl - The SH coefficients of the relative sea level
%   plmGeoid - The SH coefficients of the geoid
%   plmBedrock - The SH coefficients of the bedrock
%
% Data source
%   The SLF coefficients are from
%       Adhikari et al. (2019)
%       doi: 10.5194/essd-11-629-2019
%
% Last modified by
%   2024/08/22, williameclee@arizona.edu(@williameclee)

function varargout = slfa2plmt(varargin)
    %% Initialisation
    % Parse inputs
    [pcenter, rframe, rotation, timerange, L, ...
         forceNew, beQuiet, saveData] = parseinputs(varargin{:});

    % Get input and output files
    [time, inputPaths, outputPath] = getIOfiles(pcenter, rframe, rotation);

    %% Computing or loading data
    if isfile(outputPath) && ~forceNew
        load(outputPath, 'time', 'Ldata', ...
            'plmLandload', 'plmRsl', 'plmGeoid', 'plmSd', 'plmBedrock');

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), outputPath);
        end

    else
        [plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock, Ldata] = ...
            readdata(inputPaths, time);

        % Save data before removing data files outside the timerange
        if saveData
            save(outputPath, 'time', 'Ldata', 'plmLandload', ...
                'plmRsl', 'plmGeoid', 'plmSd', 'plmBedrock', '-v7.3');

            if ~beQuiet
                fprintf('%s saved %s\n', upper(mfilename), outputPath);
            end

        end

    end

    %% Post-processing
    % Remove data higher than L
    if isempty(L) || L == Ldata
        % Do not truncate
    elseif L < Ldata
        % Truncate if the requested degree is lower than the data files
        [plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock] = ...
            truncatedegree(L, plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock);
    elseif L > Ldata
        % Warn if the requested degree is higher than the data files
        warning('The maximum degree of the data files (%i) is lower than the requested degree (%i)', ...
            Ldata, L);
    end

    % Remove data files outside the timerange
    [time, plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock] = ...
        truncatetimerange(time, timerange, ...
        plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock);

    %% Collecting and displaying outputs
    varargout = {time, plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock};

    if nargout > 0
        return
    end

    varargout = {};
    plotslfmaps(plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock)

end

%% Subfunctions
function varargout = parseinputs(varargin)
    p = inputParser;
    addOptional(p, 'pcenter', 'CSR', ...
        @(x) (ischar(x) || isstring(x)));
    addOptional(p, 'rframe', 'CF', ...
        @(x) (ischar(x) || isstring(x)) && ismember(x, {'CF', 'CM'}));
    addOptional(p, 'rotation', true, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(p, 'L', [], ...
        @(x) isscalar(x) && isnumeric(x));
    addOptional(p, 'TimeRange', [], ...
        @(x) isdatetime(x) || isnumeric(x));
    addParameter(p, 'ForceNew', false, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeQuiet', false, ...
        @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'SaveData', true, ...
        @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    pcenter = upper(p.Results.pcenter);
    rframe = upper(p.Results.rframe);
    rotation = logical(p.Results.rotation);
    timerange = p.Results.TimeRange;
    L = p.Results.L;
    forceNew = logical(p.Results.ForceNew);
    beQuiet = logical(p.Results.BeQuiet);
    saveData = logical(p.Results.SaveData);

    if ~ismember(pcenter, {'CSR', 'JPL', 'GFZ'})
        error('Unknown pcenter: %s, must be CSR, JPL, or GFZ', pcenter);
    end

    if any(strcmpi(rframe, {'FIGURE', 'COF'}))
        rframe = 'CF';
    elseif any(strcmpi(rframe, {'MASS', 'COM'}))
        rframe = 'CM';
    elseif ~any(strcmpi(rframe, {'CF', 'CM'}))
        error('Unknown reference frame: %s, must be CF or CM', rframe);
    end

    if isnumeric(timerange)
        timerange = datetime(timerange, "ConvertFrom", 'datenum');
    end

    varargout = ...
        {pcenter, rframe, rotation, timerange, L, ...
         forceNew, beQuiet, saveData};

end

function varargout = getIOfiles(center, frame, rotation)
    %% Find folder names
    switch center
        case 'CSR'
            centerFolder = 'UTCSR';
        case 'JPL'
            centerFolder = 'JPLEM';
        case 'GFZ'
            centerFolder = 'GFZOP';
        otherwise
            error('Unknown pcenter: %s', center);
    end

    switch frame
        case 'CF'
            frameFolder = 'CF';
        case 'CM'
            frameFolder = 'CM';
        otherwise
            error('Unknown reference frame: %s', frame);
    end

    switch rotation
        case true
            rotFolder = 'WITHrotation';
        case false
            rotFolder = 'WOUTrotation';
        otherwise
            error('Unknown rotation: %s', rotation);
    end

    dataFolder = fullfile(getenv('IFILES'), ...
        'SLF', 'Adhikari2019', 'SLFsh_coefficients');

    if ~exist(dataFolder, 'dir')
        dataUrl = 'https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/8UC8IR/H9FTN4&version=3.0#';
        error('%s cannot find the data directory\nDownload the data from %s\and put it under %s', ...
            upper(mfilename), dataUrl, dataFolder);
    end

    %% Find data files
    inputFolder = fullfile( ...
        dataFolder, centerFolder, frameFolder, rotFolder);

    if ~exist(inputFolder, 'dir')
        error('%s cannot find the directory %s', ...
            upper(mfilename), inputFolder);
    end

    inputFiles = dir(inputFolder);
    inputFiles = ...
        {inputFiles(~[inputFiles.isdir] & ...
         contains({inputFiles.name}, 'SLF-2')).name};
    inputFiles = inputFiles(:);

    % Find timestamps
    timePath = fullfile(dataFolder, sprintf('time_%s.txt', center));
    time = year2date(readmatrix(timePath));

    if length(inputFiles) ~= length(time)
        error('The number of data files (%i) does not match the number of timestamps (%i)', ...
            length(inputFiles), length(time));
    end

    inputPaths = fullfile(inputFolder, inputFiles);

    % Find output file location
    outputFile = sprintf('SLF-%s-%s-%d.mat', ...
        center, frame, rotation);
    outputPath = fullfile(dataFolder, outputFile);

    varargout = {time, inputPaths, outputPath};
end

function varargout = readdata(inputPaths, time)
    Ldata = slf2plm(inputPaths{1});
    plmLandload = nan([addmup(double(Ldata)), 4, length(time)]);
    plmRsl = nan(size(plmLandload));
    plmGeoid = nan(size(plmLandload));
    plmSd = nan(size(plmLandload));
    plmBedrock = nan(size(plmLandload));

    parfor iTime = 1:length(time)
        [plmLandload(:, :, iTime), plmRsl(:, :, iTime), ...
             plmGeoid(:, :, iTime), plmBedrock(:, :, iTime)] = ...
            slf2plm(inputPaths{iTime});
        plmSd(:, :, iTime) = plm2pot(plmGeoid(:, :, iTime), [], [], [], 4);
    end

    varargout = {plmLandload, plmRsl, plmGeoid, plmSd, plmBedrock, Ldata};

end

function varargout = slf2plm(inputPath)
    inputId = fopen(inputPath, 'r');
    inputData = textscan(inputId, repmat('%f ', [1, 18]), ...
        "HeaderLines", 40);
    fclose(inputId);

    % If only one output is requested, return the degree of the input data
    if nargout == 1
        dof = length(inputData{1}(:));
        syms L y
        eqn = y == (L + 1) * (L + 2) / 2;
        eqn = subs(eqn, y, dof);
        L = max(uint16(solve(eqn, L)));

        varargout = {L};
        return
    end

    degree = inputData{1}(:);
    order = inputData{2}(:);
    degreeorder = [degree, order];
    [degreeorder, sortId] = sortrows(degreeorder, [1, 2]);

    plmLandload = [inputData{3}(:), inputData{5}(:)];
    plmLandload = [degreeorder, plmLandload(sortId, :)];

    plmRsl = [inputData{7}(:), inputData{9}(:)];
    plmRsl = [degreeorder, plmRsl(sortId, :)];

    plmGeoid = [inputData{11}(:), inputData{13}(:)];
    plmGeoid = [degreeorder, plmGeoid(sortId, :)];

    plmBedrock = [inputData{15}(:), inputData{17}(:)];
    plmBedrock = [degreeorder, plmBedrock(sortId, :)];

    varargout = {plmLandload, plmRsl, plmGeoid, plmBedrock};
end

function varargout = truncatedegree(L, varargin)
    truncation = addmup(L);
    varargin = varargin(:)';
    varargout = cell(size(varargin));

    for i = 1:length(varargin)
        Plm = varargin{i};

        if ismatrix(Plm)

            if size(Plm, 2) ~= 4
                error('The input data must have 4 columns');
            end

            Plm = Plm(1:truncation, :);
        elseif ndims(Plm) == 3

            if size(Plm, 2) == 4
                Plm = Plm(1:truncation, :, :);
            elseif size(Plm, 3) == 4
                Plm = Plm(:, 1:truncation, :);
            else
                error('The input data must have 4 columns');
            end

        end

        varargout{i} = Plm;
    end

end

function varargout = truncatetimerange(time, trange, varargin)
    maxTrange = [min(time), max(time)];

    if isempty(trange)
        % varargout
        if nargin == 2
            varargout = {time};
        else
            varargout = [{time}, varargin];
        end

        return
    elseif isscalar(trange)
    elseif all(isnat(trange))
        trange = maxTrange;
    elseif any(isnat(trange)) || ...
            trange(1) < maxTrange(1) || trange(2) > maxTrange(2)
        trange = ...
            [max(trange(1), maxTrange(1), "omitmissing"), ...
             min(trange(2), maxTrange(2), "omitmissing")];
        warning('The timerange is partially invalid or outside the range of the data files, truncated to match the data files:\n%s - %s', ...
            datetime(trange(1), "Format", 'yyyy/MM/dd'), ...
            datetime(trange(2), "Format", 'yyyy/MM/dd'));
    end

    if isscalar(trange)
        %#ok<*DATNM>
        timeId = find(min(datenum(time(:) - trange), [], 'omitnan'), 1);
    else
        timeId = time >= trange(1) & time <= trange(2);
    end

    time = time(timeId);

    if nargin == 2
        varargout = {time};
        return
    end

    varargin = varargin(:)';

    for i = 1:length(varargin)
        varargin{i} = squeeze(varargin{i}(:, :, timeId));
    end

    varargout = [{time}, varargin];
end

function plotslfmaps(plmLandload, plmRsl, plmGeoid, plmBedrock)
    meshSize = 1;
    [meshLandload, lon, lat] = ...
        plm2xyz(plmLandload, meshSize, "BeQuiet", true);
    meshRsl = plm2xyz(plmRsl, meshSize, "BeQuiet", true);
    meshGeoid = plm2xyz(plmGeoid, meshSize, "BeQuiet", true);
    meshBedrock = plm2xyz(plmBedrock, meshSize, "BeQuiet", true);

    figure(1)
    clf

    subplot(2, 2, 1)
    contourf(lon, lat, meshLandload, 100, 'LineStyle', 'none')
    colorbar
    title('Land load')

    subplot(2, 2, 2)
    contourf(lon, lat, meshRsl, 100, 'LineStyle', 'none')
    colorbar
    title('RSL')

    subplot(2, 2, 3)
    contourf(lon, lat, meshGeoid, 100, 'LineStyle', 'none')
    colorbar
    title('Geoid')

    subplot(2, 2, 4)
    contourf(lon, lat, meshBedrock, 100, 'LineStyle', 'none')
    colorbar
    title('Bedrock')
end
