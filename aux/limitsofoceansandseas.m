%% LIMITSOFOCEANSANDSEAS 
% Returns the boundary of oceans and seas
%
% Data source
%   The ocean boundaries are based on IHO's 'Limits of oceans and seas':
%       International Hydrographic Organization & Sieger, R. (2012).
%       doi: 10.1594/PANGAEA.777975
%
% Last modified by
%   2024/11/20, williameclee@arizona.edu (@williameclee)

function LimitsOfOceansAndSeas = limitsofoceansandseas(varargin)
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    p = inputParser;
    addParameter(p, 'BeQuiet', false, ...
        @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    beQuiet = logical(p.Results.BeQuiet);

    matFileName = ...
        fullfile(getenv('COASTS'), 'Limits_of_oceans_and_seas.mat');

    if isfile(matFileName)
        load(matFileName, 'LimitsOfOceansAndSeas')

        if ~beQuiet
            fprintf('%s loaded %s\n', upper(mfilename), matFileName)
        end

        return
    end

    rawFileName = ...
        fullfile(getenv('COASTS'), 'Limits_of_oceans_and_seas.tab');

    if ~isfile(rawFileName)
        error('The file %s does not exist.', rawFileName)
    end

    rawFile = ...
        fopen(rawFileName, 'r');
    oceanLimitsData = textscan(rawFile, '%s %f32 %f32 %u %u', ...
        'HeaderLines', 16, 'Delimiter', '\t');
    fclose(rawFile);

    oceanNames = oceanLimitsData{1};
    lat = oceanLimitsData{2};
    lon = oceanLimitsData{3};

    % Correct or precision issues
    lon(abs(lon) <= 2e-3) = 0;
    lon(abs(lon - 180) <= 2e-3) = 180;
    lon(abs(lon + 180) <= 2e-3) = -180;

    clear oceanLimitsData filePath rawFile

    LimitsOfOceansAndSeas = struct('Name', unique(oceanNames), 'XY', []);

    for iOcean = 1:length(LimitsOfOceansAndSeas)
        oceanName = LimitsOfOceansAndSeas(iOcean).Name;
        LimitsOfOceansAndSeas(iOcean).XY = ...
            [lon(strcmp(oceanNames, oceanName)), ...
             lat(strcmp(oceanNames, oceanName))];

        % Make sure the shape of the Arctic Ocean is properly closed
        switch oceanName
            case 'Arctic Ocean, eastern part'
                LimitsOfOceansAndSeas(iOcean).XY = ...
                    [LimitsOfOceansAndSeas(iOcean).XY(1:end - 1, :); ...
                     [180, 90]; ...
                     LimitsOfOceansAndSeas(iOcean).XY(end, :)];
            case 'Arctic Ocean, western part'
                LimitsOfOceansAndSeas(iOcean).XY = ...
                    [LimitsOfOceansAndSeas(iOcean).XY(1:end - 2, :); ...
                     [-180, 90]; ...
                     LimitsOfOceansAndSeas(iOcean).XY(end - 1:end, :)];
        end

        LimitsOfOceansAndSeas(iOcean).poly = polyshape( ...
            LimitsOfOceansAndSeas(iOcean).XY(:, 1), ...
            LimitsOfOceansAndSeas(iOcean).XY(:, 2));
    end

    save(matFileName, 'LimitsOfOceansAndSeas')

    if ~beQuiet
        fprintf('%s saved %s\n', upper(mfilename), matFileName)
    end

end
