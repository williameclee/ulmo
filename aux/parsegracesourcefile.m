%% PARSEGRACESOURCEFILE - Read spherical harmonic coefficients from a GRACE source file.
%
% Syntax
%   [gravSph, gravStdSph, epochMidDate, epochRange, gravParam, eqRad] = ...
%       parsegracesourcefile(dataPath)
%
% Input arguments
%   dataPath - Path to a GSM, GAC, or GAD text file with a YAML header.
%
% Output arguments
%   gravSph - Gravity field spherical harmonics.
%       Layout is in lmcosi order.
%   gravStdSph - Gravity field stanard deviation spherical harmonics.
%   epochMidDate - Midpoint of each epoch.
%   epochRange - Start and end times of the epoch.
%       Size: (1, 2).
%   gravParam - Earth's gravitational parameter from the header.
%       Unit: m^3/s^2.
%   eqRad - Mean equatorial radius from the header.
%       Unit: m.
%
% Notes
%   If the file starts at degree 2, the function inserts zero degrees 0 and 1.
%
% See also
%   GRACE2PLMT (GRACE2PLMT_NEW), AOD1B2PLMT
%
% Created by
%   2025/05/20, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/01, En-Chi Lee (williameclee@gmail.com)

function [gravSph, gravStdSph, epochMidDate, epochRange, gravParam, eqRad] = ...
        parsegracesourcefile(dataPath)

    arguments (Input)
        dataPath
    end

    arguments (Output)
        gravSph (:, 4) double
        gravStdSph (:, 4) double
        epochMidDate (1, 1) datetime
        epochRange (1, 2) datetime
        gravParam (1, 1) double
        eqRad (1, 1) double
    end

    %% Loading data
    if ~exist(dataPath, 'file')
        error('ULMO:parsegracesourcefile:FileNotFound', ...
            'File %s does not exist', dataPath)
    end

    % Extract the data from the file
    dataBarrier = '# End of YAML header\n';
    dataStr = fileread(dataPath);
    dataStr = strsplit(dataStr, dataBarrier);
    data = dataStr{2};
    data = textscan(data, '%s%f%f%f%f%f%f%s%s%s');
    gravSph = [data{2}, data{3}, data{4}, data{5}];

    % Sort the data
    [gravSph, rowOrder] = sortrows(gravSph, [1, 2], 'ascend');

    % Pad the degree 0 and 1 data (for JPL)
    padLowDegrees = gravSph(1, 1) == 2;

    if padLowDegrees
        gravSph = ...
            [[0, 0, 0, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravSph];
    end

    %% Quality check
    if size(gravSph, 1) ~= addmup(gravSph(end, 1))
        error('ULMO:parsegracesourcefile:DataSizeNotMatch', ...
            'Size of the gravity field (%d) does not match the expectation (%d) for degree %d', ...
            size(gravSph, 1), addmup(gravSph(end, 1)), gravSph(end, 1))
    end

    if nargout < 2
        return
    end

    gravStdSph = [data{2}, data{3}, data{6}, data{7}];
    gravStdSph = gravStdSph(rowOrder, :);

    if padLowDegrees
        gravStdSph = ...
            [[0, 0, 0, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravStdSph];
    end

    if nargout < 3
        return
    end

    %% Computing the dates
    startDate = datetime(data{8}{1}, ...
        "InputFormat", 'yyyyMMdd.HHmm', "Format", 'yyyy/MM/dd HH:mm');
    endDate = datetime(data{9}{1}, ...
        "InputFormat", 'yyyyMMdd.HHmm', "Format", 'yyyy/MM/dd HH:mm');
    epochRange = [startDate, endDate];
    epochMidDate = mean(epochRange);

    if nargout <= 4
        return
    end

    %% Fetching the parameters
    header = dataStr{1};
    lines = strtrim(strsplit(header, '\n'));

    gravParam = extractheadervalue(lines, 'earth_gravity_param');

    if nargout == 6
        eqRad = extractheadervalue(lines, 'mean_equator_radius');
    end

end

%% Subfunctions
% Extract non-standard attributes from the header
function param = extractheadervalue(header, paramName)

    if ischar(header)
        header = strtrim(strsplit(header, '\n'));
    end

    nameId = find(contains(header, paramName));

    for valueId = nameId + 1:length(header)

        if ~contains(header{valueId}, 'value')
            continue
        end

        param = strsplit(header{valueId}, ':');
        param = strtrim(param{end});
        param = str2double(param);
        break
    end

end
