%% PARSEGRACESOURCEFILE - Read spherical harmonic coefficients from a GRACE source file.
%
% Syntax
%   [gravSph, gravStdSph, epochCentres, epochRange, gravParam, eqRad] = ...
%       parsegracesourcefile(dataPath)
%
% Input arguments
%   dataPath - Path to a GSM, GAC, or GAD text file with an RL05 SHM or RL06 YAML header.
%
% Output arguments
%   gravSph - Gravity field spherical harmonics.
%       Layout is in lmcosi order.
%   gravStdSph - Gravity field stanard deviation spherical harmonics.
%   epochCentres - Midpoint of each epoch.
%   epochRange - Start and end times of the epoch.
%       Size: (1, 2).
%   gravParam - Earth's gravitational parameter from the header.
%       Unit: m^3/s^2.
%   eqRad - Mean equatorial radius from the header.
%       Unit: m.
%
% Notes
%   If the file starts at degree 2, the function inserts C00 = 1 and zeros
%   for the remaining degree 0 and 1 coefficients and their uncertainties.
%   Requested reference parameters must have finite positive values.
%
% See also
%   GRACE2PLMT (GRACE2PLMT_NEW), AOD1B2PLMT
%
% Created by
%   2025/05/20, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/01, En-Chi Lee (williameclee@gmail.com)

function [gravSph, gravStdSph, epochCentres, epochRange, gravParam, eqRad] = ...
        parsegracesourcefile(dataPath)

    arguments (Input)
        dataPath
    end

    arguments (Output)
        gravSph (:, 4) double
        gravStdSph (:, 4) double
        epochCentres (1, 1) datetime
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
    dataStr = fileread(dataPath);
    dataStart = regexp(dataStr, '(?m)^GRCOF2[ \t]', 'once');

    if isempty(dataStart)
        error('ULMO:parsegracesourcefile:InvalidFormat', ...
            'No GRCOF2 records found in %s', dataPath);
    end

    data = dataStr(dataStart:end);
    data = textscan(data, '%s%f%f%f%f%f%f%s%s%s');
    gravSph = [data{2}, data{3}, data{4}, data{5}];

    % Sort the data
    [gravSph, rowOrder] = sortrows(gravSph, [1, 2], 'ascend');

    % Pad the degree 0 and 1 data (for JPL)
    padLowDegrees = gravSph(1, 1) == 2;

    if padLowDegrees
        gravSph = ...
            [[0, 0, 1, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravSph];
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
    epochCentres = mean(epochRange);

    if nargout <= 4
        return
    end

    %% Fetching the parameters
    header = dataStr(1:dataStart - 1);
    % RL05
    % RL05 stores GM and radius, in that order, on one EARTH record
    earth = regexp(header, '(?m)^EARTH[ \t]+([^\r\n]+)', 'tokens', 'once');

    if ~isempty(earth)
        params = sscanf(regexprep(earth{1}, '[dD]', 'E'), '%f');

        if numel(params) < 2 || any(~isfinite(params(1:2))) || any(params(1:2) <= 0)
            error('ULMO:parsegracesourcefile:InvalidHeader', ...
                'Invalid EARTH parameters in %s', dataPath);
        end

        gravParam = params(1);
        eqRad = params(2);
        return
    end

    % RL06
    % Keep the leading whitespace. In RL06 YAML it identifies which
    % attribute owns each value.
    lines = regexp(header, '\r?\n', 'split');

    gravParam = extractheadervalue(lines, 'earth_gravity_param');

    if nargout == 6
        eqRad = extractheadervalue(lines, 'mean_equator_radius');
    end

end

%% Subfunctions
% Extract non-standard attributes from the header
function param = extractheadervalue(lines, paramName)
    % Match the complete attribute key and capture its indentation in match{1}
    key = ['^([ \t]*)', paramName, '[ \t]*:'];

    for iKey = 1:numel(lines)
        match = regexp(lines{iKey}, key, 'tokens', 'once');

        if isempty(match)
            continue
        end

        for iLine = iKey + 1:numel(lines)
            line = lines{iLine};

            if isempty(strtrim(line)) || startsWith(strtrim(line), '#')
                continue
            end

            indent = regexp(line, '^[ \t]*', 'match', 'once');

            % The next field at this level or above ends the attribute.
            % Otherwise a missing GM could inherit the following radius value.
            if numel(indent) <= numel(match{1})
                break
            end

            value = regexp(line, '^\s*value\s*:\s*(\S+)', 'tokens', 'once');

            if ~isempty(value)
                % Convert Fortran D exponents before checking the physical constant
                param = str2double(regexprep(value{1}, '[dD]', 'E'));

                if isreal(param) && isfinite(param) && param > 0
                    return
                end

                break
            end

        end

        break
    end

    error('ULMO:parsegracesourcefile:InvalidHeader', ...
        'Missing or invalid %s in GRACE header', paramName);
end
