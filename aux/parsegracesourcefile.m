%% PARSEGRACESOURCEFILE
% Parses RL05 legacy SHM and RL06 YAML GRCOF2 files (GSM, GSU, GAC/GAD).
% Returns coefficients, uncertainties, dates, GM and radius. Reference
% parameters are read from EARTH or YAML attributes; malformed metadata errors.
%
% See also
%	GRACE2PLMT (GRACE2PLMT_NEW), AOD1B2PLMT
%
% Notes
%	This is a helper function with limited documentation.
%
% Authored by
%   2025/05/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/07/21, williameclee@arizona.edu (@williameclee)

function varargout = parsegracesourcefile(dataPath)
    %% Loading data
    if ~exist(dataPath, 'file')
        error(sprintf('%s:FileNotFound', upper(mfilename)), ...
            'File %s does not exist', dataPath)
    end

    % Both YAML and legacy SHM headers precede GRCOF2 records. Locating
    % the records avoids depending on a particular degree or comment line.
    dataStr = fileread(dataPath);
    dataStart = regexp(dataStr, '(?m)^GRCOF2\s', 'start', 'once');
    if isempty(dataStart)
        error('ULMO:parsegracesourcefile:InvalidFormat', ...
            'No GRCOF2 records found in %s', dataPath);
    end
    header = dataStr(1:dataStart - 1);
    data = dataStr(dataStart:end);
    data = textscan(data, '%s%f%f%f%f%f%f%s%s%s');
    gravitySph = [data{2}, data{3}, data{4}, data{5}];
    gravityStdSph = [data{2}, data{3}, data{6}, data{7}];

    % Sort the data
    gravitySph = sortrows(gravitySph, [1, 2], 'ascend');
    gravityStdSph = sortrows(gravityStdSph, [1, 2], 'ascend');

    % Pad the degree 0 and 1 data (for JPL)
    if gravitySph(1, 1) == 2
        gravitySph = ...
            [[0, 0, 1, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravitySph];
    end

    if gravityStdSph(1, 1) == 2
        gravityStdSph = ...
            [[0, 0, 0, 0; 1, 0, 0, 0; 1, 1, 0, 0]; gravityStdSph];
    end

    %% Quality check
    % Make sure the sizes of the matrices are correct
    if ~isequal(size(gravitySph), size(gravityStdSph))
        error(sprintf('%s:DataSizeNotMatch', upper(mfilename)), ...
            'Size of the gravity field (%d) and the uncertainty (%d) does not match. Please check the input file %s', ...
            size(gravitySph, 1), size(gravityStdSph, 1), dataPath)
    end

    if size(gravitySph, 1) ~= addmup(gravitySph(end, 1))
        error(sprintf('%s:DataSizeNotMatch', upper(mfilename)), ...
            'Size of the gravity field (%d) does not match the expectation (%d) for degree %d', ...
            size(gravitySph, 1), addmup(gravitySph(end, 1)), gravitySph(end, 1))
    end

    if nargout <= 2
        varargout = {gravitySph, gravityStdSph};
        return
    end

    %% Computing the dates
    startDate = datetime(data{8}{1}, ...
        "InputFormat", 'yyyyMMdd.HHmm', "Format", 'yyyy/MM/dd HH:mm');
    endDate = datetime(data{9}{1}, ...
        "InputFormat", 'yyyyMMdd.HHmm', "Format", 'yyyy/MM/dd HH:mm');
    meanDate = mean([startDate, endDate]);

    varargout = {gravitySph, gravityStdSph, meanDate, [startDate, endDate]};

    if nargout <= 4
        return
    end

    %% Fetching the parameters
    earth = regexp(header, '(?m)^EARTH[ \t]+([^\r\n]+)', 'tokens', 'once');
    if ~isempty(earth)
        parameters = sscanf(regexprep(earth{1}, '[dD]', 'E'), '%f');
        if numel(parameters) < 2 || any(~isfinite(parameters(1:2))) || any(parameters(1:2) <= 0)
            error('ULMO:parsegracesourcefile:InvalidHeader', ...
                'Invalid EARTH parameters in %s', dataPath);
        end
        gravityParam = parameters(1);
        equatorRadius = parameters(2);
    else
        gravityParam = extractheadervalue(header, 'earth_gravity_param');
        equatorRadius = extractheadervalue(header, 'mean_equator_radius');
    end

    varargout = ...
        {gravitySph, gravityStdSph, meanDate, [startDate, endDate], gravityParam, equatorRadius};
end

%% Subfunctions
% Extract non-standard attributes from the header
function param = extractheadervalue(header, paramName)
    % Restrict the value search to this YAML attribute's indented block.
    lines = regexp(header, '\r?\n', 'split');
    key = ['^([ \t]*)', paramName, '[ \t]*:'];
    for k = 1:numel(lines)
        match = regexp(lines{k}, key, 'tokens', 'once');
        if isempty(match)
            continue
        end
        indent = length(match{1});
        for j = k + 1:numel(lines)
            if isempty(strtrim(lines{j}))
                continue
            end
            whitespace = regexp(lines{j}, '^[ \t]*', 'match', 'once');
            if length(whitespace) <= indent
                break
            end
            value = regexp(lines{j}, '^\s*value\s*:\s*(\S+)', 'tokens', 'once');
            if ~isempty(value)
                param = str2double(regexprep(value{1}, '[dD]', 'E'));
                if isfinite(param) && param > 0
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
