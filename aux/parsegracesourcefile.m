%% PARSEGRACESOURCEFILE
% Parses a GRACE source file and returns the gravity field and its uncertainty. Should work for both GSM and GAC/GAD files.
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

    % Extract the data from the file
    dataBarriers = ...
        ["# End of YAML header", ...
         "     90   90,", ...
         "SHM      90   90 1.00 fully normalized exclusive permanent tide", ...
         "     60   60,", ...
     "CMMNT Reported standard deviations are formal (not calibrated)"];
    dataStr = fileread(dataPath);

    for iBarrier = 1:length(dataBarriers)

        if contains(dataStr, dataBarriers(iBarrier))
            dataStr = strsplit(dataStr, dataBarriers(iBarrier));
            break
        end

        if iBarrier ~= length(dataBarriers)
            continue
        end

        error("Unrecognised file format")
    end

    data = dataStr{end};
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
    header = dataStr{1};
    lines = strtrim(strsplit(header, '\n'));

    try
        gravityParam = extractheadervalue(lines, 'earth_gravity_param');
    catch
        gravityParam = 0.3986004415E+15;
    end

    try
        equatorRadius = extractheadervalue(lines, 'mean_equator_radius');
    catch
        equatorRadius = 0.6378136300E+07;
    end

    varargout = ...
        {gravitySph, gravityStdSph, meanDate, [startDate, endDate], gravityParam, equatorRadius};
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
