%% FETCHURLDATA - downloads and unzips files from URLs in a text file
%
% Author
%   2026/02/13, williameclee@arizona.edu (@williameclee)

function fetchUrlData(urlFilePath, outputDir)

    arguments (Input)
        urlFilePath (1, 1) string {mustBeFile(urlFilePath)}
        outputDir (1, 1) string
    end

    % Create the output directory if it doesn't exist
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Read the URLs from the text file
    fileID = fopen(urlFilePath, 'r');
    urls = textscan(fileID, '%s');
    fclose(fileID);

    % Convert cell array to a string array
    urls = string(urls{1});

    % Loop through each URL, download the file, and unzip it
    if license('test', 'Distrib_Computing_Toolbox')

        parfor i = 1:length(urls)

            try
                fetchOneUrl(urls{i}, outputDir);
            catch ME
                fprintf('Error processing URL: %s\n', urls{i});
                fprintf('Error message: %s\n', ME.message);
            end

        end

    else

        for i = 1:length(urls)

            try
                fetchOneUrl(urls{i}, outputDir);
            catch ME
                fprintf('Error processing URL: %s\n', urls{i});
                fprintf('Error message: %s\n', ME.message);
            end

        end

    end

    fprintf('Download and extraction complete.\n');
end

%% Subfunction
function fetchOneUrl(url, outputDir)
    % Extract the file name from the URL
    [~, fileName, ext] = fileparts(url);
    fileName = strcat(fileName, ext);

    % Define the full path to save the file
    outputFilePath = fullfile(outputDir, fileName);

    % Download the file
    fprintf('Downloading: %s\n', url);
    websave(outputFilePath, url);

    % Unzip the file
    fprintf('Unzipping: %s\n', outputFilePath);
    unzip(outputFilePath, outputDir);

    % Optionally, delete the zip file after extraction
    delete(outputFilePath);
end
