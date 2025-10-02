%% ISOLDER
% Checks if the first directory is older than the second one.
%
% Input arguments
%	checkContent - Logical flag to check the content of the directory
%		- true: Check the content of the directory (excluding . and ..).
%		- false: Check the date of the directory only.
%
% See also
%	GRACEDEG1, SSH2XYZ
%
% Notes
%	This is a helper function with limited documentation.
%
% Authored by
%   2025/05/20, williameclee@arizona.edu (@williameclee)

function isolder = isolder(dir1, dir2, checkContent)

    dirDates = NaT([2, 1]);

    for iDir = 1:2

        if iDir == 1
            dirPath = dir1;
        else
            dirPath = dir2;
        end

        dirState = exist(dirPath, 'file');

        switch dirState
            case 2 % is a file
                dirDate = datetime(dir(dirPath).date, ...
                    "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
            case 7 % is a folder
                dirInfo = dir(dirPath);

                if ~checkContent
                    dirDate = datetime(dirInfo(strcmp({dirInfo.name}, '.')).date, ...
                        "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
                else
                    disDates = ...
                        datetime({dirInfo(~ismember({dirInfo.name}, {'.', '..', '.DS_Store'})).date}, ...
                        "InputFormat", 'dd-MMM-yyyy HH:mm:ss');
                    dirDate = max(disDates);
                end

            case 0 % does not exist
                dirDate = datetime('now');
            otherwise
                error('Unknown file type for %s, EXIST output: %d', dirPath, dirState);
        end

        dirDates(iDir) = dirDate;

    end

    isolder = dirDates(1) < dirDates(2);

end
