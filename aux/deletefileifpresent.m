%% DELETEFILEIFPRESENT - Deletes a file if it exists.
%
% Syntax
%   deletefileifpresent(path)
%
% Input arguments
%   path - Path to the file to remove.
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function deletefileifpresent(path)

    arguments (Input)
        path {mustBeTextScalar}
    end

    if isfile(path)
        % Temporary links must also be removed when MATLAB recycling is on.
        folder = java.io.File(pwd).toPath();
        java.nio.file.Files.deleteIfExists(folder.resolve(char(path)));
    end

end
