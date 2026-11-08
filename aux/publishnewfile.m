%% PUBLISHNEWFILE - Publishes a completed file under a new name.
% Create a hard link in one filesystem operation. An existing destination
% remains intact even when another writer creates it during processing.
% The caller owns the temporary file and removes its name after publication.
%
% Syntax
%   publishnewfile(tmpPath, outPath)
%
% Input arguments
%   tmpPath - Path to a completed temporary file on the destination filesystem.
%   outPath - New file path on a filesystem that supports hard links.
%
% Created by
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

function publishnewfile(tmpPath, outPath)

    arguments (Input)
        tmpPath {mustBeTextScalar}
        outPath {mustBeTextScalar}
    end

    % MATLAB's current folder may differ from the JVM working directory.
    folder = java.io.File(pwd).toPath();
    dest = folder.resolve(char(outPath));
    src = folder.resolve(char(tmpPath));

    try
        java.nio.file.Files.createLink(dest, src);
    catch exception

        if isa(exception, 'matlab.exception.JavaException') && ...
                isa(exception.ExceptionObject, 'java.nio.file.FileAlreadyExistsException')
            error('ULMO:OutputExists', 'Output already exists: %s', outPath);
        end

        rethrow(exception)
    end

end
