%% MATFILEHASVARIABLES - Checks whether a MAT file contains all requested variables.
% Inspects variable names without loading their values. Missing or
% unreadable files are treated as incomplete and return false.
%
% Syntax
%   ok = matFileHasVariables(path, names)
%
% Input arguments
%   path - MAT-file path as a character vector or string scalar.
%   names - Variable names as a cell array of character vectors or a
%       string array. Every requested name must be present.
%
% Output arguments
%   ok - Scalar logical indicating whether every requested variable exists.
%       An empty names array returns true when the MAT file is readable.
%
% Notes
%   Tests presence only: values, types, dimensions, and nonemptiness are
%   not checked. Errors raised while inspecting the file or comparing
%   names are caught and return false.
%
% Last modified
%   2026/10/09, En-Chi Lee (williameclee@gmail.com)

function ok = matFileHasVariables(path, names)
    ok = false;

    try
        info = whos('-file', path);
        ok = all(ismember(names, {info.name}));
    catch
    end

end
