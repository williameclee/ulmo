%% matFileHasVariables - Checks if a variable is saved in a .mat file.
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
