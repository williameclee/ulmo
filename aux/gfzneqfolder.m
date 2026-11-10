%% GFZNEQFOLDER - Resolves the local folder for processed GFZ covariance factors.
%
% See also
%   PROCESSGFZNORMALEQUATIONS, RANDOMGRACENOISE
%
% Created by
%   2026/11/10
function folder = gfzneqfolder()
    root = getenv('GRACEDATA');
    if isempty(root)
        root = getenv('IFILES');
        assert(~isempty(root), 'ULMO:gfzneq:NoRoot', 'Set GRACEDATA or IFILES.');
        root = fullfile(root, 'GRACE');
    end
    folder = fullfile(root, 'GFZ', 'NEQ');
end
