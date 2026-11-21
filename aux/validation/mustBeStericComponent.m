%% MUSTBESTERICCOMPONENT
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function mustBeStericComponent(comp)
    mustBeMember(comp, ...
        ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"]);
end
