%% CALLCHAINTEXT - Generates a call stack for a list of functions.
% Syntax
%   txt = callchaintext(fun)
%   txt = callchaintext([fun1, fun2, ...])
%
% Input arguments
%   fun - Function name as a character vector or string scalar
%
% Output arguments
%   txt - Call stack as a string scalar with clickable links to the
%       functions
%
% Last modified
%   2026/02/14, williameclee@arizona.edu (@williameclee)
%      - Removed consecutive identical functions from the input list

function callChainTxt = callchaintext(funs)

    arguments (Input)
        funs {mustBeA(funs, {'char', 'string', 'cell'})}
    end

    arguments (Output)
        callChainTxt (1, 1) string
    end

    if iscell(funs)

        if any(cellfun(@(x) ~(ischar(x) || isstring(x)), funs))
            error('All elements of must be a string');
        end

    else
        funs = {funs};
    end

    % Remove consecutive identical elements
    funs = funs([true, ~strcmp(funs(1:end - 1), funs(2:end))]);

    callChainTxt = cellfun(@(x) sprintf('<a href="matlab: open(''%s'')">%s</a>', which(x), x), funs, "UniformOutput", false);
    callChainTxt = strjoin(callChainTxt, '>');

end
