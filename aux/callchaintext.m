function callChainTxt = callchaintext(funs)

    arguments
        funs {mustBeA(funs, {'char', 'string', 'cell'})}
    end

    if iscell(funs)

        if any(cellfun(@(x) ~(ischar(x) || isstring(x)), funs))
            error('All elements of must be a string');
        end

    else
        funs = {funs};
    end

    callChainTxt = cellfun(@(x) sprintf('<a href="matlab: open(''%s'')">%s</a>', which(x), x), funs, "UniformOutput", false);
    callChainTxt = strjoin(callChainTxt, '>');

end
