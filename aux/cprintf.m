%% CPRINTF - FPRINTF wrapper removing HTML tags when appropriate
% Should behave exactly like FPRINTF, but if HotLinks is disabled (i.e. in 
% a terminal environment), it will remove any <a></a> tags from the input 
% string before printing.
%
% Author
%   2026/02/14, williameclee@arizona.edu (@williameclee)

function varargout = cprintf(varargin)
    txt = sprintf(varargin{:});

    if feature('HotLinks')
        cnt = fprintf('%s', txt);

        if nargout > 0
            varargout = {cnt};
        else
            varargout = {};
        end

        return
    end

    % Remove <a></a> tags if present
    txt = regexprep(txt, '<a[^>]*>', '');
    txt = strrep(txt, '</a>', '');

    cnt = fprintf('%s', txt);

    if nargout > 0
        varargout = {cnt};
    else
        varargout = {};
    end

end
