%% CPRINTF - FPRINTF wrapper removing HTML tags when appropriate
% Should behave exactly like FPRINTF when logging to the command window
% (won't work for file output), but if HotLinks is disabled (i.e. in a
% terminal environment), it will remove any <a></a> tags from the input
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

    %% Removing tags
    % Function/script files: only leave the name
    txt = regexprep(txt, '<a[^>]*href="matlab: open\(''(.*?\.m)''\)"[^>]*>(.*?)</a>', '$2');
    % Data files: show alt name and path
    txt = regexprep(txt, '<a[^>]*href="matlab: fprintf\(''(.*?)\\n''\);open\(''(.*?)''\)"[^>]*>(.*?)</a>', '$3 ($2)');
    % External links: show alt name and URL
    txt = regexprep(txt, '<a[^>]*href="(https?://[^"]*)"[^>]*>(.*?)</a>', '$2 ($1)');
    % Remaining tags
    txt = regexprep(txt, '<a[^>]*>', '');
    txt = strrep(txt, '</a>', '');

    cnt = fprintf('%s', txt);

    if nargout > 0
        varargout = {cnt};
    else
        varargout = {};
    end

end
