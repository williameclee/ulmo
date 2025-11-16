%% FILEHREF
% Generate an HTML hyperlink to a file in MATLAB command window.
%
% Last modified by
%	2025/11/11, williameclee@arizona.edu (@williameclee)

function urltxt = filehref(url, txt)

    arguments (Input)
        url {mustBeTextScalar}
        txt {mustBeTextScalar} = NaN
    end

    arguments (Output)
        urltxt (1, :) char
    end

    if isnan(txt)
        [~, name, ext] = fileparts(url);
        txt = [name, ext];
    end

    urltxt = sprintf('<a href="matlab: fprintf(''%s\\n'');open(''%s'')">%s</a>', url, url, txt);

end
