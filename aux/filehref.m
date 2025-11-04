function urltxt = filehref(url, txt)

    arguments (Input)
        url (1, :) char
        txt (1, :) char = NaN
    end

    arguments (Output)
        urltxt (1, :) char
    end

    if isnan(txt)
        [~, name, ext] = fileparts(url);
        txt = [name, ext];
    end

    urltxt = sprintf('<a href="matlab: fprintf(''%s'');open(''%s'')">%s</a>', url, url, txt);

end
