%% CONDDEFVAL
% Conditionally defines a value based on whether it is empty.
%
% Syntax
%   v = conddefval(v, v0)
%
% Input arguments
%   v - The input value to be checked
%   v0 - The default value to be returned if the input value is empty
%
% Output argument
%   v - The resulting value based on the condition.
%
% Example
%   The following example will return v = 10:
%   >>  v = [];
%   >>  v = conddefval(v, 10);
%
% Last modified by
%   2024/08/28, williameclee@arizona.edu (@williameclee)

function v = conddefval(v, v0, options)

    arguments
        v
        v0
        options.IncludeNan (1, 1) logical = false
    end

    includeNan = logical(options.IncludeNan);

    if isempty(v)
        v = v0;
        return
    end

    if includeNan && isscalar(v) && isnumeric(x) && isnan(v)
        v = v0;
    end

end
