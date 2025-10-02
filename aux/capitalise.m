%% CAPITALISE
% Capitalises the first letter of a string.
%
% Syntax
%   sc = capitalise(s)
%
% Input arguments
%   s - The string to capitalise
%
% Output arguments
%   sc - The capitalised string
%
% Last modified by
%   2025/08/15, williameclee@arizona.edu (@williameclee)

function capitalisedString = capitalise(inputString)
    isString = isstring(inputString);
    % Empty string
    if isempty(inputString)

        if isString
            capitalisedString = string.empty;
        else
            capitalisedString = char.empty;
        end

        return
    end

    if isString
        inputString = char(inputString);
    end

    % Single character
    if isscalar(inputString)
        capitalisedString = upper(inputString);
        return
    end

    % Multiple characters
    capitalisedString = ...
        [upper(extractBefore(inputString, 2)), ...
         extractAfter(inputString, 1)];
end
