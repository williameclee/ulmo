%% NUM2FILENAME - Formats a finite number for a readable filename.
% Uses the shortest significant-digit precision that round-trips through
% STR2DOUBLE. Decimal points become p and minus signs become n, including
% negative exponents. Zero has no sign; positive exponents omit the plus.
%
% Syntax
%   text = num2filename(value)
%
% Input arguments
%   value - Finite real numeric scalar, converted to double.
%
% Output arguments
%   text - Character vector, for example 'n30p72' or '1en8'.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function text = num2filename(value)
    arguments (Input)
        value (1, 1) double {mustBeReal, mustBeFinite}
    end
    if value == 0
        text = '0';
        return
    end
    for precision = 1:17
        text = sprintf('%.*g', precision, value);
        if str2double(text) == value
            break
        end
    end
    if abs(value) >= 1e-6 && abs(value) < 1e15
        decimalPlaces = max(0, precision - 1 - floor(log10(abs(value))));
        text = sprintf('%.*f', decimalPlaces, value);
    end
    text = strrep(strrep(strrep(text, '.', 'p'), '-', 'n'), '+', '');
end
