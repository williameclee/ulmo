%% NUM2FILENAME - Formats a finite number for a readable filename.
% Uses the shortest significant-digit precision that round-trips through
% STR2DOUBLE. Decimal points become p and minus signs become n, including
% exponents.
%
% Syntax
%   txt = num2filename(val)
%
% Input arguments
%   val - Finite real numeric scalar, converted to double.
%
% Output arguments
%   txt - Character vector, for example 'n30p72' or '1en8'.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function txt = num2filename(val)

    arguments (Input)
        val (1, 1) double {mustBeReal, mustBeFinite}
    end

    if val == 0
        txt = '0';
        return
    end

    for precision = 1:17
        txt = sprintf('%.*g', precision, val);

        if str2double(txt) == val
            break
        end

    end

    if abs(val) >= 1e-6 && abs(val) < 1e15
        decimalPlaces = max(0, precision - 1 - floor(log10(abs(val))));
        decimalText = sprintf('%.*f', decimalPlaces, val);

        if str2double(decimalText) == val
            txt = decimalText;
        end

    end

    txt = strrep(strrep(strrep(txt, '.', 'p'), '-', 'n'), '+', '');
end
