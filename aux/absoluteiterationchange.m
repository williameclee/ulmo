%% ABSOLUTEITERATIONCHANGE - Largest absolute update over coefficients and epochs.
% Syntax
%   maxChange = ABSOLUTEITERATIONCHANGE(crnt, prev)
%
% Input arguments
%   crnt, prev - Coefficient matrices of equal size, with epochs in columns.
%
% Output arguments
%   maxChange - Maximum absolute elementwise difference, in input units.
%
% Notes
%   Nonfinite inputs return Inf and cannot satisfy a finite tolerance.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

function maxChange = absoluteiterationchange(crnt, prev)

    arguments (Input)
        crnt (:, :) double
        prev (:, :) double
    end

    arguments (Output)
        maxChange (1, 1) double
    end

    assert(isequal(size(crnt), size(prev)), ...
        'ULMO:absoluteiterationchange:SizeMismatch', ...
    'The arrays being compared must have equal sizes.');

    if any(~isfinite(crnt), 'all') || any(~isfinite(prev), 'all')
        maxChange = Inf;
        return
    end

    maxChange = max(abs(crnt - prev), [], 'all');
end
