%% RELATIVEITERATIONCHANGE - Largest relative infinity-norm update over epochs.
% Syntax
%   maxRelChange = RELATIVEITERATIONCHANGE(crnt, prev)
%
% Input arguments
%   crnt, prev - Coefficient matrices.
%
% Output arguments
%   maxRelChange - Maximum of norm(crnt-prev, Inf) divided by
%       max(norm(crnt, Inf), norm(prev, Inf)). Dimensionless.
%
% Notes
%   An unchanged zero vector has zero relative change. A transition between
%   zero and a nonzero vector has relative change one. No dimensional floor
%   is used, so rescaling any columns does not change its stopping criterion.
%   Nonfinite inputs return Inf and cannot satisfy a finite tolerance.
%
% Created by
%   2026/10/25, En-Chi Lee (williameclee@gmail.com)

function maxRelChange = relativeiterationchange(crnt, prev)

    arguments (Input)
        crnt (:, :) double
        prev (:, :) double
    end

    arguments (Output)
        maxRelChange (1, 1) double
    end

    assert(isequal(size(crnt), size(prev)), ...
        'ULMO:relativeiterationchange:SizeMismatch', ...
    'The arrays being compared must have equal sizes.');

    if any(~isfinite(crnt), 'all') || any(~isfinite(prev), 'all')
        maxRelChange = Inf;
        return
    end

    scale = max(max(abs(crnt), [], 1), max(abs(prev), [], 1)); % Max magnitude per column
    change = max(abs(crnt - prev), [], 1); % Max change magnitude per column
    relChange = zeros(size(scale));
    active = scale > 0;
    relChange(active) = change(active) ./ scale(active);
    maxRelChange = max(relChange, [], 'all');
end
