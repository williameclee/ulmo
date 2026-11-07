%% FITGRIDDEDSLEPIAN - Fit named rectangular-grid fields to a retained basis.
%
% Syntax
%   [coefficients, fits, buffer, weighting] = FITGRIDDEDSLEPIAN( ...
%       fields, dates, lon, lat, domain, L, basis, SaveData=false)
%
% Input arguments
%   fields - Scalar structure of real single/double grids. Each field must
%       use (longitude, latitude, time) ordering, including square grids.
%       Missing observations are NaN. Field names have no physical meaning
%       here. Arrays share coordinates and dates, but may have different masks.
%   dates - Strictly increasing datetime vector, one date per epoch.
%   lon, lat - Coordinate vectors accepted by PREPAREGRIDCOORDINATES.
%       Unit: degrees east and north, respectively.
%   domain - Target domain accepted by DOMAINMASK.
%   L - Nonnegative integer harmonic degree. Positive degrees use a 90/L
%       degree spherical exterior buffer; degree zero uses only the target.
%   basis - Retained {G,V,N} prepared by PREPARESLEPIANBASIS.
%   SaveData (name-value) - Whether exterior-mask caches may be saved.
%       Default flag: true.
%   ForceNew (name-value) - Whether exterior masks must be recalculated.
%       Default flag: false.
%   RequireGlobalLongitude (name-value) - Whether full longitude coverage
%       is required by the caller. Default flag: false.
%
% Output arguments
%   coefficients - Structure with the input field names, each (nTime, J).
%       Unit: same as the corresponding input field.
%   fits - Structure of per-field tables with nObs, areaFraction, condNum
%       and residualRms. Counts include positive-area buffer observations;
%       areaFraction uses the complete target-grid area as its denominator.
%       Residual RMS has input field units; other diagnostics are dimensionless.
%   buffer - Exterior buffer width, in degrees.
%   weighting - Relative area convention from PREPAREGRIDCOORDINATES.
%
% Notes
%   Empty epochs retain NaN coefficients, conditioning and residual RMS,
%   and zero counts and coverage. Entirely empty fields and underdetermined
%   nonempty fits raise errors. The sampled design is reused across fields,
%   including points first observed in a later field. Temporal processing,
%   uncertainty propagation, output formatting and projection caches belong
%   to callers. No gridded data or projection files are written here.
%
% See also
%   preparegridcoordinates, prepareslepianbasis, xyzs2slep, domainmask
%
% Created by
%   2026/11/06, En-Chi Lee

function [coefficients, fits, buffer, weighting] = fitgriddedslepian(fields, dates, lon, lat, domain, L, basis, options)
    arguments (Input)
        fields (1, 1) struct
        dates datetime
        lon
        lat
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        basis cell
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.RequireGlobalLongitude (1, 1) logical = false
    end
    nLon = numel(lon); nLat = numel(lat);
    [lon, lat, keepLon, orderLat, weights, weighting] = preparegridcoordinates( ...
        lon, lat, RequireGlobalLongitude=options.RequireGlobalLongitude);
    assert(isvector(dates) && ~isempty(dates) && ~any(isnat(dates)) && ...
        all(diff(dates(:)) > seconds(0)), 'ULMO:fitgriddedslepian:EmptyRecord', ...
        'Supply nonempty strictly increasing epoch dates.');
    dates = dates(:);
    names = fieldnames(fields);
    assert(~isempty(names), 'ULMO:fitgriddedslepian:EmptyRecord', 'Supply at least one field.');
    observedAnywhere = false(nLon, nLat);
    for iField = 1:numel(names)
        values = fields.(names{iField});
        assert((isa(values, 'double') || isa(values, 'single')) && isreal(values) && ...
            ndims(values) <= 3 && isequal(size(values, 1:3), [nLon, nLat, numel(dates)]), ...
            'ULMO:fitgriddedslepian:InvalidGrid', ...
            'Field %s must have longitude-by-latitude-by-time dimensions.', names{iField});
        observedAnywhere = observedAnywhere | any(~isnan(values), 3);
    end
    mask = domainmask(domain, lon, lat, Format='ndgrid', SaveData=false, BeQuiet=true);
    buffer = 0;
    if L > 0, buffer = 90/L; end
    observationMask = mask;
    if buffer > 0
        observationMask = domainmask(domain, lon, lat, Format='ndgrid', ...
            exteriorBuffer=buffer, SaveData=options.SaveData, ForceNew=options.ForceNew, BeQuiet=true);
    end
    [gridLon, gridLat] = ndgrid(lon, lat);
    totalArea = sum(weights(mask));
    assert(totalArea > 0, 'ULMO:fitgriddedslepian:EmptyTarget', 'The target contains no positive-area grid cells.');
    % Retain rows used by any field, including fields fitted later. This
    % avoids sampling permanently missing cells in the single-field SSH path.
    points = find(observationMask & weights > 0 & observedAnywhere(keepLon, orderLat));
    weights = weights(points);
    target = mask(points);
    % The caller selects and validates the retained basis once. XYZS2SLEP
    % validates it again when sampling; no columns are selected here.
    J = size(basis{1}, 2);
    nTime = numel(dates);
    assert(nTime > 0, 'ULMO:fitgriddedslepian:EmptyRecord', 'Temporal selection produced no epochs.');
    design = [];
    coefficients = struct();
    fits = struct();
    for name = string(fieldnames(fields)).'
        values = fields.(name)(keepLon, orderLat, :);
        values = reshape(values, [], nTime);
        values = values(points, :);
        assert(~any(isinf(values), 'all'), 'ULMO:fitgriddedslepian:InvalidObservation', ...
            '%s contains infinite observations.', name);
        observed = isfinite(values);
        nObs = sum(observed, 1).';
        areaFraction = zeros(nTime, 1);
        if any(target)
            areaFraction = min(1, (weights(target).' * observed(target, :)).' / totalArea);
        end
        epochs = find(nObs > 0);
        assert(~isempty(epochs), 'ULMO:fitgriddedslepian:EmptyRecord', 'No observed cells for %s.', name);
        try
            if isempty(design)
                [c, ~, ~, stats, design] = xyzs2slep(values(:, epochs), gridLon(points), gridLat(points), ...
                    domain, L, GVN=basis, truncation=J, weight=weights, missingPolicy="omit");
            else
                [c, ~, ~, stats] = xyzs2slep(values(:, epochs), design = design, weight=weights, missingPolicy="omit");
            end
        catch cause
            field = regexpi(cause.message, '\bfield (\d+)\b', 'tokens', 'once');
            context = "Fit failed for " + name;
            if ~isempty(field), context = context + " at " + string(dates(epochs(str2double(field{1})))); end
            throwAsCaller(addCause(MException('ULMO:fitgriddedslepian:FitFailed', '%s: %s', context, cause.message), cause));
        end
        coeffs = nan(nTime, J);
        condNum = nan(nTime, 1);
        residualRms = nan(nTime, 1);
        coeffs(epochs, :) = c.';
        condNum(epochs) = stats.condNum;
        residualRms(epochs) = stats.residualRms;
        coefficients.(name) = coeffs;
        fits.(name) = table(nObs, areaFraction, condNum, residualRms);
    end
end
