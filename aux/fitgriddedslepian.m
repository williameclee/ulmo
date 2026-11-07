%% FITGRIDDEDSLEPIAN - Fits named rectangular-grid fields to a retained basis.
%
% Syntax
%   [coeffs, fits, buffer, weighting] = FITGRIDDEDSLEPIAN( ...
%       fields, dates, lon, lat, domain, L, GVN, SaveData=false)
%
% Input arguments
%   fields - Scalar structure of real single/double grids.
%       Each field must use (longitude, latitude, time) ordering, including square grids.
%       Missing observations are NaN. Field names have no physical meaning
%       here. Arrays share coordinates and dates, but may have different masks.
%   dates - Strictly increasing datetime vector, one date per epoch.
%   lon, lat - Coordinate vectors accepted by PREPAREGRIDCOORDINATES.
%       Unit: degrees east and north, respectively.
%   domain - Target domain accepted by DOMAINMASK.
%   L - Nonnegative integer harmonic degree. Positive degrees use a 90/L
%       degree spherical exterior buffer; degree zero uses only the target.
%   GVN - Retained {G,V,N} prepared by PREPARESLEPIANBASIS.
%   SaveData (name-value) - Whether exterior-mask caches may be saved.
%       Default flag: true.
%   ForceNew (name-value) - Whether exterior masks must be recalculated.
%       Default flag: false.
%   RequireGlobalLongitude (name-value) - Whether full longitude coverage
%       is required by the caller. Default flag: false.
%
% Output arguments
%   coeffs - Structure with the input field names, each (nTime, J).
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
%   2026/11/07, En-Chi Lee (williameclee@gmail.com)

function [coeffs, fits, buf, wghtMthd] = fitgriddedslepian( ...
        fields, dates, lon, lat, domain, L, GVN, options)

    arguments (Input)
        fields (1, 1) struct
        dates datetime {mustBeVector, mustBeNonempty, mustBeFinite}
        lon
        lat
        domain
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
        GVN cell {mustBeGVN}
        options.RequireGlobalLongitude (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
    end

    nLon = numel(lon);
    nLat = numel(lat);
    [lon, lat, keepLon, orderLat, wghts, wghtMthd] = preparegridcoordinates( ...
        lon, lat, RequireGlobalLongitude = options.RequireGlobalLongitude);
    dates = dates(:);
    names = fieldnames(fields);
    assert(~isempty(names), 'ULMO:fitgriddedslepian:EmptyRecord', 'Supply at least one field.');
    hasAnyObs = false(nLon, nLat);

    for iField = 1:numel(names)
        vals = fields.(names{iField});
        assert((isa(vals, 'double') || isa(vals, 'single')) && isreal(vals) && ...
            ndims(vals) <= 3 && isequal(size(vals, 1:3), [nLon, nLat, numel(dates)]), ...
            'ULMO:fitgriddedslepian:InvalidGrid', ...
            'Field %s must have longitude-by-latitude-by-time dimensions.', names{iField});
        hasAnyObs = hasAnyObs | any(~isnan(vals), 3);
    end

    mask = domainmask(domain, lon, lat, Format = 'ndgrid', SaveData = options.SaveData, BeQuiet = true);

    if L > 0
        buf = 90 / L;
        obsMask = domainmask(domain, lon, lat, Format = 'ndgrid', ...
            exteriorBuffer = buf, SaveData = options.SaveData, ForceNew = options.ForceNew, BeQuiet = true);
    else
        obsMask = mask;
    end

    [gridLon, gridLat] = ndgrid(lon, lat);
    totalArea = sum(wghts(mask));
    assert(totalArea > 0, 'ULMO:fitgriddedslepian:EmptyTarget', 'The target contains no positive-area grid cells.');
    % Retain rows used by any field, including fields fitted later. This
    % avoids sampling permanently missing cells in the single-field SSH path.
    pts = find(obsMask & wghts > 0 & hasAnyObs(keepLon, orderLat));
    wghts = wghts(pts);
    target = mask(pts);
    % Implicitly assumes that the supplied basis only contains the basis functions wanted
    J = size(GVN{1}, 2);
    nTime = numel(dates);
    design = [];
    coeffs = struct();
    fits = struct();

    for name = string(fieldnames(fields)).'
        vals = fields.(name)(keepLon, orderLat, :);
        vals = reshape(vals, [], nTime);
        vals = vals(pts, :);
        assert(~any(isinf(vals), 'all'), 'ULMO:fitgriddedslepian:InvalidObservation', ...
            '%s contains infinite observations.', name);
        observed = isfinite(vals);
        nObs = sum(observed, 1).';
        areaFraction = zeros(nTime, 1);

        if any(target)
            areaFraction = min(1, (wghts(target).' * observed(target, :)).' / totalArea);
        end

        epochs = find(nObs > 0);
        assert(~isempty(epochs), 'ULMO:fitgriddedslepian:EmptyRecord', ...
            'No observed cells for %s across all epochs.', name);

        try

            if isempty(design)
                [c, ~, ~, stats, design] = xyzs2slep( ...
                    vals(:, epochs), gridLon(pts), gridLat(pts), ...
                    domain, L, GVN = GVN, truncation = J, weight = wghts, missingPolicy = "omit");
            else
                [c, ~, ~, stats] = xyzs2slep(vals(:, epochs), design = design, ...
                    weight = wghts, missingPolicy = "omit");
            end

        catch cause
            field = regexpi(cause.message, '\bfield (\d+)\b', 'tokens', 'once');
            ctx = "Fit failed for " + name;

            if ~isempty(field)
                ctx = ctx + " at " + string(dates(epochs(str2double(field{1}))));
            end

            throwAsCaller(addCause(MException('ULMO:fitgriddedslepian:FitFailed', '%s: %s', ctx, cause.message), cause));
        end

        coeffs = nan(nTime, J);
        condNum = nan(nTime, 1);
        resRms = nan(nTime, 1);
        coeffs(epochs, :) = c.';
        condNum(epochs) = stats.condNum;
        resRms(epochs) = stats.residualRms;
        coeffs.(name) = coeffs;
        fits.(name) = table(nObs, areaFraction, condNum, resRms);
    end

end
