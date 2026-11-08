%% FITGRIDDEDSLEPIAN - Fits gridded fields to a retained Slepian basis.
%
% Syntax
%   [coeffs, fits, buf, wghtMthd] = FITGRIDDEDSLEPIAN(fields, dates, lon, lat, domain, L, GVN)
%   [__] = FITGRIDDEDSLEPIAN(__, "Name", value)
%
% Input arguments
%   fields - Scalar structure of real fields.
%       Fields share the same grid and dates.
%       Size of each field: (nLon, nLat, nTime), including square grids.
%   dates - Strictly increasing datetime vector with one date per epoch.
%   lon, lat - Coordinate vectors accepted by PREPAREGRIDCOORDINATES.
%       Unit: degrees east and north, respectively.
%   domain - Target domain.
%   L - Maximum harmonic degree for the basis, used to determine the buffer zone size.
%   GVN - Retained {G, V, N} from PREPARESLEPIANBASIS. The fit uses every column of G.
%   RequireGlobalLongitude (name-value) - Whether to require a complete longitude cycle.
%       Default flag: false.
%   SaveData (name-value) - Whether to save eligible target and exterior masks.
%       Default flag: true.
%   ForceNew (name-value) - Whether to recompute the exterior mask.
%       Default flag: false.
%
% Output arguments
%   coeffs - Structure of coefficient matrices with the same field names as the input.
%       Size of each matrix: (nTime, J), where J is the number of columns of G.
%       Unit: same as the corresponding input field.
%   fits - Structure of diagnostic tables with one row per epoch.
%       - nObs counts positive-weight observations in the target and buffer.
%       - areaFraction is the observed fraction of the total target area.
%       - condNum is the fit condition number.
%       - residualRms is the residual root mean square in the units of the input field.
%   buf - Exterior buffer width, 90/L for positive L and zero for L = 0.
%       Unit: degrees.
%   wghtMthd - Relative area convention returned by PREPAREGRIDCOORDINATES.
%
% Notes
%   Empty epochs have zero counts and coverage, and NaN fit outputs.
%   Entirely empty fields and underdetermined fits raise errors.
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

    dates = dates(:);
    assert(all(diff(dates) > seconds(0)), 'ULMO:fitgriddedslepian:EmptyRecord', ...
    'Supply strictly increasing epoch dates.');
    nLon = numel(lon);
    nLat = numel(lat);
    [lon, lat, keepLon, orderLat, wghts, wghtMthd] = preparegridcoordinates( ...
        lon, lat, RequireGlobalLongitude = options.RequireGlobalLongitude);
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
        buf = 0;
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
    % The caller has already selected the retained columns.
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

        fieldCoeff = nan(nTime, J);
        condNum = nan(nTime, 1);
        residualRms = nan(nTime, 1);
        fieldCoeff(epochs, :) = c.';
        condNum(epochs) = stats.condNum;
        residualRms(epochs) = stats.residualRms;
        coeffs.(name) = fieldCoeff;
        fits.(name) = table(nObs, areaFraction, condNum, residualRms);
    end

end
