%% GLDAS2PLMT - Converts monthly Noah storage to harmonics over observation windows.
%
% Syntax
%   plmt = gldas2plmt(epochs, L, unit, loveNumSrc, eqRad)
%
% Input arguments
%   epochs - Observation windows as unzoned datetimes interpreted as UTC.
%       Each window must have positive duration. Size: (nDates, 2).
%   L - Maximum harmonic degree, from 0 to 179 for the 1-degree input grid.
%   unit - SD, POT or GRAV, matching the selected GRACE record.
%   loveNumSrc - Source passed to CONVERTGRAVITY.
%   eqRad - GRACE reference equatorial radius.
%       Unit: metres.
%
% Input file content
%   GLDAS_NOAH10_M.AYYYYMM.021.nc4 - Official Noah 2.1 monthly 1-degree
%       files, under GLDASDATA (recursively), or IFILES/EARTHMODELS/GLDAS2.1.
%       Requires lon, lat and the six kg/m^2 storage fields: SoilMoi0_10cm_inst,
%       SoilMoi10_40cm_inst, SoilMoi40_100cm_inst, SoilMoi100_200cm_inst,
%       CanopInt_inst and SWE_inst. Native grid is 360 by 150, north of 60 S.
%
% Output arguments
%   plmt - Time-first lmcosi array. Size: (nDates, addmup(L), 4).
%       Unit: requested unit. Degree zero is zero, as in the GRACE loader.
%
% Notes
%   Monthly means are constant over their calendar month and weighted by
%   overlap duration. Every part of every observation window must be covered.
%   The six fields must share a time-invariant native missing-data mask;
%   that mask and latitudes south of 60 S contribute zero. No extra glacier
%   mask or groundwater term is added. No files are written.
%   XYZ2PLM uses Simpson integration with four intervals per degree in
%   latitude. The field is constant within each cell. At cell boundaries,
%   the adjacent values are averaged because both cells share the endpoint.
%   The longitude coefficients account for the cell centres and widths.
%
% Data sources
%   NASA GLDAS: https://ldas.gsfc.nasa.gov/gldas/model-output
%   https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH10_M.2.1/
%
% See also
%   convertgravity, xyz2plm
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function plmt = gldas2plmt(epochs, L, unit, loveNumSrc, eqRad)

    arguments (Input)
        epochs (:, 2) datetime
        L (1, 1) double {mustBeInteger, mustBeNonnegative, mustBeLessThan(L, 180)}
        unit (1, 1) string
        loveNumSrc (1, 1) string
        eqRad (1, 1) double {mustBeFinite, mustBePositive}
    end

    root = getenv('GLDASDATA');

    if isempty(root)
        root = fullfile(getenv('IFILES'), 'EARTHMODELS', 'GLDAS2.1');
    end

    files = dir(fullfile(root, '**', 'GLDAS_NOAH10_M.A*.021.nc4'));
    assert(~isempty(files), 'ULMO:gldas2plmt:MissingData', ...
        'No GLDAS_NOAH10_M 2.1 monthly files found under %s.', root);
    months = NaT(numel(files), 1);

    for iFile = 1:numel(files)
        token = regexp(files(iFile).name, '^GLDAS_NOAH10_M\.A(\d{6})\.021\.nc4$', 'tokens', 'once');
        assert(~isempty(token), 'ULMO:gldas2plmt:Product', 'Unexpected GLDAS filename.');
        months(iFile) = datetime(token{1}, 'InputFormat', 'yyyyMM');
    end

    [months, order] = sort(months);
    files = files(order);
    ends = months + calmonths(1);
    overlap = max(0, seconds(min(epochs(:, 2), ends.') - max(epochs(:, 1), months.')));
    used = any(overlap > 0, 1);
    months = months(used);
    files = files(used);
    overlap = overlap(:, used);
    assert(all(diff(months) > seconds(0)), 'ULMO:gldas2plmt:DuplicateMonth', ...
    'GLDAS must contain only one file per required calendar month.');
    lengths = seconds(epochs(:, 2) - epochs(:, 1));
    assert(all(abs(sum(overlap, 2) - lengths) < 1e-6), ...
        'ULMO:gldas2plmt:Coverage', 'GLDAS does not cover every observation window.');
    weights = overlap ./ lengths;
    [~, ~, ~, template] = addmon(L);
    nRows = size(template, 1);
    monthly = zeros(numel(files), nRows, 2);
    vars = {'SoilMoi0_10cm_inst', 'SoilMoi10_40cm_inst', ...
                'SoilMoi40_100cm_inst', 'SoilMoi100_200cm_inst', 'CanopInt_inst', 'SWE_inst'};
    nativeMask = [];
    intervalsPerDegree = 4;
    latSamples = 90 : -1 / intervalsPerDegree : -90;
    cellIndices = min(floor((0 : numel(latSamples) - 1) / intervalsPerDegree) + 1, 180);
    boundaryRows = 1 + intervalsPerDegree : intervalsPerDegree : numel(latSamples) - 1;
    boundaryCells = cellIndices(boundaryRows);
    orders = template(:, 2);
    lonCorrection = exp(-1i * orders * deg2rad(0.5)) .* sinc(orders / 360);

    for iMonth = 1:numel(files)
        path = fullfile(files(iMonth).folder, files(iMonth).name);
        lon = double(ncread(path, 'lon'));
        lat = double(ncread(path, 'lat'));
        assert(numel(lon) == 360 && numel(lat) == 150 && ...
            max(abs(lon(:) - (-179.5:179.5).')) < 1e-6 && ...
            max(abs(lat(:) - (-59.5:89.5).')) < 1e-6, ...
            'ULMO:gldas2plmt:Grid', 'GLDAS requires the complete native 1-degree grid.');
        field = zeros(360, 150);

        for iVar = 1:numel(vars)
            name = vars{iVar};
            info = ncinfo(path, name);
            dims = string({info.Dimensions.Name});
            assert(isequal(dims, ["lon", "lat", "time"]) && ...
                isequal(info.Size, [360, 150, 1]), ...
                'ULMO:gldas2plmt:Dimensions', '%s must have lon, lat and one monthly time.', name);
            units = regexprep(string(ncreadatt(path, name, 'units')), '\s', '');
            assert(ismember(units, ["kgm-2", "kg/m^2", "kgm^-2"]), ...
                'ULMO:gldas2plmt:Unit', '%s must use kg/m^2.', name);
            values = double(ncread(path, name));
            names = {info.Attributes.Name};

            if ismember('missing_value', names)
                missingValue = double(ncreadatt(path, name, 'missing_value'));
                values(values == missingValue) = NaN;
            end

            missing = ~isfinite(values);

            if isempty(nativeMask)
                nativeMask = missing;
            end

            assert(isequal(missing, nativeMask) && any(~missing, 'all') && ~any(isinf(values), 'all'), ...
                'ULMO:gldas2plmt:MissingValues', ...
            'GLDAS storage components must share a time-invariant native missing-data mask.');
            values(missing) = 0;
            field = field + values;
        end

        % Start at 0.5 E, repeat the first longitude, and pad south of 60 S.
        nativeCells = field([181:360, 1:181], end:-1:1).';
        globalCells = [nativeCells; zeros(30, size(nativeCells, 2))];
        samples = globalCells(cellIndices, :);
        samples(boundaryRows, :) = (globalCells(boundaryCells - 1, :) + globalCells(boundaryCells, :)) / 2;
        harmonics = xyz2plm(samples, L, 'simpson', latSamples);
        corrected = complex(harmonics(:, 3), -harmonics(:, 4)) .* lonCorrection;
        harmonics(:, 3:4) = [real(corrected), -imag(corrected)];

        harmonics(1, 3:4) = 0;
        harmonics = convertgravity(harmonics, 'SD', char(unit), ...
            'InputFormat', 'lmcosi', 'EquatorRadius', eqRad, 'LoveNumSource', char(loveNumSrc));
        monthly(iMonth, :, :) = harmonics(:, 3:4);
    end

    plmt = cat(3, repmat(reshape(template(:, 1:2), 1, nRows, 2), size(epochs, 1), 1, 1), ...
        reshape(weights * reshape(monthly, numel(files), []), size(epochs, 1), nRows, 2));
end
