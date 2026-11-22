%% SEISMICGLDAS - Averages Noah monthly storage harmonics over GRACE windows.
%
% Syntax
%   plmt = seismicgldas(epochs, L, unit, loveNumSource, equatorRadius)
%
% Input arguments
%   epochs - Actual GRACE windows as UTC datetimes. Size: (nDates, 2).
%   L - Maximum harmonic degree, from 0 to 179 for the 1-degree input grid.
%   unit - SD, POT or GRAV, matching the selected GRACE record.
%   loveNumSource - Source passed to CONVERTGRAVITY.
%   equatorRadius - GRACE reference equatorial radius in metres.
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
%   overlap duration. Every part of every GRACE window must be covered.
%   The six fields must share a time-invariant native missing-data mask;
%   that mask and latitudes south of 60 S contribute zero. No extra glacier
%   mask or groundwater term is added. No files are written.
%   Cell-constant fields are integrated with exact longitude integrals and
%   eight-point Gauss-Legendre latitude quadrature. Real harmonics have
%   squared norm 4*pi, matching the Slepian package.
%
% Data sources
%   NASA GLDAS: https://ldas.gsfc.nasa.gov/gldas/model-output
%   https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH10_M.2.1/
%
% See also
%   grace2seismic, convertgravity, gausslegendrecof
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function plmt = seismicgldas(epochs, L, unit, loveNumSource, equatorRadius)
    arguments (Input)
        epochs (:, 2) datetime
        L (1, 1) double {mustBeInteger, mustBeNonnegative, mustBeLessThan(L, 180)}
        unit {mustBeTextScalar}
        loveNumSource {mustBeTextScalar}
        equatorRadius (1, 1) double {mustBeFinite, mustBePositive}
    end
    root = getenv('GLDASDATA');
    if isempty(root), root = fullfile(getenv('IFILES'), 'EARTHMODELS', 'GLDAS2.1'); end
    files = dir(fullfile(root, '**', 'GLDAS_NOAH10_M.A*.021.nc4'));
    assert(~isempty(files), 'ULMO:seismicgldas:MissingData', ...
        'No GLDAS_NOAH10_M 2.1 monthly files found under %s.', root);
    months = NaT(numel(files), 1);
    for iFile = 1:numel(files)
        token = regexp(files(iFile).name, '^GLDAS_NOAH10_M\.A(\d{6})\.021\.nc4$', 'tokens', 'once');
        assert(~isempty(token), 'ULMO:seismicgldas:Product', 'Unexpected GLDAS filename.');
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
    assert(all(diff(months) > seconds(0)), 'ULMO:seismicgldas:DuplicateMonth', ...
        'GLDAS must contain only one file per required calendar month.');
    lengths = seconds(epochs(:, 2) - epochs(:, 1));
    assert(all(abs(sum(overlap, 2) - lengths) < 1e-6), ...
        'ULMO:seismicgldas:Coverage', 'GLDAS does not cover every actual GRACE observation window.');
    weights = overlap ./ lengths;
    [~, ~, ~, template] = addmon(L);
    nRows = size(template, 1);
    monthly = zeros(numel(files), nRows, 2);
    variables = {'SoilMoi0_10cm_inst', 'SoilMoi10_40cm_inst', ...
        'SoilMoi40_100cm_inst', 'SoilMoi100_200cm_inst', 'CanopInt_inst', 'SWE_inst'};
    nativeMask = [];
    latitudeIntegral = [];
    for iMonth = 1:numel(files)
        path = fullfile(files(iMonth).folder, files(iMonth).name);
        lon = double(ncread(path, 'lon'));
        lat = double(ncread(path, 'lat'));
        assert(numel(lon) == 360 && numel(lat) == 150 && ...
            max(abs(lon(:) - (-179.5:179.5).')) < 1e-6 && ...
            max(abs(lat(:) - (-59.5:89.5).')) < 1e-6, ...
            'ULMO:seismicgldas:Grid', 'GLDAS requires the complete native 1-degree grid.');
        field = zeros(360, 150);
        for iVariable = 1:numel(variables)
            name = variables{iVariable};
            info = ncinfo(path, name);
            dimensions = string({info.Dimensions.Name});
            assert(isequal(dimensions, ["lon", "lat", "time"]) && ...
                isequal(info.Size, [360, 150, 1]), ...
                'ULMO:seismicgldas:Dimensions', '%s must have lon, lat and one monthly time.', name);
            units = regexprep(string(ncreadatt(path, name, 'units')), '\s', '');
            assert(ismember(units, ["kgm-2", "kg/m^2", "kgm^-2"]), ...
                'ULMO:seismicgldas:Unit', '%s must use kg/m^2.', name);
            values = double(ncread(path, name));
            names = {info.Attributes.Name};
            if ismember('missing_value', names)
                missingValue = double(ncreadatt(path, name, 'missing_value'));
                values(values == missingValue) = NaN;
            end
            missing = ~isfinite(values);
            if isempty(nativeMask), nativeMask = missing; end
            assert(isequal(missing, nativeMask) && any(~missing, 'all') && ~any(isinf(values), 'all'), ...
                'ULMO:seismicgldas:MissingValues', ...
                'GLDAS storage components must share a time-invariant native missing-data mask.');
            values(missing) = 0;
            field = field + values;
        end
        if isempty(latitudeIntegral)
            latitudeIntegral = latitudeintegrals(lat, template(:, 1:2));
        end
        % The Fourier phase accounts for the native -179.5-degree cell centre.
        orders = (0:L).';
        longitudeFactor = deg2rad(1) * sinc(orders / 360) .* exp(-1i * orders * deg2rad(lon(1)));
        fourier = fft(field.', [], 2);
        fourier = fourier(:, 1:L + 1) .* longitudeFactor.';
        harmonics = template;
        for iOrder = 0:L
            selected = template(:, 2) == iOrder;
            integral = latitudeIntegral(:, selected).' * fourier(:, iOrder + 1) / (4 * pi);
            harmonics(selected, 3:4) = [real(integral), -imag(integral)];
        end
        harmonics(1, 3:4) = 0;
        harmonics = convertgravity(harmonics, 'SD', char(unit), ...
            'InputFormat', 'lmcosi', 'EquatorRadius', equatorRadius, 'LoveNumSource', char(loveNumSource));
        monthly(iMonth, :, :) = harmonics(:, 3:4);
    end
    plmt = cat(3, repmat(reshape(template(:, 1:2), 1, nRows, 2), size(epochs, 1), 1, 1), ...
        reshape(weights * reshape(monthly, numel(files), []), size(epochs, 1), nRows, 2));
end

function integrals = latitudeintegrals(lat, labels)
    [weights, nodes] = gausslegendrecof(15, [], [-1, 1]);
    phi = deg2rad(lat(:)) + deg2rad(0.5) * nodes(:).';
    measure = deg2rad(0.5) * cos(phi) .* weights(:).';
    integrals = zeros(numel(lat), size(labels, 1));
    for degree = 0:max(labels(:, 1))
        polynomials = legendre(degree, sin(phi(:)).', 'sch') * sqrt(2 * degree + 1);
        for order = 0:degree
            values = reshape(polynomials(order + 1, :), size(phi));
            integrals(:, labels(:, 1) == degree & labels(:, 2) == order) = sum(values .* measure, 2);
        end
    end
end
