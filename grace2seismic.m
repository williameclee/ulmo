%% GRACE2SEISMIC - Estimates independent earthquake corrections from GRACE.
% Fits each event separately in its own land-and-ocean cap. Static fields
% are cached independently and combined in spherical harmonic space using
% the fraction of each actual GRACE window after the rounded event day.
%
% Syntax
%   [seisPlmt, seisStdPlmt, dates, eventPlm] = ...
%       grace2seismic(product, earthquakes, radius, eventDates, L)
%   [__] = grace2seismic(__, RemoveGLDAS=true, SaveData=false)
%
% Input arguments
%   product - GRACE {centre, release, degree} cell, as in GRACE2PLMT_NEW.
%   earthquakes - Name, name vector, (nEvents, 2) [longitude, latitude]
%       matrix, or cell vector mixing scalar names and coordinate pairs.
%       Names: Sumatra (2004), Tohoku (2011), Chile (2010), and aliases
%       documented in docs/grace-seismic.md. Coordinates are in degrees.
%   radius - Positive cap radius, at most 180 degrees; scalar or event vector.
%   eventDates - Datetime or MATLAB datenum vector, one time per event.
%       Rounded to the nearest UTC midnight; noon rounds to the next day.
%       Supplied dates are authoritative even for named events.
%   L - Maximum localisation/output degree, at least 2 and at most product degree.
%       The Slepian basis uses degrees 2:L; degrees 0 and 1 are excluded.
%   Unit (name-value) - GRACE coefficient unit: SD, POT or GRAV.
%       Default unit: SD (kg/m^2).
%   OutputFormat (name-value) - timefirst or traditional.
%       Default format: timefirst.
%   TimeFormat (name-value) - datetime or datenum for returned dates.
%       Default format: datenum, as in GRACE2PLMT_NEW.
%   TimeRange (name-value) - Inclusive returned-midpoint range. Does not
%       restrict observations used to estimate an event.
%       Default range: [] (all available GRACE observations).
%   RemoveGLDAS (name-value) - Whether to subtract Noah 2.1 monthly storage
%       before fitting, using overlap averages over actual GRACE windows.
%       Default flag: false. See docs/grace-seismic.md for input data.
%   truncation (name-value) - Positive retained mode count or [].
%       Default count: [], using max(1, round(N)), bounded by basis size.
%   polynomialDegree (name-value) - Nonnegative polynomial degree.
%       Default degree: 2.
%   periods (name-value) - Positive harmonic periods in numeric years or durations.
%       Default periods: years([1, 0.5]); [] omits harmonics.
%   FitRange (name-value) - Optional fit interval for uncached experiments.
%       Default range: [] (all available GRACE observations).
%   LoveNumSource (name-value) - Load Love-number source for GRACE and GLDAS.
%       Default source: ISSM.
%   Deg1Correction, C20Correction, C30Correction (name-value) - Whether to
%       apply GRACE replacements. Default flags: [] (GRACE loader defaults).
%   RecomputeDegree1 (name-value) - GRACE low-degree recomputation options.
%       Default option: false. Forwarded to GRACE2PLMT_NEW.
%   ForceNew (name-value) - Whether to bypass event and processed GRACE caches.
%       Default flag: false. Does not refresh reusable geometric caches.
%   SaveData (name-value) - Whether to save event and supported upstream caches.
%       Default flag: true.
%   BeQuiet (name-value) - Whether to suppress progress messages.
%       Default flag: false.
%
% Output arguments
%   seisPlmt, seisStdPlmt - Combined correction and standard deviation in
%       lmcosi layout, with degree/order labels in both arrays.
%       Size: (nDates, addmup(L), 4), or (addmup(L), 4, nDates) for traditional.
%       Unit: selected Unit for coefficient columns.
%       Degrees 0 and 1 have zero coefficients and standard deviations.
%   dates - Actual GRACE interval midpoints. Size: (nDates, 1).
%   eventPlm - Static jumps in input event order and the same Unit.
%       Size: (nEvents, addmup(L), 4), or (addmup(L), 4, nEvents).
%
% Output file content
%   degree, order, cs, time, coefficients, uncertainty - Standard WRITESHNC
%       variables; time has one entry, the rounded event day. Uncertainty
%       stores standard deviations. Attributes record epicentre, radius,
%       retainedModeCount, shannonNumber, background, fit dates and GLDAS.
%
% Notes
%   The default fit is ordinary least squares over the full available
%   record. Errors assume independent Slepian modes and events; GLDAS is
%   fixed. No extra physical or basis-normalisation factor is applied.
%   Nondefault background, fit range, truncation, Love numbers or low-degree
%   corrections require ForceNew=true and SaveData=false, because the short
%   cache names deliberately omit those settings. Users clear outdated caches.
%   PLM2ROT's DLMB dependency has no save control. See docs/grace-seismic.md
%   for filename rules, catalogue sources and GLDAS coverage requirements.
%
% See also
%   slept2seismic, grace2plmt_new, graceepochs, seismiccapbasis, writeshnc, readshnc
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/13, En-Chi Lee (williameclee@arizona.edu)

function [seisPlmt, seisStdPlmt, dates, eventPlm] = ...
        grace2seismic(product, earthquakes, radius, eventDates, L, options)
    arguments (Input)
        product (1, 3) cell
        earthquakes
        radius double {mustBeReal, mustBeFinite, mustBePositive, mustBeLessThanOrEqual(radius, 180)}
        eventDates {mustBeA(eventDates, {'datetime', 'numeric'}), mustBeVector, mustBeNonempty}
        L (1, 1) double {mustBeFinite, mustBeInteger, mustBeGreaterThanOrEqual(L, 2)}
        options.Unit (1, 1) string {mustBeMember(options.Unit, ["SD", "POT", "GRAV"])} = "SD"
        options.OutputFormat (1, 1) string {mustBeMember(options.OutputFormat, ["timefirst", "traditional"])} = "timefirst"
        options.TimeFormat (1, 1) string {mustBeMember(options.TimeFormat, ["datetime", "datenum"])} = "datenum"
        options.TimeRange {mustBeTimeRange} = []
        options.RemoveGLDAS (1, 1) logical = false
        options.truncation double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive, mustBeFinite} = []
        options.polynomialDegree (1, 1) double {mustBeInteger, mustBeNonnegative, mustBeFinite} = 2
        options.periods {mustBeA(options.periods, {'numeric', 'duration'})} = years([1, 0.5])
        options.FitRange {mustBeTimeRange} = []
        options.LoveNumSource (1, 1) string = "ISSM"
        options.Deg1Correction {mustBeScalarOrEmpty, mustBeNumericOrLogical} = []
        options.C20Correction {mustBeScalarOrEmpty, mustBeNumericOrLogical} = []
        options.C30Correction {mustBeScalarOrEmpty, mustBeNumericOrLogical} = []
        options.RecomputeDegree1 = false
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) logical = false
    end
    product = normaliseproduct(product);
    assert(L <= product{3}, 'ULMO:grace2seismic:Bandwidth', 'L must not exceed the GRACE product degree.');
    [epicentres, names, expectedDates] = resolveearthquakes(earthquakes);
    nEvents = size(epicentres, 1);
    assert(numel(eventDates) == nEvents && isvector(radius) && ...
        (isscalar(radius) || numel(radius) == nEvents), ...
        'ULMO:grace2seismic:EventCount', 'Supply one date and one radius per event, or a shared scalar radius.');
    eventDates = convertdates(eventDates(:));
    assert(all(isfinite(eventDates)), 'ULMO:grace2seismic:InvalidDates', 'Event times must be finite.');
    eventDates = dateshift(eventDates + hours(12), 'start', 'day');
    if isscalar(radius), radius = repmat(radius, nEvents, 1); end
    periodYears = options.periods;
    if isduration(periodYears), periodYears = years(periodYears); end
    assert(isreal(periodYears) && (isempty(periodYears) || isvector(periodYears)) && ...
        all(isfinite(periodYears) & periodYears > 0, 'all'), ...
        'ULMO:grace2seismic:Periods', 'Periods must be finite positive years or durations.');
    defaultCorrections = strcmp(product{2}, 'RL06');
    corrections = {options.Deg1Correction, options.C20Correction, options.C30Correction};
    for iCorrection = 1:3
        value = corrections{iCorrection};
        assert(isempty(value) || ismember(value, [0, 1]), ...
            'ULMO:grace2seismic:Correction', 'Correction flags must be [], false or true.');
        if isempty(value), value = defaultCorrections; end
        corrections{iCorrection} = logical(value);
    end
    nondefault = ~isempty(options.truncation) || options.polynomialDegree ~= 2 || ...
        ~isequal(periodYears(:).', [1, 0.5]) || ~isempty(options.FitRange) || ...
        ~strcmpi(options.LoveNumSource, 'ISSM') || ...
        any([corrections{:}] ~= defaultCorrections) || ...
        ~(isempty(options.RecomputeDegree1) || isequal(options.RecomputeDegree1, false));
    assert(~nondefault || (options.ForceNew && ~options.SaveData), ...
        'ULMO:grace2seismic:NondefaultCache', ...
        'Nondefault fit, truncation or correction settings require ForceNew=true and SaveData=false.');

    [epochs, synthetic, equatorRadius] = graceepochs(product);
    epochs = epochs(~synthetic, :);
    dates = mean(epochs, 2);
    assert(all(diff(dates) > seconds(0)) && all(epochs(:, 2) > epochs(:, 1)), ...
        'ULMO:grace2seismic:InvalidWindows', 'GRACE windows must have unique midpoints and positive duration.');
    [~, ~, ~, template, ~, ~, ~, ~, ~, mapping] = addmon(L);
    [~, degrees] = addmout(L);
    nHarmonics = (L + 1) ^ 2;
    amplitudes = zeros(nEvents, nHarmonics);
    variances = zeros(nEvents, nHarmonics);
    source = [];
    for iEvent = 1:nEvents
        cache = cachepath(product, names(iEvent), radius(iEvent), eventDates(iEvent), ...
            expectedDates(iEvent), L, options.Unit, options.RemoveGLDAS);
        required = {'degree', 'order', 'cs', 'time', 'coefficients', 'uncertainty'};
        complete = false;
        if ~options.ForceNew && isfile(cache)
            info = ncinfo(cache);
            complete = all(ismember(required, {info.Variables.Name}));
        end
        if complete
            [field, sigma] = readshnc(cache);
            assert(isequal(size(field), [size(template, 1), 4]), ...
                'ULMO:grace2seismic:CacheShape', 'The cache must contain one degree-L event field.');
        else
            if isempty(source)
                redo = recomputeoptions(options.RecomputeDegree1, options.SaveData);
                [plmt, ~, loadedDates] = grace2plmt_new(product, ...
                    'Unit', char(options.Unit), 'LoveNumSource', char(options.LoveNumSource), ...
                    'OutputFormat', 'timefirst', 'TimeFormat', 'datetime', ...
                    'Deg1Correction', corrections{1}, 'C20Correction', corrections{2}, ...
                    'C30Correction', corrections{3}, 'RecomputeDegree1', redo, ...
                    'ForceNew', options.ForceNew, 'SaveData', options.SaveData, 'BeQuiet', options.BeQuiet);
                [matched, locations] = ismember(dates, loadedDates);
                assert(all(matched) && numel(loadedDates) == numel(dates), ...
                    'ULMO:grace2seismic:Dates', 'GRACE coefficients must match every actual solution window.');
                assert(size(plmt, 2) >= size(template, 1), ...
                    'ULMO:grace2seismic:Bandwidth', 'The loaded GRACE record does not cover L.');
                plmt = plmt(locations, 1:size(template, 1), :);
                source = reshape(plmt(:, :, 3:4), numel(dates), []);
                source = source(:, mapping);
                if options.RemoveGLDAS
                    storage = seismicgldas(epochs, L, options.Unit, options.LoveNumSource, equatorRadius);
                    storage = reshape(storage(:, :, 3:4), numel(dates), []);
                    source = source - storage(:, mapping);
                end
                source(:, degrees < 2) = 0;
                assert(all(isfinite(source), 'all'), 'ULMO:grace2seismic:MissingCoefficients', ...
                    'GRACE degrees 2:L must be finite at every observation.');
            end
            GVN = seismiccapbasis(epicentres(iEvent, :), radius(iEvent), L, ...
                options.truncation, SaveData = options.SaveData);
            G = GVN{1};
            fitRange = options.FitRange;
            if isempty(fitRange), fitRange = dates([1, end]).'; else, fitRange = convertdates(fitRange); end
            [step, stepStd] = slept2seismic(source * G, dates, eventDates(iEvent), [], ...
                polynomialDegree = options.polynomialDegree, periods = options.periods, FitRange = fitRange);
            field = packfield(step * G.', template, mapping);
            sigma = packfield(sqrt(stepStd .^ 2 * (G .^ 2).'), template, mapping);
            if options.SaveData
                if ~isfolder(fileparts(cache)), mkdir(fileparts(cache)); end
                attributes = {'epicentre', epicentres(iEvent, :), 'radius', radius(iEvent), ...
                    'retainedModeCount', size(G, 2), 'shannonNumber', GVN{3}, ...
                    'polynomialDegree', options.polynomialDegree, 'periodsYears', periodYears, ...
                    'fitStart', char(string(dates(1), 'yyyy-MM-dd HH:mm:ss')), ...
                    'fitEnd', char(string(dates(end), 'yyyy-MM-dd HH:mm:ss')), ...
                    'loveNumberSource', char(options.LoveNumSource), ...
                    'removeGLDAS', double(options.RemoveGLDAS), 'gldasProduct', 'none', ...
                    'description', 'Static independent event jump; uncertainty is standard deviation'};
                if options.RemoveGLDAS
                    attributes{end - 2} = 'GLDAS_NOAH10_M.2.1; four soil layers + canopy + snow; monthly overlap means';
                end
                writeshnc(cache, field, sigma, eventDates(iEvent), char(options.Unit), Attributes = attributes);
                if ~options.BeQuiet, fprintf('[ULMO>grace2seismic] Saved %s\n', cache); end
            end
        end
        coefficients = field(:, 3:4);
        uncertainties = sigma(:, 3:4);
        amplitudes(iEvent, :) = coefficients(mapping);
        variances(iEvent, :) = uncertainties(mapping) .^ 2;
    end
    weights = max(0, seconds(epochs(:, 2) - max(epochs(:, 1), eventDates.'))) ./ ...
        seconds(epochs(:, 2) - epochs(:, 1));
    keep = true(numel(dates), 1);
    if ~isempty(options.TimeRange)
        range = convertdates(options.TimeRange);
        keep = dates >= range(1) & dates <= range(2);
    end
    seisPlmt = packseries(weights(keep, :) * amplitudes, template, mapping);
    seisStdPlmt = packseries(sqrt(weights(keep, :) .^ 2 * variances), template, mapping);
    eventPlm = packseries(amplitudes, template, mapping);
    dates = dates(keep);
    if options.OutputFormat == "traditional"
        seisPlmt = permute(seisPlmt, [2, 3, 1]);
        seisStdPlmt = permute(seisStdPlmt, [2, 3, 1]);
        eventPlm = permute(eventPlm, [2, 3, 1]);
    end
    if options.TimeFormat == "datenum", dates = datenum(dates); end %#ok<DATNM>
end

function dates = convertdates(dates)
    if isnumeric(dates)
        mustBeReal(dates);
        dates = datetime(double(dates), 'ConvertFrom', 'datenum');
    elseif ~isempty(dates.TimeZone)
        dates.TimeZone = 'UTC';
        dates.TimeZone = '';
    end
end

function redo = recomputeoptions(redo, saveData)
    if isequal(redo, true) || ischar(redo) || (isstring(redo) && isscalar(redo))
        model = 'ice6gd';
        if ischar(redo) || isstring(redo), model = char(redo); end
        redo = {60, [], 180, model, GeoDomain('alloceans', 'Buffer', 0.5)};
    end
    if iscell(redo) && ~saveData
        redo = [reshape(redo, 1, []), {'SaveData', false}];
    end
end

function product = normaliseproduct(product)
    product{1} = validatestring(product{1}, {'CSR', 'GFZ', 'JPL', 'CSR mascon'});
    if isnumeric(product{2})
        validateattributes(product{2}, {'numeric'}, {'scalar', 'integer', 'finite'});
        product{2} = sprintf('RL%02d', product{2});
    end
    product{2} = validatestring(product{2}, {'RL05', 'RL06'});
    validateattributes(product{3}, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
    if strcmp(product{1}, 'CSR mascon')
        assert(strcmp(product{2}, 'RL06'), 'ULMO:grace2seismic:Product', 'CSR mascons require RL06.');
        product{3} = 720;
    elseif strcmp(product{2}, 'RL05')
        if strcmp(product{1}, 'CSR'), product{3} = 60; else, product{3} = 90; end
    else
        assert(ismember(product{3}, [60, 96]), 'ULMO:grace2seismic:Product', 'RL06 GSM supports degree 60 or 96.');
    end
end

function path = cachepath(product, name, radius, date, expected, L, unit, gldas)
    root = getenv('GRACEDATA');
    if isempty(root), root = fullfile(getenv('IFILES'), 'GRACE'); end
    tag = string(replace(product{1}, ' ', '')) + extractAfter(string(product{2}), 2) + product{3};
    stem = "grace2seismic-" + tag + "-" + name + "_" + num2filename(radius);
    if isnat(expected) || date ~= expected, stem = stem + "-T" + string(date, 'yyyyMMdd'); end
    if L ~= product{3}, stem = stem + "-L" + L; end
    stem = stem + "-" + unit;
    if gldas, stem = stem + "-GLDAS"; end
    path = fullfile(root, 'Seismic', stem + ".nc");
end

function field = packfield(values, template, mapping)
    coefficients = zeros(size(template, 1), 2);
    coefficients(mapping) = values;
    field = [template(:, 1:2), coefficients];
end

function series = packseries(values, template, mapping)
    nRows = size(template, 1);
    coefficients = zeros(size(values, 1), nRows * 2);
    coefficients(:, mapping) = values;
    series = cat(3, repmat(reshape(template(:, 1:2), 1, nRows, 2), size(values, 1), 1, 1), ...
        reshape(coefficients, size(values, 1), nRows, 2));
end
