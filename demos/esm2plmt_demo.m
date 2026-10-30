%% ESM2PLMT_DEMO - Processes and plots ESM components for one GRACE epoch.
% Saves separate H, I and O fields in native-degree GRAV and SD units, plus
% native and truncated HI/HIO sums. Generates a two-column component map.
%
% Syntax
%   products = esm2plmt_demo(graceFile, outputFolder)
%   products = esm2plmt_demo(__, "Loutput", 45, "MakeFigure", false)
%
% Input arguments
%   graceFile - GRACE file readable by PARSEGRACESOURCEFILE.
%   outputFolder - Destination for MAT files and the figure. Created as needed.
%   Loutput (name-value) - Truncated output degree, no greater than the source.
%       Default degree: 60.
%   MakeFigure (name-value) - Whether to save the component comparison figure.
%       Default option: true.
%
% Output arguments
%   products - Saved product paths for H, I, O, HI and HIO. Individual
%       components have nativeGRAV, nativeSD and truncatedSD fields; sums
%       have nativeSD and truncatedSD. Figure path is returned when generated.
%
% Output file content
%   esmPlmt, esmStdPlmt, dates - See ESM2PLMT.
%
% Notes
%   HIO is not GSM: it excludes atmosphere and solid Earth and has no DEAL
%   subtraction. One epoch demonstrates ingestion, not trend recovery.
%   Requires the standard ULMO paths, Earth constants and ISSM Love numbers.
%   Exact GRACE-window means are formed from cached raw samples. These can
%   differ from the calendar-month means returned by ESM2PLMT(component, "midmonth").
%   ESM inputs are read from $IFILES/ESM/ESAESM2/raw/COMPONENT.
%
% See also
%   ESM2PLMT, PLOT_ESM_EPOCH, PARSEGRACESOURCEFILE
%
% Created by
%   2026/10/29
%
% Last modified
%   2026/10/30

function products = esm2plmt_demo(graceFile, outputFolder, options)
    arguments (Input)
        graceFile (1, 1) string {mustBeFile}
        outputFolder (1, 1) string
        options.Loutput (1, 1) double {mustBeInteger, mustBePositive} = 60
        options.MakeFigure (1, 1) logical = true
    end
    [~, ~, graceDate, intervals] = parsegracesourcefile(graceFile);
    assert(all(mod(hour(intervals), 6) == 0 & minute(intervals) == 0 & ...
        second(intervals) == 0), 'ULMO:esm2plmt:InvalidIntervals', ...
        'The GRACE interval must align with six-hour ESM sampling.');
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    stamp = string(intervals, 'yyyyMMdd_HH');
    prefix = "ESAESM2_";
    suffix = "_" + stamp(1) + "_" + stamp(2);
    mass = struct;
    products = struct;
    for component = ["H", "I", "O"]
        [gravity, uncertainty, dates] = epochmean(intervals, graceDate, ...
            'Component', component, 'Unit', 'GRAV');
        assert(isequal(dates, datenum(graceDate)), 'ESM and GRACE epoch dates differ.');
        assert(options.Loutput <= gravity(1, end, 1), ...
            'Demo output degree must not exceed the source degree.');
        stem = fullfile(outputFolder, prefix + component + suffix);
        products.(component).nativeGRAV = savefield(stem, gravity, uncertainty, dates, "GRAV");
        [density, uncertainty, dates] = epochmean(intervals, graceDate, ...
            'Component', component, 'Unit', 'SD');
        mass.(component) = density;
        products.(component).nativeSD = savefield(stem, density, uncertainty, dates, "SD");
        products.(component).truncatedSD = savetruncated(stem, density, uncertainty, ...
            dates, options.Loutput);
    end
    for component = ["HI", "HIO"]
        density = mass.H;
        density(:, :, 3:4) = mass.H(:, :, 3:4) + mass.I(:, :, 3:4);
        if component == "HIO"
            density(:, :, 3:4) = density(:, :, 3:4) + mass.O(:, :, 3:4);
        end
        uncertainty = density;
        uncertainty(:, :, 3:4) = NaN;
        stem = fullfile(outputFolder, prefix + component + suffix);
        products.(component).nativeSD = savefield(stem, density, uncertainty, dates, "SD");
        products.(component).truncatedSD = savetruncated(stem, density, uncertainty, ...
            dates, options.Loutput);
    end
    if options.MakeFigure
        products.figure = fullfile(outputFolder, "ESAESM2" + suffix + ".png");
        plot_esm_epoch(products, products.figure);
    end
end

function file = savefield(stem, esmPlmt, esmStdPlmt, dates, unit)
    file = stem + "_L" + esmPlmt(1, end, 1) + "_" + unit + ".mat";
    save(file, 'esmPlmt', 'esmStdPlmt', 'dates');
end

function file = savetruncated(stem, field, uncertainty, dates, degree)
    nRows = (degree + 1) * (degree + 2) / 2;
    file = savefield(stem, field(:, 1:nRows, :), uncertainty(:, 1:nRows, :), dates, "SD");
end

function [plmt, uncertainty, dates] = epochmean(interval, midpoint, options)
    arguments (Input)
        interval (1, 2) datetime
        midpoint (1, 1) datetime
        options.Component (1, 1) string
        options.Unit (1, 1) string
    end
    [raw, ~, rawDates] = esm2plmt(options.Component, [], interval, ...
        'Unit', options.Unit, 'TimeFormat', 'datetime');
    expected = (interval(1):hours(6):interval(2) - hours(6))';
    assert(isequal(rawDates(rawDates < interval(2)), expected), ...
        'ULMO:esm2plmt:IncompleteSampling', ...
        'The cached raw data do not cover the complete GRACE interval.');
    [plmt, uncertainty, dates] = averageesmepochs(rawDates, raw, interval);
    assert(isequal(dates, midpoint), 'ESM and GRACE epoch dates differ.');
    dates = datenum(dates); %#ok<DATNM>
end
