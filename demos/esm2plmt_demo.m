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
%   differ from the calendar-month means returned by ESM2PLMT(component, [], "midmonth").
%   Formal errors are propagated assuming independent samples/components;
%   missing formal errors use the empirical standard error of each mean.
%   HI/HIO empirical errors use the summed samples, retaining their observed
%   cross-component covariance. No temporal correlation correction is applied.
%   Prepared inputs are $IFILES/ESM/ESAESM2/ESAESM2_COMPONENT.mat.
%   Use PREPAREESM separately to import raw files from any source folder.
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
        graceFile {mustBeTextScalar, mustBeFile}
        outputFolder {mustBeTextScalar}
        options.Loutput (1, 1) double {mustBeInteger, mustBePositive} = 60
        options.MakeFigure (1, 1) logical = true
    end
    [~, ~, graceDate, intervals] = parsegracesourcefile(graceFile);
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    stamp = string(intervals, 'yyyyMMdd_HH');
    prefix = "ESAESM2_";
    suffix = "_" + stamp(1) + "_" + stamp(2);
    products = struct;
    for component = ["H", "I", "O"]
        [gravity, uncertainty, dates] = epochmean(intervals, graceDate, ...
            component, 'GRAV');
        assert(isequal(dates, datenum(graceDate)), 'ESM and GRACE epoch dates differ.');
        assert(options.Loutput <= gravity(1, end, 1), ...
            'Demo output degree must not exceed the source degree.');
        stem = fullfile(outputFolder, prefix + component + suffix);
        products.(component).nativeGRAV = savefield(stem, gravity, uncertainty, dates, "GRAV");
        [density, uncertainty, dates] = epochmean(intervals, graceDate, ...
            component, 'SD');
        products.(component).nativeSD = savefield(stem, density, uncertainty, dates, "SD");
        products.(component).truncatedSD = savetruncated(stem, density, uncertainty, ...
            dates, options.Loutput);
    end
    for component = ["HI", "HIO"]
        components = ["H", "I"];
        if component == "HIO"
            components = ["H", "I", "O"];
        end
        [density, uncertainty, dates] = epochmean(intervals, graceDate, components, 'SD');
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

function [plmt, uncertainty, dates] = epochmean(interval, midpoint, components, unit)
    arguments (Input)
        interval (1, 2) datetime
        midpoint (1, 1) datetime
        components {mustBeText}
        unit {mustBeTextScalar}
    end
    for component = string(components)
        [field, fieldStd, fieldDates] = esm2plmt(component, interval, '', ...
            'Unit', unit, 'TimeFormat', 'datetime');
        if component == string(components(1))
            raw = field;
            rawStd = fieldStd;
            rawDates = fieldDates;
        else
            assert(isequal(rawDates, fieldDates), 'ESM component sample times differ.');
            raw(:, :, 3:4) = raw(:, :, 3:4) + field(:, :, 3:4);
            rawStd(:, :, 3:4) = hypot(rawStd(:, :, 3:4), fieldStd(:, :, 3:4));
        end
    end
    % Sum samples before estimating empirical error, retaining their observed
    % cross-component covariance. Supplied formal errors assume independence.
    [coefficients, dates, coefficientStd] = averageepochs(rawDates, ...
        permute(raw(:, :, 3:4), [2, 3, 1]), interval, hours(6), ...
        'CheckCompleteness', true, 'SamplingOrigin', datetime(2000, 1, 1), ...
        'SampleStd', permute(rawStd(:, :, 3:4), [2, 3, 1]));
    assert(all(isfinite(coefficients), 'all'), 'ULMO:esm2plmt:IncompleteSampling', ...
        'The cached raw data do not cover the complete GRACE interval.');
    plmt = raw(1, :, :);
    uncertainty = rawStd(1, :, :);
    plmt(:, :, 3:4) = permute(coefficients, [3, 1, 2]);
    uncertainty(:, :, 3:4) = permute(coefficientStd, [3, 1, 2]);
    assert(isequal(dates, midpoint), 'ESM and GRACE epoch dates differ.');
    dates = datenum(dates); %#ok<DATNM>
end
