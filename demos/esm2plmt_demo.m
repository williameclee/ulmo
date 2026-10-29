%% ESM2PLMT_DEMO - Processes and plots ESM components for one GRACE epoch.
% Saves separate H, I and O fields in native-degree GRAV and SD units, plus
% native and truncated HI/HIO sums. Generates a two-column component map.
%
% Syntax
%   products = esm2plmt_demo(graceFile, esmFolder, outputFolder)
%   products = esm2plmt_demo(__, "Loutput", 45, "MakeFigure", false)
%
% Input arguments
%   graceFile - GRACE file readable by PARSEGRACESOURCEFILE.
%   esmFolder - Folder containing the extracted six-hour H, I and O files.
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
%   esmPlmt, esmStdPlmt, dates, metadata - See ESM2PLMT.
%
% Notes
%   HIO is not GSM: it excludes atmosphere and solid Earth and has no DEAL
%   subtraction. One epoch demonstrates ingestion, not trend recovery.
%   Requires the standard ULMO paths, Earth constants and ISSM Love numbers.
%
% See also
%   ESM2PLMT, PLOT_ESM_EPOCH, PARSEGRACESOURCEFILE
%
% Created by
%   2026/10/29

function products = esm2plmt_demo(graceFile, esmFolder, outputFolder, options)
    arguments (Input)
        graceFile (1, 1) string {mustBeFile}
        esmFolder (1, 1) string {mustBeFolder}
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
    [~, graceName, extension] = fileparts(graceFile);
    mass = struct;
    info = struct;
    products = struct;
    for component = ["H", "I", "O"]
        [gravity, uncertainty, dates, metadata] = esm2plmt(esmFolder, intervals, ...
            'Component', component, 'Unit', 'GRAV');
        assert(isequal(dates, datenum(graceDate)), 'ESM and GRACE epoch dates differ.');
        assert(options.Loutput <= metadata.sourceDegree, ...
            'Demo output degree must not exceed the source degree.');
        metadata.graceSource = graceName + extension;
        stem = fullfile(outputFolder, prefix + component + suffix);
        products.(component).nativeGRAV = savefield(stem, gravity, uncertainty, dates, metadata);
        density = convertgravity(gravity, 'GRAV', 'SD', ...
            metadata.gravityParam, metadata.equatorRadius, 'LoveNumSource', 'ISSM');
        metadata.unit = "SD";
        mass.(component) = density;
        info.(component) = metadata;
        products.(component).nativeSD = savefield(stem, density, uncertainty, dates, metadata);
        products.(component).truncatedSD = savetruncated(stem, density, uncertainty, ...
            dates, metadata, options.Loutput);
    end
    for component = ["HI", "HIO"]
        density = mass.H;
        density(:, :, 3:4) = mass.H(:, :, 3:4) + mass.I(:, :, 3:4);
        inputs = [products.H.nativeSD, products.I.nativeSD];
        if component == "HIO"
            density(:, :, 3:4) = density(:, :, 3:4) + mass.O(:, :, 3:4);
            inputs(end + 1) = products.O.nativeSD;
        end
        metadata = rmfield(info.H, {'files', 'sampleDates'});
        metadata.component = component;
        metadata.componentFiles = inputs;
        uncertainty = density;
        uncertainty(:, :, 3:4) = NaN;
        stem = fullfile(outputFolder, prefix + component + suffix);
        products.(component).nativeSD = savefield(stem, density, uncertainty, dates, metadata);
        products.(component).truncatedSD = savetruncated(stem, density, uncertainty, ...
            dates, metadata, options.Loutput);
    end
    if options.MakeFigure
        products.figure = fullfile(outputFolder, "ESAESM2" + suffix + ".png");
        plot_esm_epoch(products, products.figure);
    end
end

function file = savefield(stem, esmPlmt, esmStdPlmt, dates, metadata)
    file = stem + "_L" + metadata.outputDegree + "_" + metadata.unit + ".mat";
    save(file, 'esmPlmt', 'esmStdPlmt', 'dates', 'metadata');
end

function file = savetruncated(stem, field, uncertainty, dates, metadata, degree)
    nRows = (degree + 1) * (degree + 2) / 2;
    metadata.outputDegree = degree;
    file = savefield(stem, field(:, 1:nRows, :), uncertainty(:, 1:nRows, :), dates, metadata);
end
