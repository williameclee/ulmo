%% AGGREGATESTERICSEALEVEL - Aggregates steric sea level data from multiple files
%
% Last modified
%   2026/02/17, williameclee@arizona.edu (@williameclee)
%     - Extracted from PROCESSSTERICDATAEN4 for reusability

function aggregateStericSeaLevel(inputFolder, inputFiles, outputPath, options)

    arguments (Input)
        inputFolder {mustBeTextScalar, mustBeFolder}
        inputFiles (1, :) cell
        outputPath {mustBeTextScalar}
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    [~, ~, outputExt] = fileparts(outputPath);

    matStericVars = ...
        {'stericSl', 'shallowStericSl', 'deepStericSl', ...
         'thermostericSl', 'shallowThermostericSl', 'deepThermostericSl', ...
         'halostericSl', 'shallowHalostericSl', 'deepHalostericSl'};

    switch outputExt
        case '.mat'
            vars = ...
                [{'dates', 'lat', 'lon'}, matStericVars];

            if ~options.ForceNew && exist(outputPath, 'file') && ...
                    all(ismember(vars, who('-file', outputPath)))

                if ~options.BeQuiet
                    cprintf('[ULMO>%s] Skipped computing %s, already exist.\n', ...
                        callchaintext(callChain), filehref(outputPath, 'steric sea level data'));
                end

                return
            end

        case '.nc'

            if ~options.ForceNew && exist(outputPath, 'file')

                try
                    info = ncinfo(outputPath);
                    varNamesInFile = {info.Variables.Name};

                    if all(ismember( ...
                            {'time', 'lat', 'lon', 'steric', 'steric_shallow', 'steric_deep', ...
                             'thermosteric', 'thermosteric_shallow', 'thermosteric_deep', ...
                             'halosteric', 'halosteric_shallow', 'halosteric_deep'}, ...
                            varNamesInFile))

                        if ~options.BeQuiet
                            cprintf('[ULMO>%s] Skipped computing %s, already exist.\n', ...
                                callchaintext(callChain), filehref(outputPath, 'steric sea level data'));
                        end

                        return
                    end

                catch ME
                    warning('Error checking existing NetCDF file: %s\nProceeding to create new file. Error details:\n%s', ...
                        outputPath, getReport(ME));
                end

            end

        otherwise
            error('Unsupported output file extension: %s', outputExt);
    end

    for iFile = 1:length(inputFiles)
        inputFile = inputFiles{iFile};
        inputPath = fullfile(inputFolder, inputFile);

        if ~exist(inputPath, 'file')
            error('Input file %s does not exist.', inputFile);
        end

        if any(~ismember([{'date', 'lat', 'lon'}, matStericVars], who('-file', inputPath)))
            missingVars = setdiff([{'date', 'lat', 'lon'}, matStericVars], who('-file', inputPath));
            error('Input file %s is missing variables: %s\n', inputFile, strjoin(missingVars, ', '));
        end

        load(inputPath, 'date', matStericVars{:});

        if ~exist('lat', 'var') || ~exist('lon', 'var')
            load(inputPath, 'lat', 'lon');
        end

        if ~exist('stericSls', 'var')
            stericSls = nan([size(stericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            shallowStericSls = nan([size(shallowStericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            deepStericSls = nan([size(deepStericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            thermostericSls = nan([size(thermostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            shallowThermostericSls = nan([size(shallowThermostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            deepThermostericSls = nan([size(deepThermostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            halostericSls = nan([size(halostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            shallowHalostericSls = nan([size(shallowHalostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            deepHalostericSls = nan([size(deepHalostericSl), length(inputFiles)], 'single'); %#ok<NODEF>
            dates = NaT([1, length(inputFiles)]);
        end

        stericSls(:, :, iFile) = stericSl;
        shallowStericSls(:, :, iFile) = shallowStericSl;
        deepStericSls(:, :, iFile) = deepStericSl;
        thermostericSls(:, :, iFile) = thermostericSl;
        shallowThermostericSls(:, :, iFile) = shallowThermostericSl;
        deepThermostericSls(:, :, iFile) = deepThermostericSl;
        halostericSls(:, :, iFile) = halostericSl;
        shallowHalostericSls(:, :, iFile) = shallowHalostericSl;
        deepHalostericSls(:, :, iFile) = deepHalostericSl;
        dates(iFile) = date;
    end

    stericSl = stericSls;
    shallowStericSl = shallowStericSls;
    deepStericSl = deepStericSls;
    thermostericSl = thermostericSls;
    shallowThermostericSl = shallowThermostericSls;
    deepThermostericSl = deepThermostericSls;
    halostericSl = halostericSls;
    shallowHalostericSl = shallowHalostericSls;
    deepHalostericSl = deepHalostericSls;

    switch outputExt
        case '.mat'
            save(outputPath, vars{:}, '-v7.3');
        case '.nc'

            if exist(outputPath, 'file')
                delete(outputPath);
            end

            % Coordinates
            nccreate(outputPath, 'lat', ...
                Dimensions = {'lat', length(lat)}, ...
                Datatype = 'single', Format = 'netcdf4');
            ncwrite(outputPath, 'lat', lat);
            ncwriteatt(outputPath, 'lat', "Units", 'degrees north');
            nccreate(outputPath, 'lon', ...
                Dimensions = {'lon', length(lon)}, ...
                Datatype = 'single', Format = 'netcdf4');
            ncwrite(outputPath, 'lon', lon);
            ncwriteatt(outputPath, 'lon', "Units", 'degrees east');
            nccreate(outputPath, 'time', ...
                Dimensions = {'time', length(dates)}, ...
                Datatype = 'single', Format = 'netcdf4');
            timeDays = single(days(dates - datetime(1800, 1, 1)));
            ncwrite(outputPath, 'time', timeDays);
            ncwriteatt(outputPath, 'time', "Units", 'days since 1800-01-01');

            % Variables
            write2nc(outputPath, 'steric', single(stericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Total steric sea level'});
            write2nc(outputPath, 'steric_shallow', single(shallowStericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow total steric sea level above 2000 m'});
            write2nc(outputPath, 'steric_deep', single(deepStericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep total steric sea level below 2000 m'});
            write2nc(outputPath, 'thermosteric', single(thermostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Thermosteric sea level'});
            write2nc(outputPath, 'thermosteric_shallow', single(shallowThermostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow thermosteric sea level above 2000 m'});
            write2nc(outputPath, 'thermosteric_deep', single(deepThermostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep thermosteric sea level below 2000 m'});
            write2nc(outputPath, 'halosteric', single(halostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Halosteric sea level'});
            write2nc(outputPath, 'halosteric_shallow', single(shallowHalostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Shallow halosteric sea level above 2000 m'});
            write2nc(outputPath, 'halosteric_deep', single(deepHalostericSl), ...
                Dimensions = {'lat', 'lon', 'time'}, ...
                Attributes = {'Units', 'meters', 'LongName', 'Deep halosteric sea level below 2000 m'});

            try
                ncinfo(outputPath);
            catch ME
                error('Error creating NetCDF file: %s\nDetails:\n%s', outputPath, getReport(ME));
            end

        otherwise
            error('Unsupported output file extension: %s', outputExt);
    end

    if ~options.BeQuiet
        cprintf('[ULMO>%s] Aggregated %s.\n', ...
            callchaintext(callChain), filehref(outputPath, 'steric sea level data'));
    end

end
