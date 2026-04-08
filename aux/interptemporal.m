%% INTERPTEMPORAL - Interpolates data temporally to a specified time step
%
% Syntax
%   [meshIntp, datesIntp] = INTERPTEMPORAL(dates, mesh, timeStep, intpMthd)
%   [meshIntp, datesIntp] = INTERPTEMPORAL(__, "Name", value)
%
% Input arguments
%   dates - Vector of DATETIME or DATENUM values corresponding to the time
%       dimension of the input mesh
%   mesh - 3D numeric array with the first two dimensions corresponding to
%       spatial dimensions and the third dimension corresponding to time
%       Size: [nX, nY, nDates]
%   timeStep - Desired time step for interpolation
%       Can be a duration (e.g. hours(6), days(1)) or the string "midmonth"
%       , which will interpolate to the middle of each month.
%   intpMthd (optional) - Interpolation method to use (passed to INTERP1)
%       Default method is "mean".
%   TimeRange (name-value) - Two-element vector of DATETIME or DATENUM
%       values specifying the time range to interpolate to
%       If left empty, will interpolate to the full range of input dates.
%       Default is empty.
%   BeQuiet (name-value) - If true, suppresses progress messages
%
% Last modified by
%   2026/04/08, En-Chi Lee (williameclee@arizona.edu)
%     - Added fallback to "nearest" interpolation when "mean" method is not
%       compatible with INTERP1
%     - Added option to specify a time range for interpolation

function [meshIntp, datesIntp] = ...
        interptemporal(dates, mesh, timeStep, intpMthd, options)

    arguments (Input)
        dates {mustBeA(dates, {'datetime', 'numeric'}), mustBeVector}
        mesh {mustBeNumeric}
        timeStep {mustBeTimeStep}
        intpMthd {mustBeTextScalar} = "mean"
        options.TimeRange = []
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    callChain = [options.CallChain, {mfilename}];

    if ~options.BeQuiet
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO>%s] Interpolating temporally, %s\n', ...
            callchaintext(callChain), templine);
    end

    if (ischar(timeStep) || isstring(timeStep)) && strcmpi(timeStep, 'midmonth')
        startDate = datetime(year(dates(1)), month(dates(1)), 1, 0, 0, 0);
        endDate = datetime(year(dates(end)), month(dates(end)), 1, 0, 0, 0) + calmonths(1);
        dmonths = ceil((endDate - startDate) / days(28));
        startDates = startDate + calmonths(0:dmonths - 1);
        endDates = startDate + calmonths(1:dmonths);
        datesIntp = startDates + (endDates - startDates) / 2;
        datesIntp = datesIntp(datesIntp >= dates(1) & datesIntp <= dates(end));
    else

        if mean(diff(dates)) > timeStep
            warning('Slepian:interpTemporal:InterpolationStepTooSmall', ...
                'The interpolation time step (%s) is smaller than the mean data resolution (%s)', ...
                timeStep, mean(diff(dates)));
        end

        datesIntp = dates(1):timeStep:dates(end);
    end

    if ~isempty(options.TimeRange)
        datesIntp = datesIntp(datesIntp >= options.TimeRange(1) & datesIntp <= options.TimeRange(2));
    end

    meshFlat = reshape(mesh, [prod(size(mesh, 1:2)), size(mesh, 3)])';

    if intpMthd == "linear" && mean(diff(dates)) * 5 <= mean(diff(datesIntp))
        % Take the average of all points within each interpolation bin instead
        datesDiff = abs(dates(:)' - datesIntp(:));
        [~, closestIndices] = min(datesDiff, [], 2);

        meshIntp = nan(prod(size(mesh, 1:2)), length(datesIntp));

        for i = 1:length(datesIntp)
            meshIntp(:, i) = mean(meshFlat(:, closestIndices == i), 2, "omitmissing");
        end

    else

        if strcmpi(intpMthd, "mean")
            warning('Slepian:interpTemporal:InterpolationMethodMean', ...
            'Interpolation method "mean" is not compatible with INTERP1. Using "nearest" instead.');
            intpMthd = "nearest";
        end

        meshIntp = interp1(dates, meshFlat, datesIntp, intpMthd)';
    end

    meshIntp = reshape(meshIntp, [size(mesh, 1:2), length(datesIntp)]);

    if ~options.BeQuiet
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end
