function [meshIntp, datesIntp] = ...
        interptemporal(dates, mesh, timeStep, intpMthd, options)

    arguments (Input)
        dates {mustBeA(dates, {'datetime', 'numeric'}), mustBeVector}
        mesh {mustBeNumeric}
        timeStep {mustBeTimeStep}
        intpMthd char = "mean"
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

    if ischar(timeStep) && strcmpi(timeStep, 'midmonth')
        startDate = datetime(year(dates(1)), month(dates(1)), 1, 0, 0, 0);
        endDate = datetime(year(dates(end)), month(dates(end)), 1, 0, 0, 0) + calmonths(1);
        dmonths = ceil((endDate - startDate) / days(28));
        startDates = startDate + calmonths(0:dmonths - 1);
        endDates = startDate + calmonths(1:dmonths);
        datesIntp = startDates + (endDates - startDates) / 2;
        datesIntp = datesIntp(datesIntp >= dates(1) & datesIntp <= dates(end));
    else

        if mean(diff(dates)) > timeStep
            warning(sprintf('%s:InterpolationStepTooSmall', upper(mfilename)), ...
                'The interpolation time step (%s) is smaller than the mean data resolution (%s)', timeStep, mean(diff(dates)));
        end

        datesIntp = dates(1):timeStep:dates(end);
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
        meshIntp = interp1(dates, meshFlat, datesIntp, intpMthd)';
    end

    meshIntp = reshape(meshIntp, [size(mesh, 1:2), length(datesIntp)]);

    if ~options.BeQuiet
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end
