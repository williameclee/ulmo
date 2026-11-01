%% INTERPTEMPORAL - Interpolate a gridded time series at a fixed time step.
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
%   timeStep - Positive duration or number of days between output samples.
%       The first output time is anchored to the first native sample.
%   intpMthd (optional) - Interpolation method to use (passed to INTERP1)
%       Default method: "linear".
%   TimeRange (name-value) - Two-element vector of DATETIME or DATENUM
%       values specifying the time range to interpolate to
%       If left empty, will interpolate to the full range of input dates.
%       Default is empty.
%   BeQuiet (name-value) - If true, suppresses progress messages
%
% Last modified
%   2026/11/01, En-Chi Lee (williameclee@arizona.edu)

function [meshIntp, datesIntp] = interptemporal(dates, mesh, timeStep, intpMthd, options)

    arguments (Input)
        dates {mustBeA(dates, {'datetime', 'numeric'}), mustBeVector}
        mesh {mustBeNumeric}
        timeStep {mustBeA(timeStep, {'numeric', 'duration'}), mustBeNonempty, mustBeTimeStep}
        intpMthd {mustBeTextScalar} = "linear"
        options.TimeRange = []
        options.BeQuiet (1, 1) logical = false
        options.CallChain (1, :) cell = {}
    end

    intpMthd = validatestring(intpMthd, ...
        {'linear', 'nearest', 'previous', 'next', 'pchip', 'cubic', 'spline', 'makima', 'v5cubic'});
    if size(mesh, 3) ~= numel(dates)
        error('ULMO:interpTemporal:DimensionMismatch', ...
            'The mesh time dimension must match the number of dates.');
    end

    callChain = [options.CallChain, {mfilename}];

    if ~options.BeQuiet
        t = tic;
        templine = 'this may take a while...';
        fprintf('[ULMO>%s] Interpolating temporally, %s\n', ...
            callchaintext(callChain), templine);
    end

    limits = options.TimeRange;
    if isdatetime(dates)
        if isnumeric(timeStep)
            timeStep = days(timeStep);
        end
        if isnumeric(limits)
            limits = datetime(limits, 'ConvertFrom', 'datenum');
        end
    else
        if isduration(timeStep)
            timeStep = days(timeStep);
        end
        if isdatetime(limits)
            limits = datenum(limits); %#ok<DATNM>
        end
    end

    if mean(diff(dates)) > timeStep
        warning('ULMO:interpTemporal:InterpolationStepTooSmall', ...
            'The interpolation time step (%s) is smaller than the mean data resolution (%s)', ...
            string(timeStep), string(mean(diff(dates))));
    end

    datesIntp = dates(1):timeStep:dates(end);
    if ~isempty(limits)
        datesIntp = datesIntp(datesIntp >= limits(1) & datesIntp <= limits(2));
    end

    if isempty(datesIntp)
        meshIntp = mesh(:, :, []);
    elseif isequal(dates(:), datesIntp(:))
        meshIntp = mesh;
    else
        meshFlat = reshape(mesh, [], numel(dates)).';
        values = interp1(dates, meshFlat, datesIntp, intpMthd);
        meshIntp = reshape(values.', [size(mesh, 1:2), numel(datesIntp)]);
    end

    if ~options.BeQuiet
        fprintf(repmat('\b', 1, length(templine) + 1));
        fprintf('took %.1f seconds.\n', toc(t));
    end

end
