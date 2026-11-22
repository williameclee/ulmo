%% SLEPT2SEISMIC - Estimates static earthquake fields in a supplied Slepian basis.
% Fits permanent steps jointly with polynomial and periodic backgrounds.
% Each output row is an incremental event field, in input event order.
%
% Syntax
%   [seis, seisSigma] = slept2seismic(slept, dates, eventDates)
%   [seis, seisSigma] = slept2seismic(slept, dates, eventDates, sigma)
%   [seis, seisSigma] = slept2seismic(__, "Name", Value)
%
% Input arguments
%   slept - Slepian coefficients.
%       Size: (nDates, J).
%   dates - Observation dates.
%       Type: DOUBLE (datenum) | DATETIME
%   eventDates - One or more event dates.
%       Multiple events are fitted jointly.
%       Type: DOUBLE (datenum) | DATETIME
%   sigma (optional) - Observation standard deviations for relative weighting.
%       A column is shared by all modes; a matrix supplies per-mode values.
%       Size: [], (nDates, 1) or (nDates, J).
%       Unit: same as slept.
%       Default value: [], for ordinary least squares.
%   polynomialDegree (name-value) - Degree for fitting the background signal.
%       Default degree: 2.
%   periods (name-value) - Harmonic periods for fitting the background signal.
%       Each period adds cosine and sine terms.
%       Type: DOUBLE (in years) | DURATION
%       Default periods: years([1, 0.5]), annual and semiannual.
%   FitRange (name-value) - Inclusive interval selecting observations to fit.
%       Unit and type: Same as eventDates.
%       Size: (1, 2).
%       Default interval: earliest to latest nonmissing observation date.
%
% Output arguments
%   seis - Signed incremental jump coefficients for each event.
%       Size: (nEvents, J).
%       Unit: same as slept.
%   seisSigma - Coefficient standard errors from the joint fit covariance.
%       Size: (nEvents, J).
%       Unit: same as slept.
%
% See also
%   fittimeseries
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@gmail.com)

function [seis, seisSigma] = slept2seismic(slept, dates, eventDates, sigma, options)

    arguments (Input)
        slept {mustBeNumeric, mustBeReal, mustBeMatrix, mustBeNonempty}
        dates {mustBeA(dates, {'numeric', 'datetime'}), mustBeVector, mustBeNonempty}
        eventDates {mustBeA(eventDates, {'numeric', 'datetime'}), mustBeVector, mustBeNonempty}
        sigma {mustBeNumeric, mustBeReal, mustBeMatrix} = []
        options.polynomialDegree (1, 1) double {mustBeInteger, mustBeNonnegative, mustBeFinite} = 2
        options.periods {mustBeA(options.periods, {'numeric', 'duration'})} = years([1, 0.5])
        options.FitRange (1, 2) {mustBeA(options.FitRange, {'numeric', 'datetime'})} = ...
            [min(dates(:), [], 'omitmissing'), max(dates(:), [], 'omitmissing')]
    end

    if size(slept, 1) ~= numel(dates)
        error('ULMO:slept2seismic:SizeMismatch', ...
        'slept must have one row per observation date.');
    end

    if ~isempty(sigma) && (size(sigma, 1) ~= numel(dates) ...
            || (size(sigma, 2) ~= 1 && size(sigma, 2) ~= size(slept, 2)))
        error('ULMO:slept2seismic:SizeMismatch', ...
        'sigma must be empty, a column per date, or the same size as slept.');
    end

    dates = convertdates(dates);
    eventDates = convertdates(eventDates);
    fitRange = convertdates(options.FitRange);

    [~, ~, ~, ~, ~, ~, seis, seisSigma] = fittimeseries( ...
        dates, slept, sigma, options.polynomialDegree, options.periods, ...
        JumpTimes = eventDates, FitRange = fitRange);
end

function dates = convertdates(dates)

    if isnumeric(dates)
        mustBeReal(dates);
        dates = datetime(double(dates), 'ConvertFrom', 'datenum');
    end

end
