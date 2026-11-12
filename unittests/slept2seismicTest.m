%% SLEPT2SEISMICTEST - Checks static event fields and the fitter interface.
% Uses synthetic polynomial, seasonal and step signals to check amplitudes,
% units, ordering, date conversion, relative weights and inherited errors.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

classdef slept2seismicTest < matlab.unittest.TestCase

    properties (TestParameter)
        weightMode = {'none', 'shared', 'perMode'}
        background = struct('default', 2, 'linear', 1, 'constant', 0)
        numericDates = {false, true}
        numericEvents = {false, true}
        numericRange = {false, true}
        eventCount = {1, 2}
        modeCount = {1, 3}
        scale = {1, -1000}
    end

    methods (TestClassSetup)

        function addUlmoPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'aux')));
        end

    end

    methods (Test)

        function matchesDirectFitter(testCase, weightMode, background)
            [dates, events, slept, sigma] = syntheticseries();
            slept(17, 1) = NaN;
            slept(23, 2) = Inf;
            dates(29) = NaT;
            if strcmp(weightMode, 'none')
                sigma = [];
            elseif strcmp(weightMode, 'shared')
                sigma = sigma(:, 1);
            else
                sigma(35, 3) = NaN;
            end
            periods = years([1, 0.5]);
            if background == 1, periods = 1.7; end
            if background == 0, periods = []; end
            range = dates([8, end - 8]).';
            [~, ~, ~, ~, ~, ~, expected, expectedSigma] = fittimeseries( ...
                dates, slept, sigma, background, periods, ...
                JumpTimes = events, FitRange = range);
            [actual, actualSigma] = slept2seismic(slept, dates, events, sigma, ...
                polynomialDegree = background, periods = periods, FitRange = range);
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(actualSigma, expectedSigma);
            testCase.verifyGreaterThan(actualSigma, zeros(size(actualSigma)));
        end

        function convertsDateInputs(testCase, numericDates, numericEvents, numericRange)
            [dates, events, slept, sigma] = syntheticseries();
            range = dates([9, end - 7]).';
            [expected, expectedSigma] = slept2seismic(slept, dates, events, sigma, FitRange = range);
            if numericDates, dates = datenum(dates); end %#ok<DATNM>
            if numericEvents, events = datenum(events); end %#ok<DATNM>
            if numericRange, range = datenum(range); end %#ok<DATNM>
            % Date orientation must not change event or mode ordering.
            [actual, actualSigma] = slept2seismic(slept, dates.', events.', sigma, FitRange = range);
            testCase.verifyEqual(actual, expected, AbsTol = 1e-10);
            testCase.verifyEqual(actualSigma, expectedSigma, AbsTol = 1e-10);
        end

        function recoversSignedIncrementalFields(testCase, eventCount, modeCount, scale)
            [dates, events, ~, ~, signal, amplitudes] = syntheticseries();
            events = events(1:eventCount);
            amplitudes = amplitudes(1:eventCount, 1:modeCount);
            slept = scale * (signal(:, 1:modeCount) + double(dates >= events) * amplitudes);
            [actual, actualSigma] = slept2seismic(slept, dates, events);
            testCase.verifySize(actual, [eventCount, modeCount]);
            testCase.verifySize(actualSigma, [eventCount, modeCount]);
            testCase.verifyEqual(actual, scale * amplitudes, AbsTol = 1e-8);
            testCase.verifyLessThan(max(actualSigma, [], 'all'), 1e-8);
        end

        function defaultsMatchExplicitFit(testCase)
            [dates, events, slept] = syntheticseries();
            [actual, actualSigma] = slept2seismic(slept, dates, events);
            [~, ~, ~, ~, ~, ~, expected, expectedSigma] = fittimeseries( ...
                dates, slept, [], 2, years([1, 0.5]), JumpTimes = events);
            testCase.verifyEqual(actual, expected);
            testCase.verifyEqual(actualSigma, expectedSigma);
            [serial, serialSigma] = slept2seismic(slept, datenum(dates), datenum(events)); %#ok<DATNM>
            testCase.verifyEqual(serial, expected, AbsTol = 1e-10);
            testCase.verifyEqual(serialSigma, expectedSigma, AbsTol = 1e-10);
            [emptySigma, emptySigmaErrors] = slept2seismic(slept, dates, events, []);
            testCase.verifyEqual(emptySigma, actual);
            testCase.verifyEqual(emptySigmaErrors, actualSigma);
        end

        function preservesModeOrderAndUnits(testCase)
            [dates, events, slept, sigma] = syntheticseries();
            [expected, expectedSigma] = slept2seismic(slept, dates, events, sigma);
            order = [3, 1, 2];
            scaleFactors = [-1000, 0.01, 2];
            [actual, actualSigma] = slept2seismic(slept(:, order) .* scaleFactors, ...
                dates, events, sigma(:, order) .* abs(scaleFactors));
            testCase.verifyEqual(actual, expected(:, order) .* scaleFactors, AbsTol = 1e-8);
            testCase.verifyEqual(actualSigma, expectedSigma(:, order) .* abs(scaleFactors), AbsTol = 1e-8);
        end

        function fitsEventsJointly(testCase)
            [dates, events, slept] = syntheticseries();
            [joint, jointSigma] = slept2seismic(slept, dates, events);
            [reversed, reversedSigma] = slept2seismic(slept, dates, fliplr(events));
            testCase.verifyEqual(reversed, flipud(joint), AbsTol = 1e-10);
            testCase.verifyEqual(reversedSigma, flipud(jointSigma), AbsTol = 1e-10);
            independent = slept2seismic(slept, dates, events(1));
            testCase.verifyGreaterThan(norm(independent - joint(1, :)), 0.1);
        end

        function rejectsNonTimeFirstInputs(testCase)
            [dates, events, slept, sigma] = syntheticseries();
            identifier = 'ULMO:slept2seismic:SizeMismatch';
            testCase.verifyError(@() slept2seismic(slept.', dates, events), identifier);
            testCase.verifyError(@() slept2seismic(slept, dates, events, sigma.'), identifier);
            testCase.verifyError(@() slept2seismic(slept, dates, events, sigma(:, 1).'), identifier);
            testCase.verifyError(@() slept2seismic(slept, dates, events, sigma(:, 1:2)), identifier);
            testCase.verifyError(@() slept2seismic(slept, dates(1:end - 1), events), identifier);
        end

        function inheritsEventValidation(testCase)
            [dates, events, slept] = syntheticseries();
            testCase.verifyError(@() slept2seismic(slept, dates, [events(1), events(1)]), ...
                'ULMO:fittimeseries:InvalidJumpTimes');
            testCase.verifyError(@() slept2seismic(slept, dates, NaT), ...
                'ULMO:fittimeseries:InvalidJumpTimes');
            testCase.verifyError(@() slept2seismic(slept, dates, dates(1)), ...
                'ULMO:fittimeseries:RankDeficient');
            testCase.verifyError(@() slept2seismic(slept, dates, dates(end) + days(1)), ...
                'ULMO:fittimeseries:RankDeficient');
        end

        function inheritsSelectionAndWeightValidation(testCase)
            [dates, events, slept, sigma] = syntheticseries();
            testCase.verifyError(@() slept2seismic(slept, dates, events, ...
                FitRange = [events(1), dates(end)]), 'ULMO:fittimeseries:RankDeficient');
            slept(dates < events(1), 2) = NaN;
            testCase.verifyError(@() slept2seismic(slept, dates, events), ...
                'ULMO:fittimeseries:RankDeficient');
            sigma(1) = 0;
            testCase.verifyError(@() slept2seismic(slept, dates, events, sigma), ...
                'ULMO:fittimeseries:InvalidSigma');
        end

    end

end

function [dates, events, slept, sigma, background, amplitudes] = syntheticseries()
    % Whole-day dates round-trip exactly through MATLAB serial date numbers.
    dates = datetime(2003, 1, 1) + days((0:15:2190)');
    t = years(dates - dates(1));
    events = dates([99, 50]).';
    amplitudes = [4, -3, 0.5; -2, 6, -1];
    background = [2 + 3 * t + 0.2 * t .^ 2, ...
        -1 + 0.5 * t - 0.1 * t .^ 2, 4 - 0.7 * t + 0.05 * t .^ 2] ...
        + cos(2 * pi * t) * [1, 2, -1] + sin(2 * pi * t) * [2, -1, 3] ...
        + cos(4 * pi * t) * [-1, 0.2, 0.7] + sin(4 * pi * t) * [0.5, 1, -2];
    slept = background + double(dates >= events) * amplitudes ...
        + [0.02 * sin(5 * t), 0.03 * cos(3 * t), 0.01 * sin(7 * t)];
    sigma = [1 + 0.1 * t, 2 - 0.1 * t, 0.5 + 0.05 * t];
end
