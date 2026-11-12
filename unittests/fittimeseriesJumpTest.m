%% FITTIMESERIESJUMPTEST - Checks joint jump fitting and coefficient uncertainties.
% Uses analytic signals and an independent QR covariance reference, including
% per-series data selection, time conversion and reconstruction components.
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

classdef fittimeseriesJumpTest < matlab.unittest.TestCase

    properties (TestParameter)
        weighted = {false, true}
        timeType = {'numeric', 'datetime', 'duration'}
        transposeData = {false, true}
        eventTimes = struct('single', 2.5, 'multiple', [4, 2])
        format = {'cos-sin', 'sin-cos', 'amp-phase'}
        reconstruction = {'includeharmonics', 'omitharmonics', 'onlyharmonics'}
    end

    methods (TestClassSetup)

        function addUlmoPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'aux')));
        end

    end

    methods (Test)

        function recoversJointSignal(testCase, eventTimes, weighted)
            t = (0:0.05:6)';
            polynomial = [2 + 3 * t + 0.2 * t .^ 2, -1 + 0.5 * t - 0.1 * t .^ 2];
            harmonic = [cos(2 * pi * t) + 2 * sin(4 * pi * t), ...
                3 * sin(2 * pi * t) - cos(4 * pi * t)];
            amplitudes = [(1:numel(eventTimes))', -(2:1 + numel(eventTimes))'];
            x = polynomial + harmonic + double(t >= eventTimes) * amplitudes;
            sigma = [];
            if weighted, sigma = 1 + t / 10; end
            [poly, harmon, fitted, ~, ~, ~, jumps, errors] = ...
                fittimeseries(t, x, sigma, 2, [1, 0.5], JumpTimes = eventTimes);
            testCase.verifyEqual(poly, [2, -1; 3, 0.5; 0.2, -0.1], AbsTol = 1e-9);
            testCase.verifyEqual(harmon, cat(3, [1, 0; 0, 2], [0, 3; -1, 0]), AbsTol = 1e-9);
            testCase.verifyEqual(jumps, amplitudes, AbsTol = 1e-9);
            testCase.verifyEqual(fitted, x, AbsTol = 1e-9);
            testCase.verifyLessThan(max(errors, [], 'all'), 1e-9);
        end

        function agreesWithJointQR(testCase, weighted, transposeData)
            t = (0:0.05:6)';
            events = [4, 2];
            design = [ones(size(t)), t, t .^ 2, cos(2 * pi * t), ...
                sin(2 * pi * t), double(t >= events)];
            beta = [2, -1; 3, 0.5; 0.2, -0.1; 1, 2; -2, 1; 4, -3; -1, 2];
            x = design * beta + [0.02 * sin(5 * t), 0.03 * cos(3 * t)];
            x([22, 45], 1) = NaN;
            x([55, 89], 2) = Inf;
            sigma = [];
            if weighted
                sigma = [1 + 0.1 * t, 2 - 0.1 * t];
                sigma(29, 2) = NaN;
            end
            inputX = x;
            inputSigma = sigma;
            if transposeData, inputX = x.'; inputSigma = sigma.'; end
            [poly, harmon, fitted, polyErrors, harmonErrors, polyCov, jumps, errors] = ...
                fittimeseries(t, inputX, inputSigma, 2, 1, ...
                JumpTimes = events, FitRange = [0.5, 5.5]);
            [averages, amplitudePhase, ~, averageErrors, amplitudePhaseErrors, ~, ...
                transformedJumps, transformedErrors] = ...
                fittimeseries(t, inputX, inputSigma, 2, 1, ...
                JumpTimes = events, FitRange = [0.5, 5.5], ...
                PolynomialFormat = 'average-derivatives', AverageRange = [1, 3], ...
                PeriodicFormat = 'amp-phase');
            testCase.verifyEqual(transformedJumps, jumps);
            testCase.verifyEqual(transformedErrors, errors);
            testCase.verifySize(jumps, [2, 2]);
            testCase.verifySize(errors, [2, 2]);
            testCase.verifySize(fitted, size(inputX));
            if transposeData, fitted = fitted.'; end

            for iSeries = 1:2
                valid = isfinite(x(:, iSeries)) & t >= 0.5 & t <= 5.5;
                if weighted, valid = valid & isfinite(sigma(:, iSeries)); end
                selected = design(valid, :);
                data = x(valid, iSeries);
                if weighted
                    selected = selected ./ sigma(valid, iSeries);
                    data = data ./ sigma(valid, iSeries);
                end
                [Q, R] = qr(selected, 0);
                expected = R \ (Q' * data);
                inverseR = R \ eye(size(R));
                variance = sum((data - selected * expected) .^ 2) / ...
                    (nnz(valid) - size(selected, 2));
                covariance = variance * (inverseR * inverseR.');
                standardErrors = sqrt(diag(covariance));
                testCase.verifyEqual(poly(:, iSeries), expected(1:3), AbsTol = 1e-9);
                testCase.verifyEqual(harmon(:, :, iSeries), expected(4:5).', AbsTol = 1e-9);
                testCase.verifyEqual(jumps(:, iSeries), expected(6:7), AbsTol = 1e-9);
                testCase.verifyEqual(polyErrors(:, iSeries), standardErrors(1:3), AbsTol = 1e-10);
                testCase.verifyEqual(harmonErrors(:, :, iSeries), standardErrors(4:5).', AbsTol = 1e-10);
                testCase.verifyEqual(errors(:, iSeries), standardErrors(6:7), AbsTol = 1e-10);
                testCase.verifyEqual(polyCov(:, :, iSeries), covariance(1:3, 1:3), AbsTol = 1e-10);
                testCase.verifyEqual(fitted(:, iSeries), design * expected, AbsTol = 1e-9);
                transform = [1, 2, 13 / 3; 0, 1, 4; 0, 0, 2];
                testCase.verifyEqual(averages(:, iSeries), transform * expected(1:3), AbsTol = 1e-9);
                testCase.verifyEqual(averageErrors(:, iSeries), ...
                    sqrt(diag(transform * covariance(1:3, 1:3) * transform.')), AbsTol = 1e-10);
                cosine = expected(4);
                sine = expected(5);
                amplitude = hypot(cosine, sine);
                jacobian = [cosine / amplitude, sine / amplitude; ...
                    -sine / amplitude ^ 2, cosine / amplitude ^ 2];
                testCase.verifyEqual(amplitudePhase(:, :, iSeries), ...
                    [amplitude, mod(atan2(sine, cosine), 2 * pi)], AbsTol = 1e-9);
                testCase.verifyEqual(amplitudePhaseErrors(:, :, iSeries), ...
                    sqrt(diag(jacobian * covariance(4:5, 4:5) * jacobian.')).', AbsTol = 1e-10);
            end

            if weighted
                [~, ~, ~, ~, ~, ~, scaledJumps, scaledErrors] = ...
                    fittimeseries(t, inputX, 5 * inputSigma, 2, 1, ...
                    JumpTimes = events, FitRange = [0.5, 5.5]);
                testCase.verifyEqual(scaledJumps, jumps, AbsTol = 1e-9);
                testCase.verifyEqual(scaledErrors, errors, AbsTol = 1e-10);
            end
        end

        function convertsTimesAndPreservesOrder(testCase, timeType, transposeData)
            u = (12:0.125:18)';
            events = [16, 14];
            convert = @(value) value;
            if strcmp(timeType, 'duration')
                convert = @(value) days(value * 365.2425);
            elseif strcmp(timeType, 'datetime')
                convert = @(value) datetime(2005, 7, 1, 'TimeZone', 'UTC') + years(value - 12);
            end
            t = convert(u);
            x = [2 + 0.3 * u, 4 - 0.2 * u] + double(u >= events) * [3, -2; -1, 4];
            x(20, 1) = NaN;
            % The first valid datetime is midyear, and observations run backwards.
            t = flipud(t);
            x = flipud(x);
            t = [t(end); t];
            if strcmp(timeType, 'datetime'), t(1) = NaT; else, t(1) = convert(NaN); end
            x = [NaN(1, 2); x];
            if transposeData, x = x.'; t = t.'; end
            [~, ~, fitted, ~, ~, ~, jumps, errors] = fittimeseries(t, x, [], 1, [], ...
                JumpTimes = convert(events.'), FitRange = convert([12.5, 17.5]));
            testCase.verifyEqual(jumps, [3, -2; -1, 4], AbsTol = 1e-9);
            testCase.verifySize(errors, [2, 2]);
            testCase.verifySize(fitted, size(x));
            valid = isfinite(x);
            testCase.verifyEqual(fitted(valid), x(valid), AbsTol = 1e-9);
        end

        function reconstructsSelectedComponents(testCase, reconstruction, format)
            t = (0:0.125:6)';
            polynomial = 2 + 3 * t + 0.2 * t .^ 2;
            harmonic = cos(2 * pi * t) + 2 * sin(2 * pi * t);
            step = 4 * double(t >= 2);
            expected = polynomial + harmonic + step;
            if strcmp(reconstruction, 'omitharmonics'), expected = polynomial + step; end
            if strcmp(reconstruction, 'onlyharmonics'), expected = harmonic; end
            [poly, harmon, fitted, ~, ~, ~, jumps] = ...
                fittimeseries(t, polynomial + harmonic + step, [], 2, 1, ...
                JumpTimes = 2, PeriodicFormat = format, Reconstruction = reconstruction, ...
                PolynomialFormat = 'average-derivatives', AverageRange = [1, 3]);
            testCase.verifyEqual(fitted, expected, AbsTol = 1e-9);
            testCase.verifyEqual(poly, [2 + 6 + 0.2 * 13 / 3; 3.8; 0.4], AbsTol = 1e-9);
            testCase.verifyEqual(jumps, 4, AbsTol = 1e-9);
            expectedHarmon = [1, 2];
            if strcmp(format, 'sin-cos'), expectedHarmon = [2, 1]; end
            if strcmp(format, 'amp-phase'), expectedHarmon = [sqrt(5), atan2(2, 1)]; end
            testCase.verifyEqual(harmon, expectedHarmon, AbsTol = 1e-9);
        end

        function usesPostEventValue(testCase, reconstruction)
            t = (0:8)';
            x = 2 + 4 * double(t >= 4);
            [~, harmon, fitted, ~, harmonErrors, ~, jumps] = ...
                fittimeseries(t, x, [], 0, [], JumpTimes = 4, Reconstruction = reconstruction);
            testCase.verifyEmpty(harmon);
            testCase.verifyEmpty(harmonErrors);
            testCase.verifyEqual(jumps, 4, AbsTol = 1e-12);
            if strcmp(reconstruction, 'onlyharmonics')
                testCase.verifyEqual(fitted, zeros(size(t)));
            else
                testCase.verifyEqual(fitted, x, AbsTol = 1e-12);
                testCase.verifyEqual(fitted(t == 4), 6, AbsTol = 1e-12);
            end
        end

        function noEventsReturnZeroRows(testCase, timeType, transposeData)
            t = (0:0.25:4)';
            convert = @(value) value;
            if strcmp(timeType, 'duration'), convert = @(value) years(value); end
            if strcmp(timeType, 'datetime'), convert = @(value) datetime(2005, 1, 1) + years(value); end
            x = [1 + t, 2 - t, 3 + 2 * t];
            if transposeData, x = x.'; end
            implicit = cell(1, 8);
            explicit = cell(1, 8);
            [implicit{:}] = fittimeseries(convert(t), x, [], 1);
            [explicit{:}] = fittimeseries(convert(t), x, [], 1, [], JumpTimes = convert([]));
            testCase.verifyEqual(explicit, implicit);
            testCase.verifySize(implicit{7}, [0, 3]);
            testCase.verifySize(implicit{8}, [0, 3]);
        end

        function rejectsInvalidEvents(testCase, timeType)
            t = (0:8)';
            convert = @(value) value;
            if strcmp(timeType, 'duration'), convert = @(value) years(value); end
            if strcmp(timeType, 'datetime'), convert = @(value) datetime(2005, 1, 1) + years(value); end
            for events = {[2, 2], NaN, Inf, -Inf}
                testCase.verifyError(@() fittimeseries(convert(t), t, [], 1, [], ...
                    JumpTimes = convert(events{1})), 'ULMO:fittimeseries:InvalidJumpTimes');
            end
            wrongType = years(2);
            if strcmp(timeType, 'duration'), wrongType = 2; end
            testCase.verifyError(@() fittimeseries(convert(t), t, [], 1, [], ...
                JumpTimes = wrongType), 'ULMO:fittimeseries:InvalidJumpTimes');
        end

        function rejectsShapeAndComplexEvents(testCase)
            t = (0:8)';
            for events = {[1, 2; 3, 4], 2 + 1i, 'event'}
                testCase.verifyError(@() fittimeseries(t, t, [], 1, [], JumpTimes = events{1}), ...
                    'ULMO:fittimeseries:InvalidJumpTimes');
            end
        end

        function rejectsUnidentifiableDesigns(testCase)
            t = (0:8)';
            % Constant/zero steps, distinct events with identical sampled steps,
            % and an event removed by the fitting interval.
            for events = {0, -1, 9, [2.2, 2.8]}
                testCase.verifyError(@() fittimeseries(t, t, [], 1, [], JumpTimes = events{1}), ...
                    'ULMO:fittimeseries:RankDeficient');
            end
            testCase.verifyError(@() fittimeseries(t, t, [], 1, [], ...
                JumpTimes = 3, FitRange = [4, 8]), 'ULMO:fittimeseries:RankDeficient');
            % Both sides are sampled, but a harmonic duplicates a step column.
            t = repmat([0; 1], 10, 1);
            testCase.verifyError(@() fittimeseries(t, t, [], 0, 2, JumpTimes = 0.5), ...
                'ULMO:fittimeseries:RankDeficient');
        end

        function checksEachSeriesAfterSelection(testCase)
            t = (0:10)';
            x = [t, t];
            x(t < 4, 2) = NaN;
            testCase.verifyError(@() fittimeseries(t, x, [], 1, [], JumpTimes = 4), ...
                'ULMO:fittimeseries:RankDeficient');
            x = [t, t];
            sigma = ones(size(x));
            sigma(t >= 4, 2) = NaN;
            testCase.verifyError(@() fittimeseries(t, x, sigma, 0, [], JumpTimes = 4), ...
                'ULMO:fittimeseries:RankDeficient');
        end

        function countsJumpParameters(testCase)
            t = (0:3)';
            testCase.verifyError(@() fittimeseries(t, t, [], 1, [], JumpTimes = [1, 2]), ...
                'ULMO:fittimeseries:InsufficientData');
        end

        function omitsInvalidTimesFromHarmonicReconstruction(testCase)
            t = (0:8)';
            x = 2 + 4 * double(t >= 4);
            t([2, 3, 8]) = [NaN, Inf, -Inf];
            [~, ~, fitted, ~, ~, ~, jumps] = fittimeseries(t, x, [], 0, [], ...
                JumpTimes = 4, Reconstruction = 'onlyharmonics');
            testCase.verifyEqual(jumps, 4, AbsTol = 1e-12);
            testCase.verifyTrue(all(isnan(fitted(~isfinite(t)))));
            testCase.verifyEqual(fitted(isfinite(t)), zeros(nnz(isfinite(t)), 1));
        end

    end

end
