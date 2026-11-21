%% INTERPTEMPORALTEST - Tests temporal interpolation of sparse grids.
%
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

classdef interptemporalTest < matlab.unittest.TestCase

    properties (TestParameter)
        Method = {'pchip', 'spline', 'makima'};
    end

    methods (Test)

        function nonlinearInterpolationPreservesSparseSeries(testCase, Method)
            % A fully missing series and a single observation share the grid
            % with a series that can be interpolated.
            dates = datetime(2020, 1, [1; 3; 5]);
            mesh = reshape([0, NaN, NaN; 2, NaN, 8; 4, NaN, NaN]', 1, 3, []);
            [values, labels] = interptemporal(dates, mesh, days(1), Method, BeQuiet = true);

            testCase.verifyEqual(labels(:), dates(1) + days((0:4)'));
            testCase.verifyEqual(reshape(values, 3, []), ...
                [0, 1, 2, 3, 4; NaN, NaN, NaN, NaN, NaN; NaN, NaN, 8, NaN, NaN], AbsTol = 1e-12);
        end

        function nonlinearInterpolationHandlesOnlySparseSeries(testCase, Method)
            % The grid can contain no series with enough samples to interpolate.
            dates = datetime(2020, 1, [1; 3; 5]);
            mesh = reshape([NaN; 8; NaN], 1, 1, []);
            values = interptemporal(dates, mesh, days(1), Method, BeQuiet = true);
            testCase.verifyEqual(reshape(values, 1, []), [NaN, NaN, 8, NaN, NaN]);
        end

    end

end
