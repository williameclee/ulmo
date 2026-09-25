%% oceanTopologyTest - Integration regression; requires the configured GSHHS/IHO coastline data
classdef oceanTopologyTest < matlab.unittest.TestCase

    properties
        OriginalCoasts
        ScratchDir
    end

    properties (TestParameter)
        buf = {0.5, 2}
        limit = {60, 66}
    end

    methods (TestClassSetup)

        function setupCoastsEnvironment(testCase)
            testCase.OriginalCoasts = getenv('COASTS');
            testCase.ScratchDir = tempname;
            mkdir(testCase.ScratchDir);
            copyfile(fullfile(testCase.OriginalCoasts, 'Limits_of_oceans_and_seas.mat'), ...
                testCase.ScratchDir);
            setenv('COASTS', testCase.ScratchDir);
        end

    end

    methods (TestClassTeardown)

        function restoreCoastsEnvironment(testCase)
            setenv('COASTS', testCase.OriginalCoasts);

            if exist(testCase.ScratchDir, 'dir')
                rmdir(testCase.ScratchDir, 's');
            end

        end

    end

    methods (Test)

        function testOceanTopology(testCase, buf, limit)
            parts = {'Atlantic Ocean', 'Indian Ocean', ...
                         'Pacific Ocean, eastern part', 'Pacific Ocean, western part'};

            args = {'Buffer', buf, 'Latlim', limit, ...
                        'MoreBuffers', {'earthquakes', 10}, 'BeQuiet', true};
            [xy, p] = oceans(args{:}, 'SaveData', false);
            [op, ~, lo] = oceanpoly(parts, [-limit limit], 200, 'BeQuiet', true);
            [~, land] = gshhscoastline('l', 'Buffer', buf, ...
                'LatLim', [-90 90], 'LonLim', lo, 'LonOrigin', 200, ...
                'SaveData', false, 'BeQuiet', true);

            testCase.verifyLessThan(area(intersect(p, land)), 1e-6);
            testCase.verifyLessThan(area(subtract(p, op)), 1e-6);
            testCase.verifyLessThan(area(xor(p, polyshape(xy))), 1e-6);

            % Every southern cutoff edge must cross water, not an Antarctic
            % land hole separated from it by a roundoff-sized sliver.
            edges = find(abs(xy(1:end - 1, 2) + limit) < 1e-8 & ...
                abs(xy(2:end, 2) + limit) < 1e-8);

            for edge = edges'
                x = linspace(xy(edge, 1), xy(edge + 1, 1), 101);
                testCase.verifyFalse(any(isinterior(land, x(2:end - 1), ...
                    repmat(-limit +1e-7, 1, 99))));
            end

        end

    end

end
