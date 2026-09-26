%% oceanTopologyTest - Integration regression; requires the configured GSHHS/IHO coastline data
%
% Created by
%   2026/09/25, williameclee@arizona.edu (@williameclee)
% Last modified by
%   2026/09/26, williameclee@arizona.edu (@williameclee)

classdef oceanTopologyTest < matlab.unittest.TestCase

    properties
        OriginalCoasts
        OriginalGshhs
        ScratchDir
    end

    properties (TestParameter)
        buf = {0.5, 2}
        limit = {60, 66}
    end

    methods (TestClassSetup)

        function setupCoastsEnvironment(testCase)
            % Isolate ocean and coastline caches, including cold-cache writes.
            testCase.OriginalCoasts = getenv('COASTS');
            testCase.OriginalGshhs = getenv('GSHHS');
            testCase.ScratchDir = tempname;
            mkdir(testCase.ScratchDir);
            % Register cleanup before copying inputs so setup failures also clean up.
            testCase.addTeardown(@() testCase.restoreCoastsEnvironment());
            scratchGshhs = fullfile(testCase.ScratchDir, 'gshhs');
            mkdir(scratchGshhs);
            copyfile(fullfile(testCase.OriginalCoasts, 'Limits_of_oceans_and_seas.mat'), ...
                testCase.ScratchDir);
            % Copy only raw input; all derived coastline caches belong to this run.
            copyfile(fullfile(testCase.OriginalGshhs, 'gshhs_l.b'), scratchGshhs);
            setenv('COASTS', testCase.ScratchDir);
            setenv('GSHHS', scratchGshhs);
        end

    end

    methods (Access = private)

        function restoreCoastsEnvironment(testCase)
            % Restore both environment variables before removing scratch caches.
            setenv('COASTS', testCase.OriginalCoasts);
            setenv('GSHHS', testCase.OriginalGshhs);

            if exist(testCase.ScratchDir, 'dir')
                rmdir(testCase.ScratchDir, 's');
            end

        end

    end

    methods (Test)

        function testOceanTopology(testCase, buf, limit)
            % Verify ocean polygon topology, land non-overlap, and boundary consistency
            parts = {'Atlantic Ocean', 'Indian Ocean', ...
                         'Pacific Ocean, eastern part', 'Pacific Ocean, western part'};

            args = {'Buffer', buf, 'Latlim', limit, ...
                        'MoreBuffers', {'earthquakes', 10}, 'BeQuiet', true};
            [xy, p] = oceans(args{:}, 'SaveData', false);
            [op, ~, lo] = oceanpoly(parts, [-limit limit], 200, 'BeQuiet', true);
            [~, land] = gshhscoastline('l', 'Buffer', buf, ...
                'LatLim', [-90 90], 'LonLim', lo, 'LonOrigin', 200, ...
                'SaveData', false, 'BeQuiet', true);

            % 1. Verify ocean polygon and land polygon do not overlap
            testCase.verifyLessThan(area(intersect(p, land)), 1e-6);

            % 2. Verify ocean polygon is fully contained within reference domain
            testCase.verifyLessThan(area(subtract(p, op)), 1e-6);

            % 3. Verify consistency between output boundary coordinates and polyshape
            testCase.verifyLessThan(area(xor(p, polyshape(xy))), 1e-6);

            % 4. Verify every southern cutoff edge crosses water, not an Antarctic
            % land hole separated from it by a roundoff-sized sliver
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
