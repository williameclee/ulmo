%% DOMAINMASKLONGITUDETEST - Verify masks for equivalent geographic coordinates.
% Tests use temporary mask caches and do not require steric observations.
%
% Created by
%   2026/11/05, OpenAI Codex

classdef domainmaskLongitudeTest < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function temporaryRoot(t)
            folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', folder.Folder));
        end
    end

    methods (Test)
        function shiftedSouthAtlanticGrid(t)
            % Issue #101: the same native grid must select the same ocean cells.
            domain = GeoDomain('satlantic', Buffer=0.5, ...
                DefaultParams=true, Latlim=[-60, 60]);
            lon = 20.5:379.5;
            lat = (-64.5:79.5).';
            shifted = domainmask(domain, lon, lat, ...
                Forcenew=true, SaveData=false, BeQuiet=true);
            [wrappedLon, order] = sort(mod(lon + 180, 360) - 180);
            wrapped = domainmask(domain, wrappedLon, lat, ...
                Forcenew=true, SaveData=false, BeQuiet=true);

            t.verifyEqual(nnz(xor(shifted(:, order), wrapped)), 0, ...
                'Equivalent longitude grids must select identical physical cells.');
        end
    end
end
