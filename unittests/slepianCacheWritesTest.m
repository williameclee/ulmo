%% SLEPIANCACHEWRITESTEST - Verifies geographic basis and kernel cache controls.
% Isolated data roots cover fresh calculations, cached reads and forced
% calculations without creating or replacing basis, kernel or quadrature files.
%
% Created by
%   2026/11/05, En-Chi Lee (williameclee@gmail.com)

classdef slepianCacheWritesTest < matlab.unittest.TestCase

    properties
        Folder
        Polygon = [0 -60; 120 -60; 120 60; 0 60; 0 -60]
    end

    properties (TestParameter)
        target = {'basis', 'kernel'}
    end

    methods (TestMethodSetup)

        function prepare(t)
            tmp = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Folder = tmp.Folder;
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', t.Folder));

            for folder = ["GLMALPHA", "KERNELC", "KERNELCP", "LEGENDRE", "EARTHMODELS/CONSTANTS"]
                mkdir(fullfile(t.Folder, folder));
            end

            Earth.Radius = 6371000;
            save(fullfile(t.Folder, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat'), 'Earth');
        end

    end

    methods (Test)

        function freshCalculationDoesNotWrite(t, target)
            before = snapshot(t.Folder);
            unsaved = t.calculate(target, SaveData = false);
            t.verifyEqual(snapshot(t.Folder), before);
            saved = t.calculate(target);
            t.verifyEqual(unsaved, saved, AbsTol = 1e-12);
            t.verifyGreaterThan(size(snapshot(t.Folder), 1), size(before, 1));
            t.verifyTrue(isfile(fullfile(t.Folder, 'gausslegendrecof.mat')));
        end

        function cachedAndForcedRequestsPreserveFiles(t, target)
            saved = t.calculate(target);
            before = snapshot(t.Folder);
            cached = t.calculate(target, SaveData = false);
            t.verifyEqual(cached, saved);
            t.verifyEqual(snapshot(t.Folder), before);
            rebuilt = t.calculate(target, ForceNew = true, SaveData = false);
            t.verifyEqual(rebuilt, saved, AbsTol = 1e-12);
            t.verifyEqual(snapshot(t.Folder), before);
        end

    end

    methods (Access = private)

        function result = calculate(t, target, varargin)

            if strcmp(target, 'basis')
                [g, v] = glmalpha_new(t.Polygon, 1, 'BeQuiet', 2, varargin{:});
                result = [g; v(:).'];
            else
                result = kernelcp_new(1, t.Polygon, 10, "BeQuiet", 2, varargin{:});
            end

        end

    end

end

function files = snapshot(folder)
    entries = dir(fullfile(folder, '**', '*'));
    entries = entries(~[entries.isdir]);
    files = cell(numel(entries), 2);

    for iFile = 1:numel(entries)
        path = fullfile(entries(iFile).folder, entries(iFile).name);
        files{iFile, 1} = erase(path, folder);
        files{iFile, 2} = readbinary(path);
    end

end

function bytes = readbinary(path)
    fid = fopen(path, 'rb');
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
end
