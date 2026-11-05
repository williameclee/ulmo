%% GEODOMAINIDENTITYTEST - Check NearBy identifiers and derived cache selection.
%
% Created by
%   2026/11/04, En-Chi Lee

classdef GeoDomainIdentityTest < matlab.unittest.TestCase
    properties (TestParameter)
        nearby = {false, true}
    end
    methods (Test)
        function greenlandIdentifiers(t)
            % The default keeps its existing ID, while NearBy adds a suffix.
            domain = GeoDomain('greenland', 'Buffer', 2);
            t.verifyEqual(domain.Id, 'Greenland0-2');
            domain.NearBy = true;
            t.verifyEqual(domain.Id, 'Greenland0-2-nearby');
            domain = GeoDomain('greenland', 'NearBy', true);
            t.verifyEqual(domain.Id, 'Greenland-nearby');
        end

        function otherDomainIdentifiers(t)
            % Preserve the identifiers used by existing non-Greenland caches.
            domain = GeoDomain('arctic');
            t.verifyEqual(domain.Id, 'Arctic');
            domain = GeoDomain('oceans', 'Buffer', 1);
            t.verifyEqual(domain.Id, 'Oceans0-1');
        end

        function readsMatchingCaches(t, nearby)
            % Distinct cached values expose collisions through each public reader.
            folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = folder.Folder;
            t.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('IFILES', root));
            for name = ["DOMAINMASK", "KERNELC", "GLMALPHA"]
                mkdir(fullfile(root, name));
            end
            identifiers = {'Greenland0-2', 'Greenland0-2-nearby'};
            lon = [0 1];
            lat = [0; 1];
            for iDomain = 1:2
                id = identifiers{iDomain};
                mask = false(2);
                mask(:,iDomain) = true;
                save(fullfile(root,'DOMAINMASK', ...
                    sprintf('DOMAINMASK-%s-H1-Box0_1_0_1.mat',id)), 'mask','lon','lat');

                Klmlmp = iDomain/10;
                XY = []; K1 = []; K = [];
                save(fullfile(root,'KERNELC',sprintf('WREG-%s-0.mat',id)), ...
                    'Klmlmp','XY','K1','K');

                G = 1; V = Klmlmp; N = V; EL = 0; EM = 0;
                save(fullfile(root,'GLMALPHA',sprintf('glmalpha-%s-0-1.mat',id)), ...
                    'G','V','N','EL','EM');
            end

            domain = GeoDomain('greenland', 'Buffer', 2, 'NearBy', nearby);
            expectedIndex = 1+double(nearby);
            expectedMask = false(2);
            expectedMask(:,expectedIndex) = true;
            t.verifyEqual(domainmask(domain,lon,lat,BeQuiet=true), expectedMask);
            t.verifyEqual(kernelcp_new(0,domain,BeQuiet=true), expectedIndex/10);
            [~,eigenvalue] = glmalpha_new(domain,0,BeQuiet=true);
            t.verifyEqual(eigenvalue, expectedIndex/10);
        end
    end
end
