%% NONICESHEETLANDTEST
% Tests domain subtraction and options with deterministic coastline polygons.
%
% Created by
%   2026/10/01, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/10/01, En-Chi Lee (williameclee@arizona.edu)

function tests=nonIceSheetLandTest
    tests=functiontests(localfunctions);
end

function setupOnce(tc)
    root=fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(root));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'domains')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'aux')));
    temp=tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
    % Only source polygons are mocked; GeoDomain and boolean geometry are real.
    write(temp.Folder,'alloceans',[
        "function [xy,p]=alloceans(varargin)"
        "ip=inputParser;ip.KeepUnmatched=true;addParameter(ip,'LonOrigin',180);addParameter(ip,'Buffer',0);parse(ip,varargin{:});o=ip.Results.LonOrigin;b=ip.Results.Buffer;"
        "w=polyshape(o+[-180,180,180,-180],[-90,-90,90,90]);"
        "land=polyshape([80,140,140,80]+[-b,b,b,-b],[-20,-20,40,40]+[-b,-b,b,b]);"
        "land=union(land,polyshape([300,320,320,300],[60,60,80,80]));"
        "land=union(land,polyshape([0,360,360,0],[-90,-90,-65,-65]));"
        "land=union(land,translate(land,[-360,0]));land=union(land,translate(land,[360,0]));p=subtract(w,land);[x,y]=boundary(p);xy=[x,y];end"]);
    write(temp.Folder,'greenland',[
        "function xy=greenland(varargin)"
        "xy=[300,60;320,60;320,80;300,80;300,60];xy=flipud(xy);end"]);
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(temp.Folder));
end

function testExcludesIceAndOcean(tc)
    [xy,p]=nonicesheetland('LonOrigin',180,'SaveData',false);
    tc.verifyTrue(isinterior(p,110,0));
    tc.verifyFalse(isinterior(p,310,70));
    tc.verifyFalse(isinterior(p,220,-80));
    tc.verifyFalse(any(isinterior(p,[1,180,359],[-89,-89,-89])));
    tc.verifyFalse(isinterior(p,220,0));
    tc.verifyEqual(area(p),3600,AbsTol=1e-8);
    tc.verifyEqual(area(polyshape(xy)),area(p),AbsTol=1e-8);
    tc.verifyFalse(nonicesheetland('rotated'));
end

function testBufferAndLatitudeClip(tc)
    [~,p]=nonicesheetland('Buffer',2,'Latlim',[-10,10],'LonOrigin',180,'SaveData',false);
    tc.verifyEqual(area(p),64*20,AbsTol=1e-8);
    tc.verifyTrue(isinterior(p,79,0));
    tc.verifyFalse(isinterior(p,110,20));
    tc.verifyFalse(isinterior(p,220,0));
end

function testLongitudeWindow(tc)
    [xy,p]=nonicesheetland('LonOrigin',0,'SaveData',false);
    tc.verifyEqual(area(p),3600,AbsTol=1e-8);
    tc.verifyFalse(isinterior(p,-50,70));
    tc.verifyFalse(isinterior(p,-140,-80));
    tc.verifyTrue(all(abs(xy(isfinite(xy(:,1)),1))<=180));
end

function testGeoDomainIdentity(tc)
    d=GeoDomain('nonicesheetland','Buffer',.5);
    tc.verifyTrue(contains(d.Id,'Nonicesheetland'));
    tc.verifyEqual(d.DisplayName('long'),'Land outside Greenland and Antarctica');
    xy=d.Lonlat(180);tc.verifyTrue(all(isfinite(xy(~isnan(xy)))));
end

function testPositionalAndEmptyDefaults(tc)
    [~,named]=nonicesheetland('Upscale',0,'Buffer',2,'Latlim',10,'SaveData',false);
    [~,positional]=nonicesheetland(0,2,10,{},180,false,'SaveData',false);
    tc.verifyEqual(area(positional),area(named),AbsTol=1e-8);
    [~,defaults]=nonicesheetland([],[],[],{},[],'SaveData',false);
    tc.verifyEqual(area(defaults),3600,AbsTol=1e-8);
    tc.verifyFalse(nonicesheetland("rotated"));
end

function testInvalidLatitudeLimits(tc)
    tc.verifyError(@()nonicesheetland('Latlim',[10,-10]), ...
        'ULMO:nonicesheetland:InvalidLatlim');
    tc.verifyError(@()nonicesheetland('Latlim',[-20,0,20]), ...
        'ULMO:nonicesheetland:InvalidLatlim');
end

function write(folder,name,lines)
    fid=fopen(fullfile(folder,[name,'.m']),'w');c=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',lines);
end
