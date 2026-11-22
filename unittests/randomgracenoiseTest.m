%% RANDOMGRACENOISETEST - Checks dated GRACE noise draws against known monthly factors.
% Fixtures use non-diagonal factors and independent GRACEDATA and IFILES roots.
%
% Created by
%   2026/11/10
classdef randomgracenoiseTest < matlab.unittest.TestCase
    properties
        GraceRoot
        InputRoot
    end
    properties (Constant)
        Product = {'GFZ', 'RL06', 96}
        Dates = [datetime(2003,1,16,12,0,0); datetime(2018,6,16)]
    end

    methods (TestMethodSetup)
        function isolateRoots(testCase)
            import matlab.unittest.fixtures.*
            temp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.GraceRoot = fullfile(temp.Folder, 'processed');
            testCase.InputRoot = fullfile(temp.Folder, 'inputs');
            testCase.applyFixture(EnvironmentVariableFixture('GRACEDATA', testCase.GraceRoot));
            testCase.applyFixture(EnvironmentVariableFixture('IFILES', testCase.InputRoot));
        end
    end

    methods (Test)
        function seededSeriesMatchesBothFactorsAndAllSlots(testCase)
            factors = fixture(testCase.GraceRoot, 3, 1);
            stream = RandStream('mt19937ar','Seed',17);
            reference = RandStream('mt19937ar','Seed',17);
            globalStream = RandStream.getGlobalStream();
            state = globalStream.State;
            [noise, dates] = randomgracenoise(testCase.Product, [], 3, ...
                'NumRealisations', 4, 'Stream', stream);
            testCase.verifyEqual(dates, testCase.Dates);
            testCase.verifySize(noise, [10 4 2 4]);
            % Degree 2 and 3 parameters as [lm row, C/S output column].
            rows = [4 5 5 6 6 7 8 8 9 9 10 10];
            cols = [3 3 4 3 4 3 3 4 3 4 3 4];
            slots = sub2ind([10,4], rows, cols);
            for iMonth = 1:2
                expected = factors(:,:,iMonth) * randn(reference,12,4);
                for iEns = 1:4
                    field = [0 0; 1 0; 1 1; 2 0; 2 1; 2 2; 3 0; 3 1; 3 2; 3 3];
                    field(:,3:4) = 0;
                    field(slots) = expected(:,iEns);
                    testCase.verifyEqual(noise(:,:,iMonth,iEns), field, 'AbsTol', 1e-12);
                end
            end
            testCase.verifyEqual(stream.State, reference.State);
            testCase.verifyEqual(globalStream.State, state);
        end

        function choosesSmallestSufficientFileAndRetainsLeadingBlock(testCase)
            fallbackRoot = fullfile(testCase.InputRoot,'GRACE');
            fallback = fixture(fallbackRoot, 2, 10);
            fixture(testCase.GraceRoot, 4, 4);
            larger = fixture(testCase.GraceRoot, 3, 3);
            testCase.checkDegreeTwoDraw(larger(1:5,1:5,2));
            exact = fixture(testCase.GraceRoot, 2, 2);
            testCase.checkDegreeTwoDraw(exact(:,:,2));
            setenv('GRACEDATA', '');
            testCase.checkDegreeTwoDraw(fallback(:,:,2));
        end

        function rangeSelectionAndStreamContinuation(testCase)
            factors = fixture(testCase.GraceRoot, 2, 1);
            stream = RandStream('mt19937ar','Seed',9);
            reference = RandStream('mt19937ar','Seed',9);
            serial = convertTo(testCase.Dates, 'datenum');
            ranges = {[-Inf serial(1)], [serial(2) Inf], repmat(testCase.Dates(2),1,2)};
            epochIndex = [1 2 2];
            for iRange = 1:numel(ranges)
                idx = epochIndex(iRange);
                [noise, dates] = randomgracenoise(testCase.Product, ranges{iRange}, 2, 'Stream', stream);
                testCase.verifyEqual(dates, testCase.Dates(idx));
                expected = factors(:,:,idx)*randn(reference,5,1);
                testCase.verifyEqual(noise(4,3), expected(1), 'AbsTol', 1e-13);
                testCase.verifyEqual(noise(5,4), expected(3), 'AbsTol', 1e-13);
                testCase.verifyEqual(stream.State, reference.State);
            end
        end

        function emptySelectionSkipsFactorIOAndRandomDraws(testCase)
            % A time-only degree-60 file lets factor reads fail if attempted.
            timeFile(testCase.GraceRoot, 60);
            stream = RandStream('mt19937ar','Seed',12);
            initial = stream.State;
            [noise, dates] = randomgracenoise(testCase.Product, ...
                [datetime(2010,1,1), datetime(2010,2,1)], 'NumRealisations', 2, 'Stream', stream);
            testCase.verifySize(noise, [1891,4,0,2]);
            testCase.verifySize(dates, [0,1]);
            testCase.verifyEqual(stream.State, initial);
        end

        function usesCurrentGlobalStreamByDefault(testCase)
            factors = fixture(testCase.GraceRoot, 2, 1);
            old = RandStream.getGlobalStream();
            testCase.addTeardown(@() RandStream.setGlobalStream(old));
            globalStream = RandStream('mt19937ar','Seed',37);
            RandStream.setGlobalStream(globalStream);
            reference = RandStream('mt19937ar','Seed',37);
            [noise, dates] = randomgracenoise(testCase.Product, testCase.Dates, 2);
            testCase.verifyEqual(dates, testCase.Dates);
            for iMonth = 1:2
                expected = factors(:,:,iMonth)*randn(reference,5,1);
                testCase.verifyEqual(noise(4,3,iMonth), expected(1), 'AbsTol', 1e-13);
            end
            testCase.verifyEqual(globalStream.State, reference.State);
        end

        function reportsUnsupportedProductAndMissingDegrees(testCase)
            testCase.verifyError(@() randomgracenoise(testCase.Product), ...
                'ULMO:randomgracenoise:MissingFactors');
            fixture(testCase.GraceRoot, 3, 1);
            fixture(testCase.GraceRoot, 2, 1);
            for product = {{'GFZ','RL06',60}, {'CSR','RL06',96}, {'GFZ','RL05',96}}
                testCase.verifyError(@() randomgracenoise(product{1}, [], 2), ...
                    'ULMO:randomgracenoise:UnsupportedProduct');
            end
            try
                randomgracenoise(testCase.Product, [], 5);
            catch exception
                testCase.assertEqual(exception.identifier,'ULMO:randomgracenoise:MissingFactors');
                testCase.verifySubstring(exception.message, 'degree 5');
                testCase.verifySubstring(exception.message, '[2 3]');
                testCase.verifySubstring(exception.message, fullfile(testCase.GraceRoot,'GFZ','NEQ'));
                return
            end
            testCase.assertFail('Expected missing factors for degree 5.');
        end
    end

    methods (Access = private)
        function checkDegreeTwoDraw(testCase, factor)
            stream = RandStream('mt19937ar','Seed',5);
            reference = RandStream('mt19937ar','Seed',5);
            noise = randomgracenoise(testCase.Product, repmat(testCase.Dates(2),1,2), 2, 'Stream', stream);
            testCase.verifySize(noise, [6 4]);
            testCase.verifyEqual(noise(:,1), [0;1;1;2;2;2]);
            values = [noise(4,3);noise(5,3);noise(5,4);noise(6,3);noise(6,4)];
            testCase.verifyEqual(values, factor*randn(reference,5,1), 'AbsTol', 1e-13);
        end
    end
end

function path = timeFile(root, L)
    folder = fullfile(root, 'GFZ', 'NEQ');
    if ~isfolder(folder)
        mkdir(folder);
    end
    path = fullfile(folder, sprintf('GFZ_NEQ_L%d.nc',L));
    nccreate(path,'time','Dimensions',{'time',2},'Format','netcdf4');
    ncwrite(path,'time',seconds(randomgracenoiseTest.Dates - datetime(2000,1,1)));
    ncwriteatt(path,'time','units','seconds since 2000-01-01 00:00:00 UTC');
end

function factors = fixture(root, L, multiplier)
    path = timeFile(root, L);
    n = (L+1)^2-4;
    nPacked = n*(n+1)/2;
    labels = zeros(n,3);
    last = 0;
    for degree = 2:L
        orders = [0, repelem(1:degree,2)]';
        cs = [0; repmat([0;1],degree,1)];
        idx = last + (1:numel(orders));
        labels(idx,:) = [repmat(degree,numel(orders),1), orders, cs];
        last = idx(end);
    end
    names = ["degree","order","cs"];
    for col = 1:3
        nccreate(path,names(col),'Dimensions',{'parameter',n});
        ncwrite(path,names(col),labels(:,col));
    end
    nccreate(path,'varianceFactor','Dimensions',{'time',2});
    ncwrite(path,'varianceFactor',[7;9]);
    nccreate(path,'noiseFactor','Dimensions',{'packed',nPacked,'time',2});
    factors = zeros(n,n,2);
    for iMonth = 1:2
        factor = tril(-reshape(1:n^2,n,n)/(10*n),-1) + diag(1:n);
        factor = factor * (multiplier+iMonth);
        factors(:,:,iMonth) = factor;
        ncwrite(path,'noiseFactor',factor(tril(true(n))),[1 iMonth]);
    end
end
