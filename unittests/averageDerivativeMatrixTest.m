function tests = averageDerivativeMatrixTest
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));
end

function testAnalyticCubicAndShortInterval(tc)
    tc.verifyEqual(averagederivativematrix(3, [1,3]), ...
        [1,2,13/3,10;0,1,4,13;0,0,2,12;0,0,0,6], AbsTol=1e-12);
    tc.verifyEqual(averagederivativematrix(0, [-3,5]), 1);
    tc.verifyEqual(averagederivativematrix(2, [2,2+eps(2)]), ...
        [1,2,4;0,1,4;0,0,2], AbsTol=1e-14);
    tc.verifyError(@() averagederivativematrix(2, [1,1]), ...
        'ULMO:averagederivativematrix:InvalidRange');
    tc.verifyError(@() averagederivativematrix(2, [2,1]), ...
        'ULMO:averagederivativematrix:InvalidRange');
end

function testRawCovarianceAgainstQR(tc)
    t = linspace(0,6,80)';
    G = [ones(size(t)),t,t.^2,cos(2*pi*t),sin(2*pi*t)];
    x = [G*[2;3;.2;1;2]+.1*sin(5*t), G*[-1;.5;-.1;2;1]+.2*cos(3*t)];
    sigma = [1+.1*t,2-.1*t];
    for weighted = [false,true]
        e = []; if weighted, e = sigma; end
        [raw,~,~,~,~,C] = fittimeseries(t,x,e,2,1);
        [avg,~,~,unc,~,Cavg] = fittimeseries(t,x,e,2,1, ...
            PolynomialFormat='average-derivatives', AverageRange=[1,4]);
        tc.verifyEqual(Cavg,C);
        A = averagederivativematrix(2,[1,4]);
        for j = 1:2
            X = G; y = x(:,j);
            if weighted, X = X./sigma(:,j); y = y./sigma(:,j); end
            [Q,R] = qr(X,0); b = R\(Q'*y);
            cov = (R\eye(5))*(R\eye(5))'*sum((y-X*b).^2)/(numel(t)-5);
            tc.verifyEqual(C(:,:,j),cov(1:3,1:3),AbsTol=1e-10);
            tc.verifyEqual(avg(:,j),A*raw(:,j),AbsTol=1e-10);
            tc.verifyEqual(unc(:,j),sqrt(diag(A*cov(1:3,1:3)*A')),AbsTol=1e-10);
        end
    end
end
