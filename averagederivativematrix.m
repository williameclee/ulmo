% AVERAGEDERIVATIVEMATRIX Transform ascending polynomial coefficients to
% interval-average derivatives of orders 0:p. Interval endpoints are numeric
% offsets in the polynomial's time units, measured from its origin.
function A = averagederivativematrix(p, interval)
    arguments
        p (1, 1) double {mustBeInteger, mustBeNonnegative}
        interval (1, 2) double {mustBeReal, mustBeFinite}
    end
    if interval(2) <= interval(1)
        error('ULMO:averagederivativematrix:InvalidRange', ...
            'The averaging interval must have positive length.');
    end
    % A(k+1,j+1) is the interval average of d^k(t^j)/dt^k.
    a = interval(1);
    b = interval(2);
    A = zeros(p + 1);

    for k = 0:p

        for j = k:p
            n = j - k;
            % Divided difference of t^(n+1), evaluated without subtracting
            % nearby endpoint powers or dividing by a small interval length.
            meanPower = sum(a .^ (0:n) .* b .^ (n:-1:0)) / (n + 1);
            A(k + 1, j + 1) = prod((j - k + 1):j) * meanPower;
        end

    end

end
