% Test the integration against expanding on a Fibonacci grid and summing
function plm2avgp_demo2
    dom = 'australia';
    Lmax = 20;
    pars = 10;
    fprintf("Testing integration of EGM2008 (Lmax=%d) over %s\n", Lmax, dom);
    % Get coordinates
    XY = GeoDomain(dom, pars).Lonlat;
    % Make sure the coordinates make sense
    XY(:, 1) = XY(:, 1) - 360 * (XY(:, 1) > 360);
    % Use the first Lmax data from EGM2008 as a field of interest
    v = fralmanac('EGM2008_ZeroTide', 'SHM');
    % Note that gravity does not start at zero
    % Geoid = 3 % Free-air gravity anomaly = 2
    v = plm2pot(v(1:addmup(Lmax) - addmup(v(1) - 1), :), [], [], [], 3);
    % Create a Fibonacci grid
    [lonF, latF] = Fibonacci_grid(30000);
    % Expand EGM2008 on the Fib Grid
    [rF, ~, ~, ~] = plm2xyz(v(:, 1:4), latF, lonF);
    % Now decide if we're inside or outside of the region, and set outside
    % points to zero
    rF(~inpolygon(lonF, latF, XY(:, 1), XY(:, 2))) = 0;
    % Now compute the integration which can be
    % discretized due to the equal area Fibonacci grid
    IntF = sum(rF) / length(rF);
    % Do the same integration with plm2avg
    [Int, A] = plm2avgp(v(:, 1:4), dom);
    % Now check the averages.
    indx = find(rF);
    AF = sum(rF(indx)) / length(indx);

    fprintf('\nIntegration check... PLM2AVG should equal the Fib Grid for good resolution');
    fprintf('PLM2AVG Int: %6.7f ; FIB GRID Int: %6.7f ; ERROR: %6.2f%%\n', ...
        Int, IntF, (Int - IntF) / Int * 100)
    fprintf('Avg value check... PLM2AVG should equal the Fib Grid for good resolution\n');
    fprintf('PLM2AVG Avg: %6.7f ; FIB GRID Avg: %6.7f ; ERROR: %6.2f%%', ...
        A, AF, (A - AF) / A * 100)
end
