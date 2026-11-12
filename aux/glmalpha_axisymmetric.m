%% GLMALPHA_AXISYMMETRIC
% Is an auxiliary function separated from GLMALPHA. 
% It finds G, V, and N for the given axisymmetric region.
% Uses the fixed-order polar-cap routines from the Slepian package.
%
% Authored by:
% 	En-Chi Lee <williameclee@arizona.edu>, 2024-07-12
%
% Last modified
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function [G, V, EL, EM, N, GM2AL, MTAP, IMTAP] = ...
        glmalpha_axisymmetric(TH, sord, L, lp, bp, EM, EL, blkm, ...
        blox, upco, xver, resc)

    maxL = max(L);
    ldim = (maxL + 1) ^ 2 - bp * L(1) ^ 2;
    mesg = NaN;
    MTAP = zeros(1, ldim);
    IMTAP = zeros(1, ldim);

    G = zeros((maxL + 1) ^ 2, ldim);
    V = zeros(1, ldim);

    % Find increasing column index; that's how many belong to this order
    % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
    % The middle bit is twice for every nonzero order missing
    % alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
    %   		gamini(L(2-lp)-bp*(L(1)-1),bp*2*(L(1)-1)) ...
    %   		gamini(L(2-lp)-bp*(L(1)-1):-1:1,2)]);
    % This should be the same for L and [0 L]
    alpha = cumsum( ...
        [1 L(2 - lp) - bp * L(1) + 1 ...
         gamini(L(2 - lp) - bp * (L(1) - 1), bp * 2 * L(1)) ...
         gamini(L(2 - lp) - bp * L(1):-1:1, 2)]);

    if blox ~= 0 && blox ~= 1
        error('Specify valid block-sorting option ''blox''')
    end

    % For the SINGLE or DOUBLE POLAR CAPS
    for m = 0:maxL
        % Same results for +/- m; set nth=0 thus no sign correction!
        if sord == 1

            if lp
                [~, ~, ~, C, ~, Vp] = grunbaum(TH, L, m, 0);
            elseif bp
                % Note that the small-eigenvalue eigenfunctions might be
                % numerically degenerate and thus not as using Grunbaum - if
                % you need to compare, compare where the eigenvalues are "big"
                [~, Vp, ~, ~, C] = sdwcap(TH, L, m, 0, -1);
            end

        elseif sord == 2

            if lp
                [~, ~, ~, C, ~, Vp] = grunbaum2(TH, L, m, 0);
            elseif bp
                error('Bandpass double-cap tapers not ready yet')
            end

        else
            error('Specify single or double polar cap')
        end

        if upco ~= 0

            if upco > 0
                % The upward continuation matrix
                A = diag((1 + upco) .^ (- (m:L) - 1));
            elseif upco < 0
                % The downward continuation matrix
                A = diag((1 + abs(upco)) .^ ((m:L) + 1));
            end

            % Comparisons with Grunbaum only make sense for lowpass
            if xver == 1 && lp
                % This should give the same result, more or less, less accurate
                if sord == 1
                    [~, Vs, ~, ~, Cs, ~, ~, ~, ~, ~, D] = sdwcap(TH, L, m, 0, -1);
                else
                    [~, Vs, ~, Cs, ~, ~, D] = sdwcap2(TH, L, m, 0, -1);
                end

                % This should give the eigenvalues again, which we'd had from
                % orthocheck
                warning off
                % Check difference integration and kernel eigenvalues
                difer(Vp(:) - diag((C' * D * C) ./ (C' * C)), [], [], mesg)
                % Check difference integration and diagonalization eigenvalues
                difer(Vp(:) - Vs(:), [], [], mesg)
                % Check difference between the eigenfunctions barring the sign
                % and only wherever the eigenvalues amount to anything
                difer(abs(Cs(:, Vp > 1e-6)) - abs(C(:, Vp > 1e-6)), [], [], mesg)
                warning on
                % Vc = diag((C' * A * C * diag(Vp) * C' * A * C));
                Vp0 = Vp;
            end

            % Upward continuation from 1 to 1+a or from 1+a to 1:
            % New eigenfunctions, same name
            C = A * C;
            % Calculate new eigenvalues, same name
            [~, ~, ~, Vp, nofa] = orthocheck(C, [], TH / 180 * pi, m, sord, 1);

            % Make sure these are sorted, since that's not automatically the case
            % [Vp,ind]=sort(Vp,'descend');
            % C=C(:,ind);
            % Current thinking is: do NOT resort, as you'll want to compare the
            % best at a=0 with whatever it becomes later!

            if xver == 1
                warning off
                % Check difference integration eigenvalues and those from kernel
                difer(Vp(:) - diag((C' * D * C) ./ (C' * C)), [], [], mesg)
                warning on
                % Check how many Vp>Vp0
                fprintf('%i/%i eigenvalues greater\n', sum(Vp(:) > Vp0(:)), ...
                    length(Vp0))
                fprintf('Shannon number greater: %i\n', sum(Vp) > sum(Vp0))
            end

            if resc == 1
                % Rescale these functions to have an integral to unity over the
                % sphere; note: this doesn't make the set orthonormal of course
                C = C * diag(1 ./ nofa);
            end

        end

        % Distribute this at the right point in the huge matrix
        if m > 0
            % Here you supply the negative orders
            G(EM == -m, alpha(2 * m):alpha(2 * m + 1) - 1) = C;
            V(alpha(2 * m):alpha(2 * m + 1) - 1) = Vp;
            MTAP(alpha(2 * m):alpha(2 * m + 1) - 1) = -m;
            % It's all neatly ordered here, downgoing within every order
            IMTAP(alpha(2 * m):alpha(2 * m + 1) - 1) = 1:length(Vp);
        end

        % Duplicate for the positive order in case the region is axisymmetric
        G(EM == m, alpha(2 * m + 1):alpha(2 * m + 2) - 1) = C;
        V(alpha(2 * m + 1):alpha(2 * m + 2) - 1) = Vp;
        MTAP(alpha(2 * m + 1):alpha(2 * m + 2) - 1) = m;
        % It's all neatly ordered here, downgoing within every order
        IMTAP(alpha(2 * m + 1):alpha(2 * m + 2) - 1) = 1:length(Vp);
    end

    % Calculate the Shannon number and compare it to the theory
    N = sum(V);

    if upco == 0
        difer(N - ldim * sord * (1 - cos(TH / 180 * pi)) / 2, [], [], mesg);
    end

    % Compute the sum over all orders of the squared coefficients
    % Thus works when they have not been blocksorted yet.
    GM2AL = zeros(ldim, maxL + 1);

    for l = 0:maxL
        b = (l - 1 + 1) ^ 2 + 1;
        e = (l + 1) ^ 2;
        GM2AL(:, l + 1) = sum(G(b:e, :) .^ 2, 1)';
    end

    % Make sure that the sum over all degrees is 1 - but I forgot why
    difer(sum(GM2AL, 2) - 1, [], [], mesg)

    % This is not blockdiagonal, unless you make it thus
    if blox == 1
        G = G(blkm, :);
        EM = EM(blkm);
        EL = EL(blkm);
    end

end
