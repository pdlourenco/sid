function [X_hat, A_hat, B_hat, cost] = sidLTVdiscIOInit( ...
    Y, U, H, Rinv, HtRinvH, HtRinv, lambda, N, n, py, q, L)
% SIDLTVDISCIOINIT Initialisation for Output-COSMIC (solve J|_{A=I}).
%
%   [X_hat, A_hat, B_hat, cost] = sidLTVdiscIOInit(Y, U, H, Rinv, ...
%       HtRinvH, HtRinv, lambda, N, n, py, q, L)
%
%   Jointly estimates state sequences {x_l(k)} and input matrices {B(k)}
%   by minimising the Output-COSMIC objective with A(k) = I for all k.
%   This is a convex quadratic solved in a single forward-backward pass
%   over composite blocks.
%
%   INPUTS:
%     Y        - (N+1 x py x L) output measurements.
%     U        - (N x q x L) input data.
%     H        - (py x n) observation matrix.
%     Rinv     - (py x py) inverse measurement noise covariance.
%     HtRinvH  - (n x n) precomputed H' * Rinv * H.
%     HtRinv   - (n x py) precomputed H' * Rinv.
%     lambda   - (N-1 x 1) regularisation weights.
%     N        - Number of time steps.
%     n        - State dimension.
%     py       - Output dimension.
%     q        - Input dimension.
%     L        - Number of trajectories.
%
%   OUTPUTS:
%     X_hat - (N+1 x n x L) estimated state trajectories.
%     A_hat - (n x n x N) dynamics matrices (all set to eye(n)).
%     B_hat - (n x q x N) estimated input matrices.
%     cost  - Scalar, value of J at the initialisation solution.
%
%   ALGORITHM:
%     Composite block tridiagonal solve per docs/cosmic_output.md
%     Appendix B. Block size is (Ln + nq) for k=0..N-1 and Ln for k=N.
%
%   SPECIFICATION:
%     SPEC.md section 8.12.3 -- Initialisation
%
%   See also: sidLTVblkTriSolve, sidLTVdiscIO, sidLTVStateEst
%
%   Changelog:
%   2026-04-01: First version by Pedro Lourenço.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    In = eye(n);
    Inq = eye(n * q);
    ILn = eye(L * n);

    % Number of composite blocks: K = N+1 (spec indices 0..N)
    K = N + 1;
    nq = n * q;
    Ln = L * n;

    % ---- Precompute E(k) and P(k) for each time step ----
    % E{k}: (Ln x nq), stacked kron(u_l(k)', I_n) across trajectories
    % P{k}: (nq x nq), pooled input Gram matrix
    % Spec k=0..N-1 maps to MATLAB time index k+1 for U array
    E = cell(N, 1);
    P = cell(N, 1);
    for k = 1:N
        Ek = zeros(Ln, nq);
        Pk = zeros(nq, nq);
        for l = 1:L
            ul = U(k, :, l)';  % (q x 1)
            el = kron(ul', In);  % (n x nq)
            Ek((l-1)*n+1:l*n, :) = el;
            Pk = Pk + kron(ul * ul', In);
        end
        E{k} = Ek;
        P{k} = Pk;
    end

    % ---- Precompute common diagonal sub-blocks ----
    HRH_I  = HtRinvH + In;      % for boundary blocks (k=0, k=N)
    HRH_2I = HtRinvH + 2 * In;  % for interior blocks (0 < k < N)

    % ---- Build S, U_blk, Theta cell arrays for sidLTVblkTriSolve ----
    S_blk     = cell(K, 1);
    U_blk     = cell(K - 1, 1);
    Theta_blk = cell(K, 1);

    % -- Diagonal blocks S{k} --
    % MATLAB index j=1..K corresponds to spec k=0..N
    % j=1 (spec k=0):
    S_blk{1} = [kron(eye(L), HRH_I),  E{1}
                E{1}',                  P{1} + lambda(1) * Inq];

    % j=2..N (spec k=1..N-1, interior):
    for j = 2:N
        Ediff = E{j} - E{j-1};
        % B-block gets lambda from both adjacent smoothness terms, except
        % at j=N (spec k=N-1) which has no right neighbour.
        if j < N
            lam_sum = lambda(j-1) + lambda(j);
        else
            lam_sum = lambda(j-1);
        end
        S_blk{j} = [kron(eye(L), HRH_2I),  Ediff
                     Ediff',                 P{j} + lam_sum * Inq];
    end

    % j=K=N+1 (spec k=N): no B block
    S_blk{K} = kron(eye(L), HRH_I);

    % -- Off-diagonal blocks U_blk{j} --
    % j=1..N-1 (spec k=0..N-2): size (Ln+nq) x (Ln+nq)
    for j = 1:N-1
        U_blk{j} = [-ILn,      zeros(Ln, nq)
                     -E{j}',    -lambda(j) * Inq];
    end

    % j=N (spec k=N-1): size (Ln+nq) x Ln (no B column at final step)
    U_blk{N} = [-ILn
                -E{N}'];

    % -- Right-hand side Theta{j} --
    % j=1..N (spec k=0..N-1): size (Ln+nq) x 1
    for j = 1:N
        rhs_x = zeros(Ln, 1);
        for l = 1:L
            yl = Y(j, :, l)';  % (py x 1), Y at MATLAB row j = spec time k=j-1
            rhs_x((l-1)*n+1:l*n) = HtRinv * yl;
        end
        Theta_blk{j} = [rhs_x; zeros(nq, 1)];
    end

    % j=K=N+1 (spec k=N): size Ln x 1
    rhs_x = zeros(Ln, 1);
    for l = 1:L
        yl = Y(K, :, l)';
        rhs_x((l-1)*n+1:l*n) = HtRinv * yl;
    end
    Theta_blk{K} = rhs_x;

    % ---- Solve composite block tridiagonal system ----
    [w, ~] = sidLTVblkTriSolve(S_blk, U_blk, Theta_blk);

    % ---- Extract states and B(k) from composite solution ----
    X_hat = zeros(K, n, L);  % (N+1 x n x L)
    B_hat = zeros(n, q, N);  % (n x q x N)
    A_hat = repmat(In, [1, 1, N]);  % (n x n x N), all identity

    % j=1..N (spec k=0..N-1): w{j} has size (Ln+nq)
    for j = 1:N
        wj = w{j};
        for l = 1:L
            X_hat(j, :, l) = wj((l-1)*n+1:l*n)';
        end
        B_hat(:, :, j) = reshape(wj(Ln+1:end), n, q);
    end

    % j=K=N+1 (spec k=N): w{K} has size Ln
    wK = w{K};
    for l = 1:L
        X_hat(K, :, l) = wK((l-1)*n+1:l*n)';
    end

    % ---- Evaluate cost at initialisation ----
    cost = evaluateInitCost(X_hat, B_hat, Y, U, H, Rinv, lambda, N, n, q, L);

end

function cost = evaluateInitCost(X_hat, B_hat, Y, U, H, Rinv, lambda, N, n, q, L)
% Evaluate full J at the initialisation (A = I).

    obs_fidelity = 0;
    dyn_fidelity = 0;
    smoothness = 0;

    for l = 1:L
        for k = 0:N
            j = k + 1;  % MATLAB index
            yl = Y(j, :, l)';
            xl = X_hat(j, :, l)';
            res_obs = yl - H * xl;
            obs_fidelity = obs_fidelity + res_obs' * Rinv * res_obs;
        end

        for k = 0:N-1
            j = k + 1;
            xl = X_hat(j, :, l)';
            xl1 = X_hat(j+1, :, l)';
            ul = U(j, :, l)';
            Bk = B_hat(:, :, j);
            res_dyn = xl1 - xl - Bk * ul;  % A = I
            dyn_fidelity = dyn_fidelity + res_dyn' * res_dyn;
        end
    end

    % Smoothness: lambda(k) * ||C(k+1) - C(k)||^2_F
    % With A = I, C(k) = [I; B(k)'], so C(k+1) - C(k) = [0; B(k+1)' - B(k)']
    for k = 1:N-1
        Bdiff = B_hat(:, :, k+1) - B_hat(:, :, k);
        smoothness = smoothness + lambda(k) * norm(Bdiff, 'fro')^2;
    end

    cost = obs_fidelity + dyn_fidelity + smoothness;
end
