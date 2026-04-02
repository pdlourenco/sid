function X_hat = sidLTVStateEst(Y, U, A, B, H, varargin)
% SIDLTVSTATEEST Batch LTV state estimation (RTS smoother).
%
%   X_hat = sidLTVStateEst(Y, U, A, B, H)
%   X_hat = sidLTVStateEst(Y, U, A, B, H, 'R', R, 'Q', Q)
%
%   Estimates state trajectories for a discrete-time LTV system with
%   partial observations by minimising:
%
%     J = sum_k ||y(k) - H x(k)||^2_{R^{-1}}
%       + sum_k ||x(k+1) - A(k) x(k) - B(k) u(k)||^2_{Q^{-1}}
%
%   This is equivalent to a Rauch-Tung-Striebel (RTS) fixed-interval
%   smoother. Solved via a block tridiagonal forward-backward pass in
%   O(N n^3) per trajectory.
%
%   INPUTS:
%     Y - Output data, (N+1 x py) or (N+1 x py x L).
%     U - Input data, (N x q) or (N x q x L).
%     A - Dynamics matrices, (n x n x N).
%     B - Input matrices, (n x q x N).
%     H - Observation matrix, (py x n).
%
%   NAME-VALUE OPTIONS:
%     'R' - Measurement noise covariance, (py x py) SPD. Default: eye(py).
%     'Q' - Process noise covariance, (n x n) SPD. Default: eye(n).
%
%   OUTPUTS:
%     X_hat - Estimated states, (N+1 x n x L).
%
%   EXAMPLES:
%     % State estimation with known dynamics
%     X_hat = sidLTVStateEst(Y, U, A, B, H);
%
%     % With known measurement and process noise
%     X_hat = sidLTVStateEst(Y, U, A, B, H, 'R', R_meas, 'Q', Q_proc);
%
%   ALGORITHM:
%     Block tridiagonal Gaussian elimination per docs/cosmic_output.md
%     Appendix A. Per trajectory, builds n x n diagonal blocks S{k},
%     off-diagonal U{k} = -A(k)' Q^{-1}, and RHS Theta{k}, then
%     calls sidLTVblkTriSolve.
%
%   SPECIFICATION:
%     SPEC.md section 8.12 -- Output-COSMIC state step
%     docs/cosmic_output.md -- Appendix A
%
%   See also: sidLTVblkTriSolve, sidLTVdiscIO, sidLTVdisc
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

    % ---- Parse inputs ----
    [Y, U, A, B, H, R, Q, N, n, py, q, L] = parseInputs(Y, U, A, B, H, varargin{:});

    % ---- Precompute ----
    Rinv = R \ eye(py);
    Qinv = Q \ eye(n);
    HtRinvH = H' * Rinv * H;   % (n x n)
    HtRinv  = H' * Rinv;        % (n x py)

    % Number of blocks: K = N+1 (spec indices 0..N)
    K = N + 1;

    % ---- Precompute b(k) = B(k) u_l(k) for each trajectory ----
    % b{j, l} = B(:,:,j) * U(j,:,l)', j=1..N, l=1..L
    % (store in 3D array for efficiency)
    b_all = zeros(n, L, N);
    for j = 1:N
        Bj = B(:, :, j);
        for l = 1:L
            b_all(:, l, j) = Bj * U(j, :, l)';
        end
    end

    % ---- Build block tridiagonal system (SPEC.md §8.14) ----
    % Diagonal blocks S_k and off-diagonal U_k define the RTS smoother
    % equations. Shared across trajectories — only the RHS differs.
    S_blk  = cell(K, 1);
    Uc_blk = cell(K - 1, 1);

    % Diagonal blocks (spec k=0..N, MATLAB j=1..K)
    % j=1 (spec k=0): S_0 = H'R^{-1}H + A(0)'Q^{-1}A(0)
    S_blk{1} = HtRinvH + A(:,:,1)' * Qinv * A(:,:,1);

    % j=2..N (spec k=1..N-1): S_k = H'R^{-1}H + Q^{-1} + A(k)'Q^{-1}A(k)
    for j = 2:N
        S_blk{j} = HtRinvH + Qinv + A(:,:,j)' * Qinv * A(:,:,j);
    end

    % j=K (spec k=N): S_N = H'R^{-1}H + Q^{-1}
    S_blk{K} = HtRinvH + Qinv;

    % Off-diagonal blocks (spec k=0..N-1, MATLAB j=1..N)
    % U_k = -A(k)' Q^{-1}
    for j = 1:N
        Uc_blk{j} = -A(:,:,j)' * Qinv;
    end

    % ---- Solve per trajectory ----
    X_hat = zeros(K, n, L);

    for l = 1:L
        % Build RHS Theta{k} for trajectory l
        Theta = cell(K, 1);

        % j=1 (spec k=0): Theta_0 = H'R^{-1}y(0) - A(0)'Q^{-1}b(0)
        Theta{1} = HtRinv * Y(1, :, l)' - A(:,:,1)' * Qinv * b_all(:, l, 1);

        % j=2..N (spec k=1..N-1):
        % Theta_k = H'R^{-1}y(k) + Q^{-1}b(k-1) - A(k)'Q^{-1}b(k)
        for j = 2:N
            Theta{j} = HtRinv * Y(j, :, l)' ...
                + Qinv * b_all(:, l, j-1) ...
                - A(:,:,j)' * Qinv * b_all(:, l, j);
        end

        % j=K (spec k=N): Theta_N = H'R^{-1}y(N) + Q^{-1}b(N-1)
        Theta{K} = HtRinv * Y(K, :, l)' + Qinv * b_all(:, l, N);

        % Solve
        [w, ~] = sidLTVblkTriSolve(S_blk, Uc_blk, Theta);

        % Extract states
        for j = 1:K
            X_hat(j, :, l) = w{j}';
        end
    end

end

% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [Y, U, A, B, H, R, Q, N, n, py, q, L] = parseInputs(Y, U, A, B, H, varargin)
% PARSEINPUTS Validate and parse inputs for sidLTVStateEst.

    py = size(H, 1);
    n  = size(H, 2);
    q  = size(B, 2);
    N  = size(A, 3);

    % Ensure 3D
    if ndims(Y) == 2  %#ok<ISMAT>
        Y = reshape(Y, size(Y,1), size(Y,2), 1);
    end
    if ndims(U) == 2  %#ok<ISMAT>
        U = reshape(U, size(U,1), size(U,2), 1);
    end

    L = size(Y, 3);

    % Validate dimensions
    if size(Y, 1) ~= N + 1
        error('sid:dimMismatch', 'Y must have N+1=%d rows, got %d.', N+1, size(Y,1));
    end
    if size(Y, 2) ~= py
        error('sid:dimMismatch', 'Y has %d columns but H has %d rows.', size(Y,2), py);
    end
    if size(U, 1) ~= N
        error('sid:dimMismatch', 'U must have N=%d rows, got %d.', N, size(U,1));
    end
    if size(U, 2) ~= q
        error('sid:dimMismatch', 'U has %d columns but B has %d columns.', size(U,2), q);
    end
    if size(U, 3) ~= L
        error('sid:dimMismatch', 'U has %d trajectories but Y has %d.', size(U,3), L);
    end
    if size(A, 1) ~= n || size(A, 2) ~= n
        error('sid:dimMismatch', 'A must be (%d x %d x N), got (%d x %d x %d).', ...
            n, n, size(A,1), size(A,2), size(A,3));
    end
    if size(B, 1) ~= n
        error('sid:dimMismatch', 'B must have %d rows (state dim), got %d.', n, size(B,1));
    end

    % Parse name-value options
    R = eye(py);
    Q = eye(n);

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'r'
                    R = varargin{k+1};
                    k = k + 2;
                case 'q'
                    Q = varargin{k+1};
                    k = k + 2;
                otherwise
                    error('sid:badOption', 'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badOption', 'Expected option name (string), got %s.', class(varargin{k}));
        end
    end

    % Validate R and Q
    if ~isequal(size(R), [py, py])
        error('sid:dimMismatch', 'R must be (%d x %d).', py, py);
    end
    if ~isequal(size(Q), [n, n])
        error('sid:dimMismatch', 'Q must be (%d x %d).', n, n);
    end
end
