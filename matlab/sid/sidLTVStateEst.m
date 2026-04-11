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
%     Y - Output data, (N+1 x py), (N+1 x py x L), or cell array {L x 1}
%         where Y{l} is (N_l+1 x py). Cell arrays allow variable-length
%         trajectories.
%     U - Input data, (N x q), (N x q x L), or cell array {L x 1}
%         where U{l} is (N_l x q). Must match Y format.
%     A - Dynamics matrices, (n x n x N).
%     B - Input matrices, (n x q x N).
%     H - Observation matrix, (py x n).
%
%   NAME-VALUE OPTIONS:
%     'R' - Measurement noise covariance, (py x py) SPD. Default: eye(py).
%     'Q' - Process noise covariance, (n x n) SPD. Default: eye(n).
%
%   OUTPUTS:
%     X_hat - Estimated states, (N+1 x n x L) or cell {L x 1} where
%             X_hat{l} is (N_l+1 x n). Cell output when Y is a cell.
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
%   https://github.com/pdlourenco/sid
%  -----------------------------------------------------------------------

    % ---- Parse inputs ----
    [Y, U, A, B, H, R, Q, N, n, py, q, L, isVarLen, horizons] = ...
        parseInputs(Y, U, A, B, H, varargin{:});

    % ---- Precompute ----
    Rinv = R \ eye(py);
    Qinv = Q \ eye(n);
    HtRinvH = H' * Rinv * H;   % (n x n)
    HtRinv  = H' * Rinv;        % (n x py)

    % Number of blocks: K = N+1 (spec indices 0..N)
    K = N + 1;

    % ---- Build block tridiagonal system (SPEC.md §8.14) ----
    % Diagonal blocks S_k and off-diagonal U_k define the RTS smoother
    % equations. Shared across trajectories — only the RHS differs.
    % For variable-length trajectories, precompute all K blocks; each
    % trajectory uses only its first N_l+1 blocks.
    S_blk  = cell(K, 1);
    Uc_blk = cell(K - 1, 1);

    AtQinv = zeros(n, n, N);
    for j = 1:N
        AtQinv(:, :, j) = A(:,:,j)' * Qinv;
    end

    % Diagonal blocks (spec k=0..N, MATLAB j=1..K)
    S_blk{1} = HtRinvH + AtQinv(:,:,1) * A(:,:,1);
    for j = 2:N
        S_blk{j} = HtRinvH + Qinv + AtQinv(:,:,j) * A(:,:,j);
    end
    S_blk{K} = HtRinvH + Qinv;

    % Off-diagonal blocks (spec k=0..N-1, MATLAB j=1..N)
    for j = 1:N
        Uc_blk{j} = -AtQinv(:,:,j);
    end

    % ---- Solve per trajectory ----
    if isVarLen
        X_hat = cell(L, 1);
    else
        X_hat = zeros(K, n, L);
    end

    for l = 1:L
        if isVarLen
            Nl = horizons(l);
            Kl = Nl + 1;
            Yl = Y{l};
            Ul = U{l};
        else
            Nl = N;
            Kl = K;
            Yl = Y(:, :, l);
            Ul = U(:, :, l);
        end

        % Build RHS Theta{k} for trajectory l
        Theta = cell(Kl, 1);

        % Precompute b(j) = B(j) * u_l(j) for j=1..Nl
        b_prev = [];
        b_curr = B(:,:,1) * Ul(1, :)';

        % j=1 (spec k=0): Theta_0 = H'R^{-1}y(0) - A(0)'Q^{-1}b(0)
        Theta{1} = HtRinv * Yl(1, :)' - AtQinv(:,:,1) * b_curr;

        % j=2..Nl (spec k=1..Nl-1)
        for j = 2:Nl
            b_prev = b_curr;
            b_curr = B(:,:,j) * Ul(j, :)';
            Theta{j} = HtRinv * Yl(j, :)' ...
                + Qinv * b_prev ...
                - AtQinv(:,:,j) * b_curr;
        end

        % j=Kl (spec k=Nl): Theta_Nl = H'R^{-1}y(Nl) + Q^{-1}b(Nl-1)
        Theta{Kl} = HtRinv * Yl(Kl, :)' + Qinv * b_curr;

        % Solve (slice shared blocks to trajectory horizon)
        if Kl == K
            [w, ~] = sidLTVblkTriSolve(S_blk, Uc_blk, Theta);
        else
            [w, ~] = sidLTVblkTriSolve( ...
                S_blk(1:Kl), Uc_blk(1:Nl), Theta);
        end

        % Extract states
        if isVarLen
            X_hat{l} = zeros(Kl, n);
            for j = 1:Kl
                X_hat{l}(j, :) = w{j}';
            end
        else
            for j = 1:Kl
                X_hat(j, :, l) = w{j}';
            end
        end
    end

end

% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [Y, U, A, B, H, R, Q, N, n, py, q, L, isVarLen, horizons] = ...
    parseInputs(Y, U, A, B, H, varargin)
% PARSEINPUTS Validate and parse inputs for sidLTVStateEst.
%   Supports both 3D array input (uniform horizon) and cell array input
%   (variable-length trajectories).

    py = size(H, 1);
    n  = size(H, 2);
    q  = size(B, 2);
    N  = size(A, 3);

    if size(A, 1) ~= n || size(A, 2) ~= n
        error('sid:dimMismatch', 'A must be (%d x %d x N), got (%d x %d x %d).', ...
            n, n, size(A,1), size(A,2), size(A,3));
    end
    if size(B, 1) ~= n
        error('sid:dimMismatch', 'B must have %d rows (state dim), got %d.', n, size(B,1));
    end

    isVarLen = iscell(Y);

    if isVarLen
        % ---- Variable-length trajectory mode ----
        if ~iscell(U)
            error('sid:badInput', ...
                'When Y is a cell array, U must also be a cell array.');
        end
        L = numel(Y);
        if numel(U) ~= L
            error('sid:dimMismatch', ...
                'Y has %d trajectories but U has %d.', L, numel(U));
        end
        if L == 0
            error('sid:badInput', 'Cell arrays must not be empty.');
        end

        horizons = zeros(L, 1);
        for l = 1:L
            if size(Y{l}, 2) ~= py
                error('sid:dimMismatch', ...
                    'Y{%d} has %d columns but H has %d rows.', ...
                    l, size(Y{l}, 2), py);
            end
            if size(U{l}, 2) ~= q
                error('sid:dimMismatch', ...
                    'U{%d} has %d columns but B has %d.', ...
                    l, size(U{l}, 2), q);
            end
            Nl = size(U{l}, 1);
            if size(Y{l}, 1) ~= Nl + 1
                error('sid:dimMismatch', ...
                    'Y{%d} has %d rows but U{%d} has %d (need N_l+1 and N_l).', ...
                    l, size(Y{l}, 1), l, Nl);
            end
            if Nl > N
                error('sid:dimMismatch', ...
                    'Trajectory %d has horizon %d > size(A,3)=%d.', ...
                    l, Nl, N);
            end
            horizons(l) = Nl;
        end
    else
        % ---- Uniform-horizon mode ----
        horizons = [];

        if ndims(Y) == 2  %#ok<ISMAT>
            Y = reshape(Y, size(Y,1), size(Y,2), 1);
        end
        if ndims(U) == 2  %#ok<ISMAT>
            U = reshape(U, size(U,1), size(U,2), 1);
        end

        L = size(Y, 3);

        if size(Y, 1) ~= N + 1
            error('sid:dimMismatch', ...
                'Y must have N+1=%d rows, got %d.', N+1, size(Y,1));
        end
        if size(Y, 2) ~= py
            error('sid:dimMismatch', ...
                'Y has %d columns but H has %d rows.', size(Y,2), py);
        end
        if size(U, 1) ~= N
            error('sid:dimMismatch', ...
                'U must have N=%d rows, got %d.', N, size(U,1));
        end
        if size(U, 2) ~= q
            error('sid:dimMismatch', ...
                'U has %d columns but B has %d columns.', size(U,2), q);
        end
        if size(U, 3) ~= L
            error('sid:dimMismatch', ...
                'U has %d trajectories but Y has %d.', size(U,3), L);
        end
    end

    % Parse name-value options
    defs.R = eye(py);
    defs.Q = eye(n);
    opts = sidParseOptions(defs, varargin);
    R = opts.R;
    Q = opts.Q;

    % Validate R and Q
    if ~isequal(size(R), [py, py])
        error('sid:dimMismatch', 'R must be (%d x %d).', py, py);
    end
    if ~isequal(size(Q), [n, n])
        error('sid:dimMismatch', 'Q must be (%d x %d).', n, n);
    end
end
