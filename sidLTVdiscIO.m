function result = sidLTVdiscIO(Y, U, H, varargin)
% SIDLTVDISCIO Identify discrete-time LTV system from partial observations.
%
%   result = sidLTVdiscIO(Y, U, H, 'Lambda', lambda)
%   result = sidLTVdiscIO(Y, U, H, 'Lambda', lambda, 'R', R)
%   result = sidLTVdiscIO(Y, U, H, 'Lambda', lambda, 'TrustRegion', 1)
%
%   Identifies time-varying system matrices A(k), B(k) and estimates
%   state trajectories from input-output data when only partial state
%   observations are available:
%
%       x(k+1) = A(k) x(k) + B(k) u(k)       k = 0, ..., N-1
%       y(k)   = H x(k)
%
%   Uses the Output-COSMIC algorithm: alternating minimisation between
%   a COSMIC step (dynamics estimation) and an RTS smoother (state
%   estimation), initialised via LTI realization from the I/O transfer
%   function (sidLTIfreqIO).
%
%   When H = I (full state observation), this reduces to the standard
%   COSMIC algorithm (sidLTVdisc).
%
%   INPUTS:
%     Y - Output data, (N+1 x py) or (N+1 x py x L).
%     U - Input data, (N x q) or (N x q x L).
%     H - Observation matrix, (py x n). When rank(H) = n (including
%         py >= n), states are recovered exactly via weighted least
%         squares and no EM iterations are needed.
%
%   NAME-VALUE OPTIONS:
%     'Lambda'          - Regularisation strength. Scalar or (N-1 x 1).
%                         Required.
%     'R'               - Measurement noise covariance, (py x py) SPD.
%                         Default: eye(py).
%     'MaxIter'         - Maximum alternating iterations. Default: 50.
%     'Tolerance'       - Convergence tolerance on relative cost change.
%                         Default: 1e-6.
%     'TrustRegion'     - Trust-region parameter mu_0 in [0, 1], or 'off'.
%                         Default: 'off'.
%     'TrustRegionTol'  - Minimum mu before final pass. Default: 1e-6.
%
%   OUTPUTS:
%     result - Struct with fields:
%       .A               - (n x n x N) estimated dynamics matrices
%       .B               - (n x q x N) estimated input matrices
%       .X               - (N+1 x n x L) estimated state trajectories
%       .H               - (py x n) observation matrix (copy)
%       .R               - (py x py) noise covariance used
%       .Cost            - (n_iter x 1) cost J at each iteration
%       .Iterations      - scalar, number of alternating iterations
%       .Lambda          - (N-1 x 1) regularisation used
%       .DataLength      - N
%       .StateDim        - n
%       .OutputDim       - py
%       .InputDim        - q
%       .NumTrajectories - L
%       .Algorithm       - 'cosmic'
%       .Method          - 'sidLTVdiscIO'
%
%   EXAMPLES:
%     % Basic usage
%     result = sidLTVdiscIO(Y, U, H, 'Lambda', 1e5);
%
%     % With known measurement noise
%     result = sidLTVdiscIO(Y, U, H, 'Lambda', 1e5, 'R', R_meas);
%
%     % With trust-region for difficult convergence
%     result = sidLTVdiscIO(Y, U, H, 'Lambda', 1e5, 'TrustRegion', 1);
%
%     % Multi-trajectory
%     result = sidLTVdiscIO(Y_3d, U_3d, H, 'Lambda', 1e5);
%
%   ALGORITHM:
%     Output-COSMIC alternating minimisation:
%       1. LTI initialisation: estimate A0, B0 via Ho-Kalman realization
%          of the I/O transfer function (sidLTIfreqIO)
%       2. State step: fix dynamics, RTS smoother for states
%       3. COSMIC step: fix states, solve for A(k), B(k)
%       4. Repeat 2-3 until convergence
%     Complexity: O(T * (L*N*n^3 + N*(n+q)^3)) where T is iterations.
%
%   REFERENCES:
%     Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
%     identification from large-scale data for LTV systems."
%     arXiv:2112.04355, 2022.
%
%   SPECIFICATION:
%     SPEC.md section 8.12 -- Output-COSMIC
%     docs/cosmic_output.md -- Full derivation
%
%   See also: sidLTIfreqIO, sidLTVdisc, sidLTVStateEst, sidLTVdiscFrozen
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
    [Y, U, H, lambda, R, maxIter, tol, mu, muTol, doTrustRegion, ...
     N, n, py, q, L] = parseInputs(Y, U, H, varargin{:});

    % ---- Precompute ----
    Rinv = R \ eye(py);               % (py x py) observation precision

    % ---- Full-rank fast path (SPEC.md §8.12.2) ----
    % When rank(H) = n, state is recoverable via weighted LS — no EM needed.
    if rank(H) == n
        Hpinv = (H' * Rinv * H) \ (H' * Rinv);
        X_hat = zeros(N + 1, n, L);
        for l = 1:L
            X_hat(:, :, l) = squeeze(Y(:, :, l)) * Hpinv';
        end

        [A, B] = cosmicStep(X_hat, U, lambda, N, n, q, L);

        J = evaluateFullCost( ...
            X_hat, A, B, Y, U, H, Rinv, ...
            lambda, N, n, q, L);
        result = packResult( ...
            A, B, X_hat, H, R, J, 0, lambda, ...
            N, n, py, q, L);
        return;
    end

    % ---- LTI Initialisation (SPEC.md §8.12.4) ----
    % Estimate constant dynamics (A0, B0) from the I/O transfer function
    % via Ho-Kalman realization. This gives an observable initialisation
    % for any H (including py < n).
    [A0, B0] = sidLTIfreqIO(Y, U, H);
    A = repmat(A0, [1, 1, N]);
    B = repmat(B0, [1, 1, N]);
    A0_rep = repmat(A0, [1, 1, N]);  % trust-region target

    % ---- Alternating state-COSMIC loop (SPEC.md §8.12.3) ----
    mu_current = doTrustRegion * mu;
    costHistory = [];
    nIter = 0;

    % For trust-region accept/reject
    J_converged_mu = Inf;
    X_best = [];  A_best = A;  B_best = B;

    for iter = 1:maxIter
        % -- E-step: state estimation --
        if mu_current > 0
            A_use = (1 - mu_current) * A + mu_current * A0_rep;
        else
            A_use = A;
        end
        X_hat = sidLTVStateEst(Y, U, A_use, B, H, 'R', R);

        % -- M-step: COSMIC solve --
        [A, B] = cosmicStep(X_hat, U, lambda, N, n, q, L);

        % -- Evaluate cost and check convergence --
        J = evaluateFullCost(X_hat, A, B, Y, U, H, Rinv, lambda, N, n, q, L);
        costHistory(end + 1) = J;  %#ok<AGROW>
        nIter = nIter + 1;

        if nIter >= 2
            J_prev = costHistory(end - 1);
            relChange = abs(J - J_prev) / max(abs(J_prev), 1);

            if relChange < tol
                if doTrustRegion && mu_current > muTol
                    if J <= J_converged_mu
                        J_converged_mu = J;
                        X_best = X_hat;  A_best = A;  B_best = B;
                        mu_current = mu_current / 2;
                    else
                        X_hat = X_best;  A = A_best;  B = B_best;
                        mu_current = 0;
                    end
                elseif doTrustRegion && mu_current > 0
                    J_converged_mu = J;
                    X_best = X_hat;  A_best = A;  B_best = B;
                    mu_current = 0;
                else
                    break;
                end
            end
        end
    end

    % ---- Pack result struct ----
    result = packResult( ...
        A, B, X_hat, H, R, costHistory(:), nIter, lambda, ...
        N, n, py, q, L);

end

% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function result = packResult( ...
    A, B, X, H, R, cost, nIter, lambda, ...
    N, n, py, q, L)
% PACKRESULT Build the output struct (shared by both code paths).

    result.A               = A;
    result.B               = B;
    result.X               = X;
    result.H               = H;
    result.R               = R;
    result.Cost            = cost;
    result.Iterations      = nIter;
    result.Lambda          = lambda;
    result.DataLength      = N;
    result.StateDim        = n;
    result.OutputDim       = py;
    result.InputDim        = q;
    result.NumTrajectories = L;
    result.Algorithm       = 'cosmic';
    result.Method          = 'sidLTVdiscIO';
end

function [Y, U, H, lambda, R, maxIter, tol, mu, muTol, doTrustRegion, ...
          N, n, py, q, L] = parseInputs(Y, U, H, varargin)
% PARSEINPUTS Validate and parse inputs for sidLTVdiscIO.

    py = size(H, 1);
    n  = size(H, 2);

    % Ensure 3D
    if ndims(Y) == 2  %#ok<ISMAT>
        Y = reshape(Y, size(Y,1), size(Y,2), 1);
    end
    if ndims(U) == 2  %#ok<ISMAT>
        U = reshape(U, size(U,1), size(U,2), 1);
    end

    N = size(U, 1);
    q = size(U, 2);
    L = size(Y, 3);

    % Validate dimensions
    if size(Y, 1) ~= N + 1
        error('sid:dimMismatch', 'Y must have N+1=%d rows, got %d.', N+1, size(Y,1));
    end
    if size(Y, 2) ~= py
        error('sid:dimMismatch', 'Y has %d columns but H has %d rows.', size(Y,2), py);
    end
    if size(U, 3) ~= L
        error('sid:dimMismatch', 'U has %d trajectories but Y has %d.', size(U,3), L);
    end

    % Defaults
    lambda = [];
    R = eye(py);
    maxIter = 50;
    tol = 1e-6;
    mu = 1;
    muTol = 1e-6;
    doTrustRegion = false;

    % Parse name-value options
    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'lambda'
                    lambda = varargin{k+1};
                    k = k + 2;
                case 'r'
                    R = varargin{k+1};
                    k = k + 2;
                case 'maxiter'
                    maxIter = varargin{k+1};
                    k = k + 2;
                case 'tolerance'
                    tol = varargin{k+1};
                    k = k + 2;
                case 'trustregion'
                    val = varargin{k+1};
                    if ischar(val) && strcmpi(val, 'off')
                        doTrustRegion = false;
                    else
                        doTrustRegion = true;
                        mu = val;
                    end
                    k = k + 2;
                case 'trustregiontol'
                    muTol = varargin{k+1};
                    k = k + 2;
                otherwise
                    error('sid:badOption', 'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badOption', 'Expected option name (string), got %s.', class(varargin{k}));
        end
    end

    % Validate lambda
    if isempty(lambda)
        error('sid:badInput', 'Lambda is required. Use ''Lambda'', value.');
    end
    if isscalar(lambda)
        lambda = lambda * ones(N - 1, 1);
    end
    if length(lambda) ~= N - 1
        error('sid:dimMismatch', ...
            'Lambda must be scalar or (N-1 x 1), got length %d.', ...
            length(lambda));
    end
    lambda = lambda(:);

    % Validate R
    if ~isequal(size(R), [py, py])
        error('sid:dimMismatch', 'R must be (%d x %d).', py, py);
    end
end

function [A, B] = cosmicStep(X_hat, U, lambda, N, n, q, L)
% COSMICSTEP Standard COSMIC solve on estimated states.
%
%   Treats X_hat as observed states and solves for C(k) = [A(k)'; B(k)'].

    [D, Xl] = sidLTVbuildDataMatrices(X_hat, U, N, n, q, L);
    [S, T]  = sidLTVbuildBlockTerms(D, Xl, lambda, N, n, q);
    [C, ~]  = sidLTVcosmicSolve(S, T, lambda, N, n, q);

    A = permute(C(1:n, :, :), [2 1 3]);       % (n x n x N)
    B = permute(C(n+1:end, :, :), [2 1 3]);   % (n x q x N)
end

function J = evaluateFullCost(X_hat, A, B, Y, U, H, Rinv, lambda, N, n, q, L)
% EVALUATEFULLCOST Compute full Output-COSMIC objective.
%
%   J = obs_fidelity + dyn_fidelity + smoothness

    obs_fidelity = 0;
    dyn_fidelity = 0;
    smoothness = 0;

    for l = 1:L
        for k = 0:N
            j = k + 1;
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
            Ak = A(:, :, j);
            Bk = B(:, :, j);
            res_dyn = xl1 - Ak * xl - Bk * ul;
            dyn_fidelity = dyn_fidelity + res_dyn' * res_dyn;
        end
    end

    % Smoothness: lambda(k) * ||C(k+1) - C(k)||^2_F
    for k = 1:N-1
        Ck  = [A(:,:,k)';  B(:,:,k)'];
        Ck1 = [A(:,:,k+1)'; B(:,:,k+1)'];
        smoothness = smoothness + lambda(k) * norm(Ck1 - Ck, 'fro')^2;
    end

    J = obs_fidelity + dyn_fidelity + smoothness;
end
