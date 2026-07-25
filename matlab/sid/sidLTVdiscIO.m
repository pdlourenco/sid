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
%     Y - Output data, (N+1 x py), (N+1 x py x L), or cell array {L x 1}
%         where Y{l} is (N_l+1 x py). Cell arrays allow variable-length
%         trajectories.
%     U - Input data, (N x q), (N x q x L), or cell array {L x 1}
%         where U{l} is (N_l x q). Must match Y format.
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
%     'CovarianceMode'  - How to estimate noise covariance. Options:
%                           'diagonal'  - diagonal Sigma (default)
%                           'full'      - full n x n covariance
%                           'isotropic' - scalar * I_n
%     'TrustRegion'     - Trust-region parameter mu_0 in [0, 1], or 'off'.
%                         Default: 'off'.
%     'TrustRegionTol'  - Minimum mu before final pass. Default: 1e-6.
%
%   OUTPUTS:
%     result - Struct with fields:
%       .A               - (n x n x N) estimated dynamics matrices
%       .B               - (n x q x N) estimated input matrices
%       .X               - (N+1 x n x L) or cell {L x 1} estimated states
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
%   2026-04-06: Expose CovarianceMode option (was hardcoded 'diagonal').
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
    [Y, U, H, lambda, R, maxIter, tol, mu, muTol, doTrustRegion, ...
     covMode, N, n, py, q, L, isVarLen, horizons] = ...
        parseInputs(Y, U, H, varargin{:});

    % ---- Precompute ----
    Rinv = R \ eye(py);               % (py x py) observation precision

    % ---- Full-rank fast path (SPEC.md §8.12.2) ----
    % When rank(H) = n, state is recoverable via weighted LS — no EM needed.
    if rank(H) == n
        Hpinv = (H' * Rinv * H) \ (H' * Rinv);
        if isVarLen
            X_hat = cell(L, 1);
            for l = 1:L
                X_hat{l} = Y{l} * Hpinv';
            end
        else
            X_hat = zeros(N + 1, n, L);
            for l = 1:L
                X_hat(:, :, l) = reshape(Y(:, :, l), N + 1, py) * Hpinv';
            end
        end

        [A, B, S_c, D_c, Xl_c, C_c] = cosmicStep( ...
            X_hat, U, lambda, N, n, q, L, isVarLen, horizons);

        J = evaluateFullCost( ...
            X_hat, A, B, Y, U, H, Rinv, ...
            lambda, N, n, q, L, isVarLen, horizons);
        result = packResult( ...
            A, B, X_hat, H, R, J, 0, lambda, ...
            N, n, py, q, L, isVarLen, horizons);
        result = addUncertainty( ...
            result, S_c, D_c, Xl_c, C_c, lambda, N, n, q, covMode, ...
            isVarLen, horizons);
        return;
    end

    % ---- LTI Initialisation (SPEC.md §8.12.4) ----
    % Estimate constant dynamics (A0, B0) from the I/O transfer function
    % via Ho-Kalman realization. This gives an observable initialisation
    % for any H (including py < n).
    [A0, B0] = sidLTIfreqIO(Y, U, H);
    A_init = repmat(A0, [1, 1, N]);
    B_init = repmat(B0, [1, 1, N]);
    A0_rep = repmat(A0, [1, 1, N]);  % trust-region target A0

    % ---- Two-level alternating loop (SPEC.md §8.12.3 inner / §8.12.4 outer) ----
    args = {A0_rep, Y, U, H, R, Rinv, lambda, N, n, q, L, ...
        isVarLen, horizons, maxIter, tol};
    costHistory = [];
    nIter = 0;

    if ~doTrustRegion || mu <= 0
        % Base two-block alternation at mu = 0 (SPEC §8.12.3).
        [best, iters, converged, cv] = innerLoop(0, A_init, B_init, args{:});
        costHistory = [costHistory, cv];
        nIter = nIter + iters;
        if ~converged
            warning('sid:notConverged', ...
                ['sidLTVdiscIO: alternating loop did not converge in %d ' ...
                 'iterations (relative cost change stayed above tolerance).'], ...
                maxIter);
        end
    else
        % ---- Two-level trust-region schedule (SPEC.md §8.12.4) ----
        mu_current = mu;
        % Stage at the initial mu (typically 1): trust the LTI initialisation.
        [accepted, iters, ~, cv] = innerLoop( ...
            mu_current, A_init, B_init, args{:});
        costHistory = [costHistory, cv];
        nIter = nIter + iters;
        % Outer loop: halve mu, accept only if J*(mu') <= J*(mu); reject stops.
        while mu_current > muTol
            mu_prime = mu_current / 2;
            [best_p, iters_p, ~, cv_p] = innerLoop( ...
                mu_prime, accepted.A, accepted.B, args{:});
            costHistory = [costHistory, cv_p];  %#ok<AGROW>
            nIter = nIter + iters_p;
            if best_p.J <= accepted.J
                mu_current = mu_prime;
                accepted = best_p;
            else
                break;  % reject: restore accepted, stop reducing mu (step 4)
            end
        end
        % Guarded final mu = 0 refinement (both exits): keep only if J* drops.
        [best_zero, iters_z, ~, cv_z] = innerLoop( ...
            0, accepted.A, accepted.B, args{:});
        costHistory = [costHistory, cv_z];
        nIter = nIter + iters_z;
        if best_zero.J < accepted.J
            best = best_zero;
        else
            best = accepted;
        end
    end

    A = best.A;  B = best.B;  X_hat = best.X;
    S_c = best.S;  D_c = best.D;  Xl_c = best.Xl;  C_c = best.C;
    % Report the achieved cost as the final history entry (the best iterate may
    % not be the last one evaluated when mu > 0 stages ran).
    if isempty(costHistory) || costHistory(end) ~= best.J
        costHistory(end + 1) = best.J;
    end

    % ---- Pack result struct ----
    result = packResult( ...
        A, B, X_hat, H, R, costHistory(:), nIter, lambda, ...
        N, n, py, q, L, isVarLen, horizons);

    % ---- Bayesian uncertainty from final COSMIC step (SPEC.md §8.12.9) ----
    result = addUncertainty( ...
        result, S_c, D_c, Xl_c, C_c, lambda, N, n, q, covMode, ...
        isVarLen, horizons);

end

% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [best, iters, converged, costVec] = innerLoop( ...
    muVal, A_start, B_start, A0_rep, ...
    Y, U, H, R, Rinv, lambda, N, n, q, L, isVarLen, horizons, maxIter, tol)
% INNERLOOP Alternating state-COSMIC loop at a fixed mu (SPEC §8.12.4 step 2).
%   Runs from (A_start, B_start) until the relative cost change falls below TOL
%   (SPEC §8.12.3) or the per-stage cap MAXITER is reached. Returns the BEST
%   (lowest-cost) iterate observed — not necessarily the last, since for mu>0
%   the state step uses A_use != A so the inner iteration is a homotopy
%   fixed-point, not the monotone descent of §8.12.3. COSTVEC is this stage's
%   per-iteration cost history; CONVERGED flags whether the relative test fired.

    A_k = A_start;
    B_k = B_start;
    best = struct('J', Inf);   % populated on iteration 1 (maxIter >= 1)
    Jprev = [];
    iters = 0;
    converged = false;
    costVec = zeros(1, maxIter);

    for it = 1:maxIter
        % -- E-step: state estimation (interpolated dynamics when mu > 0) --
        if muVal > 0
            A_use = (1 - muVal) * A_k + muVal * A0_rep;
        else
            A_use = A_k;
        end
        X_hat = sidLTVStateEst(Y, U, A_use, B_k, H, 'R', R);

        % -- M-step: COSMIC solve (always free A(k), B(k)) --
        [A_k, B_k, S_c, D_c, Xl_c, C_c] = cosmicStep( ...
            X_hat, U, lambda, N, n, q, L, isVarLen, horizons);

        % -- Evaluate cost --
        J = evaluateFullCost( ...
            X_hat, A_k, B_k, Y, U, H, Rinv, ...
            lambda, N, n, q, L, isVarLen, horizons);
        iters = iters + 1;
        costVec(iters) = J;

        % -- Track the best (lowest-cost) iterate of this stage. Capture on
        %    iteration 1 unconditionally: NaN <= Inf is false in MATLAB, so a
        %    stage whose first cost is NaN would otherwise leave `best` without
        %    its A/B/... fields (matching Python's `best is None or ...`). --
        if it == 1 || J <= best.J
            best.J = J;   best.X = X_hat;   best.A = A_k;   best.B = B_k;
            best.S = S_c;  best.D = D_c;  best.Xl = Xl_c;  best.C = C_c;
        end

        % -- Relative-cost convergence (SPEC §8.12.3), with the small-cost floor
        %    that keeps the test from demanding ~machine-precision absolute steps
        %    once J < 1 (pure |dJ|/|J| there fails to converge on ordinary data).
        if ~isempty(Jprev)
            relChange = abs(J - Jprev) / max(abs(Jprev), 1);
            if relChange < tol
                converged = true;
                break;
            end
        end
        Jprev = J;
    end

    costVec = costVec(1:iters);
end

function result = packResult( ...
    A, B, X, H, R, cost, nIter, lambda, ...
    N, n, py, q, L, isVarLen, horizons)
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
    if isVarLen
        result.Horizons    = horizons;
    end
end

function result = addUncertainty( ...
    result, S, D, Xl, C, lambda, N, n, q, covMode, isVarLen, horizons)
% ADDUNCERTAINTY Append Bayesian uncertainty fields (SPEC.md §8.12.9).
%
%   Computes AStd, BStd from the block-tridiagonal Hessian inverse,
%   reusing the S matrix from the final COSMIC step.

    d = n + q;

    % Diagonal blocks of the Hessian inverse (SPEC.md §8.9.2)
    P = sidLTVuncertaintyBackwardPass(S, lambda, N, d);

    % Noise covariance from COSMIC residuals
    [Sigma, dof] = sidEstimateNoiseCov(C, D, Xl, P, covMode, N, n, q);

    % Standard deviations of A(k) and B(k) entries
    [AStd, BStd] = sidExtractStd(P, Sigma, N, n, q);

    result.AStd             = AStd;
    result.BStd             = BStd;
    result.P                = P;
    result.NoiseCov         = Sigma;
    result.NoiseCovEstimated = true;
    result.NoiseVariance    = trace(Sigma) / n;
    result.DegreesOfFreedom = dof;
end

% Local functions estimateNoiseCovLocal and extractStdLocal have been
% replaced by shared private helpers: sidEstimateNoiseCov.m and
% sidExtractStd.m

function [Y, U, H, lambda, R, maxIter, tol, mu, muTol, doTrustRegion, ...
          covMode, N, n, py, q, L, isVarLen, horizons] = parseInputs(Y, U, H, varargin)
% PARSEINPUTS Validate and parse inputs for sidLTVdiscIO.
%   Supports both 3D array input (uniform horizon) and cell array input
%   (variable-length trajectories).

    py = size(H, 1);
    n  = size(H, 2);

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

        q = size(U{1}, 2);
        horizons = zeros(L, 1);
        for l = 1:L
            if size(Y{l}, 2) ~= py
                error('sid:dimMismatch', ...
                    'Y{%d} has %d columns but H has %d rows.', ...
                    l, size(Y{l}, 2), py);
            end
            if size(U{l}, 2) ~= q
                error('sid:dimMismatch', ...
                    'U{%d} has %d columns, expected %d.', ...
                    l, size(U{l}, 2), q);
            end
            Nl = size(U{l}, 1);
            if size(Y{l}, 1) ~= Nl + 1
                error('sid:dimMismatch', ...
                    'Y{%d} has %d rows but U{%d} has %d (need N_l+1 and N_l).', ...
                    l, size(Y{l}, 1), l, Nl);
            end
            if Nl < 2
                error('sid:tooShort', ...
                    'Trajectory %d has fewer than 3 measurements.', l);
            end
            horizons(l) = Nl;
        end

        N = max(horizons);
    else
        % ---- Uniform-horizon mode ----
        horizons = [];

        if ndims(Y) == 2  %#ok<ISMAT>
            Y = reshape(Y, size(Y,1), size(Y,2), 1);
        end
        if ndims(U) == 2  %#ok<ISMAT>
            U = reshape(U, size(U,1), size(U,2), 1);
        end

        N = size(U, 1);
        q = size(U, 2);
        L = size(Y, 3);

        if size(Y, 1) ~= N + 1
            error('sid:dimMismatch', ...
                'Y must have N+1=%d rows, got %d.', N+1, size(Y,1));
        end
        if size(Y, 2) ~= py
            error('sid:dimMismatch', ...
                'Y has %d columns but H has %d rows.', size(Y,2), py);
        end
        if size(U, 3) ~= L
            error('sid:dimMismatch', ...
                'U has %d trajectories but Y has %d.', size(U,3), L);
        end
    end

    % Parse name-value options
    defs.Lambda = [];
    defs.R = eye(py);
    defs.MaxIter = 50;
    defs.Tolerance = 1e-6;
    defs.TrustRegion = 'off';
    defs.TrustRegionTol = 1e-6;
    defs.CovarianceMode = 'diagonal';
    opts = sidParseOptions(defs, varargin);
    lambda = opts.Lambda;
    R = opts.R;
    maxIter = opts.MaxIter;
    if ~isscalar(maxIter) || maxIter < 1 || maxIter ~= floor(maxIter)
        error('sid:badInput', ...
            'MaxIter must be a positive integer, got %g.', maxIter);
    end
    tol = opts.Tolerance;
    muTol = opts.TrustRegionTol;
    if ischar(opts.TrustRegion) && strcmpi(opts.TrustRegion, 'off')
        doTrustRegion = false;
        mu = 1;
    else
        doTrustRegion = true;
        mu = opts.TrustRegion;
        if ~isnumeric(mu) || ~isscalar(mu) || mu < 0 || mu > 1
            error('sid:badInput', ...
                'TrustRegion must be a scalar in [0, 1] or ''off''.');
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

    % Validate CovarianceMode
    covMode = lower(opts.CovarianceMode);
    if ~ismember(covMode, {'full', 'diagonal', 'isotropic'})
        error('sid:badCovMode', ...
            ['CovarianceMode must be ''full'', ''diagonal'', ' ...
            'or ''isotropic''. Got ''%s''.'], opts.CovarianceMode);
    end
end

function [A, B, S, D, Xl, C] = cosmicStep( ...
    X_hat, U, lambda, N, n, q, L, isVarLen, horizons)
% COSMICSTEP Standard COSMIC solve on estimated states.
%
%   Treats X_hat as observed states and solves for C(k) = [A(k)'; B(k)'].
%   Returns intermediates S, D, Xl, C for optional uncertainty computation.

    if isVarLen
        [D, Xl] = sidLTVbuildDataMatricesVarLen( ...
            X_hat, U, N, n, q, L, horizons);
    else
        [D, Xl] = sidLTVbuildDataMatrices(X_hat, U, N, n, q, L);
    end
    [S, T]  = sidLTVbuildBlockTerms(D, Xl, lambda, N, n, q);
    [C, ~]  = sidLTVcosmicSolve(S, T, lambda, N, n, q);

    A = permute(C(1:n, :, :), [2 1 3]);       % (n x n x N)
    B = permute(C(n+1:end, :, :), [2 1 3]);   % (n x q x N)
end

function J = evaluateFullCost( ...
    X_hat, A, B, Y, U, H, Rinv, lambda, N, n, q, L, isVarLen, horizons)
% EVALUATEFULLCOST Compute full Output-COSMIC objective.
%
%   J = obs_fidelity + dyn_fidelity + smoothness

    obs_fidelity = 0;
    dyn_fidelity = 0;
    smoothness = 0;

    for l = 1:L
        if isVarLen
            Nl = horizons(l);
            Yl = Y{l};
            Xl = X_hat{l};
            Ul = U{l};
        else
            Nl = N;
            Yl = Y(:, :, l);
            Xl = X_hat(:, :, l);
            Ul = U(:, :, l);
        end

        for k = 0:Nl
            j = k + 1;
            res_obs = Yl(j, :)' - H * Xl(j, :)';
            obs_fidelity = obs_fidelity + res_obs' * Rinv * res_obs;
        end

        for k = 0:Nl-1
            j = k + 1;
            res_dyn = Xl(j+1, :)' ...
                - A(:, :, j) * Xl(j, :)' ...
                - B(:, :, j) * Ul(j, :)';
            dyn_fidelity = dyn_fidelity + res_dyn' * res_dyn;
        end
    end

    % Smoothness: N * lambda(k) * ||C(k+1) - C(k)||^2_F. The COSMIC step applies
    % the 1/sqrt(N) data scaling (§8.3.2), so the objective it actually minimises
    % -- hence the monotone, reported one -- uses the effective weight N*lambda,
    % not lambda (SPEC.md §8.12.2, issue #137). The user's knob stays lambda.
    for k = 1:N-1
        Ck  = [A(:,:,k)';  B(:,:,k)'];
        Ck1 = [A(:,:,k+1)'; B(:,:,k+1)'];
        smoothness = smoothness + N * lambda(k) * norm(Ck1 - Ck, 'fro')^2;
    end

    J = obs_fidelity + dyn_fidelity + smoothness;
end
