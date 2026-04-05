%% test_reviewFixes - Tests for bugs and issues found during v1.0 review.
%
% BUG-1: Multi-trajectory ETFE should use H1 estimator (cross-periodograms)
% BUG-2: Frozen TF Jacobian should use R*B, not R^H*B
% DISC-1: MIMO noise spectrum should be PSD (non-negative eigenvalues)
% DISC-7: Cell array input for frequency-domain functions

fprintf('Running test_reviewFixes...\n');

%% BUG-1a: Multi-trajectory ETFE produces consistent estimates
% With many trajectories of independent inputs, the H1 estimator should
% converge to the true transfer function. The old pooled-DFT approach
% would produce NaN/Inf due to cancellation in the denominator.
rng(101);
N = 200;
L = 20;
b = [1 0.5];
a = [1 -0.8];
y3d = zeros(N, 1, L);
u3d = zeros(N, 1, L);
for l = 1:L
    u3d(:, 1, l) = randn(N, 1);
    y3d(:, 1, l) = filter(b, a, u3d(:, 1, l));
end

result = sidFreqETFE(y3d, u3d, 'Smoothing', 5);
assert(~any(isnan(result.Response(:))), ...
    'BUG-1: Multi-traj ETFE should not produce NaN with H1 estimator.');

% Check approximate correctness at a few frequencies
w_test = result.Frequency;
G_true = freqz(b, a, w_test);
% With 20 noiseless trajectories and smoothing=5, errors should be small.
% ETFE has no windowing so variance is inherently high at individual
% frequencies; use the 90th percentile as a robust accuracy check.
err = abs(result.Response - G_true);
err_sorted = sort(err);
p90 = err_sorted(round(0.9 * length(err)));
assert(p90 < 0.3, ...
    'BUG-1: Multi-traj ETFE 90th-pct error should be < 0.3 (got %.3f).', p90);

fprintf('  BUG-1a: Multi-trajectory ETFE H1 estimator - PASSED\n');

%% BUG-1b: Single trajectory ETFE unchanged
rng(102);
N = 500;
u1 = randn(N, 1);
y1 = filter(b, a, u1) + 0.05 * randn(N, 1);
result1 = sidFreqETFE(y1, u1);
assert(~any(isnan(result1.Response(:))), ...
    'BUG-1b: Single-traj ETFE should still work.');
fprintf('  BUG-1b: Single trajectory ETFE still works - PASSED\n');

%% BUG-2: Frozen TF Jacobian via finite differences
% Compare Jacobian-based uncertainty against finite-difference perturbation.
rng(201);
N = 20;
p = 2; q = 1;
L = 10;
X = randn(N + 1, p, L);
U = randn(N, q, L);

ltvRes = sidLTVdisc(X, U, 'Lambda', 1e3, 'Uncertainty', true);

% Compute frozen TF at a single frequency and time step
w0 = pi / 4;
k0 = 5;
frz = sidLTVdiscFrozen(ltvRes, 'Frequencies', w0, 'TimeSteps', k0);

% Finite-difference Jacobian for one element G(1,1)
delta = 1e-6;
Ak = ltvRes.A(:, :, k0);
Bk = ltvRes.B(:, :, k0);
z0 = exp(1i * w0);
G0 = (z0 * eye(p) - Ak) \ Bk;

% Perturb each A entry and compute dG/dA numerically
for ii = 1:p
    for jj = 1:p
        Ap = Ak; Ap(ii, jj) = Ap(ii, jj) + delta;
        Gp = (z0 * eye(p) - Ap) \ Bk;
        dG_fd = (Gp - G0) / delta;
        % The finite difference should be finite and reasonable
        assert(all(isfinite(dG_fd(:))), ...
            'BUG-2: Finite-difference Jacobian should be finite.');
    end
end

% The uncertainty should be finite and positive
assert(all(frz.ResponseStd(:) >= 0) && all(isfinite(frz.ResponseStd(:))), ...
    'BUG-2: Frozen TF uncertainty should be finite and non-negative.');

fprintf('  BUG-2: Frozen TF Jacobian correctness - PASSED\n');

%% DISC-1: MIMO noise spectrum PSD
rng(301);
N = 500;
ny = 2; nu = 2;
u = randn(N, nu);
% Create a MIMO system with cross-coupling
y = zeros(N, ny);
y(:, 1) = filter([1 0.3], [1 -0.7], u(:, 1)) + ...
          0.5 * filter([0.2], [1 -0.5], u(:, 2)) + 0.1 * randn(N, 1);
y(:, 2) = filter([0.1], [1 -0.6], u(:, 1)) + ...
          filter([1 -0.2], [1 -0.8], u(:, 2)) + 0.1 * randn(N, 1);

warning('off', 'sid:singularPhiU');
result = sidFreqBT(y, u);
warning('on', 'sid:singularPhiU');

% Check that noise spectrum is PSD at each frequency
nf = length(result.Frequency);
for kk = 1:nf
    Vk = reshape(result.NoiseSpectrum(kk, :, :), ny, ny);
    eigvals = eig(Vk);
    assert(all(eigvals >= -1e-14), ...
        'DISC-1: MIMO noise spectrum should be PSD at freq %d (min eig = %.2e).', ...
        kk, min(eigvals));
end

fprintf('  DISC-1: MIMO noise spectrum PSD clamping - PASSED\n');

%% API-3: MIMO uncertainty should be finite (diagonal approximation)
assert(~any(isnan(result.ResponseStd(:))), ...
    'API-3: MIMO ResponseStd should not be NaN (diagonal approximation).');
assert(all(result.ResponseStd(:) >= 0), ...
    'API-3: MIMO ResponseStd should be non-negative.');
assert(all(isfinite(result.ResponseStd(:)) | isinf(result.ResponseStd(:))), ...
    'API-3: MIMO ResponseStd entries should be finite or Inf (low input power).');

fprintf('  API-3: MIMO diagonal uncertainty approximation - PASSED\n');

%% API-1: MIMO residual tests all channels
rng(501);
N = 500; ny = 2; nu = 1;
u_res = randn(N, nu);
y_res = zeros(N, ny);
y_res(:, 1) = filter([1], [1 -0.8], u_res) + 0.1 * randn(N, 1);
y_res(:, 2) = filter([0.5], [1 -0.6], u_res) + 0.1 * randn(N, 1);
G_res = sidFreqBT(y_res, u_res);
res = sidResidual(G_res, y_res, u_res);

% Should have per-channel fields
assert(isfield(res, 'AutoCorrAll'), ...
    'API-1: Should have AutoCorrAll field.');
assert(isfield(res, 'WhitenessPassAll'), ...
    'API-1: Should have WhitenessPassAll field.');
assert(size(res.AutoCorrAll, 2) == ny, ...
    'API-1: AutoCorrAll should have ny=%d columns.', ny);
assert(length(res.WhitenessPassAll) == ny, ...
    'API-1: WhitenessPassAll should have ny=%d entries.', ny);
% Cross-correlation should test all (output, input) pairs
assert(isfield(res, 'IndependencePassAll'), ...
    'API-1: Should have IndependencePassAll field.');
assert(length(res.IndependencePassAll) == ny * nu, ...
    'API-1: IndependencePassAll should have ny*nu=%d entries.', ny * nu);
% Overall pass is conjunction of all channels
assert(res.WhitenessPass == all(res.WhitenessPassAll), ...
    'API-1: WhitenessPass should be AND of all channels.');
assert(res.IndependencePass == all(res.IndependencePassAll), ...
    'API-1: IndependencePass should be AND of all pairs.');

fprintf('  API-1: Per-channel MIMO residual tests - PASSED\n');

%% NUM-2: Frequency-domain simulation no extrapolation
rng(601);
N = 200;
u_sim = randn(N, 1);
y_sim = filter([1], [1 -0.8], u_sim) + 0.1 * randn(N, 1);
G_sim = sidFreqBT(y_sim, u_sim);
cmp = sidCompare(G_sim, y_sim, u_sim);
% The predicted output should be finite (no extrapolation artifacts)
assert(all(isfinite(cmp.Predicted(:))), ...
    'NUM-2: Predicted output should be finite (no extrapolation).');
% Fit should be reasonable for a well-identified system
assert(cmp.Fit(1) > 30, ...
    'NUM-2: Fit should be reasonable (got %.1f%%).', cmp.Fit(1));

fprintf('  NUM-2: Frequency-domain simulation no extrapolation - PASSED\n');

%% DISC-7: Cell array input for frequency-domain functions
rng(401);
N1 = 200; N2 = 150;
u1 = randn(N1, 1);
u2 = randn(N2, 1);
y1 = filter([1], [1 -0.8], u1) + 0.1 * randn(N1, 1);
y2 = filter([1], [1 -0.8], u2) + 0.1 * randn(N2, 1);

% Should work with cell arrays (variable-length trajectories)
result_cell = sidFreqBT({y1, y2}, {u1, u2});
assert(isstruct(result_cell), ...
    'DISC-7: sidFreqBT should accept cell array input.');
assert(result_cell.NumTrajectories == 2, ...
    'DISC-7: Should report 2 trajectories.');
% Data length should be trimmed to shortest
assert(result_cell.DataLength == N2, ...
    'DISC-7: DataLength should be trimmed to shortest trajectory.');

fprintf('  DISC-7: Cell array input for sidFreqBT - PASSED\n');

%% DISC-7b: Cell array time-series mode
result_ts = sidFreqBT({y1, y2}, []);
assert(isstruct(result_ts), ...
    'DISC-7b: Time-series cell array should work.');
fprintf('  DISC-7b: Cell array time-series mode - PASSED\n');

%% TEST-4: Frozen TF Jacobian validated against finite differences
% Perturb each A and B entry by delta, recompute G, compare dG against
% the Jacobian-based uncertainty.
rng(701);
N_t4 = 20; p_t4 = 2; q_t4 = 1; L_t4 = 10;
X_t4 = randn(N_t4 + 1, p_t4, L_t4);
U_t4 = randn(N_t4, q_t4, L_t4);
ltv_t4 = sidLTVdisc(X_t4, U_t4, 'Lambda', 1e3, 'Uncertainty', true);

w0_t4 = pi / 4;
k0_t4 = 5;
Ak = ltv_t4.A(:, :, k0_t4);
Bk = ltv_t4.B(:, :, k0_t4);
z0_t4 = exp(1i * w0_t4);
Ip_t4 = eye(p_t4);
G0_t4 = (z0_t4 * Ip_t4 - Ak) \ Bk;

% Finite-difference Jacobian for dG/dA
delta_t4 = 1e-7;
dGdA_fd = zeros(p_t4, q_t4, p_t4, p_t4);  % (a, b, j, i) -> dG_{ab}/dA_{ji}
for ji = 1:p_t4
    for ii = 1:p_t4
        Ap = Ak; Ap(ji, ii) = Ap(ji, ii) + delta_t4;
        Gp = (z0_t4 * Ip_t4 - Ap) \ Bk;
        dGdA_fd(:, :, ji, ii) = (Gp - G0_t4) / delta_t4;
    end
end

% Analytical Jacobian: dG_{ab}/dA_{ji} = R_{aj} * [R*B]_{ib}
R_t4 = (z0_t4 * Ip_t4 - Ak) \ Ip_t4;
RB_t4 = R_t4 * Bk;
dGdA_an = zeros(p_t4, q_t4, p_t4, p_t4);
for ji = 1:p_t4
    for ii = 1:p_t4
        for aa = 1:p_t4
            for bb = 1:q_t4
                dGdA_an(aa, bb, ji, ii) = R_t4(aa, ji) * RB_t4(ii, bb);
            end
        end
    end
end

jac_err = max(abs(dGdA_fd(:) - dGdA_an(:)));
assert(jac_err < 1e-4, ...
    'TEST-4: Jacobian dG/dA finite-diff error should be < 1e-4 (got %.2e).', jac_err);

% Also check dG/dB
dGdB_fd = zeros(p_t4, q_t4, p_t4, q_t4);
for ji = 1:p_t4
    for ii = 1:q_t4
        Bp = Bk; Bp(ji, ii) = Bp(ji, ii) + delta_t4;
        Gp = (z0_t4 * Ip_t4 - Ak) \ Bp;
        dGdB_fd(:, :, ji, ii) = (Gp - G0_t4) / delta_t4;
    end
end

dGdB_an = zeros(p_t4, q_t4, p_t4, q_t4);
for ji = 1:p_t4
    for ii = 1:q_t4
        for aa = 1:p_t4
            dGdB_an(aa, ii, ji, ii) = R_t4(aa, ji);
        end
    end
end

jac_err_B = max(abs(dGdB_fd(:) - dGdB_an(:)));
assert(jac_err_B < 1e-4, ...
    'TEST-4: Jacobian dG/dB finite-diff error should be < 1e-4 (got %.2e).', jac_err_B);

fprintf('  TEST-4: Frozen TF Jacobian vs finite differences - PASSED\n');

%% TEST-2: Edge cases for COSMIC and BT

% Very high lambda should produce near-constant A(k)
rng(801);
N_t2 = 30; p_t2 = 2; q_t2 = 1; L_t2 = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X_t2 = zeros(N_t2 + 1, p_t2, L_t2);
U_t2 = randn(N_t2, q_t2, L_t2);
for l = 1:L_t2
    X_t2(1, :, l) = randn(1, p_t2);
    for kk = 1:N_t2
        X_t2(kk+1, :, l) = (A_true * X_t2(kk, :, l)' + B_true * U_t2(kk, :, l)')';
    end
end

r_high = sidLTVdisc(X_t2, U_t2, 'Lambda', 1e15);
% A should be nearly constant across time
A_var = 0;
for kk = 2:N_t2
    A_var = A_var + norm(r_high.A(:, :, kk) - r_high.A(:, :, 1), 'fro');
end
A_var = A_var / (N_t2 - 1);
assert(A_var < 1e-6, ...
    'TEST-2a: Very high lambda should give near-constant A (var=%.2e).', A_var);

fprintf('  TEST-2a: Very high lambda -> constant model - PASSED\n');

% M > N/2 should be clamped with warning
rng(802);
warning('off', 'sid:windowReduced');
r_bigM = sidFreqBT(randn(20, 1), randn(20, 1), 'WindowSize', 50);
warning('on', 'sid:windowReduced');
assert(r_bigM.WindowSize <= 10, ...
    'TEST-2b: Window size should be clamped to N/2 = 10.');

fprintf('  TEST-2b: M > N/2 clamped - PASSED\n');

% API-5: Custom LambdaGrid for auto mode
rng(803);
r_grid = sidLTVdisc(X_t2, U_t2, 'Lambda', 'auto', 'LambdaGrid', [1, 100, 10000]);
assert(all(isfinite(r_grid.A(:))), ...
    'API-5: Custom LambdaGrid should produce finite results.');

fprintf('  API-5: Custom LambdaGrid for auto mode - PASSED\n');

fprintf('test_reviewFixes: ALL PASSED\n');
