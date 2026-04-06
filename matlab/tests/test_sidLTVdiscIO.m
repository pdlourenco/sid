%% test_sidLTVdiscIO - Unit tests for output-COSMIC LTV identification
%
% Tests sidLTVdiscIO for result structure, H=I equivalence with sidLTVdisc,
% partial observation recovery, convergence, multi-trajectory, trust-region,
% state recovery, and input validation.

fprintf('Running test_sidLTVdiscIO...\n');

%% Test 1: Output struct has all required fields
rng(100);
n = 2; q = 1; py = 2; N = 30; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = eye(2);  % full observation (py == n)
R_noise = eye(py);

X = zeros(N+1, n, L);
U = randn(N, q, L);
Y = zeros(N+1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')';
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.01 * randn(1, n);
        Y(k+1, :, l) = (H_obs * X(k+1, :, l)')' + 0.01 * randn(1, py);
    end
end

result = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 1e3);

requiredFields = {'A', 'B', 'X', 'H', 'R', 'Cost', 'Iterations', ...
    'Lambda', 'DataLength', 'StateDim', 'OutputDim', 'InputDim', ...
    'NumTrajectories', 'Algorithm', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
fprintf('  Test 1 passed: all required fields present.\n');

%% Test 2: Correct metadata and dimensions
assert(result.DataLength == N, 'DataLength should be N=%d', N);
assert(result.StateDim == n, 'StateDim should be n=%d', n);
assert(result.OutputDim == py, 'OutputDim should be py=%d', py);
assert(result.InputDim == q, 'InputDim should be q=%d', q);
assert(result.NumTrajectories == L, 'NumTrajectories should be L=%d', L);
assert(strcmp(result.Algorithm, 'cosmic'), 'Algorithm should be cosmic');
assert(strcmp(result.Method, 'sidLTVdiscIO'), 'Method should be sidLTVdiscIO');
assert(isequal(size(result.A), [n, n, N]), 'A should be (n x n x N)');
assert(isequal(size(result.B), [n, q, N]), 'B should be (n x q x N)');
assert(isequal(size(result.X), [N+1, n, L]), 'X should be (N+1 x n x L)');
assert(isequal(size(result.H), [py, n]), 'H should be (py x n)');
assert(isequal(size(result.R), [py, py]), 'R should be (py x py)');
assert(result.Iterations >= 0, 'Should have >= 0 iterations');
fprintf('  Test 2 passed: metadata and dimensions correct.\n');

%% Test 3: H = I equivalence with sidLTVdisc
% When H = I and noise is zero, sidLTVdiscIO should closely match sidLTVdisc.
rng(200);
n = 2; q = 1; N = 50; L = 10;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
sigma = 0.02;

X = zeros(N+1, n, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + sigma * randn(1, n);
    end
end

lam = 1e4;
resultStd = sidLTVdisc(X, U, 'Lambda', lam);
resultIO = sidLTVdiscIO(X, U, eye(n), 'Lambda', lam);

% Compare A and B estimates
errA = norm(mean(resultIO.A, 3) - mean(resultStd.A, 3), 'fro') / ...
       max(norm(mean(resultStd.A, 3), 'fro'), eps);
errB = norm(mean(resultIO.B, 3) - mean(resultStd.B, 3), 'fro') / ...
       max(norm(mean(resultStd.B, 3), 'fro'), eps);
assert(errA < 0.05, 'H=I: A mismatch with sidLTVdisc (errA=%.4f)', errA);
assert(errB < 0.05, 'H=I: B mismatch with sidLTVdisc (errB=%.4f)', errB);
fprintf('  Test 3 passed: H=I matches sidLTVdisc (errA=%.4f, errB=%.4f).\n', errA, errB);

%% Test 4: State estimator with double integrator (position measured)
% Direct test of sidLTVStateEst on a simple system where the answer is
% obvious: a double integrator x = [position; velocity], only position
% is measured.  Deterministic — no rng needed.
n = 2; q = 1; py = 1; N = 50; L = 1;
dt = 1;
A_di = [1 dt; 0 1];       % double integrator
B_di = [0.5*dt^2; dt];    % input = acceleration
H_di = [1 0];             % measure position only

X_true = zeros(N+1, n);
U_di = zeros(N, q);
Y_di = zeros(N+1, py);

% Drive with a known deterministic input (constant acceleration then
% coast then brake) — the state trajectory is easily predictable.
for k = 1:N
    if k <= 15
        U_di(k) = 1.0;      % accelerate
    elseif k <= 35
        U_di(k) = 0.0;      % coast
    else
        U_di(k) = -1.0;     % brake
    end
end

X_true(1, :) = [0, 0];  % start at rest
for k = 1:N
    X_true(k+1, :) = (A_di * X_true(k, :)' + B_di * U_di(k))';
end
Y_di = X_true * H_di';  % noiseless position measurements

% Reshape for sidLTVStateEst (expects 3D arrays)
A_rep = repmat(A_di, [1, 1, N]);
B_rep = repmat(B_di, [1, 1, N]);

X_est = sidLTVStateEst(reshape(Y_di, [], py, 1), ...
    reshape(U_di, [], q, 1), A_rep, B_rep, H_di);

% Position (observed) should match measurements closely
pos_err = norm(X_est(:,1) - X_true(:,1)) / norm(X_true(:,1));
assert(pos_err < 0.01, ...
    'Double integrator position error too large (%.4f)', pos_err);

% Velocity (unobserved) should be recovered via dynamics coupling
vel_err = norm(X_est(:,2) - X_true(:,2)) / max(norm(X_true(:,2)), 1);
assert(vel_err < 0.05, ...
    'Double integrator velocity error too large (%.4f)', vel_err);
fprintf('  Test 4 passed: double integrator (pos=%.4f, vel=%.4f).\n', pos_err, vel_err);

%% Test 5: Partial observation via full pipeline (double integrator)
% Run sidLTVdiscIO end-to-end on the double integrator from Test 4
% (position measured, velocity hidden).  Trust-region is enabled because
% the SPEC (§8.12.4) recommends it for partial observation to smooth the
% transition from the A=I initialisation.
n = 2; q = 1; py = 1; N = 50; L = 5;
dt = 1;
A_di = [1 dt; 0 1];
B_di = [0.5 * dt^2; dt];
H_obs = [1 0];

X = zeros(N+1, n, L);
U = zeros(N, q, L);
Y = zeros(N+1, py, L);

for l = 1:L
    for k = 1:N
        if k <= 15
            U(k, 1, l) = 0.5 * l;
        elseif k <= 35
            U(k, 1, l) = 0;
        else
            U(k, 1, l) = -0.5 * l;
        end
    end
end

for l = 1:L
    X(1, :, l) = [0, 0.1 * l];
    Y(1, :, l) = H_obs * X(1, :, l)';
    for k = 1:N
        X(k+1, :, l) = (A_di * X(k, :, l)' + B_di * U(k, :, l)')';
        Y(k+1, :, l) = H_obs * X(k+1, :, l)';
    end
end

result = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 100, 'TrustRegion', 1);

assert(~any(isnan(result.A(:))), 'Partial obs pipeline produced NaN in A');
assert(~any(isnan(result.X(:))), 'Partial obs pipeline produced NaN in X');
assert(result.Iterations >= 1, 'Partial obs pipeline did not iterate');

% Verify state estimates reproduce measurements (check pipeline coherence)
for l = 1:L
    obs = squeeze(result.X(:, :, l)) * H_obs';
    obs_err = norm(obs - squeeze(Y(:, :, l)), 'fro') / ...
        norm(squeeze(Y(:, :, l)), 'fro');
    assert(obs_err < 0.5, ...
        'Partial obs pipeline: obs mismatch traj %d (%.3f)', l, obs_err);
end
fprintf('  Test 5 passed: partial obs pipeline (%d iters, no NaN).\n', ...
    result.Iterations);

%% Test 6: Monotone cost decrease (partial obs, no trust-region)
% Monotonicity is guaranteed for the alternating minimisation when
% trust-region is off. Use partial obs to exercise the EM loop.
rng(200);
n = 2; q = 1; N = 50; L = 10;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_mono = [1 0];
sigma = 0.02;
X_mono = zeros(N+1, n, L);
U_mono = randn(N, q, L);
Y_mono = zeros(N+1, 1, L);
for l = 1:L
    X_mono(1, :, l) = randn(1, n);
    Y_mono(1, :, l) = H_mono * X_mono(1, :, l)';
    for k = 1:N
        X_mono(k+1, :, l) = (A_true * X_mono(k, :, l)' + ...
            B_true * U_mono(k, :, l)')' + ...
            sigma * randn(1, n);
        Y_mono(k+1, :, l) = H_mono * X_mono(k+1, :, l)' ...
            + sigma * randn(1, 1);
    end
end
res_mono = sidLTVdiscIO( ...
    Y_mono, U_mono, H_mono, 'Lambda', 1e4);
assert(length(res_mono.Cost) >= 2, ...
    'Need at least 2 cost evaluations');
for i = 2:length(res_mono.Cost)
    assert(res_mono.Cost(i) <= res_mono.Cost(i-1) + ...
        1e-8 * abs(res_mono.Cost(i-1)), ...
        'Cost increased at iteration %d: %.6f > %.6f', ...
        i, res_mono.Cost(i), res_mono.Cost(i-1));
end
fprintf('  Test 6 passed: monotone cost decrease (%d iters).\n', ...
    res_mono.Iterations);

%% Test 7: State recovery
% Compare estimated states to true states (for observed dimensions).
% Deterministic data to avoid RNG-dependent singular blocks.
n = 2; q = 1; py = 1; N = 40; L = 10;
A_true = [0.85 0.1; -0.1 0.85];
B_true = [1; 0.3];
H_obs = [1 0];  % observe only first state

X = zeros(N+1, n, L);
U = zeros(N, q, L);
Y = zeros(N+1, py, L);

for l = 1:L
    freq = l / (4 * N);
    for k = 1:N
        U(k, 1, l) = sin(2 * pi * freq * k) + 0.5 * (-1)^(k + l);
    end
end

for l = 1:L
    X(1, :, l) = [0.5 * cos(l), 0.5 * sin(l)];
    Y(1, :, l) = H_obs * X(1, :, l)';
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k+1, :, l) = H_obs * X(k+1, :, l)';
    end
end

result = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 100, 'TrustRegion', 1);

% Check that observed states match measurements reasonably
for l = 1:min(L, 3)
    obs_states = squeeze(result.X(:, :, l)) * H_obs';  % (N+1 x py)
    obs_err = norm(obs_states - squeeze(Y(:, :, l)), 'fro') / norm(squeeze(Y(:, :, l)), 'fro');
    assert(obs_err < 0.5, ...
        'State recovery: too far from measurements (traj %d, err=%.3f)', ...
        l, obs_err);
end
fprintf('  Test 7 passed: state recovery consistent with measurements.\n');

%% Test 8: Multi-trajectory improves estimates
% Compare single-trajectory vs multi-trajectory recovery.
rng(500);
n = 2; q = 1; py = 2; N = 30; L = 10;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = eye(2);
sigma = 0.05;

X = zeros(N+1, n, L);
U = randn(N, q, L);
Y = zeros(N+1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')' + sigma * randn(1, py);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + sigma * randn(1, n);
        Y(k+1, :, l) = (H_obs * X(k+1, :, l)')' + sigma * randn(1, py);
    end
end

% Single trajectory
res1 = sidLTVdiscIO(Y(:,:,1), U(:,:,1), H_obs, 'Lambda', 1e3);
err1 = norm(mean(res1.A, 3) - A_true, 'fro');

% Multi-trajectory
resL = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 1e3);
errL = norm(mean(resL.A, 3) - A_true, 'fro');

% Multi should be at least as good (or close)
assert(errL <= err1 * 1.5, ...
    'Multi-trajectory should improve: errL=%.4f > 1.5*err1=%.4f', errL, 1.5*err1);
fprintf('  Test 8 passed: multi-trajectory (errL=%.4f vs err1=%.4f).\n', errL, err1);

%% Test 9: R weighting
% With known R, estimates should improve over R = I when noise is anisotropic.
rng(600);
n = 2; q = 1; py = 2; N = 40; L = 10;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = eye(2);
R_true = diag([0.001, 1.0]);  % channel 1 precise, channel 2 noisy

X = zeros(N+1, n, L);
U = randn(N, q, L);
Y = zeros(N+1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = X(1, :, l) + (chol(R_true, 'lower') * randn(py, 1))';
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.01 * randn(1, n);
        Y(k+1, :, l) = X(k+1, :, l) + (chol(R_true, 'lower') * randn(py, 1))';
    end
end

resI = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 1e3);
resR = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 1e3, 'R', R_true);

errI = norm(mean(resI.A, 3) - A_true, 'fro');
errR = norm(mean(resR.A, 3) - A_true, 'fro');

% With correct R, estimates should be at least as good
assert(errR <= errI * 1.5, ...
    'R weighting should help: errR=%.4f > 1.5*errI=%.4f', errR, 1.5*errI);
fprintf('  Test 9 passed: R weighting (errR=%.4f vs errI=%.4f).\n', errR, errI);

%% Test 10: Input validation - mismatched H dimensions
passed = false;
try
    H_bad = eye(3);  % 3x3 but py=2
    sidLTVdiscIO(Y(:,:,1), U(:,:,1), H_bad, 'Lambda', 1e3);
catch e
    if ~isempty(strfind(e.identifier, 'sid:'))
        passed = true;
    end
end
assert(passed, 'Should error on H dimension mismatch.');
fprintf('  Test 10 passed: input validation rejects bad H.\n');

%% Test 11: Trust-region convergence
% Deterministic data to avoid RNG-dependent singular blocks.
n = 2; q = 1; py = 1; N = 30; L = 8;
A_true = [0.9 0.2; -0.2 0.85];
B_true = [1; 0.5];
H_obs = [1 0];

X = zeros(N+1, n, L);
U = zeros(N, q, L);
Y = zeros(N+1, py, L);

for l = 1:L
    freq = l / (3 * N);
    for k = 1:N
        U(k, 1, l) = sin(2 * pi * freq * k) + (mod(k, 5) < 2);
    end
end

for l = 1:L
    X(1, :, l) = [0.4 * l / L, -0.3 * l / L];
    Y(1, :, l) = H_obs * X(1, :, l)';
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k+1, :, l) = H_obs * X(k+1, :, l)';
    end
end

result_tr = sidLTVdiscIO(Y, U, H_obs, 'Lambda', 100, 'TrustRegion', 1);

assert(isfield(result_tr, 'A'), 'Trust-region should return valid result');
assert(result_tr.Iterations >= 1, 'Trust-region should iterate');
fprintf('  Test 11 passed: trust-region converges (%d iterations).\n', result_tr.Iterations);

%% Test 12: Mass-spring-damper LTI, full observation
% Full pipeline on 6-state MSD with all states measured.
rng(1200);
[Ad, Bd] = sidTestMSD( ...
    [1.0; 1.5; 1.0], [100; 80; 60], [2; 1.5; 1], ...
    [1; 0; 0], 0.01);
n = 6; q = 1; py = 6; N = 200; L = 10;
H_full = eye(n);

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
sigma = 0.001;
for l = 1:L
    X(1, :, l) = 0.1 * randn(1, n);
    Y(1, :, l) = (H_full * X(1, :, l)')' ...
        + sigma * randn(1, py);
    for k = 1:N
        X(k + 1, :, l) = (Ad * X(k, :, l)' ...
            + Bd * U(k, :, l)')' + sigma * randn(1, n);
        Y(k + 1, :, l) = (H_full * X(k + 1, :, l)')' ...
            + sigma * randn(1, py);
    end
end

result = sidLTVdiscIO(Y, U, H_full, 'Lambda', 1e6);

assert(~any(isnan(result.A(:))), 'MSD full obs: NaN in A');
assert(result.Iterations >= 0, 'MSD full obs: unexpected iterations');

% Check eigenvalue recovery of the mean A (similarity-invariant)
eig_true = sort(abs(eig(Ad)));
eig_est = sort(abs(eig(mean(result.A, 3))));
eig_err = norm(eig_true - eig_est) / norm(eig_true);
assert(eig_err < 1.0, ...
    'MSD full obs: eigenvalue error too large (%.4f)', eig_err);
fprintf('  Test 12 passed: MSD full obs (eig_err=%.4f, %d iters).\n', ...
    eig_err, result.Iterations);

%% Test 13: Mass-spring-damper LTI, partial obs (positions only)
% py=3, n=6. Use trust-region for this harder problem.
H_pos = [eye(3), zeros(3, 3)];
Y_pos = zeros(N + 1, 3, L);
for l = 1:L
    for k = 1:N + 1
        Y_pos(k, :, l) = (H_pos * X(k, :, l)')' ...
            + sigma * randn(1, 3);
    end
end

result = sidLTVdiscIO( ...
    Y_pos, U, H_pos, 'Lambda', 1e4, 'TrustRegion', 1);

assert(~any(isnan(result.A(:))), 'MSD partial: NaN in A');
assert(~any(isnan(result.X(:))), 'MSD partial: NaN in X');
assert(result.Iterations >= 1, 'MSD partial: no iterations');

% Observed states should be consistent with measurements
for l = 1:min(L, 3)
    obs = squeeze(result.X(:, :, l)) * H_pos';
    obs_err = norm(obs - squeeze(Y_pos(:, :, l)), 'fro') / ...
        norm(squeeze(Y_pos(:, :, l)), 'fro');
    assert(obs_err < 1.0, ...
        'MSD partial: obs mismatch traj %d (%.3f)', l, obs_err);
end
fprintf('  Test 13 passed: MSD partial obs (%d iters).\n', ...
    result.Iterations);

%% Test 14: Time-varying double integrator, partial obs
% Gain varies linearly: alpha(k) = 1 + 0.5*k/N.
n = 2; q = 1; py = 1; N = 80; L = 8;
dt = 1;
A_tv = zeros(n, n, N);
B_tv = zeros(n, q, N);
for k = 1:N
    alpha = 1.0 + 0.5 * (k - 1) / (N - 1);
    A_tv(:, :, k) = [1, alpha * dt; 0, 1];
    B_tv(:, :, k) = [0.5 * alpha * dt^2; alpha * dt];
end

X = zeros(N + 1, n, L);
U = zeros(N, q, L);
Y = zeros(N + 1, py, L);
H_pos = [1 0];
for l = 1:L
    freq = l / (3 * N);
    for k = 1:N
        U(k, 1, l) = sin(2 * pi * freq * k) + 0.5 * l;
    end
    X(1, :, l) = [0, 0.1 * l / L];
    Y(1, :, l) = H_pos * X(1, :, l)';
    for k = 1:N
        X(k + 1, :, l) = (A_tv(:, :, k) * X(k, :, l)' ...
            + B_tv(:, :, k) * U(k, :, l)')';
        Y(k + 1, :, l) = H_pos * X(k + 1, :, l)';
    end
end

result = sidLTVdiscIO( ...
    Y, U, H_pos, 'Lambda', 100, 'TrustRegion', 1);

assert(~any(isnan(result.A(:))), 'TV DI: NaN in A');
assert(result.Iterations >= 1, 'TV DI: no iterations');

% Verify state estimates reproduce measurements
for l = 1:min(L, 3)
    obs = squeeze(result.X(:, :, l)) * H_pos';
    obs_err = norm(obs - squeeze(Y(:, :, l)), 'fro') / ...
        norm(squeeze(Y(:, :, l)), 'fro');
    assert(obs_err < 0.5, ...
        'TV DI: obs mismatch traj %d (%.3f)', l, obs_err);
end
fprintf('  Test 14 passed: TV DI partial obs (%d iters).\n', ...
    result.Iterations);

%% Test 15: Time-varying mass-spring-damper, partial obs
% Stiffness k1 varies sinusoidally. Check pipeline convergence.
rng(1500);
n = 6; q = 1; py = 3; N = 100; L = 5;
H_pos = [eye(3), zeros(3, 3)];
Ts = 0.01;
m_vec = [1.0; 1.5; 1.0];
c_vec = [2; 1.5; 1];
F_vec = [1; 0; 0];

A_tv = zeros(n, n, N);
B_tv = zeros(n, q, N);
for k = 1:N
    k1 = 100 * (1 + 0.3 * sin(2 * pi * k / N));
    [Adk, Bdk] = sidTestMSD( ...
        m_vec, [k1; 80; 60], c_vec, F_vec, Ts);
    A_tv(:, :, k) = Adk;
    B_tv(:, :, k) = Bdk;
end

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = 0.05 * randn(1, n);
    Y(1, :, l) = (H_pos * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_tv(:, :, k) * X(k, :, l)' ...
            + B_tv(:, :, k) * U(k, :, l)')';
        Y(k + 1, :, l) = (H_pos * X(k + 1, :, l)')';
    end
end

result = sidLTVdiscIO( ...
    Y, U, H_pos, 'Lambda', 1e3, 'TrustRegion', 1);

assert(~any(isnan(result.A(:))), 'TV MSD: NaN in A');
assert(~any(isnan(result.X(:))), 'TV MSD: NaN in X');
assert(result.Iterations >= 1, 'TV MSD: no iterations');
fprintf('  Test 15 passed: TV MSD partial obs (%d iters).\n', ...
    result.Iterations);

%% Test 16: Full-rank fast path — H=I matches sidLTVdisc exactly
% When H=I, the fast path recovers X = Y exactly and runs a single
% COSMIC pass. The result should match sidLTVdisc on the same data.
rng(1600);
n = 2; q = 1; N = 80; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];

X = zeros(N + 1, n, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' ...
            + B_true * U(k, :, l)')';
    end
end

lam = 1e4;
res_direct = sidLTVdisc(X, U, 'Lambda', lam);
res_io = sidLTVdiscIO(X, U, eye(n), 'Lambda', lam);

assert(res_io.Iterations == 0, ...
    'Fast path should have 0 iterations, got %d', ...
    res_io.Iterations);

errA = norm(mean(res_io.A, 3) - mean(res_direct.A, 3), 'fro') / ...
    norm(mean(res_direct.A, 3), 'fro');
errB = norm(mean(res_io.B, 3) - mean(res_direct.B, 3), 'fro') / ...
    norm(mean(res_direct.B, 3), 'fro');
assert(errA < 1e-10, ...
    'Fast path A should match sidLTVdisc exactly (err=%.2e)', errA);
assert(errB < 1e-10, ...
    'Fast path B should match sidLTVdisc exactly (err=%.2e)', errB);

% X should equal Y when H=I
errX = norm(res_io.X(:) - X(:)) / norm(X(:));
assert(errX < 1e-12, ...
    'Fast path X should equal Y for H=I (err=%.2e)', errX);
fprintf('  Test 16 passed: fast path H=I matches sidLTVdisc (errA=%.2e, errB=%.2e).\n', ...
    errA, errB);

%% Test 17: Full-rank fast path — non-identity square H
% H is a rotation matrix (full rank, py = n = 2).
rng(1700);
n = 2; q = 1; N = 60; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
theta = pi / 6;
H_rot = [cos(theta) -sin(theta)
         sin(theta)  cos(theta)];

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, n, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_rot * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' ...
            + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_rot * X(k + 1, :, l)')';
    end
end

lam = 1e4;
res_rot = sidLTVdiscIO(Y, U, H_rot, 'Lambda', lam);

assert(res_rot.Iterations == 0, ...
    'Rotated H fast path should have 0 iters, got %d', ...
    res_rot.Iterations);

% Recovered states should match truth (noiseless)
errX = norm(res_rot.X(:) - X(:)) / norm(X(:));
assert(errX < 1e-10, ...
    'Rotated H: state recovery error %.2e', errX);

% A should be close to true
A_mean = mean(res_rot.A, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
assert(errA < 0.01, ...
    'Rotated H: A error %.4f', errA);
fprintf('  Test 17 passed: rotated H fast path (errA=%.4f, errX=%.2e).\n', ...
    errA, errX);

%% Test 18: Full-rank fast path — tall H (py > n)
% H is 4x2 with full column rank. The fast path should trigger and
% recover states exactly via weighted least squares.
rng(1800);
n = 2; q = 1; N = 60; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_tall = [1 0; 0 1; 1 1; 1 -1];
py_tall = 4;

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py_tall, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_tall * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' ...
            + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_tall * X(k + 1, :, l)')';
    end
end

lam = 1e4;
res_tall = sidLTVdiscIO(Y, U, H_tall, 'Lambda', lam);

assert(res_tall.Iterations == 0, ...
    'Tall H fast path should have 0 iters, got %d', ...
    res_tall.Iterations);

% Recovered states should match truth (noiseless)
errX = norm(res_tall.X(:) - X(:)) / norm(X(:));
assert(errX < 1e-10, ...
    'Tall H: state recovery error %.2e', errX);

% A should be close to true
A_mean = mean(res_tall.A, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
assert(errA < 0.01, ...
    'Tall H: A error %.4f', errA);
fprintf('  Test 18 passed: tall H fast path (errA=%.4f, errX=%.2e).\n', ...
    errA, errX);

%% Test 19: Cell input matches 3D for equal-length trajectories
rng(1900);
n = 2; q = 1; N = 40; L = 3;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_19 = eye(n);

X19 = zeros(N + 1, n, L);
U19 = randn(N, q, L);
Y19 = zeros(N + 1, n, L);
for l = 1:L
    X19(1, :, l) = randn(1, n);
    Y19(1, :, l) = X19(1, :, l);
    for k = 1:N
        X19(k + 1, :, l) = (A_true * X19(k, :, l)' ...
            + B_true * U19(k, :, l)')';
        Y19(k + 1, :, l) = X19(k + 1, :, l);
    end
end

res_3d = sidLTVdiscIO(Y19, U19, H_19, 'Lambda', 1e4);

% Build cell arrays with same data
Y_cell = cell(L, 1);
U_cell = cell(L, 1);
for l = 1:L
    Y_cell{l} = Y19(:, :, l);
    U_cell{l} = U19(:, :, l);
end
res_cell = sidLTVdiscIO(Y_cell, U_cell, H_19, 'Lambda', 1e4);

% Results should be identical
errA = norm(res_cell.A(:) - res_3d.A(:));
errB = norm(res_cell.B(:) - res_3d.B(:));
assert(errA < 1e-8, 'Cell vs 3D: A mismatch %.2e', errA);
assert(errB < 1e-8, 'Cell vs 3D: B mismatch %.2e', errB);
assert(res_cell.Iterations == res_3d.Iterations, ...
    'Cell vs 3D: iteration count mismatch');

% X should be cell output matching 3D content
assert(iscell(res_cell.X), 'Cell input should produce cell X output');
for l = 1:L
    errXl = norm(res_cell.X{l} - res_3d.X(:, :, l));
    assert(errXl < 1e-8, 'Cell vs 3D: X{%d} mismatch %.2e', l, errXl);
end
fprintf('  Test 19 passed: cell input matches 3D (errA=%.2e).\n', errA);

%% Test 20: Variable-length trajectories — fast path (rank(H) = n)
rng(2000);
n = 2; q = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_20 = eye(n);
horizons_20 = [60; 40; 50; 30];
L = length(horizons_20);
N = max(horizons_20);

Y_cell = cell(L, 1);
U_cell = cell(L, 1);
X_cell = cell(L, 1);
for l = 1:L
    Nl = horizons_20(l);
    U_cell{l} = randn(Nl, q);
    X_cell{l} = zeros(Nl + 1, n);
    X_cell{l}(1, :) = randn(1, n);
    for k = 1:Nl
        X_cell{l}(k + 1, :) = (A_true * X_cell{l}(k, :)' ...
            + B_true * U_cell{l}(k, :)')';
    end
    Y_cell{l} = X_cell{l};  % H = I, noiseless
end

res_vl = sidLTVdiscIO(Y_cell, U_cell, H_20, 'Lambda', 1e4);

assert(res_vl.Iterations == 0, ...
    'VarLen fast path should have 0 iters, got %d', res_vl.Iterations);
assert(iscell(res_vl.X), 'VarLen should return cell X');
assert(isfield(res_vl, 'Horizons'), 'VarLen should have Horizons field');

% State recovery should be exact for each trajectory
for l = 1:L
    errXl = norm(res_vl.X{l}(:) - X_cell{l}(:)) ...
        / norm(X_cell{l}(:));
    assert(errXl < 1e-10, ...
        'VarLen fast path: X{%d} error %.2e', l, errXl);
end

A_mean = mean(res_vl.A, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
assert(errA < 0.01, 'VarLen fast path: A error %.4f', errA);
fprintf('  Test 20 passed: varlen fast path (errA=%.4f).\n', errA);

%% Test 21: Variable-length trajectories — EM path (rank(H) < n)
rng(2100);
n = 2; q = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_21 = [1 0];  % partial obs: py=1 < n=2
horizons_21 = [60; 40; 50; 30];
L = length(horizons_21);
N = max(horizons_21);

Y_cell = cell(L, 1);
U_cell = cell(L, 1);
X_cell = cell(L, 1);
for l = 1:L
    Nl = horizons_21(l);
    U_cell{l} = randn(Nl, q);
    X_cell{l} = zeros(Nl + 1, n);
    X_cell{l}(1, :) = randn(1, n);
    for k = 1:Nl
        X_cell{l}(k + 1, :) = (A_true * X_cell{l}(k, :)' ...
            + B_true * U_cell{l}(k, :)')';
    end
    Y_cell{l} = X_cell{l} * H_21';
end

res_vl = sidLTVdiscIO(Y_cell, U_cell, H_21, 'Lambda', 1e4);

assert(res_vl.Iterations > 0, 'EM path should iterate');
assert(iscell(res_vl.X), 'VarLen should return cell X');

% Cost should generally decrease (small numerical increases tolerated)
costs = res_vl.Cost;
assert(costs(end) < costs(1), ...
    'Final cost %.4f should be less than initial %.4f', ...
    costs(end), costs(1));

% Measurements should be consistent with estimated states
for l = 1:L
    Nl = horizons_21(l);
    Y_recon = res_vl.X{l} * H_21';
    errY = norm(Y_recon - Y_cell{l}) / norm(Y_cell{l});
    assert(errY < 0.1, ...
        'VarLen EM: Y reconstruction error %.4f for traj %d', errY, l);
end
fprintf('  Test 21 passed: varlen EM path (%d iters, %d traj).\n', ...
    res_vl.Iterations, L);

%% Test 22: Trajectory trimming in sidLTIfreqIO
rng(2200);
n = 2; q = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_22 = eye(n);
horizons_22 = [100; 90; 80; 20; 95];
L = length(horizons_22);

Y_cell = cell(L, 1);
U_cell = cell(L, 1);
for l = 1:L
    Nl = horizons_22(l);
    U_cell{l} = randn(Nl, q);
    x = zeros(Nl + 1, n);
    x(1, :) = randn(1, n);
    for k = 1:Nl
        x(k + 1, :) = (A_true * x(k, :)' + B_true * U_cell{l}(k, :)')';
    end
    Y_cell{l} = x * H_22';
end

% Traj 4 (horizon=20) is < 2/3*100=67, should be discarded
[A0, B0] = sidLTIfreqIO(Y_cell, U_cell, H_22);

assert(isequal(size(A0), [n, n]), 'A0 should be %dx%d', n, n);
assert(isequal(size(B0), [n, q]), 'B0 should be %dx%d', n, q);

% Check that the LTI estimate is reasonable
eig_err = max(abs(sort(abs(eig(A0))) - sort(abs(eig(A_true)))));
assert(eig_err < 0.3, ...
    'LTI trimmed: eigenvalue error %.4f', eig_err);
fprintf('  Test 22 passed: LTI trimming (eig_err=%.4f).\n', eig_err);

%% Test 23: Convergence criterion — early stop when |dJ/J| < tol (SPEC §8.12.3)
rng(2300);
n = 2; q = 1; N = 40; L = 5; py = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_23 = [1 0];

X23 = zeros(N + 1, n, L);
U23 = randn(N, q, L);
Y23 = zeros(N + 1, py, L);
for l = 1:L
    X23(1, :, l) = randn(1, n);
    Y23(1, :, l) = H_23 * X23(1, :, l)';
    for k = 1:N
        X23(k+1, :, l) = (A_true * X23(k, :, l)' ...
            + B_true * U23(k, :, l)')';
        Y23(k+1, :, l) = H_23 * X23(k+1, :, l)';
    end
end

% Use generous MaxIter but moderate Tolerance — should converge early
tol_val = 1e-3;
res_conv = sidLTVdiscIO(Y23, U23, H_23, 'Lambda', 1e5, ...
    'MaxIter', 200, 'Tolerance', tol_val);

assert(res_conv.Iterations < 200, ...
    'Should converge before 200 iters, got %d', res_conv.Iterations);
assert(res_conv.Iterations >= 2, ...
    'Should take at least 2 iters, got %d', res_conv.Iterations);

% Verify the stopping criterion: final relative change < tol
costs = res_conv.Cost;
if length(costs) >= 2
    final_rel = abs(costs(end) - costs(end-1)) / max(abs(costs(end-1)), 1);
    assert(final_rel < tol_val, ...
        'Final relative change %.2e should be < tol %.2e', ...
        final_rel, tol_val);
end
fprintf('  Test 23 passed: convergence criterion (%d iters).\n', ...
    res_conv.Iterations);

%% Test 24: MaxIter limit — stops at exactly MaxIter if not converged
res_lim = sidLTVdiscIO(Y23, U23, H_23, 'Lambda', 1e4, ...
    'MaxIter', 3, 'Tolerance', 1e-15);

assert(res_lim.Iterations == 3, ...
    'Should stop at MaxIter=3, got %d', res_lim.Iterations);
assert(length(res_lim.Cost) == 3, ...
    'Cost history should have 3 entries, got %d', length(res_lim.Cost));
fprintf('  Test 24 passed: MaxIter limit (3 iters).\n');

%% Test 25: Rank-deficient square H forces EM path (SPEC §8.12.3)
rng(2500);
n = 3; q = 1; N = 40; L = 3;
A_true = [0.9 0.1 0; -0.1 0.8 0.05; 0 -0.05 0.7];
B_true = [0.5; 0.3; 0.2];
% H is 3x3 but rank 2 — should NOT trigger fast path
H_rd = [1 0 0; 0 1 0; 1 1 0];  % rank 2

X25 = zeros(N + 1, n, L);
U25 = randn(N, q, L);
Y25 = zeros(N + 1, 3, L);
for l = 1:L
    X25(1, :, l) = randn(1, n);
    Y25(1, :, l) = (H_rd * X25(1, :, l)')';
    for k = 1:N
        X25(k+1, :, l) = (A_true * X25(k, :, l)' ...
            + B_true * U25(k, :, l)')';
        Y25(k+1, :, l) = (H_rd * X25(k+1, :, l)')';
    end
end

res_rd = sidLTVdiscIO(Y25, U25, H_rd, 'Lambda', 1e4, 'MaxIter', 10);
assert(res_rd.Iterations > 0, ...
    'Rank-deficient square H should use EM, got 0 iters');
fprintf('  Test 25 passed: rank-deficient square H forces EM (%d iters).\n', ...
    res_rd.Iterations);

%% Test 26: CovarianceMode option — diagonal (default)
rng(126);
n26 = 2; q26 = 1; py26 = 2; N26 = 20; L26 = 5;
A26 = [0.9 0.1; -0.1 0.8]; B26 = [0.5; 0.3]; H26 = eye(2);
X26 = zeros(N26+1, n26, L26); U26 = randn(N26, q26, L26);
Y26 = zeros(N26+1, py26, L26);
for l = 1:L26
    X26(1, :, l) = randn(1, n26);
    Y26(1, :, l) = (H26 * X26(1, :, l)')';
    for k = 1:N26
        X26(k+1, :, l) = (A26 * X26(k, :, l)' + B26 * U26(k, :, l)')' ...
            + 0.02 * randn(1, n26);
        Y26(k+1, :, l) = (H26 * X26(k+1, :, l)')';
    end
end

res26d = sidLTVdiscIO(Y26, U26, H26, 'Lambda', 1e3, ...
    'CovarianceMode', 'diagonal');
Sig26d = res26d.NoiseCov;
offDiag26d = Sig26d - diag(diag(Sig26d));
assert(max(abs(offDiag26d(:))) < 1e-15, ...
    'Diagonal mode should produce diagonal NoiseCov');
fprintf('  Test 26 passed: CovarianceMode diagonal.\n');

%% Test 27: CovarianceMode option — full
res26f = sidLTVdiscIO(Y26, U26, H26, 'Lambda', 1e3, ...
    'CovarianceMode', 'full');
Sig26f = res26f.NoiseCov;
assert(isequal(size(Sig26f), [n26, n26]), 'Full mode NoiseCov dims');
% Full mode may have off-diagonals (no assertion on zeros)
assert(all(isfinite(Sig26f(:))), 'Full mode NoiseCov should be finite');
fprintf('  Test 27 passed: CovarianceMode full.\n');

%% Test 28: CovarianceMode option — isotropic
res26i = sidLTVdiscIO(Y26, U26, H26, 'Lambda', 1e3, ...
    'CovarianceMode', 'isotropic');
Sig26i = res26i.NoiseCov;
assert(abs(Sig26i(1,1) - Sig26i(2,2)) < 1e-15, ...
    'Isotropic mode should have equal diagonal entries');
offDiag26i = Sig26i - diag(diag(Sig26i));
assert(max(abs(offDiag26i(:))) < 1e-15, ...
    'Isotropic mode should be scalar * I');
fprintf('  Test 28 passed: CovarianceMode isotropic.\n');

%% Test 29: Invalid CovarianceMode errors
threw = false;
try
    sidLTVdiscIO(Y26, U26, H26, 'Lambda', 1e3, 'CovarianceMode', 'bad');
catch e
    threw = true;
    assert(~isempty(strfind(e.identifier, 'badCovMode')), ...
        'Expected sid:badCovMode error, got %s', e.identifier);
end
assert(threw, 'Invalid CovarianceMode should throw an error');
fprintf('  Test 29 passed: invalid CovarianceMode throws error.\n');

fprintf('test_sidLTVdiscIO: all tests passed.\n');
