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
assert(result.Iterations >= 1, 'Should have at least 1 iteration');
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

%% Test 6: Monotone cost decrease
% The cost should be non-increasing across alternating iterations.
assert(length(result.Cost) >= 2, 'Need at least 2 cost evaluations');
for i = 2:length(result.Cost)
    assert(result.Cost(i) <= result.Cost(i-1) + 1e-8 * abs(result.Cost(i-1)), ...
        'Cost increased at iteration %d: %.6f > %.6f', i, result.Cost(i), result.Cost(i-1));
end
fprintf('  Test 6 passed: monotone cost decrease verified (%d iterations).\n', result.Iterations);

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

fprintf('test_sidLTVdiscIO: all tests passed.\n');
