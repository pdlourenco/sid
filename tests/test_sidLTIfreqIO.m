%% test_sidLTIfreqIO - Unit tests for LTI frequency-domain I/O identification
%
% Tests sidLTIfreqIO for output dimensions, known system recovery,
% stability enforcement, H-basis consistency, multi-trajectory support,
% and input validation.

fprintf('Running test_sidLTIfreqIO...\n');

%% Test 1: Output dimensions
rng(100);
n = 3; q = 2; py = 2; N = 200; L = 5;
A_true = [0.8 0.1 0; -0.1 0.9 0.05; 0 -0.05 0.85];
B_true = [0.5 0; 0 0.3; 0.1 0.2];
H_obs = [1 0 0; 0 1 0];

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = 0.1 * randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_obs * X(k + 1, :, l)')';
    end
end

[A0, B0] = sidLTIfreqIO(Y, U, H_obs);

assert(isequal(size(A0), [n, n]), 'A0 should be (%d x %d)', n, n);
assert(isequal(size(B0), [n, q]), 'B0 should be (%d x %d)', n, q);
assert(isreal(A0), 'A0 should be real');
assert(isreal(B0), 'B0 should be real');
fprintf('  Test 1 passed: output dimensions correct.\n');

%% Test 2: Known system recovery (noiseless, full obs)
rng(200);
n = 2; q = 1; py = 2; N = 500; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = eye(2);

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_obs * X(k + 1, :, l)')';
    end
end

[A0, B0] = sidLTIfreqIO(Y, U, H_obs);

% With noiseless data and full observation, should recover A_true closely.
% Similarity transform may differ, so compare eigenvalues.
eig_true = sort(eig(A_true));
eig_est  = sort(eig(A0));
eig_err = norm(eig_true - eig_est) / norm(eig_true);
assert(eig_err < 0.1, ...
    'Known system: eigenvalue error too large (%.4f)', eig_err);

% With H=I, basis should match, so compare A directly
A_err = norm(A0 - A_true, 'fro') / norm(A_true, 'fro');
assert(A_err < 0.15, ...
    'Known system (H=I): A error too large (%.4f)', A_err);
fprintf('  Test 2 passed: known system recovery (eig_err=%.4f, A_err=%.4f).\n', ...
    eig_err, A_err);

%% Test 3: Stability enforcement
rng(300);
n = 2; q = 1; py = 2; N = 200; L = 3;
% Use a marginally stable system (eigenvalue near 1)
A_true = [0.99 0.05; -0.05 0.98];
B_true = [0.5; 0.3];
H_obs = eye(2);

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = 0.1 * randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_obs * X(k + 1, :, l)')';
    end
end

[A0, ~] = sidLTIfreqIO(Y, U, H_obs, 'MaxStability', 0.95);

max_eig = max(abs(eig(A0)));
assert(max_eig <= 0.95 + 1e-10, ...
    'Stability: max eigenvalue %.4f exceeds 0.95', max_eig);
fprintf('  Test 3 passed: stability enforced (max |eig|=%.4f).\n', max_eig);

%% Test 4: H-basis consistency (impulse response match)
% The realization should satisfy H * A0^k * B0 ≈ g(k+1)
rng(400);
n = 2; q = 1; py = 1; N = 300; L = 5;
A_true = [0.85 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = [1 0];  % observe only first state

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_obs * X(k + 1, :, l)')';
    end
end

[A0, B0] = sidLTIfreqIO(Y, U, H_obs);

% Compare Markov parameters: g(k) = H * A0^{k-1} * B0 vs H * A_true^{k-1} * B_true
% Use normalized error over the first 10 lags (high lags have tiny values)
g_est_vec = zeros(10, py);
g_true_vec = zeros(10, py);
for k = 1:10
    g_est_vec(k, :)  = (H_obs * A0^(k - 1) * B0)';
    g_true_vec(k, :) = (H_obs * A_true^(k - 1) * B_true)';
end
max_err = norm(g_est_vec - g_true_vec, 'fro') / norm(g_true_vec, 'fro');
assert(max_err < 0.15, ...
    'H-basis: Markov parameter mismatch (norm err=%.4f)', max_err);
fprintf('  Test 4 passed: H-basis impulse response match (err=%.4f).\n', ...
    max_err);

%% Test 5: Multi-trajectory vs single trajectory
rng(500);
n = 2; q = 1; py = 2; N = 200; L = 10;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = eye(2);
sigma = 0.05;

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')' + sigma * randn(1, py);
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' ...
            + 0.01 * randn(1, n);
        Y(k + 1, :, l) = (H_obs * X(k + 1, :, l)')' + sigma * randn(1, py);
    end
end

% Single trajectory
[A1, ~] = sidLTIfreqIO(Y(:, :, 1), U(:, :, 1), H_obs);
err1 = norm(sort(eig(A1)) - sort(eig(A_true)));

% Multi-trajectory
[AL, ~] = sidLTIfreqIO(Y, U, H_obs);
errL = norm(sort(eig(AL)) - sort(eig(A_true)));

% Multi should be at least as good (or close)
assert(errL <= err1 * 2.0, ...
    'Multi-traj should help: errL=%.4f > 2*err1=%.4f', errL, 2 * err1);
fprintf('  Test 5 passed: multi-trajectory (errL=%.4f vs err1=%.4f).\n', ...
    errL, err1);

%% Test 6: Input validation - mismatched H
passed = false;
try
    H_bad = eye(5);  % 5x5 but py=2
    sidLTIfreqIO(Y(:, :, 1), U(:, :, 1), H_bad);
catch e
    if ~isempty(strfind(e.identifier, 'sid:'))
        passed = true;
    end
end
assert(passed, 'Should error on H dimension mismatch.');
fprintf('  Test 6 passed: input validation rejects bad H.\n');

%% Test 7: Horizon parameter
rng(700);
n = 2; q = 1; py = 2; N = 200; L = 3;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
H_obs = eye(2);

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
        Y(k + 1, :, l) = (H_obs * X(k + 1, :, l)')';
    end
end

[A0_h20, ~] = sidLTIfreqIO(Y, U, H_obs, 'Horizon', 20);
[A0_h40, ~] = sidLTIfreqIO(Y, U, H_obs, 'Horizon', 40);

% Both should produce valid results
assert(all(isfinite(A0_h20(:))), 'Horizon=20 produced non-finite A');
assert(all(isfinite(A0_h40(:))), 'Horizon=40 produced non-finite A');
fprintf('  Test 7 passed: custom horizon works.\n');

%% Test 8: Mass-spring-damper, full observation, LTI
% n=6 system with all states measured. The LTI estimator should
% recover eigenvalues accurately from noiseless data.
rng(800);
[Ad, Bd] = sidTestMSD( ...
    [1.0; 1.5; 1.0], [100; 80; 60], [2; 1.5; 1], ...
    [1; 0; 0], 0.01);
n = 6; q = 1; py = 6; N = 500; L = 5;
H_full = eye(n);

X = zeros(N + 1, n, L);
U = randn(N, q, L);
Y = zeros(N + 1, py, L);
for l = 1:L
    X(1, :, l) = 0.1 * randn(1, n);
    Y(1, :, l) = (H_full * X(1, :, l)')';
    for k = 1:N
        X(k + 1, :, l) = (Ad * X(k, :, l)' ...
            + Bd * U(k, :, l)')';
        Y(k + 1, :, l) = (H_full * X(k + 1, :, l)')';
    end
end

[A0, B0] = sidLTIfreqIO(Y, U, H_full);

eig_true = sort(abs(eig(Ad)));
eig_est = sort(abs(eig(A0)));
eig_err = norm(eig_true - eig_est) / norm(eig_true);
assert(eig_err < 0.15, ...
    'MSD full obs: eigenvalue error %.4f', eig_err);
fprintf('  Test 8 passed: MSD full obs LTI (eig_err=%.4f).\n', ...
    eig_err);

%% Test 9: Mass-spring-damper, partial obs (positions only), LTI
% n=6, py=3. Only positions measured, velocities hidden.
% With 3 outputs and 1 input, the Ho-Kalman realization must
% recover 6 modes from the 3x1 transfer function. Use long data
% for accurate spectral estimation.
rng(900);
N9 = 2000; L9 = 10;
X9 = zeros(N9 + 1, n, L9);
U9 = randn(N9, q, L9);
H_pos = [eye(3), zeros(3, 3)];
Y9 = zeros(N9 + 1, 3, L9);
for l = 1:L9
    X9(1, :, l) = 0.1 * randn(1, n);
    Y9(1, :, l) = (H_pos * X9(1, :, l)')';
    for k = 1:N9
        X9(k + 1, :, l) = (Ad * X9(k, :, l)' ...
            + Bd * U9(k, :, l)')';
        Y9(k + 1, :, l) = (H_pos * X9(k + 1, :, l)')';
    end
end

[A0p, B0p] = sidLTIfreqIO(Y9, U9, H_pos);

% Compare first 10 Markov parameters (short horizon
% avoids amplifying small eigenvalue errors).
g_est_vec = zeros(10, 3);
g_true_vec = zeros(10, 3);
for k = 1:10
    g_est_vec(k, :) = (H_pos * A0p^(k - 1) * B0p)';
    g_true_vec(k, :) = (H_pos * Ad^(k - 1) * Bd)';
end
mp_err = norm(g_est_vec - g_true_vec, 'fro') / ...
    norm(g_true_vec, 'fro');
assert(mp_err < 0.5, ...
    'MSD partial obs LTI: Markov param error %.4f', mp_err);
fprintf('  Test 9 passed: MSD partial obs LTI (mp_err=%.4f).\n', ...
    mp_err);

fprintf('test_sidLTIfreqIO: all tests passed.\n');
