%% test_sidLTVStateEst - Unit tests for batch LTV state estimation
%
% Tests sidLTVStateEst on physically meaningful systems:
% double integrator and 3-mass spring-damper chain, with
% full and partial observation, constant and time-varying dynamics.

fprintf('Running test_sidLTVStateEst...\n');

%% Test 1: Double integrator, full observation, noiseless
% With H = I and exact dynamics, the estimator should recover
% the true state perfectly.
n = 2; q = 1; py = 2; N = 50; L = 1;
dt = 1;
A_di = [1 dt; 0 1];
B_di = [0.5 * dt^2; dt];
H_full = eye(2);

X_true = zeros(N + 1, n);
U_di = zeros(N, q);
for k = 1:N
    if k <= 20,     U_di(k) = 1.0;
    elseif k <= 40, U_di(k) = 0.0;
    else,           U_di(k) = -1.0;
    end
end
X_true(1, :) = [0, 0];
for k = 1:N
    X_true(k + 1, :) = (A_di * X_true(k, :)' ...
        + B_di * U_di(k))';
end
Y_di = X_true * H_full';

A_rep = repmat(A_di, [1, 1, N]);
B_rep = repmat(B_di, [1, 1, N]);

X_hat = sidLTVStateEst(Y_di, U_di, A_rep, B_rep, H_full);

err = norm(X_hat - X_true, 'fro') / norm(X_true, 'fro');
assert(err < 1e-6, ...
    'Double integrator full obs: error %.2e', err);
fprintf('  Test 1 passed: DI full obs (err=%.2e).\n', err);

%% Test 2: Double integrator, partial obs (position only), noiseless
H_pos = [1 0];
Y_pos = X_true * H_pos';

X_hat = sidLTVStateEst(Y_pos, U_di, A_rep, B_rep, H_pos);

pos_err = norm(X_hat(:, 1) - X_true(:, 1)) / ...
    norm(X_true(:, 1));
vel_err = norm(X_hat(:, 2) - X_true(:, 2)) / ...
    max(norm(X_true(:, 2)), 1);
assert(pos_err < 1e-3, ...
    'DI partial obs: position error %.4f', pos_err);
assert(vel_err < 0.05, ...
    'DI partial obs: velocity error %.4f', vel_err);
fprintf('  Test 2 passed: DI partial obs (pos=%.4f, vel=%.4f).\n', ...
    pos_err, vel_err);

%% Test 3: Mass-spring-damper, full obs, noiseless
[Ad, Bd] = sidTestMSD( ...
    [1.0; 1.5; 1.0], [100; 80; 60], [2; 1.5; 1], ...
    [1; 0; 0], 0.01);
n = 6; q = 1; py = 6; N = 200; L = 1;
H_full = eye(n);

X_true = zeros(N + 1, n);
U_msd = zeros(N, q);
for k = 1:N
    U_msd(k) = 5 * sin(2 * pi * k / 50);
end
X_true(1, :) = [0.1, 0, -0.05, 0, 0, 0];
for k = 1:N
    X_true(k + 1, :) = (Ad * X_true(k, :)' ...
        + Bd * U_msd(k))';
end
Y_msd = X_true * H_full';

A_rep = repmat(Ad, [1, 1, N]);
B_rep = repmat(Bd, [1, 1, N]);

X_hat = sidLTVStateEst( ...
    Y_msd, U_msd, A_rep, B_rep, H_full);

err = norm(X_hat - X_true, 'fro') / norm(X_true, 'fro');
assert(err < 1e-4, ...
    'MSD full obs: error %.2e', err);
fprintf('  Test 3 passed: MSD full obs (err=%.2e).\n', err);

%% Test 4: Mass-spring-damper, partial obs (positions only)
H_pos = [eye(3), zeros(3, 3)];
py = 3;
Y_pos = X_true * H_pos';

X_hat = sidLTVStateEst( ...
    Y_pos, U_msd, A_rep, B_rep, H_pos);

pos_err = norm(X_hat(:, 1:3) - X_true(:, 1:3), 'fro') / ...
    norm(X_true(:, 1:3), 'fro');
vel_err = norm(X_hat(:, 4:6) - X_true(:, 4:6), 'fro') / ...
    max(norm(X_true(:, 4:6), 'fro'), 1e-6);
assert(pos_err < 0.01, ...
    'MSD partial obs: position error %.4f', pos_err);
assert(vel_err < 0.1, ...
    'MSD partial obs: velocity error %.4f', vel_err);
fprintf('  Test 4 passed: MSD partial obs (pos=%.4f, vel=%.4f).\n', ...
    pos_err, vel_err);

%% Test 5: Mass-spring-damper, partial obs, noisy measurements
% Dynamics are exact (no process noise), so set Q small to
% tell the estimator to trust dynamics over measurements.
rng(500);
sigma_meas = 0.01;
Y_noisy = Y_pos + sigma_meas * randn(size(Y_pos));
R_meas = sigma_meas^2 * eye(3);
Q_small = 1e-6 * eye(n);

X_hat = sidLTVStateEst( ...
    Y_noisy, U_msd, A_rep, B_rep, H_pos, ...
    'R', R_meas, 'Q', Q_small);

pos_err = norm(X_hat(:, 1:3) - X_true(:, 1:3), 'fro') / ...
    norm(X_true(:, 1:3), 'fro');
vel_err = norm(X_hat(:, 4:6) - X_true(:, 4:6), 'fro') / ...
    max(norm(X_true(:, 4:6), 'fro'), 1e-6);
assert(pos_err < 0.1, ...
    'MSD noisy: position error %.4f', pos_err);
assert(vel_err < 0.3, ...
    'MSD noisy: velocity error %.4f', vel_err);
fprintf('  Test 5 passed: MSD noisy (pos=%.4f, vel=%.4f).\n', ...
    pos_err, vel_err);

%% Test 6: Time-varying double integrator, partial obs
% Gain alpha(k) increases linearly over time.
n = 2; q = 1; py = 1; N = 80; L = 1;
dt = 1;
A_tv = zeros(n, n, N);
B_tv = zeros(n, q, N);
for k = 1:N
    alpha = 1.0 + 0.5 * (k - 1) / (N - 1);
    A_tv(:, :, k) = [1, alpha * dt; 0, 1];
    B_tv(:, :, k) = [0.5 * alpha * dt^2; alpha * dt];
end

X_true = zeros(N + 1, n);
U_tv = zeros(N, q);
for k = 1:N
    U_tv(k) = sin(2 * pi * k / 40);
end
X_true(1, :) = [0, 0];
for k = 1:N
    X_true(k + 1, :) = (A_tv(:, :, k) * X_true(k, :)' ...
        + B_tv(:, :, k) * U_tv(k))';
end
H_pos = [1 0];
Y_tv = X_true * H_pos';

X_hat = sidLTVStateEst(Y_tv, U_tv, A_tv, B_tv, H_pos);

pos_err = norm(X_hat(:, 1) - X_true(:, 1)) / ...
    norm(X_true(:, 1));
vel_err = norm(X_hat(:, 2) - X_true(:, 2)) / ...
    max(norm(X_true(:, 2)), 1);
assert(pos_err < 0.01, ...
    'TV DI partial obs: position error %.4f', pos_err);
assert(vel_err < 0.1, ...
    'TV DI partial obs: velocity error %.4f', vel_err);
fprintf('  Test 6 passed: TV DI partial obs (pos=%.4f, vel=%.4f).\n', ...
    pos_err, vel_err);

%% Test 7: Multi-trajectory mass-spring-damper
rng(700);
[Ad, Bd] = sidTestMSD( ...
    [1.0; 1.5; 1.0], [100; 80; 60], [2; 1.5; 1], ...
    [1; 0; 0], 0.01);
n = 6; q = 1; py = 3; N = 100; L = 5;
H_pos = [eye(3), zeros(3, 3)];

X_true = zeros(N + 1, n, L);
U_mt = randn(N, q, L);
Y_mt = zeros(N + 1, py, L);
for l = 1:L
    X_true(1, :, l) = 0.1 * randn(1, n);
    Y_mt(1, :, l) = (H_pos * X_true(1, :, l)')';
    for k = 1:N
        X_true(k + 1, :, l) = (Ad * X_true(k, :, l)' ...
            + Bd * U_mt(k, :, l)')';
        Y_mt(k + 1, :, l) = ...
            (H_pos * X_true(k + 1, :, l)')';
    end
end

A_rep = repmat(Ad, [1, 1, N]);
B_rep = repmat(Bd, [1, 1, N]);

X_hat = sidLTVStateEst(Y_mt, U_mt, A_rep, B_rep, H_pos);

assert(isequal(size(X_hat), [N + 1, n, L]), ...
    'Multi-traj: wrong output size');
pos_diff = X_hat(:, 1:3, :) - X_true(:, 1:3, :);
pos_ref = X_true(:, 1:3, :);
pos_err = norm(pos_diff(:)) / norm(pos_ref(:));
assert(pos_err < 0.01, ...
    'Multi-traj: position error %.4f', pos_err);
fprintf('  Test 7 passed: multi-traj MSD (pos=%.4f).\n', ...
    pos_err);

fprintf('test_sidLTVStateEst: all tests passed.\n');
