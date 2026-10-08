%% exampleOutputCOSMIC - Output-COSMIC LTV identification from partial observations
%
% Demonstrates sidLTVdiscIO on a natural partial-observation
% scenario: a 2-mass SMD chain (Plant B) with position sensors on
% both masses but no velocity sensors. State dimension n = 4, output
% dimension py = 2.
%
% Output-COSMIC recovers (A, B, x) up to an unobservable similarity
% transform. Element-wise comparison of recovered A/B against
% Ad/Bd is NOT meaningful; the validation metric is the gauge-
% invariant observation reconstruction error ||H*x_hat - Y|| / ||Y||.
%
% See spec/EXAMPLES.md section 3.11 for the binding specification.

runner__nCompleted = 0;

%% System setup
% Plant B: 2-mass SMD with force only at mass 1. Observation matrix
% H selects the two position channels; the two velocities are hidden.

rng(42);

m  = [1; 1];
k  = [100; 80];
c  = [2; 2];
F  = [1; 0];
Ts = 0.01;

[Ad, Bd] = util_msd(m, k, c, F, Ts);

n = 4;  q = 1;  py = 2;
H_obs = [1 0 0 0; 0 1 0 0];

fprintf('Discrete dynamics matrix Ad:\n');  disp(Ad);
fprintf('Observation matrix H (measure positions, velocities hidden):\n');
disp(H_obs);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: System setup.\n', runner__nCompleted);

%% Simulate trajectories
% L = 10 trajectories, N = 80 samples each. Input is scaled up so
% displacements reach a few centimetres.

N = 80;  L = 10;
sigma_proc = 1e-3;
sigma_meas = 1e-4;

X = zeros(N + 1, n, L);
U = 5.0 * randn(N, q, L);
Y = zeros(N + 1, py, L);

for l = 1:L
    X(1, :, l) = 0.05 * randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)' + sigma_meas * randn(py, 1))';
    for step = 1:N
        X(step + 1, :, l) = (Ad * X(step, :, l)' ...
            + Bd * U(step, :, l)' + sigma_proc * randn(n, 1))';
        Y(step + 1, :, l) = (H_obs * X(step + 1, :, l)' ...
            + sigma_meas * randn(py, 1))';
    end
end

fprintf('Max |Y|:           %.3e m (measured positions)\n', ...
    max(abs(Y(:))));
fprintf('Max |hidden v|:    %.3e m/s\n', ...
    max(max(max(abs(X(:, 3:4, :))))));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Simulate trajectories.\n', runner__nCompleted);

%% Step 1: estimate frequency response
% Use the first trajectory. Trim Y to match U length.

G = sidFreqBT(Y(1:N, :, 1), U(:, :, 1), 'WindowSize', 20, 'SampleTime', Ts);
fprintf('G.Response shape: %s\n', mat2str(size(G.Response)));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Step 1: estimate frequency response.\n', ...
    runner__nCompleted);

%% Step 2: model-order determination
% sidModelOrder estimates n from Hankel singular values. For short
% lightly-damped records the estimate can overshoot; we pin n = 4.

[n_est, svInfo] = sidModelOrder(G);
fprintf('Hankel-SVD estimate: n = %d (true = %d)\n', n_est, n);
fprintf('(For short lightly-damped records the estimate can overshoot;\n');
fprintf(' we pin n = 4 for the identification step below.)\n');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Step 2: model-order determination.\n', ...
    runner__nCompleted);

%% Step 3: construct observation matrix
H_use = H_obs;
fprintf('Observation matrix H (%d x %d):\n', size(H_use));
disp(H_use);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Step 3: construct observation matrix.\n', ...
    runner__nCompleted);

%% Step 4: identify the LTV model via sidLTVdiscIO
fprintf('Running Output-COSMIC identification...\n');
result = sidLTVdiscIO(Y, U, H_use, 'Lambda', 1e5);

fprintf('Converged in %d iterations.\n', result.Iterations);
fprintf('Final cost: %.4f\n', result.Cost(1));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Step 4: identify the LTV model.\n', ...
    runner__nCompleted);

%% Convergence history
figure;
semilogy(1:length(result.CostHistory), result.CostHistory, 'b-o', 'MarkerSize', 4);
xlabel('Iteration');
ylabel('Cost J');
title('Output-COSMIC: convergence');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Convergence history.\n', runner__nCompleted);

%% State recovery: observed channels vs hidden channels
% 2x2 grid: top row positions (measured), bottom row velocities (hidden).

figure;
t_axis = 0:N;

subplot(2, 2, 1);
plot(t_axis, squeeze(X(:, 1, 1)), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'True x_1');
hold on;
plot(t_axis, squeeze(result.X(:, 1, 1)), 'b--', 'DisplayName', 'Estimated x_1');
plot(t_axis, squeeze(Y(:, 1, 1)), 'r.', 'MarkerSize', 4, ...
    'DisplayName', 'Measurement y_1');
hold off;
xlabel('k');  ylabel('x_1 (m)');
title('Observed: position of mass 1');
legend('Location', 'southwest');  grid on;

subplot(2, 2, 2);
plot(t_axis, squeeze(X(:, 2, 1)), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'True x_2');
hold on;
plot(t_axis, squeeze(result.X(:, 2, 1)), 'b--', 'DisplayName', 'Estimated x_2');
plot(t_axis, squeeze(Y(:, 2, 1)), 'r.', 'MarkerSize', 4, ...
    'DisplayName', 'Measurement y_2');
hold off;
xlabel('k');  ylabel('x_2 (m)');
title('Observed: position of mass 2');
legend('Location', 'southwest');  grid on;

subplot(2, 2, 3);
plot(t_axis, squeeze(X(:, 3, 1)), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'True v_1');
hold on;
plot(t_axis, squeeze(result.X(:, 3, 1)), 'b--', 'DisplayName', 'Estimated v_1');
hold off;
xlabel('k');  ylabel('v_1 (m/s)');
title('Hidden: velocity of mass 1');
legend('Location', 'southwest');  grid on;

subplot(2, 2, 4);
plot(t_axis, squeeze(X(:, 4, 1)), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'True v_2');
hold on;
plot(t_axis, squeeze(result.X(:, 4, 1)), 'b--', 'DisplayName', 'Estimated v_2');
hold off;
xlabel('k');  ylabel('v_2 (m/s)');
title('Hidden: velocity of mass 2');
legend('Location', 'southwest');  grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: State recovery.\n', runner__nCompleted);

%% Validation: observation reconstruction error
% Gauge-invariant metric. Compute H*x_hat for every trajectory and
% compare against Y.

Y_recon = zeros(N + 1, py, L);
for l = 1:L
    Y_recon(:, :, l) = (H_use * squeeze(result.X(:, :, l))')';
end

obs_err = norm(Y_recon(:) - Y(:)) / norm(Y(:));
fprintf('Observation reconstruction error: %.4f\n', obs_err);
fprintf('(Relative Frobenius over %d time steps x %d channels x %d trajectories.)\n', ...
    N + 1, py, L);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Validation: observation reconstruction error.\n', ...
    runner__nCompleted);

%% Frozen-time inspection of the recovered A and B
mid_k = round(N / 2);
fprintf('A at midpoint (k = %d):\n', mid_k);
disp(result.A(:, :, mid_k));
fprintf('B at midpoint (k = %d):\n', mid_k);
disp(result.B(:, :, mid_k));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Frozen-time inspection.\n', runner__nCompleted);

fprintf('exampleOutputCOSMIC: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
