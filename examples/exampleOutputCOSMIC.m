%% exampleOutputCOSMIC - LTV identification from partial observations
%
% This example demonstrates the Output-COSMIC workflow for identifying
% time-varying systems when only partial state measurements are available:
%
%   x(k+1) = A(k) x(k) + B(k) u(k)     (unknown dynamics)
%   y(k)   = H x(k)                      (partial observation)
%
% The workflow is:
%   1. Estimate frequency response from input-output data
%   2. Determine model order (state dimension) via Hankel SVD
%   3. Construct observation matrix H
%   4. Identify LTV system matrices using sidLTVdiscIO

%% System setup
% A 4th-order LTI system with 2 measured outputs (py = 2, n = 4).

rng(42);
n = 4;   % state dimension
q = 1;   % single input
py = 2;  % two measured outputs

% Stable dynamics: two pairs of complex poles
A_true = [0.8  0.15  0    0
         -0.15 0.8   0    0
          0    0     0.6  0.2
          0    0    -0.2  0.6];

B_true = [1; 0.3; 0.5; 0.1];

% Observation matrix: measure states 1 and 3
H_obs = [1 0 0 0
         0 0 1 0];

%% Simulate trajectories
N = 80;   % time steps per trajectory
L = 10;   % number of trajectories
sigma_proc = 0.02;  % process noise
sigma_meas = 0.05;  % measurement noise

X = zeros(N+1, n, L);
U = randn(N, q, L);
Y = zeros(N+1, py, L);

for l = 1:L
    X(1, :, l) = 0.5 * randn(1, n);
    Y(1, :, l) = (H_obs * X(1, :, l)')' + sigma_meas * randn(1, py);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' ...
            + sigma_proc * randn(1, n);
        Y(k+1, :, l) = (H_obs * X(k+1, :, l)')' + sigma_meas * randn(1, py);
    end
end

fprintf('Simulated %d trajectories of length %d.\n', L, N);
fprintf('True state dimension: n = %d, observed: py = %d.\n', n, py);

%% Step 1: Estimate frequency response
% Use the first trajectory for frequency response estimation.
G = sidFreqBT(Y(:,:,1), U(:,:,1), 'WindowSize', 20);

%% Step 2: Model order determination
% sidModelOrder estimates n from the Hankel singular values.
[n_est, sv] = sidModelOrder(G, 'Plot', true);
fprintf('Estimated model order: n = %d (true = %d).\n', n_est, n);

%% Step 3: Construct observation matrix
% In practice, H encodes which states are measured.
% Here we use the true H since we know which channels correspond to which states.
H_use = H_obs;

%% Step 4: Identify LTV system from partial observations
fprintf('Running Output-COSMIC identification...\n');
result = sidLTVdiscIO(Y, U, H_use, 'Lambda', 1e5);

fprintf('Converged in %d iterations.\n', result.Iterations);
fprintf('Final cost: %.4f\n', result.Cost(end));

%% Compare recovered A to true A
A_mean = mean(result.A, 3);
B_mean = mean(result.B, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
errB = norm(B_mean - B_true, 'fro') / norm(B_true, 'fro');
fprintf('A recovery error (relative Frobenius): %.4f\n', errA);
fprintf('B recovery error (relative Frobenius): %.4f\n', errB);

%% Plot: recovered A(1,1) over time vs true value
figure;
plot(1:N, squeeze(result.A(1,1,:)), 'b-', 'DisplayName', 'Recovered A_{11}(k)');
hold on;
yline(A_true(1,1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'True A_{11}');
xlabel('Time step k');
ylabel('A_{11}(k)');
title('Output-COSMIC: Recovered Dynamics (LTI Case)');
legend('show');
grid on;
hold off;

%% Plot: estimated vs true states (trajectory 1)
figure;
subplot(2,1,1);
plot(0:N, X(:,1,1), 'k-', 'LineWidth', 1.5, 'DisplayName', 'True x_1');
hold on;
plot(0:N, result.X(:,1,1), 'b--', 'DisplayName', 'Estimated x_1');
plot(0:N, Y(:,1,1), 'r.', 'MarkerSize', 4, 'DisplayName', 'Measurement y_1');
legend('show'); xlabel('k'); ylabel('x_1(k)');
title('State Recovery: Observed State');
grid on; hold off;

subplot(2,1,2);
plot(0:N, X(:,2,1), 'k-', 'LineWidth', 1.5, 'DisplayName', 'True x_2');
hold on;
plot(0:N, result.X(:,2,1), 'b--', 'DisplayName', 'Estimated x_2');
legend('show'); xlabel('k'); ylabel('x_2(k)');
title('State Recovery: Hidden State');
grid on; hold off;

%% Plot: convergence history
figure;
semilogy(1:length(result.Cost), result.Cost, 'b-o', 'MarkerSize', 4);
xlabel('Iteration');
ylabel('Cost J');
title('Output-COSMIC: Convergence');
grid on;

%% Frozen transfer function (if available)
% Compute the frozen transfer function at the midpoint and plot.
midResult = struct();
midResult.A = result.A(:,:,round(N/2));
midResult.B = result.B(:,:,round(N/2));
midResult.Lambda = result.Lambda;
midResult.DataLength = 1;
midResult.StateDim = n;
midResult.InputDim = q;
midResult.NumTrajectories = 1;
midResult.Algorithm = 'cosmic';
midResult.Preconditioned = false;
midResult.Method = 'sidLTVdisc';
midResult.Cost = [0, 0, 0];

fprintf('\nExample completed successfully.\n');
