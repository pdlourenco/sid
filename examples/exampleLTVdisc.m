%% exampleLTVdisc - LTV state-space identification with sidLTVdisc
%
% This example demonstrates sidLTVdisc (COSMIC algorithm) for identifying
% time-varying discrete linear systems of the form:
%
%   x(k+1) = A(k) x(k) + B(k) u(k)
%
% It also demonstrates sidLTVdiscTune for automatic regularization tuning.

%% LTI system recovery
% First, verify that sidLTVdisc correctly identifies a known LTI system
% where A and B are constant.

rng(100);
p = 2;  q = 1;  N = 50;  L = 10;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [1; 0.5];
sigma = 0.02;

% Generate L trajectories
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)' ...
            + sigma * randn(p, 1))';
    end
end

result_lti = sidLTVdisc(X, U, 'Lambda', 1e5);

% A(:,:,k) should be approximately constant
A_mean = mean(result_lti.A, 3);
fprintf('True A:\n');  disp(A_true);
fprintf('Mean recovered A:\n');  disp(A_mean);
fprintf('Recovery error: %.4f\n', norm(A_mean - A_true, 'fro'));

%% LTV system: time-varying pole
% The (1,1) entry of A ramps from 0.5 to 0.9 over time.

rng(200);
N = 80;  L = 15;
a_ramp = linspace(0.5, 0.9, N)';
A_tv = zeros(p, p, N);
for k = 1:N
    A_tv(:,:,k) = [a_ramp(k), 0.1; -0.1, 0.8];
end

X_tv = zeros(N+1, p, L);
U_tv = randn(N, q, L);
for l = 1:L
    X_tv(1, :, l) = randn(1, p);
    for k = 1:N
        X_tv(k+1, :, l) = (A_tv(:,:,k) * X_tv(k, :, l)' ...
            + B_true * U_tv(k, :, l)' + 0.05 * randn(p, 1))';
    end
end

result_tv = sidLTVdisc(X_tv, U_tv, 'Lambda', 1e4);

% Plot recovered A(1,1,k) vs true ramp
figure;
plot(1:N, squeeze(result_tv.A(1,1,:)), 'b', 'DisplayName', 'Recovered A_{11}(k)');
hold on;
plot(1:N, a_ramp, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True A_{11}(k)');
xlabel('Time step k');
ylabel('A_{11}(k)');
title('LTV Identification: Time-Varying Pole Recovery');
legend('show');
grid on;
hold off;

%% Automatic lambda selection (L-curve)
% With 'Lambda', 'auto', sidLTVdisc uses the L-curve method to find the
% regularization that best balances data fidelity and smoothness.

result_auto = sidLTVdisc(X_tv, U_tv, 'Lambda', 'auto');
fprintf('\nAutomatic lambda: %.2e\n', result_auto.Lambda(1));

figure;
plot(1:N, squeeze(result_auto.A(1,1,:)), 'r', 'DisplayName', 'Auto \lambda');
hold on;
plot(1:N, squeeze(result_tv.A(1,1,:)), 'b', 'DisplayName', '\lambda = 10^4');
plot(1:N, a_ramp, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
xlabel('Time step k');
ylabel('A_{11}(k)');
title('Manual vs Automatic Lambda');
legend('show');
grid on;
hold off;

%% Multi-trajectory benefit
% More trajectories provide more information and reduce estimation error.

rng(300);
L_few  = 3;
L_many = 20;

X_few  = X_tv(:, :, 1:L_few);
U_few  = U_tv(:, :, 1:L_few);

X_many = zeros(N+1, p, L_many);
U_many = randn(N, q, L_many);
for l = 1:L_many
    X_many(1, :, l) = randn(1, p);
    for k = 1:N
        X_many(k+1, :, l) = (A_tv(:,:,k) * X_many(k, :, l)' ...
            + B_true * U_many(k, :, l)' + 0.05 * randn(p, 1))';
    end
end

r_few  = sidLTVdisc(X_few,  U_few,  'Lambda', 1e4);
r_many = sidLTVdisc(X_many, U_many, 'Lambda', 1e4);

figure;
plot(1:N, squeeze(r_few.A(1,1,:)),  'r', 'DisplayName', sprintf('L = %d', L_few));
hold on;
plot(1:N, squeeze(r_many.A(1,1,:)), 'b', 'DisplayName', sprintf('L = %d', L_many));
plot(1:N, a_ramp, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
xlabel('Time step k');
ylabel('A_{11}(k)');
title('Effect of Number of Trajectories');
legend('show');
grid on;
hold off;

%% Validation-based lambda tuning with sidLTVdiscTune
% Split trajectories into training and validation sets, then search over
% a grid of lambda values to minimize validation prediction error.

rng(400);
L_total = 20;  L_train = 14;  L_val = L_total - L_train;

X_all = zeros(N+1, p, L_total);
U_all = randn(N, q, L_total);
for l = 1:L_total
    X_all(1, :, l) = randn(1, p);
    for k = 1:N
        X_all(k+1, :, l) = (A_tv(:,:,k) * X_all(k, :, l)' ...
            + B_true * U_all(k, :, l)' + 0.05 * randn(p, 1))';
    end
end

X_train = X_all(:, :, 1:L_train);
U_train = U_all(:, :, 1:L_train);
X_val   = X_all(:, :, L_train+1:end);
U_val   = U_all(:, :, L_train+1:end);

grid_lam = logspace(0, 10, 30);
[bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, ...
    X_val, U_val, 'LambdaGrid', grid_lam);

fprintf('\nValidation-tuned lambda: %.2e\n', bestLambda);

% Plot validation RMSE vs lambda
figure;
semilogx(grid_lam, allLosses, 'b.-');
hold on;
semilogx(bestLambda, min(allLosses), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlabel('\lambda');
ylabel('Validation RMSE');
title('Lambda Tuning: Validation Loss Curve');
grid on;

%% Preconditioning for numerical stability
% Block-diagonal preconditioning can improve conditioning of the solve.

result_pre = sidLTVdisc(X_tv, U_tv, 'Lambda', 1e4, 'Precondition', true);
fprintf('\nPreconditioned: %d\n', result_pre.Preconditioned);

%% Cost decomposition
% result.Cost = [total, data_fidelity, regularization]

fprintf('\nCost decomposition:\n');
fprintf('  Total:          %.4f\n', result_tv.Cost(1));
fprintf('  Data fidelity:  %.4f\n', result_tv.Cost(2));
fprintf('  Regularization: %.4f\n', result_tv.Cost(3));
fprintf('  Check: %.4e (should be ~0)\n', ...
    result_tv.Cost(1) - result_tv.Cost(2) - result_tv.Cost(3));

%% Uncertainty quantification
% Enable Bayesian posterior uncertainty to get standard deviations for each
% A(k) and B(k) entry. The noise covariance is estimated from residuals.

result_unc = sidLTVdisc(X_tv, U_tv, 'Lambda', 1e4, 'Uncertainty', true);

fprintf('\nUncertainty results:\n');
fprintf('  Noise covariance estimated: %d\n', result_unc.NoiseCovEstimated);
fprintf('  Noise variance (trace/p):   %.6f\n', result_unc.NoiseVariance);
fprintf('  Degrees of freedom:         %.1f\n', result_unc.DegreesOfFreedom);
fprintf('  NoiseCov:\n');  disp(result_unc.NoiseCov);

% Plot A(1,1,k) with +/-2 sigma uncertainty band
a11     = squeeze(result_unc.A(1,1,:));
a11_std = squeeze(result_unc.AStd(1,1,:));
kk = (1:N)';

figure;
fill([kk; flipud(kk)], [a11 - 2*a11_std; flipud(a11 + 2*a11_std)], ...
    [0.8 0.8 1], 'EdgeColor', 'none', 'DisplayName', '\pm 2\sigma');
hold on;
plot(kk, a11, 'b', 'LineWidth', 1.5, 'DisplayName', 'Recovered A_{11}(k)');
plot(kk, a_ramp, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True A_{11}(k)');
xlabel('Time step k');
ylabel('A_{11}(k)');
title('LTV Identification with Uncertainty Bands');
legend('show', 'Location', 'southeast');
grid on;
hold off;

%% Frozen transfer function with sidLTVdiscFrozen
% Compute the instantaneous frequency response G(w,k) = (e^{jw}I - A(k))^{-1} B(k)
% at selected time steps. When the LTV result includes uncertainty,
% ResponseStd is propagated via first-order linearization.

kSteps = [1, round(N/2), N];
frz = sidLTVdiscFrozen(result_unc, 'TimeSteps', kSteps);

fprintf('\nFrozen transfer function:\n');
fprintf('  Method: %s\n', frz.Method);
fprintf('  Response size: %s\n', mat2str(size(frz.Response)));
fprintf('  ResponseStd available: %d\n', ~isempty(frz.ResponseStd));

% Plot Bode magnitude at the three time steps
figure;
w = frz.Frequency;
colors = {'b', 'r', [0 0.6 0]};
for i = 1:length(kSteps)
    G_k = squeeze(frz.Response(:, 1, 1, i));
    mag_dB = 20*log10(abs(G_k));
    semilogx(w, mag_dB, 'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('k = %d', kSteps(i)));
    hold on;
end

% Add +/-2 sigma band for the last time step
G_last = squeeze(frz.Response(:, 1, 1, length(kSteps)));
G_std  = squeeze(frz.ResponseStd(:, 1, 1, length(kSteps)));
mag_upper = 20*log10(abs(G_last) + 2*G_std);
mag_lower = 20*log10(max(abs(G_last) - 2*G_std, eps));
fill([w; flipud(w)], [mag_lower; flipud(mag_upper)], ...
    [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'HandleVisibility', 'off');
% Re-plot the line on top
semilogx(w, 20*log10(abs(G_last)), 'Color', colors{end}, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');

xlabel('Frequency (rad/sample)');
ylabel('|G(w,k)| (dB)');
title('Frozen Transfer Function at Selected Time Steps');
legend('show');
grid on;
hold off;


%% 6. Frequency-Based Lambda Tuning (no validation data needed)
%
% sidLTVdiscTune with 'Method','frequency' selects lambda by comparing
% the COSMIC frozen transfer function against a non-parametric sidFreqMap
% estimate. No validation data is required — only training data.

fprintf('\n--- Frequency-based lambda tuning ---\n');

% Use a smaller lambda grid for speed in this demo
grid_freq = logspace(1, 8, 12);
[bestResult_freq, bestLambda_freq, info_freq] = sidLTVdiscTune(X, U, ...
    'Method', 'frequency', 'LambdaGrid', grid_freq, 'SegmentLength', 20);

fprintf('Frequency-tuned lambda: %.2e (consistency: %.1f%%)\n', ...
    bestLambda_freq, info_freq.bestFraction * 100);
fprintf('Validation-tuned lambda: %.2e\n', bestLambda);

% Plot consistency fraction vs lambda
figure;
semilogx(info_freq.lambdaGrid, info_freq.fractions * 100, 'b.-', 'MarkerSize', 10);
hold on;
semilogx(bestLambda_freq, info_freq.bestFraction * 100, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(xlim, [90 90], 'k--');  % 90% threshold line
hold off;
xlabel('\lambda');
ylabel('Consistency fraction (%)');
title('Frequency-Response Consistency Scoring');
legend('Consistency', 'Selected \lambda', '90% threshold');
grid on;

fprintf('\nFrequency-based tuning complete.\n');

%% Model validation with sidCompare and sidResidual
% Compare COSMIC model predictions against the training data.

comp = sidCompare(result_auto, X_tv, U_tv);
fprintf('\nCOSMIC model fit (per state component):\n');
for ch = 1:p
    fprintf('  x_%d: %.1f%%\n', ch, comp.Fit(ch));
end

% Residual diagnostics
resid = sidResidual(result_auto, X_tv, U_tv);
if resid.WhitenessPass
    fprintf('Residual whiteness test: PASS\n');
else
    fprintf('Residual whiteness test: FAIL (model may need refinement)\n');
end
