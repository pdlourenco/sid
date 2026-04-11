%% exampleLTVdisc - COSMIC LTV state-space identification on a 1-DoF SMD
%
% Full walk-through of the COSMIC workflow on a physical 1-DoF
% spring-mass-damper plant: LTI recovery, LTV recovery with ramping
% stiffness, multi-trajectory benefit, validation-based lambda
% tuning, cost decomposition, uncertainty quantification, frozen
% transfer function, and compare/residual validation. Ends with a
% weakly-nonlinear Duffing section that recovers the amplitude-
% dependent linearization.
%
% Plant A: m = 1 kg, k_baseline = 100 N/m, c = 2 N.s/m, Ts = 0.01 s.
% Plant E: same m, k, c with cubic stiffness k_cub = 1e5.
%
% For the 1-DoF plant with small Ts, Ad(2,1) ~ -k*Ts, so as k ramps
% 200 -> 50 the recovered A(2,1,k) should drift from -2.0 to -0.5.
%
% See spec/EXAMPLES.md section 3.9 for the binding specification.

runner__nCompleted = 0;

%% 1. LTI system recovery
% Plant A with constant stiffness. Simulate L trajectories with
% process noise and verify sidLTVdisc recovers Ad at high lambda.

rng(100);

m  = 1.0;    k  = 100.0;    c  = 2.0;    F  = 1.0;
Ts = 0.01;
[Ad_true, Bd_true] = util_msd(m, k, c, F, Ts);
p = 2;  q = 1;

N = 50;  L = 10;
sigma = 0.01;

X = zeros(N + 1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for step = 1:N
        X(step + 1, :, l) = (Ad_true * X(step, :, l)' ...
            + Bd_true * U(step, :, l)' + sigma * randn(p, 1))';
    end
end

result_lti = sidLTVdisc(X, U, 'Lambda', 1e5);
A_mean = mean(result_lti.A, 3);
fprintf('True Ad:\n');  disp(Ad_true);
fprintf('Mean recovered A(k):\n');  disp(A_mean);
fprintf('Recovery error (Frobenius): %.4e\n', ...
    norm(A_mean - Ad_true, 'fro'));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: LTI system recovery.\n', runner__nCompleted);

%% 2. LTV system: time-varying stiffness
% k ramps from 200 to 50 N/m over N = 250 steps. We build the per-
% step Ad_tv stack with util_msd_ltv and simulate the LTV recursion.
% Auto-lambda recovers A(2,1,k) ~ -k(t)*Ts. The long record (N=250)
% and large ensemble (L=30) are chosen so that the downstream
% residual whiteness test (Section 8) passes cleanly.

rng(200);

N_tv = 250;  L_tv = 30;
k_tv = linspace(200.0, 50.0, N_tv);
[Ad_tv, Bd_tv] = util_msd_ltv(m, k_tv, c, F, Ts);

X_tv = zeros(N_tv + 1, p, L_tv);
U_tv = randn(N_tv, q, L_tv);
for l = 1:L_tv
    X_tv(1, :, l) = randn(1, p);
    for step = 1:N_tv
        X_tv(step + 1, :, l) = (Ad_tv(:, :, step) * X_tv(step, :, l)' ...
            + Bd_tv(:, :, step) * U_tv(step, :, l)' ...
            + 0.01 * randn(p, 1))';
    end
end

result_tv = sidLTVdisc(X_tv, U_tv, 'Lambda', 'auto');
fprintf('auto lambda: %.3e\n', result_tv.Lambda(1));

figure;
kk = (1:N_tv)';
A21_true = squeeze(Ad_tv(2, 1, :));
A21_rec  = squeeze(result_tv.A(2, 1, :));
plot(kk, A21_true, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True A(2,1,k) \approx -k(t) T_s');
hold on;
plot(kk, A21_rec, 'b', 'DisplayName', 'COSMIC (auto \lambda)');
hold off;
xlabel('Time step k');
ylabel('A(2,1,k)');
title('LTV identification: recovered A(2,1) tracks the ramping stiffness');
legend('Location', 'southeast');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: LTV system: time-varying stiffness.\n', ...
    runner__nCompleted);

%% 3. Multi-trajectory benefit
% Compare L = 5 against L = 30 trajectories on the same LTV plant.

rng(300);
L_few = 5;  L_many = 30;

X_few = X_tv(:, :, 1:L_few);
U_few = U_tv(:, :, 1:L_few);

X_many = zeros(N_tv + 1, p, L_many);
U_many = randn(N_tv, q, L_many);
for l = 1:L_many
    X_many(1, :, l) = randn(1, p);
    for step = 1:N_tv
        X_many(step + 1, :, l) = (Ad_tv(:, :, step) * X_many(step, :, l)' ...
            + Bd_tv(:, :, step) * U_many(step, :, l)' ...
            + 0.01 * randn(p, 1))';
    end
end

r_few  = sidLTVdisc(X_few,  U_few,  'Lambda', 'auto');
r_many = sidLTVdisc(X_many, U_many, 'Lambda', 'auto');

figure;
plot(kk, A21_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
hold on;
plot(kk, squeeze(r_few.A(2, 1, :)),  'r', ...
    'DisplayName', sprintf('L = %d', L_few));
plot(kk, squeeze(r_many.A(2, 1, :)), 'b', ...
    'DisplayName', sprintf('L = %d', L_many));
hold off;
xlabel('Time step k');
ylabel('A(2,1,k)');
title('Effect of number of trajectories');
legend('Location', 'southeast');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Multi-trajectory benefit.\n', ...
    runner__nCompleted);

%% 4. Validation-based lambda tuning with sidLTVdiscTune
% Split L = 30 into 20 train + 10 validation, grid-search lambda.

rng(400);
L_total = 30;  L_train = 20;

X_all = zeros(N_tv + 1, p, L_total);
U_all = randn(N_tv, q, L_total);
for l = 1:L_total
    X_all(1, :, l) = randn(1, p);
    for step = 1:N_tv
        X_all(step + 1, :, l) = (Ad_tv(:, :, step) * X_all(step, :, l)' ...
            + Bd_tv(:, :, step) * U_all(step, :, l)' ...
            + 0.01 * randn(p, 1))';
    end
end
X_train = X_all(:, :, 1:L_train);
U_train = U_all(:, :, 1:L_train);
X_val   = X_all(:, :, L_train + 1:end);
U_val   = U_all(:, :, L_train + 1:end);

grid_lam = logspace(-3, 6, 30);
[bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, ...
    X_val, U_val, 'LambdaGrid', grid_lam);
fprintf('Validation-tuned lambda: %.2e\n', bestLambda);

figure;
semilogx(grid_lam, allLosses, 'b.-');
hold on;
semilogx(bestLambda, min(allLosses), 'ro', 'MarkerSize', 10);
hold off;
xlabel('\lambda');
ylabel('Validation RMSE');
title('Lambda tuning: validation loss curve');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Validation-based lambda tuning.\n', ...
    runner__nCompleted);

%% 5. Cost decomposition
fprintf('Cost decomposition:\n');
fprintf('  Total:          %.4f\n', result_tv.Cost(1));
fprintf('  Data fidelity:  %.4f\n', result_tv.Cost(2));
fprintf('  Regularisation: %.4f\n', result_tv.Cost(3));
fprintf('  Check (should be ~0): %.4e\n', ...
    result_tv.Cost(1) - result_tv.Cost(2) - result_tv.Cost(3));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Cost decomposition.\n', runner__nCompleted);

%% 6. Uncertainty quantification
% Compute Bayesian posterior uncertainty. Plot A(2,1,k) +- 2sigma
% against the true curve.

result_unc = sidLTVdisc(X_tv, U_tv, 'Lambda', 'auto', 'Uncertainty', true);
fprintf('Noise variance: %.4e\n', result_unc.NoiseVariance);
fprintf('Degrees of freedom: %.1f\n', result_unc.DegreesOfFreedom);

a21 = squeeze(result_unc.A(2, 1, :));
a21_std = squeeze(result_unc.AStd(2, 1, :));

figure;
fill([kk; flipud(kk)], ...
    [a21 - 2*a21_std; flipud(a21 + 2*a21_std)], ...
    [0.85 0.85 1], 'EdgeColor', 'none', ...
    'DisplayName', '\pm 2\sigma');
hold on;
plot(kk, a21, 'b', 'LineWidth', 1.5, 'DisplayName', 'Recovered A(2,1,k)');
plot(kk, A21_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
hold off;
xlabel('Time step k');
ylabel('A(2,1,k)');
title('LTV identification with uncertainty bands');
legend('Location', 'southeast');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Uncertainty quantification.\n', ...
    runner__nCompleted);

%% 7. Frozen transfer function with sidLTVdiscFrozen
% Instantaneous frequency response at three time steps.

kSteps = [1, round(N_tv / 2), N_tv];
frz = sidLTVdiscFrozen(result_unc, 'TimeSteps', kSteps);

fprintf('Frozen response shape: %s\n', mat2str(size(frz.Response)));
fprintf('ResponseStd available: %d\n', ~isempty(frz.ResponseStd));

figure;
w = frz.Frequency;
colors = {'b', 'r', [0 0.6 0]};
for i = 1:length(kSteps)
    G_k = squeeze(frz.Response(:, 1, 1, i));
    semilogx(w, 20*log10(abs(G_k)), 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('k = %d', kSteps(i)));
    hold on;
end
hold off;
xlabel('Frequency (rad/sample)');
ylabel('|G(\omega, k)| (dB)');
title('Frozen transfer function at selected time steps');
legend;
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Frozen transfer function.\n', ...
    runner__nCompleted);

%% 8. Model validation with sidCompare and sidResidual
% With the longer record (N = 250) and larger ensemble (L = 30)
% from Section 2, the residual whiteness test passes cleanly --
% COSMIC recovers a well-identified LTV model with residuals that
% look white at the 99% confidence level.

comp = sidCompare(result_tv, X_tv, U_tv);
fprintf('COSMIC model fit (per state component):\n');
for ch = 1:p
    fprintf('  x_%d: %.1f%%\n', ch, comp.Fit(ch));
end

resid = sidResidual(result_tv, X_tv, U_tv);
if resid.WhitenessPass
    fprintf('Residual whiteness test: PASS\n');
else
    fprintf('Residual whiteness test: FAIL\n');
end

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Model validation.\n', runner__nCompleted);

%% 9. Weakly-nonlinear Duffing oscillator
% The cubic stiffness means the effective stiffness depends on
% displacement amplitude. We drive the plant with a ramped-amplitude
% white force (0.5 -> 8 over the record) so typical displacement
% grows through the record, and fit an LTV model. With a small
% manual lambda (0.1), COSMIC recovers an A(2,1) that drifts toward
% more negative values (stiffer apparent spring) as amplitude grows.

rng(500);

m_nl     = 1.0;
k_lin    = 100.0;
k_cubic  = 1e5;
c_nl     = 2.0;
F_nl     = 1.0;

N_nl = 400;  L_nl = 12;
amp_profile = linspace(0.5, 8.0, N_nl)';

X_nl = zeros(N_nl + 1, p, L_nl);
U_nl = zeros(N_nl, q, L_nl);
for l = 1:L_nl
    U_nl(:, 1, l) = amp_profile .* randn(N_nl, 1);
    x_nl = util_msd_nl(m_nl, k_lin, k_cubic, c_nl, F_nl, Ts, ...
                       U_nl(:, 1, l), [], 4);
    X_nl(:, :, l) = x_nl;
end

% Auto-lambda is over-regularising on this dataset; use a small manual
% lambda to let COSMIC track the amplitude-dependent drift.
result_nl = sidLTVdisc(X_nl, U_nl, 'Lambda', 0.1);

[Ad_linear, ~] = util_msd(m_nl, k_lin, c_nl, F_nl, Ts);
fprintf('Linear Ad(2,1) (small amplitude): %.4f\n', Ad_linear(2, 1));
fprintf('COSMIC A(2,1) early (small x):    %.4f\n', ...
    mean(result_nl.A(2, 1, 1:N_nl/4)));
fprintf('COSMIC A(2,1) late  (large x):    %.4f\n', ...
    mean(result_nl.A(2, 1, end - N_nl/4 + 1:end)));

figure;
subplot(2, 1, 1);
plot((1:N_nl)', amp_profile, 'r');
ylabel('Excitation amplitude (N)');
title('Duffing oscillator: amplitude grows over the record');
grid on;

subplot(2, 1, 2);
plot((1:N_nl)', squeeze(result_nl.A(2, 1, :)), 'b', ...
    'DisplayName', 'COSMIC A(2,1,k)');
hold on;
plot([1 N_nl], [Ad_linear(2,1) Ad_linear(2,1)], 'k--', ...
    'DisplayName', 'Small-amplitude linearisation');
hold off;
xlabel('Time step k');
ylabel('A(2,1,k)');
title('Recovered local linearisation stiffens with amplitude');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Weakly-nonlinear Duffing oscillator.\n', ...
    runner__nCompleted);

fprintf('exampleLTVdisc: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
