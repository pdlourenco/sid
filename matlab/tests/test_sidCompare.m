% test_sidCompare.m - Test model output comparison function

fprintf('Running test_sidCompare...\n');
runner__nPassed = 0;

%% Test 1: Perfect state-space model → ~100% fit
rng(7001);
p = 2; q = 1; N = 50; L = 3;
A_true = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
    end
end

% Build a "perfect" model struct
perfect_model.A = repmat(A_true, [1, 1, N]);
perfect_model.B = repmat(B_true, [1, 1, N]);
perfect_model.DataLength = N;
perfect_model.StateDim = p;
perfect_model.InputDim = q;
perfect_model.Method = 'sidLTVdisc';

comp = sidCompare(perfect_model, X, U);

assert(all(comp.Fit > 99), ...
    'Perfect model should give >99%% fit, got [%.1f, %.1f]', comp.Fit(1), comp.Fit(2));
% Multi-trajectory Predicted/Residual are per-trajectory (N x p x L), not the
% ensemble mean — SPEC §15.3, issue #140.
assert(isequal(size(comp.Predicted), [N, p, L]), 'Predicted should be N x p x L');
assert(isequal(size(comp.Residual), [N, p, L]), 'Residual should be N x p x L');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: perfect state-space model gives %.1f%% fit.\n', min(comp.Fit));

%% Test 1b: mirror-image trajectories are perfectly identifiable, not NaN (#140)
% X2 = [X, -X], U2 = [U, -U] sum to zero; the old signal-averaging path gave
% measured == 0 and fit = NaN. The per-trajectory fit is ~100% for each.
rng(7011);
p = 2; q = 1; N = 60;
A_true = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
Xs = zeros(N+1, p); Us = randn(N, q); Xs(1, :) = randn(1, p);
for k = 1:N
    Xs(k+1, :) = (A_true * Xs(k, :)' + B_true * Us(k, :)')';
end
mirror_model.A = repmat(A_true, [1, 1, N]);
mirror_model.B = repmat(B_true, [1, 1, N]);
mirror_model.DataLength = N;
mirror_model.StateDim = p;
mirror_model.InputDim = q;
mirror_model.Method = 'sidLTVdisc';
Xmir = cat(3, Xs, -Xs); Umir = cat(3, Us, -Us);
compMir = sidCompare(mirror_model, Xmir, Umir);
assert(all(isfinite(compMir.Fit)) && all(compMir.Fit > 99), ...
    'Mirror-image trajectories should give ~100%% fit, got [%.2f, %.2f]', ...
    compMir.Fit(1), compMir.Fit(2));
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1b passed: mirror-image trajectories fit %.1f%% (#140).\n', ...
    min(compMir.Fit));

%% Test 2: Frequency-domain model fit
rng(7002);
N = 2000;
u = randn(N, 1);
y = filter(1, [1 -0.85], u) + 0.05 * randn(N, 1);

result = sidFreqBT(y, u, 'WindowSize', 40);
comp = sidCompare(result, y, u);

assert(comp.Fit > 50, 'Good freq-domain model should give >50%% fit, got %.1f%%', comp.Fit);
assert(isequal(size(comp.Predicted), [N, 1]), 'Predicted should be N x 1');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: frequency-domain model fit = %.1f%%.\n', comp.Fit);

%% Test 3: State-space model on noisy validation data
rng(7003);
p = 2; q = 1; N = 60; L = 5;
A_true = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
sigma = 0.03;
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + sigma * randn(1, p);
    end
end

ltv = sidLTVdisc(X, U, 'Lambda', 1e4);
comp_ss = sidCompare(ltv, X, U);

assert(all(comp_ss.Fit > 50), 'COSMIC model fit should be decent, got %.1f%%', min(comp_ss.Fit));
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: COSMIC model on noisy data, fit = [%.1f%%, %.1f%%].\n', ...
    comp_ss.Fit(1), comp_ss.Fit(2));

%% Test 4: Multi-channel fit percentages
rng(7004);
N = 2000;
u = randn(N, 1);
y1 = filter(1, [1 -0.5], u) + 0.1 * randn(N, 1);
y2 = filter(0.3, [1 -0.7], u) + 0.5 * randn(N, 1);  % noisier
y = [y1, y2];

result_mimo = sidFreqBT(y, u, 'WindowSize', 30);
comp_mimo = sidCompare(result_mimo, y, u);

assert(length(comp_mimo.Fit) == 2, 'Should have per-channel fit');
% Channel 1 (cleaner) should have better fit than channel 2 (noisier)
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: multi-channel fit = [%.1f%%, %.1f%%].\n', ...
    comp_mimo.Fit(1), comp_mimo.Fit(2));

%% Test 5: Plot option (headless check)
rng(7005);
N = 500;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.1 * randn(N, 1);
result_plot = sidFreqBT(y, u);

try
    sidCompare(result_plot, y, u, 'Plot', true);
    close all;
    plotOk = true;
catch
    plotOk = false;
end
assert(plotOk, 'Plot option should not error');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: plot option works.\n');

%% Test 6: Output struct fields
rng(7006);
N = 300;
u = randn(N, 1);
y = filter(1, [1 -0.7], u) + 0.1 * randn(N, 1);
result = sidFreqBT(y, u);
comp = sidCompare(result, y, u);

assert(isfield(comp, 'Predicted'), 'Should have Predicted field');
assert(isfield(comp, 'Measured'), 'Should have Measured field');
assert(isfield(comp, 'Fit'), 'Should have Fit field');
assert(isfield(comp, 'Residual'), 'Should have Residual field');
assert(isfield(comp, 'Method'), 'Should have Method field');
assert(strcmp(comp.Method, 'sidFreqBT'), 'Method should match source model');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: output struct has correct fields.\n');

%% Test 7: NRMSE formula verification (SPEC §15.3)
% fit = 100 * (1 - ||y - y_pred|| / ||y - mean(y)||)
rng(7007);
N = 500;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.05 * randn(N, 1);
result_f = sidFreqBT(y, u, 'WindowSize', 30);
comp_f = sidCompare(result_f, y, u);

% Manually compute NRMSE from predicted and measured
y_pred = comp_f.Predicted;
y_meas = comp_f.Measured;
manual_fit = 100 * (1 - norm(y_meas - y_pred) / norm(y_meas - mean(y_meas)));
assert(abs(comp_f.Fit - manual_fit) < 1e-8, ...
    'Fit=%.4f should match manual NRMSE=%.4f', comp_f.Fit, manual_fit);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: NRMSE formula verified (fit=%.1f%%).\n', ...
    comp_f.Fit);

%% Test 8: Negative fit — model worse than mean (SPEC §15.3)
rng(7008);
N = 500;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.05 * randn(N, 1);

% Build a deliberately terrible model by using wrong system
y_wrong = filter(1, [1 0.99], u);  % completely wrong dynamics
result_wrong = sidFreqBT(y_wrong, u, 'WindowSize', 10);
comp_wrong = sidCompare(result_wrong, y, u);

% Fit can be negative when model is worse than predicting the mean
assert(comp_wrong.Fit < 50, ...
    'Wrong model should give poor fit, got %.1f%%', comp_wrong.Fit);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: poor model fit = %.1f%%.\n', comp_wrong.Fit);

%% Test 9: NRMSE per channel for MIMO (SPEC §15.3)
rng(7009);
N = 1000;
u = randn(N, 1);
y1 = filter(1, [1 -0.9], u);             % very clean
y2 = filter(0.3, [1 -0.5], u) + randn(N, 1);  % very noisy
y = [y1, y2];

result_m = sidFreqBT(y, u, 'WindowSize', 40);
comp_m = sidCompare(result_m, y, u);

assert(length(comp_m.Fit) == 2, 'Should have 2 fit values');
% Channel 1 (clean) should have much better fit than channel 2 (noisy)
assert(comp_m.Fit(1) > comp_m.Fit(2), ...
    'Clean channel fit %.1f%% should exceed noisy %.1f%%', ...
    comp_m.Fit(1), comp_m.Fit(2));

% Verify per-channel formula
for ch = 1:2
    ym = comp_m.Measured(:, ch);
    yp = comp_m.Predicted(:, ch);
    manual = 100 * (1 - norm(ym - yp) / norm(ym - mean(ym)));
    assert(abs(comp_m.Fit(ch) - manual) < 1e-8, ...
        'Channel %d: Fit=%.4f != manual=%.4f', ch, comp_m.Fit(ch), manual);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: per-channel NRMSE (%.1f%%, %.1f%%).\n', ...
    comp_m.Fit(1), comp_m.Fit(2));

%% Test 10: InitialState parameter (SPEC §15.4)
rng(7010);
p = 2; q = 1; N = 50; L = 1;
A_true = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
x0_true = [3; -2];  % known initial state

X = zeros(N+1, p);
U = randn(N, q);
X(1, :) = x0_true';
for k = 1:N
    X(k+1, :) = (A_true * X(k, :)' + B_true * U(k, :)')';
end

model.A = repmat(A_true, [1 1 N]);
model.B = repmat(B_true, [1 1 N]);
model.DataLength = N;
model.StateDim = p;
model.InputDim = q;
model.Method = 'sidLTVdisc';

% With correct initial state
comp_x0 = sidCompare(model, X, U, 'InitialState', x0_true);
assert(all(comp_x0.Fit > 99), ...
    'Correct x0 should give >99%% fit, got %.1f%%', min(comp_x0.Fit));

% With wrong initial state
comp_bad = sidCompare(model, X, U, 'InitialState', [0; 0]);
% Wrong x0 should give worse fit (transient error)
assert(all(comp_bad.Fit < comp_x0.Fit), ...
    'Wrong x0 should give worse fit');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: InitialState parameter (%.1f%% vs %.1f%%).\n', ...
    min(comp_x0.Fit), min(comp_bad.Fit));

fprintf('test_sidCompare: %d/%d passed\n', runner__nPassed, runner__nPassed);
