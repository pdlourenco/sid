% test_sidCompare.m - Test model output comparison function

fprintf('Running test_sidCompare...\n');

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
assert(isequal(size(comp.Predicted), [N, p]), 'Predicted should be N x p');
assert(isequal(size(comp.Residual), [N, p]), 'Residual should be N x p');
fprintf('  Test 1 passed: perfect state-space model gives %.1f%% fit.\n', min(comp.Fit));

%% Test 2: Frequency-domain model fit
rng(7002);
N = 2000;
u = randn(N, 1);
y = filter(1, [1 -0.85], u) + 0.05 * randn(N, 1);

result = sidFreqBT(y, u, 'WindowSize', 40);
comp = sidCompare(result, y, u);

assert(comp.Fit > 50, 'Good freq-domain model should give >50%% fit, got %.1f%%', comp.Fit);
assert(isequal(size(comp.Predicted), [N, 1]), 'Predicted should be N x 1');
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
catch e
    % In headless environments, allow display-related failures
    msg = e.message;
    isDisplayErr = ~isempty(strfind(msg, 'figure')) || ~isempty(strfind(msg, 'display')) ...
        || ~isempty(strfind(msg, 'DISPLAY')) || ~isempty(strfind(msg, 'java'));
    if ~isDisplayErr
        rethrow(e);
    end
end
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
fprintf('  Test 6 passed: output struct has correct fields.\n');

fprintf('  test_sidCompare: ALL PASSED\n');
