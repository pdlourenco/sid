% test_sidResidual.m - Test residual analysis function

fprintf('Running test_sidResidual...\n');
runner__nPassed = 0;

%% Test 1: Good model — white residuals
rng(6001);
N = 2000;
u = randn(N, 1);
y = filter(1, [1 -0.85], u) + 0.1 * randn(N, 1);

result = sidFreqBT(y, u, 'WindowSize', 40);
resid = sidResidual(result, y, u);

assert(isfield(resid, 'Residual'), 'Should have Residual field');
assert(isfield(resid, 'WhitenessPass'), 'Should have WhitenessPass field');
assert(isfield(resid, 'IndependencePass'), 'Should have IndependencePass field');
assert(isequal(size(resid.Residual), [N, 1]), 'Residual should be N x 1');
assert(resid.ConfidenceBound > 0, 'Confidence bound should be positive');
% With a good model and enough data, whiteness test should usually pass
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: residual struct valid (whiteness=%d, independence=%d).\n', ...
    resid.WhitenessPass, resid.IndependencePass);

%% Test 2: Bad model — coloured residuals
rng(6002);
N = 2000;
u = randn(N, 1);
% True system: second order
y = filter([1 0.5], [1 -0.9 0.4], u) + 0.05 * randn(N, 1);

% Estimate with very coarse model (tiny window → biased)
result_bad = sidFreqBT(y, u, 'WindowSize', 2);
resid_bad = sidResidual(result_bad, y, u);

% With a deliberately poor model, residuals should show structure
assert(~isempty(resid_bad.AutoCorr), 'AutoCorr should be computed');
assert(length(resid_bad.AutoCorr) > 1, 'AutoCorr should have multiple lags');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: bad model residual analysis works.\n');

%% Test 3: State-space path (sidLTVdisc)
rng(6003);
p = 2; q = 1; N = 50; L = 5;
A_true = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end

ltv = sidLTVdisc(X, U, 'Lambda', 1e4);
resid_ss = sidResidual(ltv, X, U);

% Multi-trajectory residual is per-trajectory (N x p x L), and the whiteness
% bound uses the pooled sample count N_eff = L*N (issue #140).
assert(isequal(size(resid_ss.Residual), [N, p, L]), ...
    'SS multi-traj residual should be N x p x L');
assert(isscalar(resid_ss.ConfidenceBound), 'Confidence bound should be scalar');
assert(abs(resid_ss.ConfidenceBound - 2.58 / sqrt(L * N)) < 1e-12, ...
    'Multi-traj bound should be 2.58/sqrt(L*N), got %.5f', resid_ss.ConfidenceBound);
assert(resid_ss.DataLength == L * N, ...
    'Multi-traj DataLength should be L*N = %d, got %d', L * N, resid_ss.DataLength);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: state-space multi-trajectory residual (#140).\n');

%% Test 4: Time-series mode (u=[])
rng(6004);
N = 1000;
y = filter(1, [1 -0.6], randn(N, 1));

result_ts = sidFreqBT(y, []);
resid_ts = sidResidual(result_ts, y, []);

assert(isempty(resid_ts.CrossCorr), 'CrossCorr should be empty for time series');
assert(resid_ts.IndependencePass == true, 'Independence should be vacuously true');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: time-series residual analysis.\n');

%% Test 5: Multi-channel (MIMO)
rng(6005);
N = 2000;
u = randn(N, 1);
y1 = filter(1, [1 -0.5], u) + 0.1 * randn(N, 1);
y2 = filter(0.3, [1 -0.7], u) + 0.1 * randn(N, 1);
y = [y1, y2];

result_mimo = sidFreqBT(y, u, 'WindowSize', 30);
resid_mimo = sidResidual(result_mimo, y, u);

assert(size(resid_mimo.Residual, 2) == 2, 'MIMO residual should have 2 columns');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: MIMO residual analysis.\n');

%% Test 6: Plot option (headless check)
rng(6006);
N = 500;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.1 * randn(N, 1);
result_plot = sidFreqBT(y, u);

% Should not error even in headless environment
try
    sidResidual(result_plot, y, u, 'Plot', true);
    close all;
    plotOk = true;
catch
    plotOk = false;
end
assert(plotOk, 'Plot option should not error');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: plot option works.\n');

%% Test 7: Confidence bound formula verification (SPEC §14.3)
% ConfidenceBound should be 2.58/sqrt(N)
rng(6007);
N = 1600;
u = randn(N, 1);
y = filter(1, [1 -0.7], u) + 0.1 * randn(N, 1);
result_cb = sidFreqBT(y, u, 'WindowSize', 30);
resid_cb = sidResidual(result_cb, y, u);

expected_bound = 2.58 / sqrt(N);
assert(abs(resid_cb.ConfidenceBound - expected_bound) < 1e-10, ...
    'ConfidenceBound should be 2.58/sqrt(N)=%.6f, got %.6f', ...
    expected_bound, resid_cb.ConfidenceBound);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: confidence bound = 2.58/sqrt(N).\n');

%% Test 8: Normalised autocorrelation formula (SPEC §14.2)
% r_ee(0) should be exactly 1 (normalised by R_ee(0))
rng(6008);
N = 1000;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.1 * randn(N, 1);
result_ac = sidFreqBT(y, u, 'WindowSize', 30);
resid_ac = sidResidual(result_ac, y, u);

assert(abs(resid_ac.AutoCorr(1) - 1.0) < 1e-10, ...
    'r_ee(0) should be 1.0, got %.6f', resid_ac.AutoCorr(1));

% WhitenessPass should be: all |r_ee(tau)| < 2.58/sqrt(N) for tau >= 1
bound = resid_ac.ConfidenceBound;
autoLags = resid_ac.AutoCorr(2:end);  % exclude lag 0
manual_pass = all(abs(autoLags) < bound);
assert(resid_ac.WhitenessPass == manual_pass, ...
    'WhitenessPass should match manual check');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: autocorrelation normalised, r_ee(0)=1.\n');

%% Test 9: Cross-correlation normalisation (SPEC §14.2)
% r_eu(tau) = R_eu(tau) / sqrt(R_ee(0) * R_uu(0))
% IndependencePass: all |r_eu(tau)| < 2.58/sqrt(N)
rng(6009);
N = 1000;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.1 * randn(N, 1);
result_cc = sidFreqBT(y, u, 'WindowSize', 30);
resid_cc = sidResidual(result_cc, y, u);

assert(~isempty(resid_cc.CrossCorr), 'CrossCorr should not be empty');
bound_cc = resid_cc.ConfidenceBound;
manual_indep = all(abs(resid_cc.CrossCorr) < bound_cc);
assert(resid_cc.IndependencePass == manual_indep, ...
    'IndependencePass should match manual check');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: cross-correlation independence check.\n');

fprintf('test_sidResidual: %d/%d passed\n', runner__nPassed, runner__nPassed);
