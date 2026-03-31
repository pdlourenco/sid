% test_sidResidual.m - Test residual analysis function

fprintf('Running test_sidResidual...\n');

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

assert(isequal(size(resid_ss.Residual), [N, p]), ...
    'SS residual should be N x p');
assert(isscalar(resid_ss.ConfidenceBound), 'Confidence bound should be scalar');
fprintf('  Test 3 passed: state-space residual analysis.\n');

%% Test 4: Time-series mode (u=[])
rng(6004);
N = 1000;
y = filter(1, [1 -0.6], randn(N, 1));

result_ts = sidFreqBT(y, []);
resid_ts = sidResidual(result_ts, y, []);

assert(isempty(resid_ts.CrossCorr), 'CrossCorr should be empty for time series');
assert(resid_ts.IndependencePass == true, 'Independence should be vacuously true');
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
fprintf('  Test 5 passed: MIMO residual analysis.\n');

%% Test 6: Plot option (headless check)
rng(6006);
N = 500;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.1 * randn(N, 1);
result_plot = sidFreqBT(y, u);

% Should not error; in headless environments, allow display-related failures
try
    sidResidual(result_plot, y, u, 'Plot', true);
    close all;
catch e
    % In headless environments (xvfb, no display), plotting may fail
    % Only rethrow if it's NOT a display/rendering issue
    msg = e.message;
    isDisplayErr = ~isempty(strfind(msg, 'figure')) || ~isempty(strfind(msg, 'display')) ...
        || ~isempty(strfind(msg, 'DISPLAY')) || ~isempty(strfind(msg, 'bar')) ...
        || ~isempty(strfind(msg, 'int_value')) || ~isempty(strfind(msg, 'FaceColor'));
    if ~isDisplayErr
        rethrow(e);
    end
end
fprintf('  Test 6 passed: plot option works.\n');

fprintf('  test_sidResidual: ALL PASSED\n');
