%% test_sidFreqBTFDR - Unit tests for Blackman-Tukey with frequency-dependent resolution
%
% Tests sidFreqBTFDR for result structure, resolution parameter behavior,
% time-series mode, and SISO/MIMO support.

fprintf('Running test_sidFreqBTFDR...\n');
runner__nPassed = 0;

%% Test 1: Result struct has all required fields
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1], [1 -0.8], u) + 0.1 * randn(N, 1);
result = sidFreqBTFDR(y, u);

requiredFields = {'Frequency', 'FrequencyHz', 'Response', 'ResponseStd', ...
    'NoiseSpectrum', 'NoiseSpectrumStd', 'Coherence', 'SampleTime', ...
    'WindowSize', 'DataLength', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: result struct has all required fields.\n');

%% Test 2: Method identifier
assert(strcmp(result.Method, 'sidFreqBTFDR'), 'Method should be sidFreqBTFDR');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: method identifier.\n');

%% Test 3: WindowSize is a vector (per-frequency)
nf = length(result.Frequency);
assert(isequal(size(result.WindowSize), [nf, 1]), 'WindowSize should be (nf x 1) vector');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: WindowSize is a vector (per-frequency).\n');

%% Test 4: Default resolution gives uniform window size
% Default R = 2*pi/M_default, so all M_k should be the same
Mk = result.WindowSize;
assert(all(Mk == Mk(1)), 'Default resolution should give uniform window size');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: default resolution gives uniform window size.\n');

%% Test 5: Custom scalar resolution
result_fine = sidFreqBTFDR(y, u, 'Resolution', 0.1);
result_coarse = sidFreqBTFDR(y, u, 'Resolution', 1.0);
% Finer resolution => larger window
assert(result_fine.WindowSize(1) > result_coarse.WindowSize(1), ...
    'Finer resolution should use larger window');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: custom scalar resolution.\n');

%% Test 6: Per-frequency resolution vector
nf = 128;
R_vec = linspace(0.2, 2.0, nf)';
result_vec = sidFreqBTFDR(y, u, 'Resolution', R_vec);
Mk = result_vec.WindowSize;
% Window size should generally decrease as resolution increases
assert(Mk(1) >= Mk(end), 'Larger resolution should give smaller window');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: per-frequency resolution vector.\n');

%% Test 7: Error on negative resolution
try
    sidFreqBTFDR(y, u, 'Resolution', -0.5);
    error('Should have thrown sid:badResolution');
catch e
    assert(strcmp(e.identifier, 'sid:badResolution'), 'Expected sid:badResolution error');
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: error on negative resolution.\n');

%% Test 8: Error on mismatched resolution vector length
try
    sidFreqBTFDR(y, u, 'Resolution', [0.1; 0.2; 0.3]);
    error('Should have thrown sid:badResolution');
catch e
    assert(strcmp(e.identifier, 'sid:badResolution'), 'Expected sid:badResolution error');
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: error on mismatched resolution vector length.\n');

%% Test 9: Time series mode
y_ts = randn(300, 1);
result_ts = sidFreqBTFDR(y_ts, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(isempty(result_ts.ResponseStd), 'Time series: ResponseStd should be empty');
assert(isempty(result_ts.Coherence), 'Time series: Coherence should be empty');
assert(length(result_ts.NoiseSpectrum) == 128, 'Time series should have 128 spectral values');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: time series mode.\n');

%% Test 10: SISO coherence is in [0, 1]
result = sidFreqBTFDR(y, u);
assert(all(result.Coherence >= 0) && all(result.Coherence <= 1), ...
    'Coherence should be in [0, 1]');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: SISO coherence is in [0, 1].\n');

%% Test 11: Noise spectrum is non-negative (SISO)
assert(all(result.NoiseSpectrum >= 0), 'Noise spectrum should be non-negative');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: noise spectrum is non-negative (SISO).\n');

%% Test 12: MIMO mode
rng(7);
N = 500;
u = randn(N, 2);
y = [u(:,1) + 0.5*u(:,2), 0.3*u(:,1) + u(:,2)] + 0.1*randn(N, 2);
result_mimo = sidFreqBTFDR(y, u);
nf = length(result_mimo.Frequency);
assert(isequal(size(result_mimo.Response), [nf, 2, 2]), 'MIMO Response should be (nf x 2 x 2)');
assert(isempty(result_mimo.Coherence), 'MIMO Coherence should be empty');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 12 passed: MIMO mode.\n');

%% Test 13: Custom frequencies
w = linspace(0.1, pi, 50)';
result = sidFreqBTFDR(y(:,1), u(:,1), 'Frequencies', w);
assert(length(result.Frequency) == 50, 'Should have 50 frequencies');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 13 passed: custom frequencies.\n');

%% Test 14: Window sizes are clamped to [2, N/2]
N = 20;
y = randn(N, 1);
u = randn(N, 1);
result = sidFreqBTFDR(y, u, 'Resolution', 0.01);  % Very fine => large M
assert(all(result.WindowSize <= floor(N/2)), 'M should be clamped to N/2');
assert(all(result.WindowSize >= 2), 'M should be at least 2');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 14 passed: window sizes are clamped to [2, N/2].\n');

%% Test 15: Known first-order system
rng(1);
N = 5000;
u = randn(N, 1);
y = filter(1, [1 -0.9], u);
result = sidFreqBTFDR(y, u, 'Resolution', 0.2);
w = result.Frequency;
G_true = 1 ./ (1 - 0.9 * exp(-1j * w));
idx = [10, 30, 60, 100];
for i = 1:length(idx)
    k = idx(i);
    relErr = abs(abs(result.Response(k)) - abs(G_true(k))) / abs(G_true(k));
    assert(relErr < 0.25, 'Magnitude at freq %d should match (relErr=%.3f)', k, relErr);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 15 passed: known first-order system.\n');

%% Test 16: Multi-trajectory — NumTrajectories and variance reduction
rng(16);
N16 = 2000; L16 = 8;
u16 = randn(N16, 1, L16);
y16 = zeros(N16, 1, L16);
for l = 1:L16
    y16(:, :, l) = filter(1, [1 -0.9], u16(:, :, l)) ...
        + 0.3 * randn(N16, 1);
end

res_mt = sidFreqBTFDR(y16, u16, 'Resolution', 0.2);
assert(res_mt.NumTrajectories == L16, ...
    'NumTrajectories should be %d, got %d', L16, res_mt.NumTrajectories);

res_st = sidFreqBTFDR(y16(:, :, 1), u16(:, :, 1), 'Resolution', 0.2);
assert(res_st.NumTrajectories == 1, 'Single traj should have NumTrajectories=1');

ratio = median(res_mt.ResponseStd) / median(res_st.ResponseStd);
expected = 1 / sqrt(L16);
assert(ratio < expected * 2.5 && ratio > expected * 0.2, ...
    'Variance reduction ratio %.3f should be near %.3f', ...
    ratio, expected);

runner__nPassed = runner__nPassed + 1;
fprintf('  Test 16 passed: multi-trajectory (L=%d, ratio=%.3f).\n', ...
    L16, ratio);

%% Test 17: Collinear MIMO inputs -> NaN + Inf sigma, no crash (#141)
rng(2);
N17 = 400;
u0 = randn(N17, 1);
u17 = [u0, 2 * u0];                     % rank-deficient Phi_u
y17 = [u0 + 0.1 * randn(N17, 1), u0 + 0.1 * randn(N17, 1)];
ws = warning('off', 'all');
r17 = sidFreqBTFDR(y17, u17, 'Resolution', 0.3);
warning(ws);
assert(all(isnan(r17.Response(:))), 'Collinear MIMO -> Response all NaN');
assert(all(isinf(r17.ResponseStd(:))), 'Collinear MIMO -> sigma_G all Inf');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 17 passed: collinear MIMO -> NaN + Inf sigma (#141).\n');

%% Test 18: Constant input -> all NaN with warning (SPEC 10.3)
rng(3);
N18 = 300;
y18 = filter(1, [1 -0.8], randn(N18, 1));
% Warnings ENABLED so lastwarn captures the id on CI's Octave 8.4.
lastwarn('');
r18 = sidFreqBTFDR(y18, ones(N18, 1), 'Resolution', 0.3);
[~, id18] = lastwarn();
assert(all(isnan(r18.Response)), 'Constant input -> Response all NaN');
assert(strcmp(id18, 'sid:constantInput'), 'Constant input should warn sid:constantInput');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 18 passed: constant input -> all NaN + warning.\n');

%% Test 19: Resolution too coarse (M_k < 2) is an error, not a silent clamp
rng(4);
N19 = 400;
threw19 = false;
try
    sidFreqBTFDR(randn(N19, 1), randn(N19, 1), 'Resolution', 10);
catch err19
    threw19 = strcmp(err19.identifier, 'sid:badResolution');
end
assert(threw19, 'Coarse resolution (M_k < 2) should raise sid:badResolution');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 19 passed: coarse resolution errors (SPEC 5.2).\n');

%% Test 20: Oracle -- BTFDR at constant resolution equals BT with M = ceil(2*pi/R)
rng(6);
N20 = 800;
u20 = randn(N20, 1);
y20 = filter(1, [1 -0.8], u20) + 0.05 * randn(N20, 1);
w20 = linspace(0.05, pi, 100)';
M20 = 20;
rBT  = sidFreqBT(y20, u20, 'WindowSize', M20, 'Frequencies', w20);
rFDR = sidFreqBTFDR(y20, u20, 'Resolution', 2 * pi / M20, 'Frequencies', w20);
assert(all(rFDR.WindowSize == M20), 'Constant resolution should give uniform M');
assert(max(abs(rFDR.Response - rBT.Response)) < 1e-8, 'BTFDR Response must match BT');
assert(max(abs(rFDR.NoiseSpectrum - rBT.NoiseSpectrum)) < 1e-8, ...
    'BTFDR NoiseSpectrum must match BT');
assert(max(abs(rFDR.Coherence - rBT.Coherence)) < 1e-8, 'BTFDR Coherence must match BT');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 20 passed: BTFDR == BT at constant resolution (oracle).\n');

fprintf('test_sidFreqBTFDR: %d/%d passed\n', runner__nPassed, runner__nPassed);
