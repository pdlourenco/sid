%% test_sidFreqETFE - Unit tests for Empirical Transfer Function Estimate
%
% Tests sidFreqETFE for result structure, smoothing, time-series mode,
% and correctness against known systems.

fprintf('Running test_sidFreqETFE...\n');
runner__nPassed = 0;

%% Test 1: Result struct has all required fields
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1], [1 -0.8], u) + 0.1 * randn(N, 1);
result = sidFreqETFE(y, u);

requiredFields = {'Frequency', 'FrequencyHz', 'Response', 'ResponseStd', ...
    'NoiseSpectrum', 'NoiseSpectrumStd', 'Coherence', 'SampleTime', ...
    'WindowSize', 'DataLength', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: result struct has all required fields.\n');

%% Test 2: Method and metadata
assert(strcmp(result.Method, 'sidFreqETFE'), 'Method should be sidFreqETFE');
assert(result.DataLength == N, 'DataLength should be N');
assert(result.WindowSize == N, 'WindowSize should be N for ETFE');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: method and metadata.\n');

%% Test 3: SISO dimensions
nf = length(result.Frequency);
assert(isequal(size(result.Response), [nf, 1]), 'Response should be (nf x 1)');
assert(isempty(result.Coherence), 'ETFE Coherence should be empty');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: SISO dimensions.\n');

%% Test 4: ResponseStd is NaN (no asymptotic formula for ETFE)
assert(all(isnan(result.ResponseStd)), 'ETFE ResponseStd should be NaN');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: ResponseStd is NaN (no asymptotic formula for ETFE).\n');

%% Test 5: Noise spectrum is non-negative
assert(all(result.NoiseSpectrum >= 0), 'Noise spectrum should be non-negative');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: noise spectrum is non-negative.\n');

%% Test 6: Time series mode (periodogram)
rng(10);
y = randn(200, 1);
result_ts = sidFreqETFE(y, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(all(result_ts.NoiseSpectrum >= 0), 'Periodogram should be non-negative');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: time series mode (periodogram).\n');

%% Test 7: Smoothing parameter
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1], [1 -0.8], u);
result_s1 = sidFreqETFE(y, u, 'Smoothing', 1);
result_s11 = sidFreqETFE(y, u, 'Smoothing', 11);
% Smoothed version should have less variation
var_s1 = var(abs(result_s1.Response));
var_s11 = var(abs(result_s11.Response));
assert(var_s11 < var_s1, 'Smoothing should reduce variance of response');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: smoothing parameter.\n');

%% Test 8: Error on even smoothing parameter
try
    sidFreqETFE(y, u, 'Smoothing', 4);
    error('Should have thrown sid:badSmoothing');
catch e
    assert(strcmp(e.identifier, 'sid:badSmoothing'), 'Expected sid:badSmoothing error');
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: error on even smoothing parameter.\n');

%% Test 9: Error on non-integer smoothing
try
    sidFreqETFE(y, u, 'Smoothing', 3.5);
    error('Should have thrown sid:badSmoothing');
catch e
    assert(strcmp(e.identifier, 'sid:badSmoothing'), 'Expected sid:badSmoothing error');
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: error on non-integer smoothing.\n');

%% Test 10: Custom frequencies
w = linspace(0.1, pi, 50)';
result_cust = sidFreqETFE(y, u, 'Frequencies', w);
assert(length(result_cust.Frequency) == 50, 'Should have 50 custom frequencies');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: custom frequencies.\n');

%% Test 11: Known system - pure gain
% y = 2*u (no noise, no dynamics)
rng(5);
N = 1024;
u = randn(N, 1);
y = 2 * u;
result = sidFreqETFE(y, u);
% G should be approximately 2 at all frequencies
G_mag = abs(result.Response);
assert(max(abs(G_mag - 2)) < 0.01, 'ETFE of pure gain y=2u should give |G|~2');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: known system - pure gain.\n');

%% Test 12: Known system - pure delay
% y(t) = u(t-1), so G(z) = z^{-1}, |G| = 1, phase = -w
rng(5);
N = 1024;
u = randn(N, 1);
y = [0; u(1:end-1)];
result = sidFreqETFE(y, u);
G_mag = abs(result.Response);
G_phase = angle(result.Response);
expected_phase = -result.Frequency;
% Magnitude should be ~1
assert(median(abs(G_mag - 1)) < 0.05, 'ETFE of pure delay: |G| should be ~1');
% Phase should be ~-w (allow some tolerance)
phase_err = abs(G_phase - expected_phase);
% Wrap phase errors
phase_err = min(phase_err, 2*pi - phase_err);
assert(median(phase_err) < 0.1, 'ETFE of pure delay: phase should be ~-w');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 12 passed: known system - pure delay.\n');

%% Test 13: MIMO mode
rng(7);
N = 500;
u = randn(N, 1);
y = [filter([1], [1 -0.5], u), filter([0.3], [1 -0.7], u)];
result_mimo = sidFreqETFE(y, u);
nf = length(result_mimo.Frequency);
assert(size(result_mimo.Response, 1) == nf && size(result_mimo.Response, 2) == 2, ...
    'MIMO Response should be (nf x 2 x 1)');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 13 passed: MIMO mode.\n');

%% Test 14: Custom sample time
result = sidFreqETFE(y, u, 'SampleTime', 0.001);
assert(result.SampleTime == 0.001, 'SampleTime should be 0.001');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 14 passed: custom sample time.\n');

%% Test 15: Multi-trajectory — NumTrajectories and ensemble averaging
rng(15);
N15 = 2000; L15 = 8;
u15 = randn(N15, 1, L15);
y15 = zeros(N15, 1, L15);
for l = 1:L15
    y15(:, :, l) = filter(1, [1 -0.5], u15(:, :, l)) ...
        + 0.3 * randn(N15, 1);
end

res_mt = sidFreqETFE(y15, u15);
assert(res_mt.NumTrajectories == L15, ...
    'NumTrajectories should be %d, got %d', L15, res_mt.NumTrajectories);

res_st = sidFreqETFE(y15(:, :, 1), u15(:, :, 1));
assert(res_st.NumTrajectories == 1, 'Single traj should have NumTrajectories=1');

% Ensemble-averaged magnitude should be closer to truth
w = res_mt.Frequency;
G_true = 1 ./ (1 - 0.5 * exp(-1j * w));
err_mt = median(abs(abs(res_mt.Response) - abs(G_true)));
err_st = median(abs(abs(res_st.Response) - abs(G_true)));
assert(err_mt < err_st * 1.5, ...
    'Multi-traj error %.4f should improve on single %.4f', ...
    err_mt, err_st);

runner__nPassed = runner__nPassed + 1;
fprintf('  Test 15 passed: multi-trajectory ETFE (L=%d).\n', L15);

%% Test 16: Constant input -> all NaN + warning (SPEC 10.3)
% A constant signal's DFT is a nonzero Dirichlet kernel off DC, so the
% relative floor cannot see it; the whole-signal 10.3 check forces NaN.
rng(16);
N16 = 300;
y16 = filter(1, [1 -0.8], randn(N16, 1));
% Warnings ENABLED so lastwarn captures the id (portable across Octave releases).
lastwarn('');
r16 = sidFreqETFE(y16, ones(N16, 1));
[~, id16] = lastwarn();
assert(all(isnan(r16.Response)), 'Constant input -> Response all NaN');
assert(strcmp(id16, 'sid:constantInput'), 'Constant input should warn sid:constantInput');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 16 passed: constant input -> all NaN + warning (SPEC 10.3).\n');

%% Test 17: Collinear MIMO inputs warn (previously silent) and NaN the response
rng(17);
N17 = 400;
u0 = randn(N17, 1);
u17 = [u0, 2 * u0];
y17 = [u0 + 0.1 * randn(N17, 1), u0 + 0.1 * randn(N17, 1)];
lastwarn('');
r17 = sidFreqETFE(y17, u17);
[~, id17] = lastwarn();
assert(any(isnan(r17.Response(:))), 'Collinear MIMO -> some Response NaN');
assert(strcmp(id17, 'sid:singularPhiU'), 'Collinear MIMO should warn sid:singularPhiU');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 17 passed: collinear MIMO warns + NaN (SPEC 10.2).\n');

fprintf('test_sidFreqETFE: %d/%d passed\n', runner__nPassed, runner__nPassed);
