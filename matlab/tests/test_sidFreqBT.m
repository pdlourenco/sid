%% test_sidFreqBT - Unit tests for Blackman-Tukey spectral analysis
%
% Tests sidFreqBT for result structure, SISO/MIMO/time-series modes,
% and correctness of spectral estimates.

fprintf('Running test_sidFreqBT...\n');
runner__nPassed = 0;

%% Test 1: Result struct has all required fields
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1 0.5], [1 -0.8], u) + 0.1 * randn(N, 1);
result = sidFreqBT(y, u);

requiredFields = {'Frequency', 'FrequencyHz', 'Response', 'ResponseStd', ...
    'NoiseSpectrum', 'NoiseSpectrumStd', 'Coherence', 'SampleTime', ...
    'WindowSize', 'DataLength', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: result struct has all required fields.\n');

%% Test 2: Correct metadata
assert(result.DataLength == N, 'DataLength should be N');
assert(strcmp(result.Method, 'sidFreqBT'), 'Method should be sidFreqBT');
assert(result.SampleTime == 1.0, 'Default SampleTime should be 1.0');
assert(length(result.Frequency) == 128, 'Default should be 128 frequencies');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: correct metadata.\n');

%% Test 3: SISO dimensions
nf = length(result.Frequency);
assert(isequal(size(result.Response), [nf, 1]), 'SISO Response should be (nf x 1)');
assert(isequal(size(result.ResponseStd), [nf, 1]), 'SISO ResponseStd should be (nf x 1)');
assert(isequal(size(result.NoiseSpectrum), [nf, 1]), 'SISO NoiseSpectrum should be (nf x 1)');
assert(isequal(size(result.Coherence), [nf, 1]), 'SISO Coherence should be (nf x 1)');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: SISO dimensions.\n');

%% Test 4: Coherence is in [0, 1]
assert(all(result.Coherence >= 0) && all(result.Coherence <= 1), ...
    'Coherence should be in [0, 1]');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: coherence is in [0, 1].\n');

%% Test 5: Noise spectrum is non-negative
assert(all(result.NoiseSpectrum >= 0), 'Noise spectrum should be non-negative');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: noise spectrum is non-negative.\n');

%% Test 6: Time series mode
y_ts = randn(300, 1);
result_ts = sidFreqBT(y_ts, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(isempty(result_ts.ResponseStd), 'Time series: ResponseStd should be empty');
assert(isempty(result_ts.Coherence), 'Time series: Coherence should be empty');
assert(all(result_ts.NoiseSpectrum >= -1e-10), 'Time series spectrum should be non-negative');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: time series mode.\n');

%% Test 7: Custom window size
result_cust = sidFreqBT(y, u, 'WindowSize', 20);
assert(result_cust.WindowSize == 20, 'Custom WindowSize should be 20');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: custom window size.\n');

%% Test 8: Custom frequencies
w = linspace(0.1, pi, 50)';
result_cust = sidFreqBT(y, u, 'WindowSize', 20, 'Frequencies', w);
assert(length(result_cust.Frequency) == 50, 'Should have 50 custom frequencies');
assert(max(abs(result_cust.Frequency - w)) < 1e-12, 'Frequencies should match input');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: custom frequencies.\n');

%% Test 9: Custom sample time affects FrequencyHz
result_st = sidFreqBT(y, u, 'SampleTime', 0.01);
expected_hz = result_st.Frequency / (2 * pi * 0.01);
assert(max(abs(result_st.FrequencyHz - expected_hz)) < 1e-12, ...
    'FrequencyHz should be freq / (2*pi*Ts)');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: custom sample time affects FrequencyHz.\n');

%% Test 10: Positional syntax works
result_pos = sidFreqBT(y, u, 15);
assert(result_pos.WindowSize == 15, 'Positional M should work');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: positional syntax works.\n');

%% Test 11: MIMO mode (2 outputs, 1 input)
rng(7);
N = 500;
u = randn(N, 1);
y = [filter([1], [1 -0.8], u), filter([0.5], [1 -0.5], u)] + 0.1*randn(N, 2);
result_mimo = sidFreqBT(y, u);
nf = length(result_mimo.Frequency);
assert(size(result_mimo.Response, 1) == nf && size(result_mimo.Response, 2) == 2, ...
    'MIMO Response size should be (nf x 2)');
assert(isempty(result_mimo.Coherence), 'MIMO Coherence should be empty');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: MIMO mode (2 outputs, 1 input).\n');

%% Test 12: Uncertainty std is finite and non-negative for SISO
rng(42);
u = randn(500, 1);
y = filter([1], [1 -0.8], u) + 0.1 * randn(500, 1);
result = sidFreqBT(y, u);
assert(all(isfinite(result.ResponseStd)), 'SISO ResponseStd should be finite');
assert(all(result.ResponseStd >= 0), 'SISO ResponseStd should be non-negative');
assert(all(isfinite(result.NoiseSpectrumStd)), 'NoiseSpectrumStd should be finite');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 12 passed: uncertainty std is finite and non-negative for SISO.\n');

%% Test 13: Known first-order system identification
% y(t) = 0.9*y(t-1) + u(t) => G(z) = 1/(1-0.9*z^{-1})
% |G(e^{jw})| = 1/|1 - 0.9*e^{-jw}|
rng(1);
N = 5000;
u = randn(N, 1);
y = filter(1, [1 -0.9], u);
result = sidFreqBT(y, u, 'WindowSize', 50);
w = result.Frequency;
G_true = 1 ./ (1 - 0.9 * exp(-1j * w));
% Check magnitude match at a few frequencies
idx = [10, 30, 60, 100];
for i = 1:length(idx)
    k = idx(i);
    relErr = abs(abs(result.Response(k)) - abs(G_true(k))) / abs(G_true(k));
    assert(relErr < 0.15, 'Magnitude at freq %d should match true system (relErr=%.3f)', k, relErr);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 13 passed: known first-order system identification.\n');

%% Test 14: Multi-trajectory — NumTrajectories field and variance reduction
rng(14);
N14 = 2000; L14 = 8;
u14 = randn(N14, 1, L14);
y14 = zeros(N14, 1, L14);
for l = 1:L14
    y14(:, :, l) = filter(1, [1 -0.9], u14(:, :, l)) ...
        + 0.3 * randn(N14, 1);
end

% Multi-trajectory result
res_mt = sidFreqBT(y14, u14, 'WindowSize', 50);
assert(res_mt.NumTrajectories == L14, ...
    'NumTrajectories should be %d, got %d', L14, res_mt.NumTrajectories);

% Single-trajectory result
res_st = sidFreqBT(y14(:, :, 1), u14(:, :, 1), 'WindowSize', 50);
assert(res_st.NumTrajectories == 1, 'Single traj should have NumTrajectories=1');

% Variance reduction: multi-traj std should be roughly std_1 / sqrt(L)
ratio = median(res_mt.ResponseStd) / median(res_st.ResponseStd);
expected = 1 / sqrt(L14);
assert(ratio < expected * 2 && ratio > expected * 0.3, ...
    'Variance reduction ratio %.3f should be near 1/sqrt(%d)=%.3f', ...
    ratio, L14, expected);

% Ensemble-averaged estimate should be closer to truth than single
w = res_mt.Frequency;
G_true = 1 ./ (1 - 0.9 * exp(-1j * w));
err_mt = median(abs(abs(res_mt.Response) - abs(G_true)));
err_st = median(abs(abs(res_st.Response) - abs(G_true)));
assert(err_mt < err_st, ...
    'Multi-traj error %.4f should be less than single %.4f', ...
    err_mt, err_st);

runner__nPassed = runner__nPassed + 1;
fprintf('  Test 14 passed: multi-trajectory (L=%d, ratio=%.3f).\n', ...
    L14, ratio);

%% Test 15: WindowSize M < 2 should error (SPEC §2.2)
rng(15);
N = 100;
u = randn(N, 1);
y = filter(1, [1 -0.5], u);
passed = false;
try
    sidFreqBT(y, u, 'WindowSize', 1);
catch e
    if ~isempty(strfind(e.identifier, 'sid:'))
        passed = true;
    end
end
assert(passed, 'WindowSize=1 should raise an error');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 15 passed: M < 2 rejected.\n');

%% Test 16: Near-zero input — should not crash
rng(16);
N = 500;
u = 1e-15 * randn(N, 1);  % near-zero input
y = randn(N, 1);
result = sidFreqBT(y, u, 'WindowSize', 20);
% Response may be huge/NaN but should not crash
assert(isfield(result, 'Response'), 'Should return a result struct');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 16 passed: near-zero input handled.\n');

%% Test 17: Large window size (M in the old FFT-path bug envelope)
% Regression for the sidWindowedDFT L = 2*nf hardcode (SPEC.md §2.5.1).
% With the default 128-point grid, M in [128, 255] silently produced
% the wrong spectrum (positive/negative lag overlap) and M >= 256
% silently truncated (MATLAB auto-extended s past L, then fft(s, L)
% dropped the extra entries). This end-to-end test exercises sidFreqBT
% on the default grid across that range and cross-checks against the
% direct DFT path (forced by offsetting the frequency vector slightly).
for M_i = [128, 200, 256, 500]
    rng(123);
    N_i = 2 * M_i + 100;  % ensures sidFreqBT does not clamp M
    u_i = randn(N_i, 1);
    y_i = filter([1 0.5], [1 -0.8], u_i) + 0.05 * randn(N_i, 1);

    res_fft = sidFreqBT(y_i, u_i, 'WindowSize', M_i);
    assert(res_fft.WindowSize == M_i, 'WindowSize should be %d', M_i);
    assert(isequal(size(res_fft.Response), [128, 1]), ...
        'Default grid gives (128 x 1) Response');
    assert(all(isfinite(res_fft.Response)), ...
        'M=%d: Response should be finite', M_i);
    assert(all(isfinite(res_fft.NoiseSpectrum)), ...
        'M=%d: NoiseSpectrum should be finite', M_i);
    assert(all(res_fft.NoiseSpectrum >= -1e-10), ...
        'M=%d: NoiseSpectrum should be non-negative', M_i);

    % Cross-check against the direct DFT path. sidIsDefaultFreqs returns
    % true for vectors that match the default grid to within 1e-12, so
    % we offset by -1e-11 to escape that check and force the direct
    % path. The shift is negative so the last element stays inside the
    % required (0, pi] range. DO NOT "clean this up" — the offset is
    % what makes the cross-check exercise both code paths.
    w_direct = (1:128)' * pi / 128 - 1e-11;
    res_dir = sidFreqBT(y_i, u_i, 'WindowSize', M_i, 'Frequencies', w_direct);
    relErr = max(abs(res_fft.Response - res_dir.Response)) / ...
        max(abs(res_dir.Response));
    assert(relErr < 1e-8, ...
        'M=%d: FFT path vs direct path relErr=%.2e', M_i, relErr);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 17 passed: large WindowSize (FFT-path regression).\n');

fprintf('test_sidFreqBT: %d/%d passed\n', runner__nPassed, runner__nPassed);
