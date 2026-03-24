%% test_sidFreqBT - Unit tests for Blackman-Tukey spectral analysis
%
% Tests sidFreqBT for result structure, SISO/MIMO/time-series modes,
% and correctness of spectral estimates.

fprintf('Running test_sidFreqBT...\n');

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

%% Test 2: Correct metadata
assert(result.DataLength == N, 'DataLength should be N');
assert(strcmp(result.Method, 'sidFreqBT'), 'Method should be sidFreqBT');
assert(result.SampleTime == 1.0, 'Default SampleTime should be 1.0');
assert(length(result.Frequency) == 128, 'Default should be 128 frequencies');

%% Test 3: SISO dimensions
nf = length(result.Frequency);
assert(isequal(size(result.Response), [nf, 1]), 'SISO Response should be (nf x 1)');
assert(isequal(size(result.ResponseStd), [nf, 1]), 'SISO ResponseStd should be (nf x 1)');
assert(isequal(size(result.NoiseSpectrum), [nf, 1]), 'SISO NoiseSpectrum should be (nf x 1)');
assert(isequal(size(result.Coherence), [nf, 1]), 'SISO Coherence should be (nf x 1)');

%% Test 4: Coherence is in [0, 1]
assert(all(result.Coherence >= 0) && all(result.Coherence <= 1), ...
    'Coherence should be in [0, 1]');

%% Test 5: Noise spectrum is non-negative
assert(all(result.NoiseSpectrum >= 0), 'Noise spectrum should be non-negative');

%% Test 6: Time series mode
y_ts = randn(300, 1);
result_ts = sidFreqBT(y_ts, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(isempty(result_ts.ResponseStd), 'Time series: ResponseStd should be empty');
assert(isempty(result_ts.Coherence), 'Time series: Coherence should be empty');
assert(all(result_ts.NoiseSpectrum >= -1e-10), 'Time series spectrum should be non-negative');

%% Test 7: Custom window size
result_cust = sidFreqBT(y, u, 'WindowSize', 20);
assert(result_cust.WindowSize == 20, 'Custom WindowSize should be 20');

%% Test 8: Custom frequencies
w = linspace(0.1, pi, 50)';
result_cust = sidFreqBT(y, u, 'WindowSize', 20, 'Frequencies', w);
assert(length(result_cust.Frequency) == 50, 'Should have 50 custom frequencies');
assert(max(abs(result_cust.Frequency - w)) < 1e-12, 'Frequencies should match input');

%% Test 9: Custom sample time affects FrequencyHz
result_st = sidFreqBT(y, u, 'SampleTime', 0.01);
expected_hz = result_st.Frequency / (2 * pi * 0.01);
assert(max(abs(result_st.FrequencyHz - expected_hz)) < 1e-12, ...
    'FrequencyHz should be freq / (2*pi*Ts)');

%% Test 10: Positional syntax works
result_pos = sidFreqBT(y, u, 15);
assert(result_pos.WindowSize == 15, 'Positional M should work');

%% Test 11: MIMO mode (2 outputs, 1 input)
rng(7);
N = 500;
u = randn(N, 1);
y = [filter([1], [1 -0.8], u), filter([0.5], [1 -0.5], u)] + 0.1*randn(N, 2);
result_mimo = sidFreqBT(y, u);
nf = length(result_mimo.Frequency);
assert(isequal(size(result_mimo.Response), [nf, 2, 1]), 'MIMO Response size should be (nf x ny x nu)');
assert(isempty(result_mimo.Coherence), 'MIMO Coherence should be empty');

%% Test 12: Uncertainty std is finite and non-negative for SISO
rng(42);
u = randn(500, 1);
y = filter([1], [1 -0.8], u) + 0.1 * randn(500, 1);
result = sidFreqBT(y, u);
assert(all(isfinite(result.ResponseStd)), 'SISO ResponseStd should be finite');
assert(all(result.ResponseStd >= 0), 'SISO ResponseStd should be non-negative');
assert(all(isfinite(result.NoiseSpectrumStd)), 'NoiseSpectrumStd should be finite');

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
    assert(relErr < 0.1, 'Magnitude at freq %d should match true system (relErr=%.3f)', k, relErr);
end

fprintf('  test_sidFreqBT: ALL PASSED\n');
