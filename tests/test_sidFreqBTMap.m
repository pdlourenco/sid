%% test_sidFreqBTMap - Unit tests for time-varying frequency response map
%
% Tests sidFreqBTMap for result structure, segmentation, time series mode,
% LTI constancy, and alignment with sidSpectrogram.

fprintf('Running test_sidFreqBTMap...\n');

%% Test 1: Result struct has all required fields
rng(42);
N = 2000;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);
result = sidFreqBTMap(y, u);

requiredFields = {'Time', 'Frequency', 'FrequencyHz', 'Response', 'ResponseStd', ...
    'NoiseSpectrum', 'NoiseSpectrumStd', 'Coherence', 'SampleTime', ...
    'SegmentLength', 'Overlap', 'WindowSize', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
assert(strcmp(result.Method, 'sidFreqBTMap'), 'Method should be sidFreqBTMap');
fprintf('  Test 1 passed: result struct has all fields\n');

%% Test 2: Default parameters
L_exp = min(floor(N / 4), 256);
P_exp = floor(L_exp / 2);
M_exp = min(floor(L_exp / 10), 30);
assert(result.SegmentLength == L_exp, 'Default SegmentLength mismatch');
assert(result.Overlap == P_exp, 'Default Overlap mismatch');
assert(result.WindowSize == M_exp, 'Default WindowSize mismatch');
fprintf('  Test 2 passed: default parameters correct\n');

%% Test 3: Segment count and dimensions
L = result.SegmentLength;
P = result.Overlap;
step = L - P;
K = floor((N - L) / step) + 1;
nf = length(result.Frequency);
assert(length(result.Time) == K, 'Time vector length should be K=%d', K);
assert(isequal(size(result.Response), [nf, K]), 'SISO Response should be (nf x K)');
assert(isequal(size(result.ResponseStd), [nf, K]), 'SISO ResponseStd should be (nf x K)');
assert(isequal(size(result.NoiseSpectrum), [nf, K]), 'NoiseSpectrum should be (nf x K)');
assert(isequal(size(result.Coherence), [nf, K]), 'Coherence should be (nf x K)');
fprintf('  Test 3 passed: segment count and dimensions correct\n');

%% Test 4: Time vector correctness
expectedTime = ((0:K-1)' * step + L / 2) * result.SampleTime;
assert(max(abs(result.Time - expectedTime)) < 1e-12, 'Time vector mismatch');
fprintf('  Test 4 passed: time vector correct\n');

%% Test 5: Coherence in [0, 1]
assert(all(result.Coherence(:) >= -1e-10) && all(result.Coherence(:) <= 1 + 1e-10), ...
    'Coherence should be in [0, 1]');
fprintf('  Test 5 passed: coherence in valid range\n');

%% Test 6: Noise spectrum is non-negative
assert(all(result.NoiseSpectrum(:) >= -1e-10), 'Noise spectrum should be non-negative');
fprintf('  Test 6 passed: noise spectrum non-negative\n');

%% Test 7: LTI system - map should be approximately constant along time
% For an LTI system, the response should not vary much across segments
magMap = abs(result.Response);
meanMag = mean(magMap, 2);
stdMag = std(magMap, 0, 2);
% Coefficient of variation should be small for most frequencies
cv = stdMag ./ max(meanMag, eps);
% Allow some frequencies to be noisy, but median CV should be < 0.5
assert(median(cv) < 0.5, 'LTI map should be roughly constant (median CV = %.2f)', median(cv));
fprintf('  Test 7 passed: LTI constancy check\n');

%% Test 8: Time series mode
y_ts = randn(1000, 1);
result_ts = sidFreqBTMap(y_ts, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(isempty(result_ts.ResponseStd), 'Time series: ResponseStd should be empty');
assert(isempty(result_ts.Coherence), 'Time series: Coherence should be empty');
assert(~isempty(result_ts.NoiseSpectrum), 'Time series: NoiseSpectrum should exist');
fprintf('  Test 8 passed: time series mode\n');

%% Test 9: Custom parameters
L9 = 128; P9 = 64; M9 = 10;
result9 = sidFreqBTMap(y, u, 'SegmentLength', L9, 'Overlap', P9, 'WindowSize', M9);
assert(result9.SegmentLength == L9, 'Custom SegmentLength');
assert(result9.Overlap == P9, 'Custom Overlap');
assert(result9.WindowSize == M9, 'Custom WindowSize');
fprintf('  Test 9 passed: custom parameters\n');

%% Test 10: Custom sample time
Ts10 = 0.001;
result10 = sidFreqBTMap(y, u, 'SampleTime', Ts10, 'SegmentLength', 128);
assert(result10.SampleTime == Ts10, 'Custom SampleTime');
fprintf('  Test 10 passed: custom sample time\n');

%% Test 11: Time axis alignment with sidSpectrogram
L11 = 200; P11 = 100; Ts11 = 0.01;
N11 = 2000;
rng(7);
u11 = randn(N11, 1);
y11 = filter([1], [1 -0.8], u11) + 0.1 * randn(N11, 1);
mapResult = sidFreqBTMap(y11, u11, 'SegmentLength', L11, 'Overlap', P11, 'SampleTime', Ts11);
specResult = sidSpectrogram(y11, 'WindowLength', L11, 'Overlap', P11, 'SampleTime', Ts11);
assert(max(abs(mapResult.Time - specResult.Time)) < 1e-12, ...
    'sidFreqBTMap and sidSpectrogram time axes should align');
fprintf('  Test 11 passed: time axis alignment with sidSpectrogram\n');

%% Test 12: Error on segment too long
try
    sidFreqBTMap(randn(50, 1), randn(50, 1), 'SegmentLength', 100);
    error('Should have thrown an error');
catch e
    assert(strcmp(e.identifier, 'sid:segmentTooLong'), 'Expected sid:segmentTooLong');
end
fprintf('  Test 12 passed: error on segment too long\n');

%% Test 13: Error on L <= 2*M
try
    sidFreqBTMap(y, u, 'SegmentLength', 20, 'WindowSize', 15);
    error('Should have thrown an error');
catch e
    assert(strcmp(e.identifier, 'sid:segmentTooShort'), 'Expected sid:segmentTooShort');
end
fprintf('  Test 13 passed: error on L <= 2*M\n');

%% Test 14: Custom frequency vector
freqs14 = linspace(0.1, pi, 64)';
result14 = sidFreqBTMap(y, u, 'SegmentLength', 128, 'Frequencies', freqs14);
assert(length(result14.Frequency) == 64, 'Custom frequency vector length');
assert(max(abs(result14.Frequency - freqs14)) < 1e-12, 'Custom frequencies preserved');
fprintf('  Test 14 passed: custom frequency vector\n');

fprintf('test_sidFreqBTMap: ALL TESTS PASSED\n');
