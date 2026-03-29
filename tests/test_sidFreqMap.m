%% test_sidFreqMap - Unit tests for time-varying frequency response map
%
% Tests sidFreqMap for result structure, segmentation, BT and Welch
% algorithms, time series mode, LTI constancy, and alignment with
% sidSpectrogram.

fprintf('Running test_sidFreqMap...\n');

%% Test 1: BT - Result struct has all required fields
rng(42);
N = 2000;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);
result = sidFreqMap(y, u);

requiredFields = {'Time', 'Frequency', 'FrequencyHz', 'Response', 'ResponseStd', ...
    'NoiseSpectrum', 'NoiseSpectrumStd', 'Coherence', 'SampleTime', ...
    'SegmentLength', 'Overlap', 'WindowSize', 'Algorithm', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
assert(strcmp(result.Method, 'sidFreqMap'), 'Method should be sidFreqMap');
assert(strcmp(result.Algorithm, 'bt'), 'Default Algorithm should be bt');
fprintf('  Test 1 passed: BT result struct has all fields\n');

%% Test 2: BT - Default parameters
L_exp = min(floor(N / 4), 256);
P_exp = floor(L_exp / 2);
M_exp = min(floor(L_exp / 10), 30);
assert(result.SegmentLength == L_exp, 'Default SegmentLength mismatch');
assert(result.Overlap == P_exp, 'Default Overlap mismatch');
assert(result.WindowSize == M_exp, 'Default WindowSize mismatch');
fprintf('  Test 2 passed: BT default parameters correct\n');

%% Test 3: BT - Segment count and dimensions
L = result.SegmentLength;
P = result.Overlap;
step = L - P;
K = floor((N - L) / step) + 1;
nf = length(result.Frequency);
assert(length(result.Time) == K, 'Time vector length should be K=%d', K);
assert(isequal(size(result.Response), [nf, K]), 'SISO Response should be (nf x K)');
assert(isequal(size(result.Coherence), [nf, K]), 'Coherence should be (nf x K)');
fprintf('  Test 3 passed: BT segment count and dimensions correct\n');

%% Test 4: BT - Time vector correctness
expectedTime = ((0:K-1)' * step + L / 2) * result.SampleTime;
assert(max(abs(result.Time - expectedTime)) < 1e-12, 'Time vector mismatch');
fprintf('  Test 4 passed: BT time vector correct\n');

%% Test 5: BT - Coherence in [0, 1]
assert(all(result.Coherence(:) >= -1e-10) && all(result.Coherence(:) <= 1 + 1e-10), ...
    'Coherence should be in [0, 1]');
fprintf('  Test 5 passed: BT coherence in valid range\n');

%% Test 6: BT - Noise spectrum non-negative
assert(all(result.NoiseSpectrum(:) >= -1e-10), 'Noise spectrum should be non-negative');
fprintf('  Test 6 passed: BT noise spectrum non-negative\n');

%% Test 7: BT - LTI system map approximately constant along time
magMap = abs(result.Response);
meanMag = mean(magMap, 2);
stdMag = std(magMap, 0, 2);
cv = stdMag ./ max(meanMag, eps);
assert(median(cv) < 0.5, 'LTI map should be roughly constant (median CV = %.2f)', median(cv));
fprintf('  Test 7 passed: BT LTI constancy check\n');

%% Test 8: BT - Time series mode
y_ts = randn(1000, 1);
result_ts = sidFreqMap(y_ts, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(isempty(result_ts.Coherence), 'Time series: Coherence should be empty');
assert(~isempty(result_ts.NoiseSpectrum), 'Time series: NoiseSpectrum should exist');
fprintf('  Test 8 passed: BT time series mode\n');

%% Test 9: BT - Custom parameters
L9 = 128; P9 = 64; M9 = 10;
result9 = sidFreqMap(y, u, 'SegmentLength', L9, 'Overlap', P9, 'WindowSize', M9);
assert(result9.SegmentLength == L9, 'Custom SegmentLength');
assert(result9.Overlap == P9, 'Custom Overlap');
assert(result9.WindowSize == M9, 'Custom WindowSize');
fprintf('  Test 9 passed: BT custom parameters\n');

%% Test 10: BT - Time axis alignment with sidSpectrogram
L10 = 200; P10 = 100; Ts10 = 0.01;
rng(7);
u10 = randn(N, 1);
y10 = filter([1], [1 -0.8], u10) + 0.1 * randn(N, 1);
mapResult = sidFreqMap(y10, u10, 'SegmentLength', L10, 'Overlap', P10, 'SampleTime', Ts10);
specResult = sidSpectrogram(y10, 'WindowLength', L10, 'Overlap', P10, 'SampleTime', Ts10);
assert(max(abs(mapResult.Time - specResult.Time)) < 1e-12, ...
    'sidFreqMap and sidSpectrogram time axes should align');
fprintf('  Test 10 passed: BT time axis alignment with sidSpectrogram\n');

%% Test 11: BT - Custom frequency vector
freqs11 = linspace(0.1, pi, 64)';
result11 = sidFreqMap(y, u, 'SegmentLength', 128, 'Frequencies', freqs11);
assert(length(result11.Frequency) == 64, 'Custom frequency vector length');
assert(max(abs(result11.Frequency - freqs11)) < 1e-12, 'Custom frequencies preserved');
fprintf('  Test 11 passed: BT custom frequency vector\n');

%% Test 12: BT - Error on segment too long
try
    sidFreqMap(randn(50, 1), randn(50, 1), 'SegmentLength', 100);
    error('Should have thrown an error');
catch e
    assert(strcmp(e.identifier, 'sid:segmentTooLong'), 'Expected sid:segmentTooLong');
end
fprintf('  Test 12 passed: error on segment too long\n');

%% Test 13: BT - Error on L <= 2*M
try
    sidFreqMap(y, u, 'SegmentLength', 20, 'WindowSize', 15);
    error('Should have thrown an error');
catch e
    assert(strcmp(e.identifier, 'sid:segmentTooShort'), 'Expected sid:segmentTooShort');
end
fprintf('  Test 13 passed: error on L <= 2*M\n');

%% Test 14: Welch - Runs without error and returns correct struct
rng(50);
N14 = 4000;
u14 = randn(N14, 1);
y14 = filter([1], [1 -0.9], u14) + 0.1 * randn(N14, 1);
result_w = sidFreqMap(y14, u14, 'Algorithm', 'welch', 'SegmentLength', 512);
assert(strcmp(result_w.Method, 'sidFreqMap'), 'Welch Method should be sidFreqMap');
assert(strcmp(result_w.Algorithm, 'welch'), 'Algorithm should be welch');
assert(isempty(result_w.WindowSize), 'Welch should have empty WindowSize');
assert(~isempty(result_w.Response), 'Welch should have Response');
assert(~isempty(result_w.Coherence), 'Welch SISO should have Coherence');
fprintf('  Test 14 passed: Welch runs and returns correct struct\n');

%% Test 15: Welch - Frequency grid is FFT bins
nf_w = length(result_w.Frequency);
Lsub_default = floor(512 / 4.5);
nfft_default = max(256, 2^nextpow2(Lsub_default));
expected_nf = floor(nfft_default / 2);
assert(nf_w == expected_nf, 'Welch should have %d frequency bins, got %d', expected_nf, nf_w);
assert(result_w.Frequency(1) > 0, 'Welch frequencies should skip DC');
assert(result_w.Frequency(end) <= pi + 1e-10, 'Welch frequencies should be <= pi');
fprintf('  Test 15 passed: Welch frequency grid correct\n');

%% Test 16: Welch - Coherence in [0, 1]
assert(all(result_w.Coherence(:) >= -1e-10) && all(result_w.Coherence(:) <= 1 + 1e-10), ...
    'Welch coherence should be in [0, 1]');
fprintf('  Test 16 passed: Welch coherence in valid range\n');

%% Test 17: Welch - Noise spectrum non-negative
assert(all(result_w.NoiseSpectrum(:) >= -1e-10), 'Welch noise spectrum non-negative');
fprintf('  Test 17 passed: Welch noise spectrum non-negative\n');

%% Test 18: Welch - Time series mode
result_wts = sidFreqMap(randn(2000, 1), [], 'Algorithm', 'welch', 'SegmentLength', 256);
assert(isempty(result_wts.Response), 'Welch time series: Response should be empty');
assert(~isempty(result_wts.NoiseSpectrum), 'Welch time series: NoiseSpectrum should exist');
fprintf('  Test 18 passed: Welch time series mode\n');

%% Test 19: Welch - Custom sub-segment parameters
result_wc = sidFreqMap(y14, u14, 'Algorithm', 'welch', 'SegmentLength', 512, ...
    'SubSegmentLength', 128, 'SubOverlap', 64, 'NFFT', 512, 'Window', 'hamming');
assert(length(result_wc.Frequency) == 256, 'Custom NFFT=512 should give 256 bins');
fprintf('  Test 19 passed: Welch custom sub-segment parameters\n');

%% Test 20: Welch - LTI system map roughly constant
magMap_w = abs(result_w.Response);
meanMag_w = mean(magMap_w, 2);
stdMag_w = std(magMap_w, 0, 2);
cv_w = stdMag_w ./ max(meanMag_w, eps);
assert(median(cv_w) < 0.5, 'Welch LTI map should be roughly constant (median CV = %.2f)', median(cv_w));
fprintf('  Test 20 passed: Welch LTI constancy check\n');

%% Test 21: Welch - Time axis alignment with sidSpectrogram
L21 = 256; P21 = 128; Ts21 = 0.01;
mapW = sidFreqMap(y14, u14, 'Algorithm', 'welch', 'SegmentLength', L21, ...
    'Overlap', P21, 'SampleTime', Ts21);
specW = sidSpectrogram(y14, 'WindowLength', L21, 'Overlap', P21, 'SampleTime', Ts21);
assert(max(abs(mapW.Time - specW.Time)) < 1e-12, ...
    'Welch sidFreqMap and sidSpectrogram time axes should align');
fprintf('  Test 21 passed: Welch time axis alignment\n');

%% Test 22: Invalid algorithm
try
    sidFreqMap(y, u, 'Algorithm', 'foobar');
    error('Should have thrown an error');
catch e
    assert(strcmp(e.identifier, 'sid:invalidAlgorithm'), 'Expected sid:invalidAlgorithm');
end
fprintf('  Test 22 passed: error on invalid algorithm\n');

%% Test 23: BT and Welch produce correlated results on same LTI data
rng(99);
N23 = 4000;
u23 = randn(N23, 1);
y23 = filter([1 0.5], [1 -0.8], u23) + 0.1 * randn(N23, 1);
res_bt = sidFreqMap(y23, u23, 'Algorithm', 'bt', 'SegmentLength', 512);
res_wl = sidFreqMap(y23, u23, 'Algorithm', 'welch', 'SegmentLength', 512);
% Average across time segments (LTI so should be roughly constant)
mag_bt = mean(abs(res_bt.Response), 2);
mag_wl = mean(abs(res_wl.Response), 2);
% Interpolate Welch onto BT frequency grid for comparison
mag_wl_interp = interp1(res_wl.Frequency, mag_wl, res_bt.Frequency, 'linear', 'extrap');
corr_mat = corrcoef(mag_bt, mag_wl_interp);
assert(corr_mat(1,2) > 0.8, 'BT and Welch magnitude shapes should correlate (r=%.2f)', corr_mat(1,2));
fprintf('  Test 23 passed: BT and Welch correlated (r=%.4f)\n', corr_mat(1,2));

fprintf('test_sidFreqMap: ALL TESTS PASSED\n');
