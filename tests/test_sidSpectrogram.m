%% test_sidSpectrogram - Unit tests for sidSpectrogram
%
% Tests sidSpectrogram for result structure, correctness, edge cases,
% and window types.

fprintf('Running test_sidSpectrogram...\n');

%% Test 1: Result struct has all required fields
rng(42);
N = 1000;
x = randn(N, 1);
result = sidSpectrogram(x);

requiredFields = {'Time', 'Frequency', 'FrequencyRad', 'Power', 'PowerDB', ...
    'Complex', 'SampleTime', 'WindowLength', 'Overlap', 'NFFT', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
assert(strcmp(result.Method, 'sidSpectrogram'), 'Method should be sidSpectrogram');
fprintf('  Test 1 passed: result struct has all fields\n');

%% Test 2: Default parameters
assert(result.WindowLength == 256, 'Default WindowLength should be 256');
assert(result.Overlap == 128, 'Default Overlap should be floor(256/2)');
assert(result.NFFT == 256, 'Default NFFT should be max(256, 2^nextpow2(256))');
assert(result.SampleTime == 1.0, 'Default SampleTime should be 1.0');
fprintf('  Test 2 passed: default parameters correct\n');

%% Test 3: Output dimensions
L = 256; P = 128; nfft = 256;
step = L - P;
K = floor((N - L) / step) + 1;
nBins = floor(nfft / 2) + 1;
assert(length(result.Time) == K, 'Time vector length should be K=%d', K);
assert(length(result.Frequency) == nBins, 'Frequency vector length should be %d', nBins);
assert(isequal(size(result.Power), [nBins, K]), 'Power dimensions wrong');
assert(isequal(size(result.PowerDB), [nBins, K]), 'PowerDB dimensions wrong');
assert(isequal(size(result.Complex), [nBins, K]), 'Complex dimensions wrong');
fprintf('  Test 3 passed: output dimensions correct\n');

%% Test 4: Frequency vector starts at 0 Hz
assert(result.Frequency(1) == 0, 'First frequency should be 0 Hz (DC)');
fprintf('  Test 4 passed: frequency vector starts at DC\n');

%% Test 5: Power is non-negative
assert(all(result.Power(:) >= 0), 'Power should be non-negative');
fprintf('  Test 5 passed: power is non-negative\n');

%% Test 6: PowerDB consistent with Power
expectedDB = 10 * log10(max(result.Power, eps));
assert(max(abs(result.PowerDB(:) - expectedDB(:))) < 1e-10, ...
    'PowerDB should be 10*log10(Power)');
fprintf('  Test 6 passed: PowerDB consistent with Power\n');

%% Test 7: Known sinusoid - peak at correct frequency
Fs = 1000; Ts = 1/Fs;
N2 = 4000;
f0 = 100;  % Hz
t = (0:N2-1)' * Ts;
x2 = sin(2 * pi * f0 * t);
L2 = 256;
result2 = sidSpectrogram(x2, 'WindowLength', L2, 'SampleTime', Ts);

% For each segment, the peak should be near f0
for k = 1:length(result2.Time)
    [~, peakIdx] = max(result2.Power(:, k, 1));
    peakFreq = result2.Frequency(peakIdx);
    assert(abs(peakFreq - f0) < Fs / L2, ...
        'Peak at segment %d should be near %g Hz, got %g Hz', k, f0, peakFreq);
end
fprintf('  Test 7 passed: sinusoid peak at correct frequency\n');

%% Test 8: Custom window size and overlap
result3 = sidSpectrogram(randn(500, 1), 'WindowLength', 64, 'Overlap', 32, ...
    'NFFT', 128);
assert(result3.WindowLength == 64, 'Custom WindowLength');
assert(result3.Overlap == 32, 'Custom Overlap');
assert(result3.NFFT == 128, 'Custom NFFT');
nBins3 = floor(128 / 2) + 1;
assert(size(result3.Power, 1) == nBins3, 'Frequency bins should match NFFT');
fprintf('  Test 8 passed: custom parameters\n');

%% Test 9: Multi-channel data
x_mc = randn(1000, 3);
result_mc = sidSpectrogram(x_mc, 'WindowLength', 128);
assert(size(result_mc.Power, 3) == 3, 'Multi-channel: 3 channels expected');
assert(size(result_mc.Complex, 3) == 3, 'Multi-channel: 3 channels in Complex');
fprintf('  Test 9 passed: multi-channel data\n');

%% Test 10: Hamming window
result_ham = sidSpectrogram(randn(500, 1), 'WindowLength', 64, 'Window', 'hamming');
assert(all(result_ham.Power(:) >= 0), 'Hamming: power should be non-negative');
fprintf('  Test 10 passed: hamming window\n');

%% Test 11: Rectangular window
result_rect = sidSpectrogram(randn(500, 1), 'WindowLength', 64, 'Window', 'rect');
assert(all(result_rect.Power(:) >= 0), 'Rect: power should be non-negative');
fprintf('  Test 11 passed: rectangular window\n');

%% Test 12: Custom window vector
w = ones(64, 1) * 0.5;
result_cust = sidSpectrogram(randn(500, 1), 'WindowLength', 64, 'Window', w);
assert(all(result_cust.Power(:) >= 0), 'Custom window: power should be non-negative');
fprintf('  Test 12 passed: custom window vector\n');

%% Test 13: Time vector correctness
L4 = 100; P4 = 50; Ts4 = 0.01;
N4 = 500;
result4 = sidSpectrogram(randn(N4, 1), 'WindowLength', L4, 'Overlap', P4, ...
    'SampleTime', Ts4);
step4 = L4 - P4;
K4 = floor((N4 - L4) / step4) + 1;
expectedTime = ((0:K4-1)' * step4 + L4 / 2) * Ts4;
assert(max(abs(result4.Time - expectedTime)) < 1e-12, 'Time vector mismatch');
fprintf('  Test 13 passed: time vector correct\n');

%% Test 14: Error on short data
try
    sidSpectrogram(randn(10, 1), 'WindowLength', 256);
    error('Should have thrown an error for short data');
catch e
    assert(strcmp(e.identifier, 'sid:tooShort'), 'Expected sid:tooShort error');
end
fprintf('  Test 14 passed: error on short data\n');

%% Test 15: Error on invalid overlap
try
    sidSpectrogram(randn(500, 1), 'WindowLength', 64, 'Overlap', 64);
    error('Should have thrown an error for P >= L');
catch e
    assert(strcmp(e.identifier, 'sid:invalidOverlap'), 'Expected sid:invalidOverlap error');
end
fprintf('  Test 15 passed: error on invalid overlap\n');

%% Test 16: FrequencyRad consistent with Frequency
assert(max(abs(result.FrequencyRad - 2 * pi * result.Frequency)) < 1e-12, ...
    'FrequencyRad should be 2*pi*Frequency');
fprintf('  Test 16 passed: FrequencyRad consistent\n');

fprintf('test_sidSpectrogram: ALL TESTS PASSED\n');
