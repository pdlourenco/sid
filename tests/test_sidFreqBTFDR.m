%% test_sidFreqBTFDR - Unit tests for Blackman-Tukey with frequency-dependent resolution
%
% Tests sidFreqBTFDR for result structure, resolution parameter behavior,
% time-series mode, and SISO/MIMO support.

fprintf('Running test_sidFreqBTFDR...\n');

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

%% Test 2: Method identifier
assert(strcmp(result.Method, 'sidFreqBTFDR'), 'Method should be sidFreqBTFDR');

%% Test 3: WindowSize is a vector (per-frequency)
nf = length(result.Frequency);
assert(isequal(size(result.WindowSize), [nf, 1]), 'WindowSize should be (nf x 1) vector');

%% Test 4: Default resolution gives uniform window size
% Default R = 2*pi/M_default, so all M_k should be the same
Mk = result.WindowSize;
assert(all(Mk == Mk(1)), 'Default resolution should give uniform window size');

%% Test 5: Custom scalar resolution
result_fine = sidFreqBTFDR(y, u, 'Resolution', 0.1);
result_coarse = sidFreqBTFDR(y, u, 'Resolution', 1.0);
% Finer resolution => larger window
assert(result_fine.WindowSize(1) > result_coarse.WindowSize(1), ...
    'Finer resolution should use larger window');

%% Test 6: Per-frequency resolution vector
nf = 128;
R_vec = linspace(0.2, 2.0, nf)';
result_vec = sidFreqBTFDR(y, u, 'Resolution', R_vec);
Mk = result_vec.WindowSize;
% Window size should generally decrease as resolution increases
assert(Mk(1) >= Mk(end), 'Larger resolution should give smaller window');

%% Test 7: Error on negative resolution
try
    sidFreqBTFDR(y, u, 'Resolution', -0.5);
    error('Should have thrown sid:badResolution');
catch e
    assert(strcmp(e.identifier, 'sid:badResolution'), 'Expected sid:badResolution error');
end

%% Test 8: Error on mismatched resolution vector length
try
    sidFreqBTFDR(y, u, 'Resolution', [0.1; 0.2; 0.3]);
    error('Should have thrown sid:badResolution');
catch e
    assert(strcmp(e.identifier, 'sid:badResolution'), 'Expected sid:badResolution error');
end

%% Test 9: Time series mode
y_ts = randn(300, 1);
result_ts = sidFreqBTFDR(y_ts, []);
assert(isempty(result_ts.Response), 'Time series: Response should be empty');
assert(isempty(result_ts.ResponseStd), 'Time series: ResponseStd should be empty');
assert(isempty(result_ts.Coherence), 'Time series: Coherence should be empty');
assert(length(result_ts.NoiseSpectrum) == 128, 'Time series should have 128 spectral values');

%% Test 10: SISO coherence is in [0, 1]
result = sidFreqBTFDR(y, u);
assert(all(result.Coherence >= 0) && all(result.Coherence <= 1), ...
    'Coherence should be in [0, 1]');

%% Test 11: Noise spectrum is non-negative (SISO)
assert(all(result.NoiseSpectrum >= 0), 'Noise spectrum should be non-negative');

%% Test 12: MIMO mode
rng(7);
N = 500;
u = randn(N, 2);
y = [u(:,1) + 0.5*u(:,2), 0.3*u(:,1) + u(:,2)] + 0.1*randn(N, 2);
result_mimo = sidFreqBTFDR(y, u);
nf = length(result_mimo.Frequency);
assert(isequal(size(result_mimo.Response), [nf, 2, 2]), 'MIMO Response should be (nf x 2 x 2)');
assert(isempty(result_mimo.Coherence), 'MIMO Coherence should be empty');

%% Test 13: Custom frequencies
w = linspace(0.1, pi, 50)';
result = sidFreqBTFDR(y(:,1), u(:,1), 'Frequencies', w);
assert(length(result.Frequency) == 50, 'Should have 50 frequencies');

%% Test 14: Window sizes are clamped to [2, N/2]
N = 20;
y = randn(N, 1);
u = randn(N, 1);
result = sidFreqBTFDR(y, u, 'Resolution', 0.01);  % Very fine => large M
assert(all(result.WindowSize <= floor(N/2)), 'M should be clamped to N/2');
assert(all(result.WindowSize >= 2), 'M should be at least 2');

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
    assert(relErr < 0.15, 'Magnitude at freq %d should match (relErr=%.3f)', k, relErr);
end

fprintf('  test_sidFreqBTFDR: ALL PASSED\n');
