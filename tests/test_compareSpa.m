% test_compareSpa.m - Compare sidFreqBT against MATLAB's spa function
%
% Requires: System Identification Toolbox
%
% Verifies that sidFreqBT produces results matching MATLAB's spa for
% various system configurations: SISO, MIMO, time-series, custom
% parameters, and different noise levels.

fprintf('Running test_compareSpa...\n');

%% Toolbox check
if ~exist('spa', 'file')
    fprintf('  test_compareSpa: SKIPPED (requires System Identification Toolbox)\n');
    return;
end

%% Test 1: SISO first-order system, default window size
rng(10);
N = 2000;
Ts = 1;
a = 0.85;
sigma_v = 0.5;  % moderate noise to avoid catastrophic cancellation in PhiV
u = randn(N, 1);
y = filter(1, [1 -a], u) + sigma_v * randn(N, 1);

result_sid = sidFreqBT(y, u);
w = result_sid.Frequency;
M = result_sid.WindowSize;

data = iddata(y, u, Ts);
G_spa = spa(data, M, w);
resp_spa = squeeze(G_spa.ResponseData);
spec_spa = real(squeeze(G_spa.SpectrumData));

% Magnitude comparison
relErr_mag = max(abs(abs(result_sid.Response) - abs(resp_spa)) ./ abs(resp_spa));
assert(relErr_mag < 0.01, ...
    'Test 1: magnitude relErr=%.6f should be <1%%', relErr_mag);

% Phase comparison
phaseErr = max(abs(angle(result_sid.Response) - angle(resp_spa)));
assert(phaseErr < 0.01, ...
    'Test 1: max phase error=%.6f rad should be <0.01', phaseErr);

% Noise spectrum comparison (looser tolerance: PhiV = PhiY - |PhiYU|^2/PhiU
% involves subtraction of similar magnitudes, amplifying small differences
% in window normalization between sid and spa)
relErr_noise = max(abs(real(result_sid.NoiseSpectrum) - spec_spa) ./ max(abs(spec_spa), 1e-10));
assert(relErr_noise < 0.10, ...
    'Test 1: noise spectrum relErr=%.6f should be <10%%', relErr_noise);

%% Test 2: SISO second-order system, custom window size
rng(20);
N = 5000;
M = 60;
u = randn(N, 1);
% Second-order system: poles at 0.9*exp(±j*pi/4)
b = 1;
a_coeff = [1, -2*0.9*cos(pi/4), 0.9^2];
y = filter(b, a_coeff, u) + 0.05 * randn(N, 1);

result_sid = sidFreqBT(y, u, 'WindowSize', M);
w = result_sid.Frequency;

data = iddata(y, u, Ts);
G_spa = spa(data, M, w);
resp_spa = squeeze(G_spa.ResponseData);

relErr = max(abs(result_sid.Response - resp_spa) ./ max(abs(resp_spa), 1e-10));
assert(relErr < 0.01, ...
    'Test 2: complex response relErr=%.6f should be <1%%', relErr);

%% Test 3: Custom frequency vector (non-default grid)
rng(30);
N = 1000;
M = 25;
u = randn(N, 1);
y = filter(1, [1 -0.7], u) + 0.2 * randn(N, 1);

w_custom = logspace(log10(0.05), log10(pi), 64)';

result_sid = sidFreqBT(y, u, 'WindowSize', M, 'Frequencies', w_custom);

data = iddata(y, u, Ts);
G_spa = spa(data, M, w_custom);
resp_spa = squeeze(G_spa.ResponseData);

relErr = max(abs(result_sid.Response - resp_spa) ./ max(abs(resp_spa), 1e-10));
assert(relErr < 0.01, ...
    'Test 3: custom freq response relErr=%.6f should be <1%%', relErr);

%% Test 4: Time-series (output spectrum only)
rng(40);
N = 2000;
M = 30;
% Colored noise: filter white noise through AR(1)
e = randn(N, 1);
y = filter(1, [1 -0.6], e);

w = (1:128)' * pi / 128;
result_sid = sidFreqBT(y, [], 'WindowSize', M);

data = iddata(y, [], Ts);
G_spa = spa(data, M, w);
spec_spa = real(squeeze(G_spa.SpectrumData));

relErr = max(abs(real(result_sid.NoiseSpectrum) - spec_spa) ./ max(abs(spec_spa), 1e-10));
assert(relErr < 0.01, ...
    'Test 4: time-series spectrum relErr=%.6f should be <1%%', relErr);

assert(isempty(result_sid.Response), 'Test 4: time-series Response should be empty');

%% Test 5: SISO with significant noise - check coherence and noise spectrum
rng(50);
N = 3000;
M = 40;
sigma_v = 0.5;
u = randn(N, 1);
y = filter([1 0.3], [1 -0.8], u) + sigma_v * randn(N, 1);

result_sid = sidFreqBT(y, u, 'WindowSize', M);
w = result_sid.Frequency;

data = iddata(y, u, Ts);
G_spa = spa(data, M, w);
resp_spa = squeeze(G_spa.ResponseData);
spec_spa = real(squeeze(G_spa.SpectrumData));

% Response comparison
relErr_resp = max(abs(result_sid.Response - resp_spa) ./ max(abs(resp_spa), 1e-10));
assert(relErr_resp < 0.01, ...
    'Test 5: noisy response relErr=%.6f should be <1%%', relErr_resp);

% Noise spectrum comparison (looser tolerance for implementation differences)
relErr_noise = max(abs(real(result_sid.NoiseSpectrum) - spec_spa) ./ max(abs(spec_spa), 1e-10));
assert(relErr_noise < 0.10, ...
    'Test 5: noise spectrum relErr=%.6f should be <10%%', relErr_noise);

%% Test 6: MIMO system (2 outputs, 1 input)
rng(60);
N = 3000;
M = 30;
u = randn(N, 1);
y1 = filter(1, [1 -0.5], u) + 0.1 * randn(N, 1);
y2 = filter([0.3], [1 -0.7], u) + 0.1 * randn(N, 1);
y = [y1, y2];

result_sid = sidFreqBT(y, u, 'WindowSize', M);
w = result_sid.Frequency;

data = iddata(y, u, Ts);
G_spa = spa(data, M, w);
resp_spa = permute(G_spa.ResponseData, [3, 1, 2]);  % (ny x nu x nf) -> (nf x ny x nu)
spec_spa = real(permute(G_spa.SpectrumData, [3, 1, 2]));

% Response comparison (handle trailing singleton collapse)
resp_sid = result_sid.Response;
nf = length(w);
for ch = 1:2
    resp_sid_ch = resp_sid(:, ch);
    resp_spa_ch = resp_spa(:, ch);
    relErr = max(abs(resp_sid_ch - resp_spa_ch) ./ max(abs(resp_spa_ch), 1e-10));
    assert(relErr < 0.01, ...
        'Test 6: MIMO channel %d response relErr=%.6f should be <1%%', ch, relErr);
end

% Noise spectrum comparison (2x2 per frequency)
noise_sid = result_sid.NoiseSpectrum;
for ii = 1:2
    for jj = 1:2
        ns_sid = real(squeeze(noise_sid(:, ii, jj)));
        ns_spa = spec_spa(:, ii, jj);
        relErr = max(abs(ns_sid - ns_spa) ./ max(abs(ns_spa), 1e-10));
        assert(relErr < 0.10, ...
            'Test 6: MIMO noise(%d,%d) relErr=%.6f should be <10%%', ii, jj, relErr);
    end
end

%% Test 7: Non-unit sample time
rng(70);
N = 2000;
Ts = 0.01;
M = 30;
u = randn(N, 1);
y = filter(1, [1 -0.9], u) + 0.1 * randn(N, 1);

% sid uses frequencies in rad/sample; convert to rad/s for spa
w_sid = (1:128)' * pi / 128;           % rad/sample
w_spa = w_sid / Ts;                     % rad/s (= rad/TimeUnit with Ts)

result_sid = sidFreqBT(y, u, 'WindowSize', M, 'SampleTime', Ts);

data = iddata(y, u, Ts);
G_spa = spa(data, M, w_spa);
resp_spa = squeeze(G_spa.ResponseData);

relErr = max(abs(result_sid.Response - resp_spa) ./ max(abs(resp_spa), 1e-10));
assert(relErr < 0.01, ...
    'Test 7: non-unit Ts response relErr=%.6f should be <1%%', relErr);

% Verify FrequencyHz matches
assert(max(abs(result_sid.FrequencyHz - G_spa.Frequency / (2*pi))) < 1e-10, ...
    'Test 7: FrequencyHz should match spa frequency in Hz');

fprintf('  test_compareSpa: ALL PASSED\n');
