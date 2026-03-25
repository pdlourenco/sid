% test_compareSpafdr.m - Compare sidFreqBTFDR against MATLAB's spafdr function
%
% Requires: System Identification Toolbox
%
% Verifies that sidFreqBTFDR produces results matching MATLAB's spafdr for
% various resolution settings, system configurations, and MIMO.

fprintf('Running test_compareSpafdr...\n');

%% Toolbox check
if ~exist('spafdr', 'file')
    fprintf('  test_compareSpafdr: SKIPPED (requires System Identification Toolbox)\n');
    return;
end

%% Test 1: SISO first-order system, default resolution
rng(21);
N = 2000;
Ts = 1;
u = randn(N, 1);
y = filter(1, [1 -0.85], u) + 0.1 * randn(N, 1);

w = (1:128)' * pi / 128;

result_sid = sidFreqBTFDR(y, u, 'Frequencies', w);

data = iddata(y, u, Ts);
G_spafdr = spafdr(data, [], w);
resp_spafdr = squeeze(G_spafdr.ResponseData);
spec_spafdr = real(squeeze(G_spafdr.SpectrumData));

% Response comparison
relErr = max(abs(result_sid.Response - resp_spafdr) ./ max(abs(resp_spafdr), 1e-10));
assert(relErr < 0.02, ...
    'Test 1: default response relErr=%.6f should be <2%%', relErr);

% Noise spectrum comparison
relErr_noise = max(abs(real(result_sid.NoiseSpectrum) - spec_spafdr) ./ max(abs(spec_spafdr), 1e-10));
assert(relErr_noise < 0.02, ...
    'Test 1: default noise spectrum relErr=%.6f should be <2%%', relErr_noise);

%% Test 2: Fine resolution (large window, low variance)
rng(22);
N = 5000;
u = randn(N, 1);
y = filter(1, [1 -0.9], u) + 0.05 * randn(N, 1);

R_fine = 0.1;  % rad/sample — corresponds to M = ceil(2*pi/0.1) = 63
w = (1:128)' * pi / 128;

result_sid = sidFreqBTFDR(y, u, 'Resolution', R_fine, 'Frequencies', w);

data = iddata(y, u, Ts);
G_spafdr = spafdr(data, R_fine, w);
resp_spafdr = squeeze(G_spafdr.ResponseData);

relErr = max(abs(result_sid.Response - resp_spafdr) ./ max(abs(resp_spafdr), 1e-10));
assert(relErr < 0.02, ...
    'Test 2: fine resolution response relErr=%.6f should be <2%%', relErr);

%% Test 3: Coarse resolution (small window, high smoothing)
rng(23);
N = 2000;
u = randn(N, 1);
y = filter([1 0.5], [1 -0.8], u) + 0.2 * randn(N, 1);

R_coarse = 1.0;  % rad/sample — corresponds to M = ceil(2*pi/1) = 7
w = (1:128)' * pi / 128;

result_sid = sidFreqBTFDR(y, u, 'Resolution', R_coarse, 'Frequencies', w);

data = iddata(y, u, Ts);
G_spafdr = spafdr(data, R_coarse, w);
resp_spafdr = squeeze(G_spafdr.ResponseData);

relErr = max(abs(result_sid.Response - resp_spafdr) ./ max(abs(resp_spafdr), 1e-10));
assert(relErr < 0.02, ...
    'Test 3: coarse resolution response relErr=%.6f should be <2%%', relErr);

%% Test 4: Time-series spectrum
rng(24);
N = 2000;
e = randn(N, 1);
y = filter(1, [1 -0.6 0.3], e);

w = (1:128)' * pi / 128;

result_sid = sidFreqBTFDR(y, [], 'Frequencies', w);

data = iddata(y, [], Ts);
G_spafdr = spafdr(data, [], w);
spec_spafdr = real(squeeze(G_spafdr.SpectrumData));

relErr = max(abs(real(result_sid.NoiseSpectrum) - spec_spafdr) ./ max(abs(spec_spafdr), 1e-10));
assert(relErr < 0.02, ...
    'Test 4: time-series spectrum relErr=%.6f should be <2%%', relErr);

%% Test 5: Per-frequency resolution vector
rng(25);
N = 3000;
u = randn(N, 1);
y = filter(1, [1 -0.85], u) + 0.1 * randn(N, 1);

w = (1:64)' * pi / 64;
% Fine resolution at low frequencies, coarse at high
R_vec = linspace(0.1, 1.0, 64)';

result_sid = sidFreqBTFDR(y, u, 'Resolution', R_vec, 'Frequencies', w);

data = iddata(y, u, Ts);
G_spafdr = spafdr(data, R_vec, w);
resp_spafdr = squeeze(G_spafdr.ResponseData);

relErr = max(abs(result_sid.Response - resp_spafdr) ./ max(abs(resp_spafdr), 1e-10));
assert(relErr < 0.02, ...
    'Test 5: per-freq resolution response relErr=%.6f should be <2%%', relErr);

%% Test 6: MIMO system (2 outputs, 2 inputs)
rng(26);
N = 3000;
u = randn(N, 2);
y1 = filter(1, [1 -0.5], u(:,1)) + filter(0.3, [1 -0.3], u(:,2)) + 0.1*randn(N,1);
y2 = filter(0.5, [1 -0.7], u(:,1)) + filter(1, [1 -0.4], u(:,2)) + 0.1*randn(N,1);
y = [y1, y2];

w = (1:64)' * pi / 64;
R_res = 0.3;

result_sid = sidFreqBTFDR(y, u, 'Resolution', R_res, 'Frequencies', w);

data = iddata(y, u, Ts);
G_spafdr = spafdr(data, R_res, w);
resp_spafdr = permute(G_spafdr.ResponseData, [3, 1, 2]);  % (ny x nu x nf) -> (nf x ny x nu)

for ii = 1:2
    for jj = 1:2
        resp_sid_ij = squeeze(result_sid.Response(:, ii, jj));
        resp_spafdr_ij = squeeze(resp_spafdr(:, ii, jj));
        relErr = max(abs(resp_sid_ij - resp_spafdr_ij) ./ max(abs(resp_spafdr_ij), 1e-10));
        assert(relErr < 0.02, ...
            'Test 6: MIMO(%d,%d) response relErr=%.6f should be <2%%', ii, jj, relErr);
    end
end

%% Test 7: Non-unit sample time
rng(27);
N = 2000;
Ts = 0.001;
u = randn(N, 1);
y = filter(1, [1 -0.9], u) + 0.1 * randn(N, 1);

w_sid = (1:64)' * pi / 64;   % rad/sample
w_spa = w_sid / Ts;            % rad/s

R_res = 0.5;  % rad/sample for sid
R_spa = R_res / Ts;  % rad/s for spafdr

result_sid = sidFreqBTFDR(y, u, 'Resolution', R_res, 'Frequencies', w_sid, 'SampleTime', Ts);

data = iddata(y, u, Ts);
G_spafdr = spafdr(data, R_spa, w_spa);
resp_spafdr = squeeze(G_spafdr.ResponseData);

relErr = max(abs(result_sid.Response - resp_spafdr) ./ max(abs(resp_spafdr), 1e-10));
assert(relErr < 0.02, ...
    'Test 7: non-unit Ts response relErr=%.6f should be <2%%', relErr);

fprintf('  test_compareSpafdr: ALL PASSED\n');
