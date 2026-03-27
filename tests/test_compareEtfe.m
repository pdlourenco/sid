% test_compareEtfe.m - Compare sidFreqETFE against MATLAB's etfe function
%
% Requires: System Identification Toolbox
%
% Verifies that sidFreqETFE produces results matching MATLAB's etfe for
% SISO, MIMO, and time-series (periodogram) cases. Only unsmoothed ETFE
% is compared since sid and MATLAB use different smoothing windows.

fprintf('Running test_compareEtfe...\n');

%% Toolbox check
if ~exist('etfe', 'file')
    fprintf('  test_compareEtfe: SKIPPED (requires System Identification Toolbox)\n');
    return;
end

%% Test 1: SISO first-order system
rng(11);
N = 1024;
Ts = 1;
u = randn(N, 1);
y = filter(1, [1 -0.85], u) + 0.1 * randn(N, 1);

% MATLAB etfe: compute at its default frequencies
data = iddata(y, u, Ts);
G_etfe = etfe(data);
w_etfe = G_etfe.Frequency;
resp_etfe = squeeze(G_etfe.ResponseData);

% sid ETFE: compute at etfe's frequency grid
result_sid = sidFreqETFE(y, u, 'Frequencies', w_etfe);

relErr = max(abs(result_sid.Response - resp_etfe) ./ max(abs(resp_etfe), 1e-10));
assert(relErr < 0.01, ...
    'Test 1: SISO response relErr=%.6f should be <1%%', relErr);

%% Test 2: Noiseless first-order system (should match very closely)
rng(12);
N = 512;
u = randn(N, 1);
y = filter(1, [1 -0.7], u);  % no noise

data = iddata(y, u, Ts);
G_etfe = etfe(data);
w_etfe = G_etfe.Frequency;
resp_etfe = squeeze(G_etfe.ResponseData);

result_sid = sidFreqETFE(y, u, 'Frequencies', w_etfe);

relErr = max(abs(result_sid.Response - resp_etfe) ./ max(abs(resp_etfe), 1e-10));
assert(relErr < 0.01, ...
    'Test 2: noiseless response relErr=%.6f should be <1%%', relErr);

%% Test 3: Time-series periodogram
rng(13);
N = 1024;
e = randn(N, 1);
y = filter(1, [1 -0.6], e);

data = iddata(y, [], Ts);
G_etfe = etfe(data);
w_etfe = G_etfe.Frequency;
spec_etfe = real(squeeze(G_etfe.SpectrumData));

result_sid = sidFreqETFE(y, [], 'Frequencies', w_etfe);

% Compare periodogram (output spectrum)
relErr = max(abs(real(result_sid.NoiseSpectrum) - spec_etfe) ./ max(abs(spec_etfe), 1e-10));
assert(relErr < 0.01, ...
    'Test 3: periodogram relErr=%.6f should be <1%%', relErr);

assert(isempty(result_sid.Response), 'Test 3: time-series Response should be empty');

%% Test 4: FIR system
rng(14);
N = 2048;
u = randn(N, 1);
b = [1, 0.5, -0.3, 0.1];  % FIR filter
y = filter(b, 1, u);

data = iddata(y, u, Ts);
G_etfe = etfe(data);
w_etfe = G_etfe.Frequency;
resp_etfe = squeeze(G_etfe.ResponseData);

result_sid = sidFreqETFE(y, u, 'Frequencies', w_etfe);

relErr = max(abs(result_sid.Response - resp_etfe) ./ max(abs(resp_etfe), 1e-10));
assert(relErr < 0.01, ...
    'Test 4: FIR response relErr=%.6f should be <1%%', relErr);

%% Test 5: MIMO system (2 outputs, 1 input)
rng(15);
N = 1024;
u = randn(N, 1);
y1 = filter(1, [1 -0.5], u);
y2 = filter([0.3], [1 -0.7], u);
y = [y1, y2];

data = iddata(y, u, Ts);
G_etfe = etfe(data);
w_etfe = G_etfe.Frequency;
resp_etfe = permute(G_etfe.ResponseData, [3, 1, 2]);  % (ny x nu x nf) -> (nf x ny x nu)

result_sid = sidFreqETFE(y, u, 'Frequencies', w_etfe);

for ch = 1:2
    resp_sid_ch = result_sid.Response(:, ch);
    resp_etfe_ch = resp_etfe(:, ch);
    relErr = max(abs(resp_sid_ch - resp_etfe_ch) ./ max(abs(resp_etfe_ch), 1e-10));
    assert(relErr < 0.01, ...
        'Test 5: MIMO channel %d response relErr=%.6f should be <1%%', ch, relErr);
end

%% Test 6: Non-unit sample time
rng(16);
N = 512;
Ts = 0.01;
u = randn(N, 1);
y = filter(1, [1 -0.8], u);

w_sid = (1:128)' * pi / 128;
w_spa = w_sid / Ts;

data = iddata(y, u, Ts);
G_etfe = etfe(data);
w_etfe = G_etfe.Frequency;

% sid at etfe's frequencies (convert from rad/s to rad/sample)
w_etfe_radsamp = w_etfe * Ts;
result_sid = sidFreqETFE(y, u, 'Frequencies', w_etfe_radsamp, 'SampleTime', Ts);

resp_etfe = squeeze(G_etfe.ResponseData);
relErr = max(abs(result_sid.Response - resp_etfe) ./ max(abs(resp_etfe), 1e-10));
assert(relErr < 0.01, ...
    'Test 6: non-unit Ts response relErr=%.6f should be <1%%', relErr);

fprintf('  test_compareEtfe: ALL PASSED\n');
