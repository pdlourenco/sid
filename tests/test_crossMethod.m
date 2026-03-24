%% test_crossMethod - Cross-method consistency and comparison tests
%
% Verifies that different estimation methods produce consistent results
% when applied to the same data and known systems.

fprintf('Running test_crossMethod...\n');

%% Test 1: BT and BTFDR agree with equivalent parameters
% When BTFDR uses a resolution that maps to the same window size as BT,
% the results should be very close (not identical due to FFT vs direct DFT).
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1], [1 -0.8], u) + 0.1 * randn(N, 1);

M = 30;
R_equivalent = 2 * pi / M;
w = linspace(0.1, pi, 50)';  % Custom frequencies to force direct DFT in both

result_bt = sidFreqBT(y, u, 'WindowSize', M, 'Frequencies', w);
result_btfdr = sidFreqBTFDR(y, u, 'Resolution', R_equivalent, 'Frequencies', w);

% Response should be very close
relErr_G = max(abs(result_bt.Response - result_btfdr.Response)) / max(abs(result_bt.Response));
assert(relErr_G < 0.01, 'BT and BTFDR should agree when resolution=2pi/M (relErr=%.4f)', relErr_G);

% Noise spectrum should be very close
relErr_Phi = max(abs(result_bt.NoiseSpectrum - result_btfdr.NoiseSpectrum)) / max(abs(result_bt.NoiseSpectrum));
assert(relErr_Phi < 0.01, 'BT and BTFDR noise spectra should agree (relErr=%.4f)', relErr_Phi);

%% Test 2: All methods identify gain of pure-gain system
rng(1);
N = 1024;
u = randn(N, 1);
y = 2.5 * u;

result_bt = sidFreqBT(y, u);
result_etfe = sidFreqETFE(y, u);
result_btfdr = sidFreqBTFDR(y, u);

% All should give |G| ~ 2.5
assert(max(abs(abs(result_bt.Response) - 2.5)) < 0.05, 'BT: pure gain should be ~2.5');
assert(max(abs(abs(result_etfe.Response) - 2.5)) < 0.01, 'ETFE: pure gain should be ~2.5');
assert(max(abs(abs(result_btfdr.Response) - 2.5)) < 0.05, 'BTFDR: pure gain should be ~2.5');

%% Test 3: All methods give consistent phase for first-order system
rng(2);
N = 5000;
u = randn(N, 1);
a = 0.8;
y = filter(1, [1 -a], u);

w = (1:128)' * pi / 128;
G_true = 1 ./ (1 - a * exp(-1j * w));

result_bt = sidFreqBT(y, u);
result_btfdr = sidFreqBTFDR(y, u);

% Compare phases (excluding very low/high frequencies where precision is lower)
idx = 10:110;
phase_true = angle(G_true(idx));

phase_err_bt = abs(angle(result_bt.Response(idx)) - phase_true);
phase_err_bt = min(phase_err_bt, 2*pi - phase_err_bt);
assert(median(phase_err_bt) < 0.10, 'BT phase should match true system');

phase_err_btfdr = abs(angle(result_btfdr.Response(idx)) - phase_true);
phase_err_btfdr = min(phase_err_btfdr, 2*pi - phase_err_btfdr);
assert(median(phase_err_btfdr) < 0.10, 'BTFDR phase should match true system');

%% Test 4: Time series spectra are consistent across methods
rng(3);
N = 2000;
% AR(1) process: x(t) = 0.9*x(t-1) + e(t), spectrum = 1/|1-0.9*e^{-jw}|^2
x = filter(1, [1 -0.9], randn(N, 1));

result_bt = sidFreqBT(x, []);
result_etfe = sidFreqETFE(x, []);
result_btfdr = sidFreqBTFDR(x, []);

w = result_bt.Frequency;
Phi_true = 1 ./ abs(1 - 0.9 * exp(-1j * w)).^2;

% BT should be smooth and close to true
relErr_bt = median(abs(result_bt.NoiseSpectrum - Phi_true) ./ Phi_true);
assert(relErr_bt < 0.15, 'BT time series spectrum should be close to true (relErr=%.2f)', relErr_bt);

% BTFDR should also be close
relErr_btfdr = median(abs(result_btfdr.NoiseSpectrum - Phi_true) ./ Phi_true);
assert(relErr_btfdr < 0.15, 'BTFDR time series spectrum should be close (relErr=%.2f)', relErr_btfdr);

% ETFE (periodogram) has high variance but correct on average
% Check that median is in right ballpark
ratio_etfe = median(result_etfe.NoiseSpectrum ./ Phi_true);
assert(ratio_etfe > 0.5 && ratio_etfe < 2.0, ...
    'ETFE periodogram should be roughly correct on average (ratio=%.2f)', ratio_etfe);

%% Test 5: All methods produce same output struct fields
fields_bt = sort(fieldnames(result_bt));
fields_etfe = sort(fieldnames(result_etfe));
fields_btfdr = sort(fieldnames(result_btfdr));
assert(isequal(fields_bt, fields_etfe), 'BT and ETFE should have same fields');
assert(isequal(fields_bt, fields_btfdr), 'BT and BTFDR should have same fields');

%% Test 6: Frequency vectors are consistent
assert(max(abs(result_bt.Frequency - result_etfe.Frequency)) < 1e-12, ...
    'Default frequencies should be the same across methods');
assert(max(abs(result_bt.Frequency - result_btfdr.Frequency)) < 1e-12, ...
    'Default frequencies should be the same across methods');

%% Test 7: BT with ETFE-like large window approaches ETFE
% When BT uses M = N/2 (maximum), it's similar to the periodogram
rng(4);
N = 200;
u = randn(N, 1);
y = filter(1, [1 -0.5], u);
M_max = floor(N/2);
w = linspace(0.1, pi, 30)';
result_bt_large = sidFreqBT(y, u, M_max, w);
result_etfe = sidFreqETFE(y, u, 'Frequencies', w);

% They won't be identical but responses should be correlated
corr_mag = corrcoef(abs(result_bt_large.Response), abs(result_etfe.Response));
assert(corr_mag(1,2) > 0.7, 'BT with large M should correlate with ETFE (corr=%.2f)', corr_mag(1,2));

%% Test 8: MIMO results have consistent dimensions across methods
rng(5);
N = 500;
u = randn(N, 2);
y = [u(:,1) + 0.3*u(:,2) + 0.1*randn(N,1), ...
     0.5*u(:,1) + u(:,2) + 0.1*randn(N,1)];

result_bt = sidFreqBT(y, u);
result_etfe = sidFreqETFE(y, u);
result_btfdr = sidFreqBTFDR(y, u);

nf = 128;
assert(isequal(size(result_bt.Response), [nf 2 2]), 'BT MIMO: wrong Response size');
assert(isequal(size(result_etfe.Response), [nf 2 2]), 'ETFE MIMO: wrong Response size');
assert(isequal(size(result_btfdr.Response), [nf 2 2]), 'BTFDR MIMO: wrong Response size');

%% Test 9: Smoothed ETFE approaches BT-like smoothness
rng(6);
N = 1000;
u = randn(N, 1);
y = filter(1, [1 -0.7], u) + 0.2*randn(N, 1);
result_etfe_smooth = sidFreqETFE(y, u, 'Smoothing', 21);
result_bt = sidFreqBT(y, u, 'WindowSize', 30);

% Both smoothed estimates should have similar variance
var_etfe = var(abs(result_etfe_smooth.Response));
var_bt = var(abs(result_bt.Response));
% They should be in the same order of magnitude
ratio = var_etfe / var_bt;
assert(ratio > 0.1 && ratio < 10, ...
    'Smoothed ETFE and BT should have comparable variance (ratio=%.2f)', ratio);

fprintf('  test_crossMethod: ALL PASSED\n');
