%% exampleMethodComparison - Comparing frequency response estimation methods
%
% This example compares the three frequency-domain estimators in sid:
%   sidFreqBT    - Blackman-Tukey (fixed window, replaces spa)
%   sidFreqBTFDR - Blackman-Tukey with frequency-dependent resolution (replaces spafdr)
%   sidFreqETFE  - Empirical Transfer Function Estimate (replaces etfe)

%% Generate test data
% First-order system with moderate noise.

rng(4);
N = 2000;
u = randn(N, 1);
y_clean = filter(1, [1 -0.85], u);
y = y_clean + 0.3 * randn(N, 1);

%% Estimate with all three methods

r_bt   = sidFreqBT(y, u, 'WindowSize', 30);
r_etfe = sidFreqETFE(y, u);
r_etfe_s = sidFreqETFE(y, u, 'Smoothing', 15);
r_fdr  = sidFreqBTFDR(y, u, 'Resolution', 0.3);

%% Compare Bode magnitude plots

w = r_bt.Frequency;
G_true = 1 ./ (1 - 0.85 * exp(-1j * w));

figure;
semilogx(w, 20*log10(abs(r_etfe.Response)), ...
    'Color', [0.8 0.8 0.8], 'DisplayName', 'ETFE (raw)');
hold on;
semilogx(w, 20*log10(abs(r_etfe_s.Response)), ...
    'g', 'DisplayName', 'ETFE (S=15)');
semilogx(w, 20*log10(abs(r_bt.Response)), ...
    'b', 'DisplayName', 'BT (M=30)');
semilogx(w, 20*log10(abs(r_fdr.Response)), ...
    'r', 'DisplayName', 'BTFDR (R=0.3)');
semilogx(w, 20*log10(abs(G_true)), ...
    'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Method Comparison: Bode Magnitude');
legend('show', 'Location', 'southwest');
grid on;
hold off;

%% Compare noise spectra
% BT and BTFDR compute the noise spectrum from covariance estimates.
% ETFE computes it from residuals.

figure;
semilogx(w, 10*log10(abs(r_bt.NoiseSpectrum)), 'b', 'DisplayName', 'BT');
hold on;
semilogx(w, 10*log10(abs(r_etfe_s.NoiseSpectrum)), 'g', 'DisplayName', 'ETFE (S=15)');
semilogx(w, 10*log10(abs(r_fdr.NoiseSpectrum)), 'r', 'DisplayName', 'BTFDR');
xlabel('Frequency (rad/sample)');
ylabel('Noise Spectrum (dB)');
title('Noise Spectrum Comparison');
legend('show');
grid on;
hold off;

%% Custom logarithmic frequency grid
% A log-spaced grid provides better low-frequency coverage.

w_log = logspace(log10(0.05), log10(pi), 200)';

r_bt_log  = sidFreqBT(y, u, 'WindowSize', 30, 'Frequencies', w_log);
r_etfe_log = sidFreqETFE(y, u, 'Smoothing', 15, 'Frequencies', w_log);
r_fdr_log = sidFreqBTFDR(y, u, 'Resolution', 0.3, 'Frequencies', w_log);

G_true_log = 1 ./ (1 - 0.85 * exp(-1j * w_log));

figure;
semilogx(w_log, 20*log10(abs(r_bt_log.Response)), 'b', 'DisplayName', 'BT');
hold on;
semilogx(w_log, 20*log10(abs(r_etfe_log.Response)), 'g', 'DisplayName', 'ETFE');
semilogx(w_log, 20*log10(abs(r_fdr_log.Response)), 'r', 'DisplayName', 'BTFDR');
semilogx(w_log, 20*log10(abs(G_true_log)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Log Frequency Grid (200 Points)');
legend('show', 'Location', 'southwest');
grid on;
hold off;

%% Time-series comparison: periodogram vs smoothed spectrum
% With no input (u=[]), all three methods estimate the output power spectrum.

rng(5);
y_ts = filter(1, [1 -0.8], randn(1000, 1));

r_ts_bt   = sidFreqBT(y_ts, [], 'WindowSize', 30);
r_ts_etfe = sidFreqETFE(y_ts, []);

w_ts = r_ts_bt.Frequency;
% True spectrum: |1/(1-0.8*e^{-jw})|^2
Phi_true = abs(1 ./ (1 - 0.8 * exp(-1j * w_ts))).^2;

figure;
semilogx(w_ts, 10*log10(abs(r_ts_etfe.NoiseSpectrum)), ...
    'Color', [0.7 0.7 0.7], 'DisplayName', 'ETFE (periodogram)');
hold on;
semilogx(w_ts, 10*log10(abs(r_ts_bt.NoiseSpectrum)), ...
    'b', 'DisplayName', 'BT (M=30)');
semilogx(w_ts, 10*log10(Phi_true), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Power Spectrum (dB)');
title('Time-Series: Periodogram vs Blackman-Tukey');
legend('show', 'Location', 'southwest');
grid on;
hold off;

%% Model output comparison using sidCompare
% Compare how well each method predicts the measured output.
comp_bt   = sidCompare(r_bt, y, u);
comp_fdr  = sidCompare(r_fdr, y, u);
comp_etfe = sidCompare(r_etfe, y, u);

fprintf('\n--- Prediction Fit (NRMSE %%) ---\n');
fprintf('  sidFreqBT:    %.1f%%\n', comp_bt.Fit);
fprintf('  sidFreqBTFDR: %.1f%%\n', comp_fdr.Fit);
fprintf('  sidFreqETFE:  %.1f%%\n', comp_etfe.Fit);

%% Summary of method trade-offs
fprintf('\n--- Method Comparison Summary ---\n');
fprintf('%-12s  %-12s  %-11s  %-9s\n', 'Method', 'WindowSize', 'Uncertainty', 'Coherence');
fprintf('%-12s  %-12s  %-11s  %-9s\n', '------', '----------', '-----------', '---------');
fprintf('%-12s  %-12s  %-11s  %-9s\n', 'sidFreqBT',   'Fixed M',      'Yes', 'Yes');
fprintf('%-12s  %-12s  %-11s  %-9s\n', 'sidFreqBTFDR','Per-freq M_k', 'Yes', 'Yes');
fprintf('%-12s  %-12s  %-11s  %-9s\n', 'sidFreqETFE', 'N (full)',     'No',  'No');
