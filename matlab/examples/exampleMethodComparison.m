%% exampleMethodComparison - Comparing frequency response estimators on Plant A
%
% Applies the three frequency-domain estimators (sidFreqBT,
% sidFreqBTFDR, sidFreqETFE) to the same physical 1-DoF SDOF plant
% used by exampleSISO and exampleETFE. On a narrow resonance the raw
% ETFE (full DFT resolution) often beats the smoothed BT estimate on
% NRMSE fit -- resolution wins over smoothness for sharp peaks.
%
% See spec/EXAMPLES.md section 3.5 for the binding specification.

runner__nCompleted = 0;

%% Generate test data
% Plant A: m = 1, k = 100, c = 2 (same as exampleSISO). N = 2048
% samples of white-force excitation and sensor-noise corrupted
% position output.

rng(4);

m  = 1.0;    k  = 100.0;    c  = 2.0;    F  = 1.0;
Ts = 0.01;   N  = 2048;

[Ad, Bd] = util_msd(m, k, c, F, Ts);
C_out = [1 0];

u = randn(N, 1);
x = zeros(N + 1, 2);
for step = 1:N
    x(step + 1, :) = (Ad * x(step, :)' + Bd * u(step))';
end
y = x(2:end, 1) + 2e-4 * randn(N, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Estimate with all three methods
% Shared custom frequency grid so the estimators compare on identical
% bins.

w_grid = linspace(0.005, pi, 512)';

r_bt     = sidFreqBT(y, u, 'WindowSize', 200, 'Frequencies', w_grid, ...
                      'SampleTime', Ts);
r_etfe   = sidFreqETFE(y, u, 'Frequencies', w_grid, 'SampleTime', Ts);
r_etfe_s = sidFreqETFE(y, u, 'Smoothing', 15, 'Frequencies', w_grid, ...
                        'SampleTime', Ts);
r_fdr    = sidFreqBTFDR(y, u, 'Resolution', 0.3, 'Frequencies', w_grid, ...
                         'SampleTime', Ts);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Estimate with all three methods.\n', ...
    runner__nCompleted);

%% Compare Bode magnitude plots
% Overlay the four curves plus the true discrete transfer function as
% a dashed reference.

w = r_bt.Frequency;
nf = length(w);
G_true = zeros(nf, 1);
I2 = eye(2);
for i = 1:nf
    G_true(i) = C_out * ((exp(1j * w(i)) * I2 - Ad) \ Bd);
end

figure;
semilogx(w, 20*log10(abs(r_etfe.Response)), 'Color', [0.8 0.8 0.8], ...
    'DisplayName', 'ETFE (raw)');
hold on;
semilogx(w, 20*log10(abs(r_etfe_s.Response)), 'g', ...
    'DisplayName', 'ETFE (S = 15)');
semilogx(w, 20*log10(abs(r_bt.Response)), 'b', ...
    'DisplayName', 'BT (M = 200)');
semilogx(w, 20*log10(abs(r_fdr.Response)), 'r', ...
    'DisplayName', 'BTFDR (R = 0.3)');
semilogx(w, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Method comparison: Bode magnitude');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Compare Bode magnitude plots.\n', ...
    runner__nCompleted);

%% Compare noise spectra
% BT and BTFDR compute the noise spectrum from covariance estimates,
% ETFE from residuals. They should all trace out roughly the same
% noise-floor shape.

figure;
% Clamp the noise spectrum with a small floor before log10: per
% SPEC.md 2.7 the noise spectrum is non-negative and may hit exact
% zero at frequencies where the estimator variance has been clamped
% from a slightly-negative value, which would otherwise produce
% log10(0) warnings / "non-positive data" log-plot warnings.
noise_floor = eps;
semilogx(w, 10*log10(max(abs(r_bt.NoiseSpectrum), noise_floor)), 'b', ...
    'DisplayName', 'BT');
hold on;
semilogx(w, 10*log10(max(abs(r_etfe_s.NoiseSpectrum), noise_floor)), 'g', ...
    'DisplayName', 'ETFE (S = 15)');
semilogx(w, 10*log10(max(abs(r_fdr.NoiseSpectrum), noise_floor)), 'r', ...
    'DisplayName', 'BTFDR');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Noise spectrum (dB)');
title('Noise spectrum comparison');
legend;
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Compare noise spectra.\n', runner__nCompleted);

%% Custom logarithmic frequency grid
% Log-spaced grid covers the low-frequency region more densely.

w_log = logspace(log10(0.005), log10(pi), 200)';

r_bt_log   = sidFreqBT(y, u, 'WindowSize', 200, 'Frequencies', w_log, ...
                       'SampleTime', Ts);
r_etfe_log = sidFreqETFE(y, u, 'Smoothing', 15, 'Frequencies', w_log, ...
                          'SampleTime', Ts);
r_fdr_log  = sidFreqBTFDR(y, u, 'Resolution', 0.3, 'Frequencies', w_log, ...
                           'SampleTime', Ts);

G_true_log = zeros(length(w_log), 1);
for i = 1:length(w_log)
    G_true_log(i) = C_out * ((exp(1j * w_log(i)) * I2 - Ad) \ Bd);
end

figure;
semilogx(w_log, 20*log10(abs(r_bt_log.Response)), 'b', 'DisplayName', 'BT');
hold on;
semilogx(w_log, 20*log10(abs(r_etfe_log.Response)), 'g', 'DisplayName', 'ETFE');
semilogx(w_log, 20*log10(abs(r_fdr_log.Response)), 'r', 'DisplayName', 'BTFDR');
semilogx(w_log, 20*log10(abs(G_true_log)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Log frequency grid (200 points)');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Custom logarithmic frequency grid.\n', ...
    runner__nCompleted);

%% Time-series comparison: periodogram vs smoothed spectrum
% With no input, all methods estimate the output power spectrum.

rng(5);
N_ts = 1000;
u_ts = randn(N_ts, 1);
x_ts = zeros(N_ts + 1, 2);
for step = 1:N_ts
    x_ts(step + 1, :) = (Ad * x_ts(step, :)' + Bd * u_ts(step))';
end
y_ts = x_ts(2:end, 1);

r_ts_bt   = sidFreqBT(y_ts, [], 'WindowSize', 100);
r_ts_etfe = sidFreqETFE(y_ts, []);

w_ts = r_ts_bt.Frequency;

figure;
semilogx(w_ts, 10*log10(max(abs(r_ts_etfe.NoiseSpectrum), noise_floor)), ...
    'Color', [0.7 0.7 0.7], 'DisplayName', 'ETFE (periodogram)');
hold on;
semilogx(w_ts, 10*log10(max(abs(r_ts_bt.NoiseSpectrum), noise_floor)), 'b', ...
    'DisplayName', 'BT (M = 100)');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Power spectrum (dB)');
title('Time-series: periodogram vs Blackman-Tukey');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Time-series comparison.\n', ...
    runner__nCompleted);

%% Model output comparison using sidCompare
% On this lightly-damped plant, the raw ETFE wins on NRMSE because it
% preserves the full resolution of the narrow resonance.

comp_bt     = sidCompare(r_bt,     y, u);
comp_fdr    = sidCompare(r_fdr,    y, u);
comp_etfe   = sidCompare(r_etfe,   y, u);
comp_etfe_s = sidCompare(r_etfe_s, y, u);

fprintf('--- Prediction fit (NRMSE %%) ---\n');
fprintf('  sidFreqBT    (M=200):  %5.1f%%\n', comp_bt.Fit(1));
fprintf('  sidFreqBTFDR (R=0.3):  %5.1f%%\n', comp_fdr.Fit(1));
fprintf('  sidFreqETFE  (raw):    %5.1f%%\n', comp_etfe.Fit(1));
fprintf('  sidFreqETFE  (S=15):   %5.1f%%\n', comp_etfe_s.Fit(1));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Model output comparison using sidCompare.\n', ...
    runner__nCompleted);

%% Summary of method trade-offs
% (Printed as a small text table since MATLAB has no native markdown.)

fprintf('\nMethod trade-offs:\n');
fprintf('  %-14s  %-15s  %-12s  %-10s\n', ...
    'Method', 'Window size', 'Uncertainty', 'Coherence');
fprintf('  %-14s  %-15s  %-12s  %-10s\n', ...
    '------', '-----------', '-----------', '---------');
fprintf('  %-14s  %-15s  %-12s  %-10s\n', ...
    'sidFreqBT',    'Fixed M',         'Yes', 'Yes');
fprintf('  %-14s  %-15s  %-12s  %-10s\n', ...
    'sidFreqBTFDR', 'Per-freq M(k)',   'Yes', 'Yes');
fprintf('  %-14s  %-15s  %-12s  %-10s\n', ...
    'sidFreqETFE',  'N (full DFT)',    'No',  'No');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Summary of method trade-offs.\n', ...
    runner__nCompleted);

fprintf('exampleMethodComparison: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
