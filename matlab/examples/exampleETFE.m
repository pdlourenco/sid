%% exampleETFE - Empirical Transfer Function Estimate with sidFreqETFE
%
% This example demonstrates sidFreqETFE, which estimates the frequency
% response as the ratio of output and input DFTs. It provides maximum
% frequency resolution but high variance. Optional smoothing reduces
% variance at the cost of resolution.

runner__nCompleted = 0;

%% Generate test data
% True system: G(z) = 1 / (1 - 0.8 z^{-1})  (first-order, pole at 0.8)

rng(1);
N = 1024;
u = randn(N, 1);
y_clean = filter(1, [1 -0.8], u);
y = y_clean + 0.3 * randn(N, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Basic ETFE (no smoothing)
% The raw ETFE has maximum resolution but is very noisy.
% Note: ResponseStd is NaN because ETFE has no asymptotic uncertainty formula.

result = sidFreqETFE(y, u);

figure;
sidBodePlot(result);
title('ETFE - Raw (No Smoothing)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Basic ETFE (no smoothing).\n', runner__nCompleted);

%% Effect of smoothing
% Smoothing averages nearby frequency bins using a Hann window of odd length S.
% Larger S = smoother estimate but coarser resolution.

r1  = sidFreqETFE(y, u, 'Smoothing', 1);   % no smoothing (raw)
r11 = sidFreqETFE(y, u, 'Smoothing', 11);  % moderate smoothing
r21 = sidFreqETFE(y, u, 'Smoothing', 21);  % heavy smoothing

w = r1.Frequency;
G_true = 1 ./ (1 - 0.8 * exp(-1j * w));

figure;
semilogx(w, 20*log10(abs(r1.Response)),  'Color', [0.7 0.7 0.7], 'DisplayName', 'S = 1 (raw)');
hold on;
semilogx(w, 20*log10(abs(r11.Response)), 'b', 'DisplayName', 'S = 11');
semilogx(w, 20*log10(abs(r21.Response)), 'r', 'DisplayName', 'S = 21');
semilogx(w, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('ETFE Smoothing Comparison');
legend('show', 'Location', 'southwest');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Effect of smoothing.\n', runner__nCompleted);

%% Known FIR system: pure delay
% For y(t) = u(t-1), the true transfer function is G(z) = z^{-1}.
% ETFE recovers this exactly when there is no noise.

rng(2);
N = 1024;
u_fir = randn(N, 1);
y_fir = [0; u_fir(1:end-1)];   % one-sample delay

result_fir = sidFreqETFE(y_fir, u_fir);
w_fir = result_fir.Frequency;

figure;
subplot(2,1,1);
plot(w_fir, abs(result_fir.Response), 'b');
ylabel('|G|');
title('ETFE of Pure Delay: |G| should be 1');
grid on;

subplot(2,1,2);
plot(w_fir, angle(result_fir.Response), 'b');
hold on;
plot(w_fir, -w_fir, 'k--', 'DisplayName', 'True phase = -\omega');
ylabel('Phase (rad)');
xlabel('Frequency (rad/sample)');
legend('ETFE', 'True', 'Location', 'southwest');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Known FIR system: pure delay.\n', runner__nCompleted);

%% Time-series mode: periodogram
% With no input signal, sidFreqETFE computes the periodogram of the output.

rng(3);
y_ts = filter(1, [1 -0.8], randn(500, 1));
result_ts = sidFreqETFE(y_ts, []);

figure;
sidSpectrumPlot(result_ts);
title('Periodogram of AR(1) Process');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Time-series mode: periodogram.\n', runner__nCompleted);

%% Custom frequency grid and Hz display
% Use a logarithmic frequency grid and plot in Hz.

Ts = 0.001;                                           % 1 kHz sampling
w_log = logspace(log10(0.05), log10(pi), 200)';       % log-spaced in rad/sample

rng(4);
N = 2048;
u_hz = randn(N, 1);
y_hz = filter(1, [1 -0.9], u_hz) + 0.1 * randn(N, 1);

result_hz = sidFreqETFE(y_hz, u_hz, 'Smoothing', 11, ...
    'Frequencies', w_log, 'SampleTime', Ts);

figure;
sidBodePlot(result_hz, 'FrequencyUnit', 'Hz');
title('ETFE with Log Frequency Grid (Hz)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Custom frequency grid and Hz display.\n', runner__nCompleted);

fprintf('exampleETFE: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
