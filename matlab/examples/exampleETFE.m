%% exampleETFE - ETFE on a 1-DoF spring-mass-damper
%
% This example estimates the frequency response of Plant A (the same
% SDOF used by exampleSISO: m = 1, k = 100, c = 2) via sidFreqETFE,
% then explores the smoothing-vs-resolution trade-off. Because the
% plant is the same as exampleSISO, this example also serves as a
% direct BT-vs-ETFE comparison.
%
% See spec/EXAMPLES.md section 3.2 for the binding specification.

runner__nCompleted = 0;

%% Generate test data
% Build Plant A with util_msd and simulate under white-force excitation.
% Also precompute the exact discrete transfer function
%    G(e^{jw}) = C (e^{jw} I - Ad)^-1 Bd
% for later overlays.

rng(1);

m  = 1.0;    k  = 100.0;    c  = 2.0;    F  = 1.0;
Ts = 0.01;   N  = 2048;

[Ad, Bd] = util_msd(m, k, c, F, Ts);
C_out = [1 0];

u = randn(N, 1);
x = zeros(N + 1, 2);
for step = 1:N
    x(step + 1, :) = (Ad * x(step, :)' + Bd * u(step))';
end
y_clean = x(2:end, 1);
y = y_clean + 2e-4 * randn(N, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Basic ETFE (no smoothing)
% The raw ETFE is just Y(w) / U(w) -- maximum resolution, extremely
% noisy. ResponseStd is NaN because ETFE has no asymptotic uncertainty
% formula.

result = sidFreqETFE(y, u, 'SampleTime', Ts);

figure;
sidBodePlot(result);
title('ETFE raw (no smoothing)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Basic ETFE (no smoothing).\n', ...
    runner__nCompleted);

%% Effect of smoothing
% Smoothing averages S adjacent frequency bins (S must be odd). Larger
% S suppresses periodogram noise but smears the resonance peak. We
% overlay the true discrete transfer function as a dashed reference.

r1  = sidFreqETFE(y, u, 'Smoothing',  1, 'SampleTime', Ts);
r11 = sidFreqETFE(y, u, 'Smoothing', 11, 'SampleTime', Ts);
r21 = sidFreqETFE(y, u, 'Smoothing', 21, 'SampleTime', Ts);

w = r1.Frequency;
nf = length(w);
G_true = zeros(nf, 1);
I2 = eye(2);
for i = 1:nf
    G_true(i) = C_out * ((exp(1j * w(i)) * I2 - Ad) \ Bd);
end

figure;
semilogx(w, 20*log10(abs(r1.Response)),  'Color', [0.7 0.7 0.7], ...
    'DisplayName', 'S = 1 (raw)');
hold on;
semilogx(w, 20*log10(abs(r11.Response)), 'b', 'DisplayName', 'S = 11');
semilogx(w, 20*log10(abs(r21.Response)), 'r', 'DisplayName', 'S = 21');
semilogx(w, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('ETFE smoothing comparison');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Effect of smoothing.\n', runner__nCompleted);

%% Known FIR system: pure delay
% For y(t) = u(t-1) the true transfer function is G(z) = z^{-1}: unit
% magnitude and phase -w at every frequency. ETFE recovers this exactly
% in the noise-free case.

rng(2);
N_fir = 1024;
u_fir = randn(N_fir, 1);
y_fir = [0; u_fir(1:end-1)];

result_fir = sidFreqETFE(y_fir, u_fir);
w_fir = result_fir.Frequency;

figure;
subplot(2, 1, 1);
plot(w_fir, abs(result_fir.Response), 'b');
ylabel('|G|');
title('ETFE of pure delay: |G| should be 1');
grid on;

subplot(2, 1, 2);
plot(w_fir, angle(result_fir.Response), 'b', 'DisplayName', 'ETFE');
hold on;
plot(w_fir, -w_fir, 'k--', 'DisplayName', 'True phase = -\omega');
hold off;
ylabel('Phase (rad)');
xlabel('Frequency (rad/sample)');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Known FIR system: pure delay.\n', ...
    runner__nCompleted);

%% Time-series mode: periodogram
% With no input, sidFreqETFE computes the periodogram of the output.
% Re-simulate the same SDOF under fresh white-force excitation and
% pass only the position record.

rng(3);
N_ts = 500;
u_ts = randn(N_ts, 1);
x_ts = zeros(N_ts + 1, 2);
for step = 1:N_ts
    x_ts(step + 1, :) = (Ad * x_ts(step, :)' + Bd * u_ts(step))';
end
y_ts = x_ts(2:end, 1);

result_ts = sidFreqETFE(y_ts, []);

figure;
sidSpectrumPlot(result_ts);
title('SDOF output periodogram');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Time-series mode: periodogram.\n', ...
    runner__nCompleted);

%% Custom frequency grid and Hz display
% For narrow resonances the uniform DFT grid can undersample the peak.
% Supply a log-spaced grid and plot in Hz.

rng(4);
w_log = logspace(log10(0.005), log10(pi), 200)';

N_hz = 2048;
u_hz = randn(N_hz, 1);
x_hz = zeros(N_hz + 1, 2);
for step = 1:N_hz
    x_hz(step + 1, :) = (Ad * x_hz(step, :)' + Bd * u_hz(step))';
end
y_hz = x_hz(2:end, 1) + 2e-4 * randn(N_hz, 1);

result_hz = sidFreqETFE(y_hz, u_hz, 'Smoothing', 11, ...
    'Frequencies', w_log, 'SampleTime', Ts);

figure;
sidBodePlot(result_hz, 'FrequencyUnit', 'Hz');
title('ETFE with log frequency grid (Hz)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Custom frequency grid and Hz display.\n', ...
    runner__nCompleted);

fprintf('exampleETFE: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
