%% exampleFreqDepRes - Frequency-dependent resolution with sidFreqBTFDR
%
% This example demonstrates sidFreqBTFDR, which uses a different window
% size at each frequency. This is valuable when a system has sharp features
% (e.g., resonances) that need fine resolution at some frequencies but not
% others. Replaces MATLAB's spafdr.

runner__nCompleted = 0;

%% Generate test data: second-order resonant system
% Poles at 0.9*exp(+/-j*pi/4) create a resonance peak near w = pi/4.

rng(2);
N = 5000;
Ts = 0.01;
b = 1;
a_coeff = [1, -2*0.9*cos(pi/4), 0.9^2];
u = randn(N, 1);
y = filter(b, a_coeff, u) + 0.1 * randn(N, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Generate test data: second-order resonant system');

%% Fixed-window sidFreqBT: the resolution-variance trade-off
% Small M (=15): smooth but misses the resonance peak.
% Large M (=80): captures the peak but has high variance.

r_small = sidFreqBT(y, u, 'WindowSize', 15, 'SampleTime', Ts);
r_large = sidFreqBT(y, u, 'WindowSize', 80, 'SampleTime', Ts);

w = r_small.Frequency;
G_true = 1 ./ (1 - 2*0.9*cos(pi/4)*exp(-1j*w) + 0.81*exp(-2j*w));

figure;
semilogx(w, 20*log10(abs(r_small.Response)), 'b', 'DisplayName', 'BT M=15');
hold on;
semilogx(w, 20*log10(abs(r_large.Response)), 'r', 'DisplayName', 'BT M=80');
semilogx(w, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Fixed Window: Resolution vs Variance Trade-off');
legend('Location', 'southwest');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, ...
    'Fixed-window sidFreqBT: the resolution-variance trade-off');

%% Scalar resolution with sidFreqBTFDR
% Resolution R sets the window size as M = round(2*pi / R).
% Smaller R = finer resolution (larger window).

result_fdr = sidFreqBTFDR(y, u, 'Resolution', 0.2, 'SampleTime', Ts);

figure;
sidBodePlot(result_fdr);
title('sidFreqBTFDR with Scalar Resolution R = 0.2');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Scalar resolution with sidFreqBTFDR.\n', runner__nCompleted);

%% Per-frequency resolution vector
% Use fine resolution near the resonance (low frequencies) and coarse
% resolution at high frequencies where the response is flat.

nf = length(result_fdr.Frequency);
R_vec = linspace(0.1, 1.5, nf)';  % fine at low freq, coarse at high freq

result_vec = sidFreqBTFDR(y, u, 'Resolution', R_vec, 'SampleTime', Ts);

% The resulting WindowSize varies across frequency
figure;
subplot(2,1,1);
plot(result_vec.Frequency, result_vec.WindowSize, 'b');
xlabel('Frequency (rad/sample)');
ylabel('Window Size M');
title('Per-Frequency Window Size');
grid on;

subplot(2,1,2);
semilogx(result_vec.Frequency, 20*log10(abs(result_vec.Response)), 'b', ...
    'DisplayName', 'BTFDR (variable R)');
hold on;
semilogx(w, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('BTFDR with Per-Frequency Resolution');
legend('Location', 'southwest');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Per-frequency resolution vector.\n', runner__nCompleted);

%% Compare BT vs BTFDR side by side
% BTFDR adapts to capture the peak while keeping variance low elsewhere.

r_bt = sidFreqBT(y, u, 'WindowSize', 30, 'SampleTime', Ts);
r_fdr = sidFreqBTFDR(y, u, 'Resolution', 0.3, 'SampleTime', Ts);

figure;
semilogx(r_bt.Frequency, 20*log10(abs(r_bt.Response)), 'b', ...
    'DisplayName', 'BT (M=30)');
hold on;
semilogx(r_fdr.Frequency, 20*log10(abs(r_fdr.Response)), 'r', ...
    'DisplayName', 'BTFDR (R=0.3)');
semilogx(r_bt.Frequency, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Blackman-Tukey: Fixed vs Frequency-Dependent Resolution');
legend('Location', 'southwest');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Compare BT vs BTFDR side by side.\n', runner__nCompleted);

fprintf('exampleFreqDepRes: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
