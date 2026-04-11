%% exampleFreqDepRes - Frequency-dependent resolution on a 3-mass chain
%
% This example demonstrates sidFreqBTFDR on a lightly-damped 3-mass
% SMD chain (Plant C) with three well-separated modes in the bottom
% decade of the spectrum. A short fixed BT window smears the modes,
% a long fixed window resolves them at the cost of high-frequency
% variance, and a per-frequency resolution vector combines the best
% of both.
%
% See spec/EXAMPLES.md section 3.3 for the binding specification.

runner__nCompleted = 0;

%% Generate test data
% Plant C: m = [1 1 1], k = [300 200 100], c = [8 8 8]. Three modes at
% approximately 6.4, 15.1, 25.1 rad/s. Force at mass 1, measure x1.

rng(2);

m  = [1; 1; 1];
k  = [300; 200; 100];
c  = [8; 8; 8];
F  = [1; 0; 0];
Ts = 0.01;
N  = 6000;

[Ad, Bd] = util_msd(m, k, c, F, Ts);
C_out = [1 0 0 0 0 0];

u = randn(N, 1);
x = zeros(N + 1, 6);
for step = 1:N
    x(step + 1, :) = (Ad * x(step, :)' + Bd * u(step))';
end
y = x(2:end, 1) + 5e-4 * randn(N, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Fixed-window sidFreqBT: the resolution-variance trade-off
% Short window (M = 15): smooth but smears the three modes.
% Long window (M = 80): captures the peaks but has high variance.

r_small = sidFreqBT(y, u, 'WindowSize', 15, 'SampleTime', Ts);
r_large = sidFreqBT(y, u, 'WindowSize', 80, 'SampleTime', Ts);

w = r_small.Frequency;
nf = length(w);
G_true = zeros(nf, 1);
I6 = eye(6);
for i = 1:nf
    G_true(i) = C_out * ((exp(1j * w(i)) * I6 - Ad) \ Bd);
end

figure;
semilogx(w, 20*log10(abs(r_small.Response)), 'b', 'DisplayName', 'BT M = 15');
hold on;
semilogx(w, 20*log10(abs(r_large.Response)), 'r', 'DisplayName', 'BT M = 80');
semilogx(w, 20*log10(abs(G_true)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Fixed window: resolution vs variance trade-off');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', runner__nCompleted, ...
    'Fixed-window sidFreqBT: the resolution-variance trade-off');

%% Scalar resolution with sidFreqBTFDR
% Resolution R sets the window size as M = ceil(2*pi / R). Smaller R
% = finer resolution (larger window).

result_fdr = sidFreqBTFDR(y, u, 'Resolution', 0.2, 'SampleTime', Ts);

figure;
sidBodePlot(result_fdr);
title('sidFreqBTFDR with scalar resolution R = 0.2');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Scalar resolution with sidFreqBTFDR.\n', ...
    runner__nCompleted);

%% Per-frequency resolution vector
% Use fine resolution (R small) at low frequencies where the modes are,
% coarse resolution (R large) at high frequencies where the response
% is smooth.

nf2 = length(result_fdr.Frequency);
R_vec = linspace(0.1, 1.5, nf2)';

result_vec = sidFreqBTFDR(y, u, 'Resolution', R_vec, 'SampleTime', Ts);

figure;
subplot(2, 1, 1);
plot(result_vec.Frequency, result_vec.WindowSize, 'b');
xlabel('Frequency (rad/sample)');
ylabel('Window size M');
title('Per-frequency window size');
grid on;

w2 = result_vec.Frequency;
G_true2 = zeros(length(w2), 1);
for i = 1:length(w2)
    G_true2(i) = C_out * ((exp(1j * w2(i)) * I6 - Ad) \ Bd);
end

subplot(2, 1, 2);
semilogx(w2, 20*log10(abs(result_vec.Response)), 'b', ...
    'DisplayName', 'BTFDR (variable R)');
hold on;
semilogx(w2, 20*log10(abs(G_true2)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('BTFDR with per-frequency resolution');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Per-frequency resolution vector.\n', ...
    runner__nCompleted);

%% Compare BT vs BTFDR side by side
% BTFDR adapts: big window near the modes, small window at high
% frequencies where the response is featureless.

r_bt  = sidFreqBT(y, u, 'WindowSize', 30, 'SampleTime', Ts);
r_fdr = sidFreqBTFDR(y, u, 'Resolution', 0.3, 'SampleTime', Ts);

w3 = r_bt.Frequency;
G_true3 = zeros(length(w3), 1);
for i = 1:length(w3)
    G_true3(i) = C_out * ((exp(1j * w3(i)) * I6 - Ad) \ Bd);
end

figure;
semilogx(r_bt.Frequency, 20*log10(abs(r_bt.Response)), 'b', ...
    'DisplayName', 'BT (M = 30)');
hold on;
semilogx(r_fdr.Frequency, 20*log10(abs(r_fdr.Response)), 'r', ...
    'DisplayName', 'BTFDR (R = 0.3)');
semilogx(w3, 20*log10(abs(G_true3)), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'True');
hold off;
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('BT: fixed vs frequency-dependent resolution');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Compare BT vs BTFDR side by side.\n', ...
    runner__nCompleted);

fprintf('exampleFreqDepRes: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
