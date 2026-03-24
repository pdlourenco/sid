%% exampleSISO - Basic frequency response estimation with sidFreqBT
%
% This example demonstrates how to estimate the frequency response of a
% simple SISO system and plot the results with confidence bands.

%% Generate test data
% True system: G(z) = 1 / (1 - 0.9 z^{-1})
% This is a stable first-order system with a pole at z = 0.9.

N = 1000;                          % Number of samples
Ts = 0.01;                         % Sample time (seconds)
u = randn(N, 1);                   % White noise input
y_clean = filter(1, [1 -0.9], u);  % Noiseless output
noise = 0.1 * randn(N, 1);         % Measurement noise
y = y_clean + noise;                % Noisy output

%% Estimate frequency response using Blackman-Tukey
result = sidFreqBT(y, u, 'SampleTime', Ts);

%% Plot Bode diagram
figure;
sidBodePlot(result);

%% Plot noise spectrum
figure;
sidSpectrumPlot(result);

%% Compare different window sizes
% Larger window = finer resolution but more variance.
r10  = sidFreqBT(y, u, 'WindowSize', 10,  'SampleTime', Ts);
r30  = sidFreqBT(y, u, 'WindowSize', 30,  'SampleTime', Ts);
r100 = sidFreqBT(y, u, 'WindowSize', 100, 'SampleTime', Ts);

figure;
freq = r30.Frequency / Ts;
semilogx(freq, 20*log10(abs(r10.Response)),  'b', 'DisplayName', 'M = 10');
hold on;
semilogx(freq, 20*log10(abs(r30.Response)),  'r', 'DisplayName', 'M = 30');
semilogx(freq, 20*log10(abs(r100.Response)), 'g', 'DisplayName', 'M = 100');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Effect of Window Size on Frequency Resolution');
legend('show');
grid on;
hold off;

%% Time series mode (no input)
% Estimate the output power spectrum of an AR(1) process.
y_ts = filter(1, [1 -0.8], randn(500, 1));
result_ts = sidFreqBT(y_ts, []);

figure;
sidSpectrumPlot(result_ts);
title('AR(1) Output Spectrum');
