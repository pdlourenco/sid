%% exampleSISO - Basic frequency response estimation with sidFreqBT
%
% This example demonstrates how to estimate the frequency response of a
% simple SISO system and plot the results with confidence bands.

runner__nCompleted = 0;

%% Generate test data
% True system: G(z) = 1 / (1 - 0.9 z^{-1})
% This is a stable first-order system with a pole at z = 0.9.

N = 1000;                          % Number of samples
Ts = 0.01;                         % Sample time (seconds)
u = randn(N, 1);                   % White noise input
y_clean = filter(1, [1 -0.9], u);  % Noiseless output
noise = 0.1 * randn(N, 1);         % Measurement noise
y = y_clean + noise;                % Noisy output

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Estimate frequency response using Blackman-Tukey
result = sidFreqBT(y, u, 'SampleTime', Ts);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Estimate frequency response using Blackman-Tukey');

%% Plot Bode diagram
figure;
sidBodePlot(result);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Plot Bode diagram.\n', runner__nCompleted);

%% Plot noise spectrum
figure;
sidSpectrumPlot(result);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Plot noise spectrum.\n', runner__nCompleted);

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
legend;
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Compare different window sizes.\n', runner__nCompleted);

%% Preprocessing: detrend data before estimation
% Add a drift to the data and show that detrending improves results.
y_drift = y + 0.01 * (1:N)';  % add linear drift
u_drift = u + 5;              % add DC offset to input

% Without detrending: drift biases the low-frequency estimate
result_raw = sidFreqBT(y_drift, u_drift, 'SampleTime', Ts);

% With detrending
y_dt = sidDetrend(y_drift);
u_dt = sidDetrend(u_drift);
result_dt = sidFreqBT(y_dt, u_dt, 'SampleTime', Ts);

fprintf('Without detrend: max |G| at low freq = %.2f\n', max(abs(result_raw.Response)));
fprintf('With detrend:    max |G| at low freq = %.2f\n', max(abs(result_dt.Response)));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Preprocessing: detrend data before estimation');

%% Model validation: residual analysis
resid = sidResidual(result, y, u);
if resid.WhitenessPass
    fprintf('Whiteness test: PASS\n');
else
    fprintf('Whiteness test: FAIL\n');
end
if resid.IndependencePass
    fprintf('Independence test: PASS\n');
else
    fprintf('Independence test: FAIL\n');
end

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Model validation: residual analysis.\n', runner__nCompleted);

%% Model validation: compare predicted vs measured
comp = sidCompare(result, y, u);
fprintf('NRMSE fit: %.1f%%\n', comp.Fit);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Model validation: compare predicted vs measured');

%% Time series mode (no input)
% Estimate the output power spectrum of an AR(1) process.
y_ts = filter(1, [1 -0.8], randn(500, 1));
result_ts = sidFreqBT(y_ts, []);

figure;
sidSpectrumPlot(result_ts);
title('AR(1) Output Spectrum');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Time series mode (no input).\n', runner__nCompleted);

fprintf('exampleSISO: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
