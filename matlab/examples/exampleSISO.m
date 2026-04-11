%% exampleSISO - SISO frequency response on a 1-DoF spring-mass-damper
%
% This example identifies the frequency response of a physical
% single-degree-of-freedom spring-mass-damper oscillator from a noisy
% input/output record using sidFreqBT (Blackman-Tukey estimator).
%
% Plant: m = 1 kg, k = 100 N/m, c = 2 N.s/m (Plant A in spec/EXAMPLES.md).
% Natural frequency omega_n = sqrt(k/m) = 10 rad/s (~ 1.59 Hz),
% damping ratio zeta = c / (2 sqrt(k m)) = 0.1 (Q = 5).
%
% Input is a force; output is the mass position x1. See
% spec/EXAMPLES.md section 3.1 for the binding specification.

runner__nCompleted = 0;

%% Generate test data
% Build the SMD state-space model with util_msd and simulate under
% white-force excitation. Measured output is the position x1 plus
% additive sensor noise.

rng(42);

% ---- Physical plant: 1-DoF SMD (omega_n = 10 rad/s, zeta = 0.1) ----
m  = 1.0;    % kg
k  = 100.0;  % N/m
c  = 2.0;    % N.s/m
F  = 1.0;    % force on the single mass
Ts = 0.01;   % s (fs = 100 Hz, Nyquist = 50 Hz)
N  = 2048;   % number of samples

[Ad, Bd] = util_msd(m, k, c, F, Ts);   % Ad: 2x2, Bd: 2x1

% ---- Simulate the plant: x[k+1] = Ad x[k] + Bd u[k] ----
u = randn(N, 1);                       % white-force excitation
x = zeros(N + 1, 2);                   % state trajectory [pos; vel]
for step = 1:N
    x(step + 1, :) = (Ad * x(step, :)' + Bd * u(step))';
end
y_clean = x(2:end, 1);                 % measured position x1 (m)
y       = y_clean + 2e-4 * randn(N, 1);  % + sensor noise (~0.2 mm)

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Estimate frequency response using Blackman-Tukey
% Lightly damped (Q = 5) -> narrow resonance bandwidth of 2 rad/s.
% BT frequency resolution is ~ pi/M rad/sample, so M must be large enough
% to see the peak. We also pass a dense custom frequency grid so the peak
% is not undersampled by the default bin spacing.

w_grid = linspace(0.005, pi, 512)';
result = sidFreqBT(y, u, 'WindowSize', 200, 'Frequencies', w_grid, ...
                   'SampleTime', Ts);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Estimate frequency response using Blackman-Tukey');

%% Plot Bode diagram
% Magnitude should peak near omega_n = 10 rad/s and phase should drop
% by about pi through the resonance.

figure;
sidBodePlot(result);
title('Bode diagram (freq\_bt)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Plot Bode diagram.\n', runner__nCompleted);

%% Plot noise spectrum
figure;
sidSpectrumPlot(result);
title('Noise spectrum');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Plot noise spectrum.\n', runner__nCompleted);

%% Compare different window sizes
% For a narrow resonance the short window (M=50) smears the peak toward
% DC while the long window (M=300) resolves it within half a bin of the
% true omega_n.

r50  = sidFreqBT(y, u, 'WindowSize',  50, 'Frequencies', w_grid, ...
                 'SampleTime', Ts);
r100 = sidFreqBT(y, u, 'WindowSize', 100, 'Frequencies', w_grid, ...
                 'SampleTime', Ts);
r300 = sidFreqBT(y, u, 'WindowSize', 300, 'Frequencies', w_grid, ...
                 'SampleTime', Ts);

figure;
freq = r300.Frequency / Ts;   % rad/s
semilogx(freq, 20*log10(abs(r50.Response)),  'b', 'DisplayName', 'M = 50');
hold on;
semilogx(freq, 20*log10(abs(r100.Response)), 'r', 'DisplayName', 'M = 100');
semilogx(freq, 20*log10(abs(r300.Response)), 'g', 'DisplayName', 'M = 300');
yl_wc = ylim;
plot([10.0 10.0], yl_wc, 'k:', 'HandleVisibility', 'off');
hold off;
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Effect of window size on resonance resolution');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Compare different window sizes.\n', ...
    runner__nCompleted);

%% Preprocessing: detrend data before estimation
% Overlay a linear drift on the output and a DC offset on the input.
% sidDetrend removes the bias before BT estimation.

y_drift = y + 2e-4 * (1:N)';       % linear drift (~0.4 m total)
u_drift = u + 5.0;                  % DC offset on input

% Without detrending: drift biases the low-frequency estimate
result_raw = sidFreqBT(y_drift, u_drift, 'WindowSize', 200, ...
                        'Frequencies', w_grid, 'SampleTime', Ts);

% With detrending
y_dt = sidDetrend(y_drift);
u_dt = sidDetrend(u_drift);
result_dt = sidFreqBT(y_dt, u_dt, 'WindowSize', 200, ...
                       'Frequencies', w_grid, 'SampleTime', Ts);

fprintf('Without detrend: max |G| at low freq = %.4f\n', ...
    max(abs(result_raw.Response)));
fprintf('With detrend:    max |G| at low freq = %.4f\n', ...
    max(abs(result_dt.Response)));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', runner__nCompleted, ...
    'Preprocessing: detrend data before estimation');

%% Model validation: residual analysis
resid = sidResidual(result, y, u);
if resid.WhitenessPass
    fprintf('Whiteness test:    PASS\n');
else
    fprintf('Whiteness test:    FAIL\n');
end
if resid.IndependencePass
    fprintf('Independence test: PASS\n');
else
    fprintf('Independence test: FAIL\n');
end

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', runner__nCompleted, ...
    'Model validation: residual analysis');

%% Model validation: compare predicted vs measured
comp = sidCompare(result, y, u);
fprintf('NRMSE fit: %.1f%%\n', comp.Fit(1));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', runner__nCompleted, ...
    'Model validation: compare predicted vs measured');

%% Time-series mode (no input)
% Re-simulate the plant and hand only the position record to sidFreqBT.
% The resonance still appears in the output spectrum because the plant
% colours the white-force excitation.

u_ts = randn(1000, 1);
x_ts = zeros(1001, 2);
for step = 1:1000
    x_ts(step + 1, :) = (Ad * x_ts(step, :)' + Bd * u_ts(step))';
end
y_ts = x_ts(2:end, 1);

result_ts = sidFreqBT(y_ts, [], 'WindowSize', 200, 'Frequencies', w_grid);

figure;
sidSpectrumPlot(result_ts);
title('SDOF output spectrum (time-series mode)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Time-series mode (no input).\n', ...
    runner__nCompleted);

fprintf('exampleSISO: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
