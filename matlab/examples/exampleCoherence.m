%% exampleCoherence - Coherence analysis on a 2-mass SMD with colored disturbance
%
% Uses the squared coherence gamma^2(w) to assess frequency-local
% estimate quality. A commanded force is applied at mass 1 and an
% unmeasured colored disturbance enters at mass 2. Coherence drops
% near frequencies where the disturbance dominates.
%
% Plant B: m = [1 1], k = [100 80], c = [2 2]. Two input channels:
% commanded force on channel 0 (at mass 1) and colored disturbance
% on channel 1 (at mass 2).
%
% See spec/EXAMPLES.md section 3.4 for the binding specification.

runner__nCompleted = 0;

%% Generate test data
% F matrix has two columns: column 1 injects u at mass 1, column 2
% injects d at mass 2. util_msd returns a single (Ad, Bd) pair with
% Bd of shape (4, 2) handling both channels.

rng(3);

m  = [1; 1];
k  = [100; 80];
c  = [2; 2];
F  = [1 0; 0 1];     % col 1: u at m1; col 2: d at m2
Ts = 0.01;
N  = 4000;

[Ad, Bd] = util_msd(m, k, c, F, Ts);

% Excitation and disturbance
u = 10.0 * randn(N, 1);                      % commanded force
e = randn(N, 1);
d = 0.5 * filter(1, [1 -0.9], e);             % DC-heavy colored disturbance

x = zeros(N + 1, 4);
for step = 1:N
    x(step + 1, :) = (Ad * x(step, :)' ...
        + Bd(:, 1) * u(step) + Bd(:, 2) * d(step))';
end
y = x(2:end, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Generate test data.\n', runner__nCompleted);

%% Estimate with sidFreqBT
% Dense custom frequency grid to resolve the narrow modes. We only
% feed the commanded input u -- the disturbance is unmodeled.

w_grid = linspace(0.01, pi, 512)';
result = sidFreqBT(y, u, 'WindowSize', 200, 'Frequencies', w_grid, ...
                   'SampleTime', Ts);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Estimate with sidFreqBT.\n', ...
    runner__nCompleted);

%% Plot Bode magnitude and coherence together
% Coherence drops where unmodeled content masks the cause-effect
% relationship: near the low-frequency disturbance band and in the
% high-frequency tail where the plant gain is below the noise floor.

w = result.Frequency;
coh = result.Coherence;

figure;
subplot(2, 1, 1);
semilogx(w, 20*log10(abs(result.Response)), 'b');
ylabel('Magnitude (dB)');
title('Frequency response estimate');
grid on;

subplot(2, 1, 2);
semilogx(w, coh, 'b');
hold on;
semilogx(w, 0.9 * ones(size(w)), 'g--', 'DisplayName', '\gamma^2 = 0.9');
semilogx(w, 0.5 * ones(size(w)), 'r--', 'DisplayName', '\gamma^2 = 0.5');
hold off;
ylabel('Coherence \gamma^2');
xlabel('Frequency (rad/sample)');
ylim([0 1]);
title('Squared coherence');
legend('Location', 'southwest');
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Plot Bode magnitude and coherence together.\n', ...
    runner__nCompleted);

%% Confidence bands reflect coherence
% Where coherence is high the uncertainty is low (narrow bands).
% sidBodePlot needs [mag_ax, phase_ax] so we build a 2-column figure.

figure;
ax_m2 = subplot(2, 2, 1);
ax_p2 = subplot(2, 2, 3);
sidBodePlot(result, 'Confidence', 2, 'Axes', [ax_m2, ax_p2]);
title(ax_m2, '2\sigma confidence bands');

ax_m3 = subplot(2, 2, 2);
ax_p3 = subplot(2, 2, 4);
sidBodePlot(result, 'Confidence', 3, 'Axes', [ax_m3, ax_p3]);
title(ax_m3, '3\sigma confidence bands');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Confidence bands reflect coherence.\n', ...
    runner__nCompleted);

%% High-disturbance vs low-disturbance comparison
% Scale the disturbance 0.1x and 2.0x and overlay the two coherence
% curves.

rng(3);
e2 = randn(N, 1);
d_low  = 0.1 * filter(1, [1 -0.9], e2);
d_high = 2.0 * filter(1, [1 -0.9], e2);

x_low = zeros(N + 1, 4);
x_high = zeros(N + 1, 4);
for step = 1:N
    x_low (step + 1, :) = (Ad * x_low (step, :)' ...
        + Bd(:, 1) * u(step) + Bd(:, 2) * d_low(step))';
    x_high(step + 1, :) = (Ad * x_high(step, :)' ...
        + Bd(:, 1) * u(step) + Bd(:, 2) * d_high(step))';
end
y_low  = x_low (2:end, 1);
y_high = x_high(2:end, 1);

r_low  = sidFreqBT(y_low,  u, 'WindowSize', 200, 'Frequencies', w_grid, ...
                   'SampleTime', Ts);
r_high = sidFreqBT(y_high, u, 'WindowSize', 200, 'Frequencies', w_grid, ...
                   'SampleTime', Ts);

figure;
semilogx(w, r_low.Coherence,  'b', 'DisplayName', 'Low disturbance (0.1x)');
hold on;
semilogx(w, r_high.Coherence, 'r', 'DisplayName', 'High disturbance (2.0x)');
hold off;
ylabel('Coherence \gamma^2');
xlabel('Frequency (rad/sample)');
title('Effect of disturbance level on coherence');
legend('Location', 'southwest');
ylim([0 1]);
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: High-disturbance vs low-disturbance comparison.\n', ...
    runner__nCompleted);

%% Note: ETFE does not provide coherence
% sidFreqETFE returns Coherence = [] because there is no windowed
% cross-spectral estimate in the FFT-ratio approach.

result_etfe = sidFreqETFE(y, u);
fprintf('ETFE coherence is empty: %d\n', isempty(result_etfe.Coherence));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Note: ETFE does not provide coherence.\n', ...
    runner__nCompleted);

fprintf('exampleCoherence: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
