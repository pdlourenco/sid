%% exampleCoherence - Coherence analysis and signal quality assessment
%
% This example demonstrates how to use the squared coherence gamma^2(w) to
% assess the quality of a frequency response estimate. Coherence near 1
% indicates a reliable estimate; near 0 means noise dominates.
% Coherence is only available for SISO systems estimated with sidFreqBT
% or sidFreqBTFDR.

runner__nCompleted = 0;

%% Generate data: ARMA system with colored noise
% G(z) = (1 + 0.5 z^{-1}) / (1 - 0.8 z^{-1})
% Noise is colored: v(t) = e(t) / (1 - 0.6 z^{-1}), so coherence varies
% across frequency.

rng(3);
N = 2000;
u = randn(N, 1);
y_clean = filter([1 0.5], [1 -0.8], u);
e = randn(N, 1);
v = 0.5 * filter(1, [1 -0.6], e);   % colored noise
y = y_clean + v;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Generate data: ARMA system with colored noise');

%% Estimate with sidFreqBT
result = sidFreqBT(y, u, 'WindowSize', 40);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Estimate with sidFreqBT.\n', runner__nCompleted);

%% Plot Bode magnitude and coherence together
w = result.Frequency;

figure;
subplot(2,1,1);
semilogx(w, 20*log10(abs(result.Response)), 'b');
ylabel('Magnitude (dB)');
title('Frequency Response Estimate');
grid on;

subplot(2,1,2);
semilogx(w, result.Coherence, 'b');
hold on;
plot(xlim, [0.9 0.9], 'g--', 'DisplayName', '\gamma^2 = 0.9');
plot(xlim, [0.5 0.5], 'r--', 'DisplayName', '\gamma^2 = 0.5');
hold off;
ylabel('Coherence \gamma^2');
xlabel('Frequency (rad/sample)');
title('Squared Coherence');
legend('show', 'Location', 'southwest');
ylim([0 1]);
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Plot Bode magnitude and coherence together');

%% Confidence bands reflect coherence
% Where coherence is high, uncertainty is low (narrow bands).
% The 'Confidence' option controls the number of standard deviations.

figure;
subplot(1,2,1);
sidBodePlot(result, 'Confidence', 2);
title('2\sigma Confidence Bands');

subplot(1,2,2);
sidBodePlot(result, 'Confidence', 3);
title('3\sigma Confidence Bands');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Confidence bands reflect coherence.\n', runner__nCompleted);

%% High-noise vs low-noise comparison
% More noise means lower coherence across all frequencies.

rng(3);
y_low  = y_clean + 0.1 * filter(1, [1 -0.6], randn(N, 1));
y_high = y_clean + 2.0 * filter(1, [1 -0.6], randn(N, 1));

r_low  = sidFreqBT(y_low,  u, 'WindowSize', 40);
r_high = sidFreqBT(y_high, u, 'WindowSize', 40);

figure;
semilogx(w, r_low.Coherence,  'b', 'DisplayName', 'Low noise (\sigma=0.1)');
hold on;
semilogx(w, r_high.Coherence, 'r', 'DisplayName', 'High noise (\sigma=2.0)');
hold off;
ylabel('Coherence \gamma^2');
xlabel('Frequency (rad/sample)');
title('Effect of Noise Level on Coherence');
legend('show', 'Location', 'southwest');
ylim([0 1]);
grid on;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: High-noise vs low-noise comparison.\n', runner__nCompleted);

%% Note: ETFE does not provide coherence
% sidFreqETFE returns Coherence = [] because there is no windowed
% cross-spectral estimate in the FFT-ratio approach.

result_etfe = sidFreqETFE(y, u);
fprintf('ETFE Coherence is empty: %d\n', isempty(result_etfe.Coherence));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Note: ETFE does not provide coherence.\n', runner__nCompleted);

fprintf('exampleCoherence: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
