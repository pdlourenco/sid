%% exampleFreqMap - Time-varying frequency response maps with sidFreqMap
%
% This example demonstrates sidFreqMap, which estimates how the frequency
% response evolves over time by applying spectral analysis to overlapping
% segments. Useful for detecting time-varying dynamics (e.g., drifting
% poles, changing gains).

runner__nCompleted = 0;

%% LTI baseline: constant system
% For a time-invariant system, the frequency map should be constant along
% the time axis.

rng(10);
N = 4000;
u = randn(N, 1);
y = filter(1, [1 -0.9], u) + 0.1 * randn(N, 1);

result_lti = sidFreqMap(y, u, 'SegmentLength', 512);

figure;
sidMapPlot(result_lti, 'PlotType', 'magnitude');
title('LTI System: Magnitude Should Be Constant Along Time');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: LTI baseline: constant system.\n', runner__nCompleted);

%% Time-varying system: drifting pole
% Simulate a first-order system x(k+1) = a(k)*x(k) + u(k) where the
% pole a(k) ramps from 0.5 to 0.95 over time. This creates a system that
% becomes more resonant (higher gain at low frequencies) as time progresses.

rng(20);
N = 6000;
u_tv = randn(N, 1);
a_k = linspace(0.5, 0.95, N)';    % time-varying pole

y_tv = zeros(N, 1);
for k = 2:N
    y_tv(k) = a_k(k) * y_tv(k-1) + u_tv(k);
end
y_tv = y_tv + 0.05 * randn(N, 1);  % small measurement noise

result_tv = sidFreqMap(y_tv, u_tv, 'SegmentLength', 512, 'Overlap', 384);

figure;
sidMapPlot(result_tv, 'PlotType', 'magnitude');
title('Time-Varying System: Pole Drifts from 0.5 to 0.95');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Time-varying system: drifting pole.\n', runner__nCompleted);

%% Coherence map
% Coherence shows how the signal-to-noise ratio evolves over time.

figure;
sidMapPlot(result_tv, 'PlotType', 'coherence');
title('Coherence Map');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Coherence map.\n', runner__nCompleted);

%% BT vs Welch algorithm
% sidFreqMap supports two algorithms:
%   'bt'    - Blackman-Tukey (default): correlogram-based, configurable lag window
%   'welch' - Welch averaged periodogram: sub-segment FFT averaging

result_bt    = sidFreqMap(y_tv, u_tv, 'SegmentLength', 512, 'Algorithm', 'bt');
result_welch = sidFreqMap(y_tv, u_tv, 'SegmentLength', 512, 'Algorithm', 'welch');

figure;
subplot(1,2,1);
sidMapPlot(result_bt, 'PlotType', 'magnitude');
title('Blackman-Tukey Algorithm');

subplot(1,2,2);
sidMapPlot(result_welch, 'PlotType', 'magnitude');
title('Welch Algorithm');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: BT vs Welch algorithm.\n', runner__nCompleted);

%% Segment length and overlap tuning
% Shorter segments = better time resolution but worse frequency resolution.
% More overlap = smoother time axis but more computation.

result_short = sidFreqMap(y_tv, u_tv, 'SegmentLength', 256, 'Overlap', 192);
result_long  = sidFreqMap(y_tv, u_tv, 'SegmentLength', 1024, 'Overlap', 768);

figure;
subplot(1,2,1);
sidMapPlot(result_short, 'PlotType', 'magnitude');
title('Short Segments (L=256): Good Time Resolution');

subplot(1,2,2);
sidMapPlot(result_long, 'PlotType', 'magnitude');
title('Long Segments (L=1024): Good Freq Resolution');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Segment length and overlap tuning.\n', runner__nCompleted);

%% Time-series mode: evolving output spectrum
% With u=[], sidFreqMap estimates how the output spectrum changes over time.
% Here we simulate a signal whose spectral content drifts.

rng(30);
N = 4000;
a_ts = linspace(0.3, 0.9, N)';
y_ts = zeros(N, 1);
for k = 2:N
    y_ts(k) = a_ts(k) * y_ts(k-1) + randn(1);
end

result_ts = sidFreqMap(y_ts, [], 'SegmentLength', 512);

figure;
sidMapPlot(result_ts, 'PlotType', 'spectrum');
title('Time-Series: Output Spectrum Evolves as AR(1) Pole Drifts');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Time-series mode: evolving output spectrum');

fprintf('exampleFreqMap: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
