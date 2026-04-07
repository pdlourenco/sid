%% exampleMultiTrajectory - Multi-Trajectory Ensemble Averaging
%
% Demonstrates how multiple independent trajectories improve spectral
% estimation by ensemble averaging. The same system is measured L times
% with different inputs and noise — averaging reduces variance by factor L
% without sacrificing frequency resolution.
%
% Functions demonstrated: sidFreqBT, sidFreqMap, sidSpectrogram, sidLTVdisc

runner__nCompleted = 0;

%% 1. LTI Ensemble Averaging — Tighter Confidence Bands
% Generate L=10 trajectories of a known second-order system.
% Compare sidFreqBT with 1 vs 10 trajectories.

rng(5001);
N = 2000; L = 10; Ts = 1;
b_true = 1; a_true = [1, -2*0.9*cos(pi/4), 0.9^2];

y_all = zeros(N, 1, L);
u_all = zeros(N, 1, L);
for l = 1:L
    u_all(:, 1, l) = randn(N, 1);
    y_all(:, 1, l) = filter(b_true, a_true, u_all(:, 1, l)) + 0.3*randn(N, 1);
end

% Single trajectory
r1 = sidFreqBT(y_all(:, :, 1), u_all(:, :, 1), 'WindowSize', 40);

% Multi-trajectory ensemble (all L=10)
rL = sidFreqBT(y_all, u_all, 'WindowSize', 40);

fprintf('Single trajectory: max ResponseStd = %.4f\n', max(r1.ResponseStd));
fprintf('10 trajectories:   max ResponseStd = %.4f\n', max(rL.ResponseStd));
fprintf('Ratio: %.2f (expected ~%.2f = 1/sqrt(%d))\n', ...
    max(rL.ResponseStd)/max(r1.ResponseStd), 1/sqrt(L), L);

% Plot comparison
figure;
subplot(2, 1, 1);
sidBodePlot(r1, 'Confidence', 3);
title('Single Trajectory');
subplot(2, 1, 2);
sidBodePlot(rL, 'Confidence', 3);
title(sprintf('Ensemble of %d Trajectories', L));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, ...
    '1. LTI Ensemble Averaging — Tighter Confidence Bands');

%% 2. LTV Time-Varying Map — Sharper Transition Detection
% System with a step change at t=N/2: pole shifts from 0.85 to 0.5.

rng(5002);
N = 4000; L = 5;
y_ltv = zeros(N, 1, L);
u_ltv = zeros(N, 1, L);
for l = 1:L
    ul = randn(N, 1);
    u_ltv(:, 1, l) = ul;
    yl = zeros(N, 1);
    for t = 2:N
        if t <= N/2
            a_t = 0.85;
        else
            a_t = 0.5;
        end
        yl(t) = a_t * yl(t-1) + ul(t) + 0.2*randn;
    end
    y_ltv(:, 1, l) = yl;
end

% Single trajectory map
r1_map = sidFreqMap(y_ltv(:, :, 1), u_ltv(:, :, 1), 'SegmentLength', 500);

% Multi-trajectory map
rL_map = sidFreqMap(y_ltv, u_ltv, 'SegmentLength', 500);

figure;
subplot(1, 2, 1);
sidMapPlot(r1_map, 'PlotType', 'magnitude');
title('Single Trajectory');
subplot(1, 2, 2);
sidMapPlot(rL_map, 'PlotType', 'magnitude');
title(sprintf('Ensemble of %d Trajectories', L));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, ...
    '2. LTV Time-Varying Map — Sharper Transition Detection');

%% 3. Spectrogram Averaging — Chirp in Noise
% A chirp signal buried in noise. Single trajectory: barely visible.
% 8 trajectories averaged: chirp track emerges clearly.

rng(5003);
N = 2000; L = 8; Fs = 1000; Ts = 1/Fs;
t = (0:N-1)'/Fs;
f0 = 50; f1 = 400;
chirp_sig = sin(2*pi * (f0*t + (f1-f0)/(2*t(end))*t.^2));

x_all = zeros(N, 1, L);
for l = 1:L
    x_all(:, 1, l) = 0.3*chirp_sig + randn(N, 1);
end

r1_spec = sidSpectrogram(x_all(:, :, 1), 'WindowLength', 128, 'SampleTime', Ts);
rL_spec = sidSpectrogram(x_all, 'WindowLength', 128, 'SampleTime', Ts);

figure;
subplot(1, 2, 1);
sidSpectrogramPlot(r1_spec);
title('Single Trajectory');
subplot(1, 2, 2);
sidSpectrogramPlot(rL_spec);
title(sprintf('Ensemble of %d Trajectories', L));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, ...
    '3. Spectrogram Averaging — Chirp in Noise');

%% 4. COSMIC + sidFreqMap Consistency
% Use the same multi-trajectory dataset for both COSMIC (state-space)
% and sidFreqMap (non-parametric frequency). Compare frozen transfer
% function from COSMIC against ensemble-averaged spectral map.

rng(5004);
p = 2; q = 1; N = 50; L_traj = 10;
A_true = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
sigma = 0.02;

X = zeros(N+1, p, L_traj); U = randn(N, q, L_traj);
for l = 1:L_traj
    X(1, :, l) = 0.1*randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + sigma*randn(1, p);
    end
end

% COSMIC identification
ltv = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true);
fprintf('\nCOSMIC identified A(1) =\n');
disp(ltv.A(:, :, 1));

% Frozen transfer function
frz = sidLTVdiscFrozen(ltv);

% Ensemble-averaged sidFreqMap (treating state x1 as output, u as input)
y_freq = reshape(X(1:N, 1, :), N, 1, L_traj);
u_freq = reshape(U, N, 1, L_traj);
fmap = sidFreqMap(y_freq, u_freq, 'SegmentLength', min(N, 30));

fprintf('sidFreqMap ensemble with %d trajectories: %d segments.\n', ...
    fmap.NumTrajectories, length(fmap.Time));
fprintf('Both COSMIC and sidFreqMap use the same %d-trajectory dataset.\n', L_traj);

fprintf('\nDone. Multi-trajectory examples completed.\n');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: 4. COSMIC + sidFreqMap Consistency.\n', runner__nCompleted);

fprintf('exampleMultiTrajectory: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
