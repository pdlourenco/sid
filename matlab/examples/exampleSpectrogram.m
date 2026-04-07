%% exampleSpectrogram - Short-time FFT spectrogram with sidSpectrogram
%
% This example demonstrates sidSpectrogram for time-frequency analysis of
% signals. It computes the STFT and visualizes how frequency content
% evolves over time.

runner__nCompleted = 0;

%% Chirp signal: frequency sweep
% A signal whose instantaneous frequency increases linearly from 50 Hz
% to 150 Hz over 5 seconds.

Fs = 1000;
Ts = 1/Fs;
N = 5000;
t = (0:N-1)' * Ts;
f0 = 50;  f1 = 150;
x_chirp = cos(2*pi * (f0 + (f1-f0)/(2*max(t)) * t) .* t);

result = sidSpectrogram(x_chirp, 'WindowLength', 256, 'SampleTime', Ts);

figure;
sidSpectrogramPlot(result);
title('Chirp Signal: 50 Hz to 150 Hz');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Chirp signal: frequency sweep.\n', runner__nCompleted);

%% Window length trade-off
% Short window: good time resolution, poor frequency resolution.
% Long window: poor time resolution, good frequency resolution.

r_short = sidSpectrogram(x_chirp, 'WindowLength', 64,  'SampleTime', Ts);
r_long  = sidSpectrogram(x_chirp, 'WindowLength', 512, 'SampleTime', Ts);

figure;
subplot(1,2,1);
sidSpectrogramPlot(r_short);
title('Short Window (L=64)');

subplot(1,2,2);
sidSpectrogramPlot(r_long);
title('Long Window (L=512)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Window length trade-off.\n', runner__nCompleted);

%% Window types
% Different windows affect the trade-off between main lobe width and
% side lobe suppression.

r_hann = sidSpectrogram(x_chirp, 'WindowLength', 256, 'Window', 'hann', ...
    'SampleTime', Ts);
r_hamm = sidSpectrogram(x_chirp, 'WindowLength', 256, 'Window', 'hamming', ...
    'SampleTime', Ts);
r_rect = sidSpectrogram(x_chirp, 'WindowLength', 256, 'Window', 'rect', ...
    'SampleTime', Ts);

figure;
subplot(1,3,1);
sidSpectrogramPlot(r_hann);
title('Hann Window');

subplot(1,3,2);
sidSpectrogramPlot(r_hamm);
title('Hamming Window');

subplot(1,3,3);
sidSpectrogramPlot(r_rect);
title('Rectangular Window');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Window types.\n', runner__nCompleted);

%% Multi-channel signal
% Two channels: channel 1 has a chirp, channel 2 has a fixed sinusoid.

x_fixed = cos(2*pi * 100 * t);   % 100 Hz constant tone
x_mc = [x_chirp, x_fixed];

result_mc = sidSpectrogram(x_mc, 'WindowLength', 256, 'SampleTime', Ts);

figure;
subplot(1,2,1);
sidSpectrogramPlot(result_mc, 'Channel', 1);
title('Channel 1: Chirp');

subplot(1,2,2);
sidSpectrogramPlot(result_mc, 'Channel', 2);
title('Channel 2: Fixed 100 Hz');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Multi-channel signal.\n', runner__nCompleted);

%% Log frequency scale and NFFT zero-padding

r_zp = sidSpectrogram(x_chirp, 'WindowLength', 128, 'NFFT', 1024, ...
    'SampleTime', Ts);

figure;
sidSpectrogramPlot(r_zp, 'FrequencyScale', 'log');
title('Log Frequency Scale with Zero-Padding (NFFT=1024)');

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Log frequency scale and NFFT zero-padding.\n', runner__nCompleted);

%% Accessing raw STFT data
% result.Power contains the one-sided PSD, result.PowerDB is in dB,
% and result.Complex holds the STFT coefficients.

fprintf('Spectrogram dimensions:\n');
fprintf('  Time points:      %d\n', length(result.Time));
fprintf('  Frequency bins:   %d\n', length(result.Frequency));
fprintf('  Power size:       [%s]\n', num2str(size(result.Power)));
fprintf('  Complex size:     [%s]\n', num2str(size(result.Complex)));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Accessing raw STFT data.\n', runner__nCompleted);

fprintf('exampleSpectrogram: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
