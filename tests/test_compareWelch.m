%% test_compareWelch - Compare sidFreqMap Welch against MathWorks functions
%
% Validates sidFreqMap with 'Algorithm', 'welch' against the Signal
% Processing Toolbox functions tfestimate, mscohere, cpsd, and pwelch.
%
% Requires: Signal Processing Toolbox.
% Skips gracefully if not available.

fprintf('Running test_compareWelch...\n');

if ~exist('tfestimate', 'file')
    fprintf('  test_compareWelch: SKIPPED (requires Signal Processing Toolbox)\n');
    return;
end

%% Test 1: SISO first-order system - tfestimate response comparison
rng(100);
N = 4000; Fs = 1; Ts = 1 / Fs;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);

Lsub = 256; Psub = 128; nfft = 512;
win = hann(Lsub);

% MathWorks: tfestimate(input, output, ...)
[Txy_mw, F_mw] = tfestimate(u, y, win, Psub, nfft, Fs);

% sid: single outer segment covering all data
result = sidFreqMap(y, u, 'Algorithm', 'welch', 'SegmentLength', N, ...
    'SubSegmentLength', Lsub, 'SubOverlap', Psub, 'NFFT', nfft, ...
    'SampleTime', Ts);

% Align frequencies: MathWorks includes DC (F_mw(1)=0), sid skips DC
% MathWorks returns nfft/2+1 bins [0, df, ..., Fs/2]
% sid returns nfft/2 bins [df, 2*df, ..., Fs/2] in rad/sample
Txy_mw_noDC = Txy_mw(2:end);  % skip DC
G_sid = result.Response(:, 1);  % single segment, column 1

% Both should have nfft/2 entries
assert(length(G_sid) == length(Txy_mw_noDC), ...
    'Frequency bin count mismatch: sid=%d, mw=%d', length(G_sid), length(Txy_mw_noDC));

relErr = abs(G_sid - Txy_mw_noDC) ./ max(abs(Txy_mw_noDC), eps);
assert(median(relErr) < 0.02, ...
    'SISO tfestimate median relative error too large: %.4f', median(relErr));
fprintf('  Test 1 passed: tfestimate response (median relErr=%.6f)\n', median(relErr));

%% Test 2: SISO mscohere coherence comparison
[Cxy_mw, ~] = mscohere(u, y, win, Psub, nfft, Fs);
Cxy_mw_noDC = Cxy_mw(2:end);
Coh_sid = result.Coherence(:, 1);

absErr = abs(Coh_sid - Cxy_mw_noDC);
assert(median(absErr) < 0.02, ...
    'SISO mscohere median absolute error too large: %.4f', median(absErr));
fprintf('  Test 2 passed: mscohere coherence (median absErr=%.6f)\n', median(absErr));

%% Test 3: SISO cpsd cross-spectrum comparison
[Pyu_mw, ~] = cpsd(y, u, win, Psub, nfft, Fs);
Pyu_mw_noDC = Pyu_mw(2:end);

% sid uses raw periodogram normalization (no Ts factor)
% MathWorks cpsd includes Ts (PSD convention: V^2/Hz)
% To compare: sid_cross * Ts should ≈ MathWorks cpsd
% But sidFreqMap doesn't expose the raw cross-spectrum directly.
% Instead, reconstruct it: Pyu_sid = G_sid .* Puu_sid
% and Puu_sid can be obtained from a time-series run or from cpsd(u,u)
[Puu_mw, ~] = cpsd(u, u, win, Psub, nfft, Fs);
Puu_mw_noDC = Puu_mw(2:end);

% Reconstruct cross-spectrum from sid: Pyu = G * Puu
% Since G_sid ≈ Txy_mw = Pyu_mw/Puu_mw, we have G_sid * Puu_mw ≈ Pyu_mw
% This is a consistency check rather than direct PSD comparison
Pyu_reconstructed = G_sid .* Puu_mw_noDC;
relErr_cross = abs(Pyu_reconstructed - Pyu_mw_noDC) ./ max(abs(Pyu_mw_noDC), eps);
assert(median(relErr_cross) < 0.05, ...
    'Cross-spectrum reconstruction median relErr too large: %.4f', median(relErr_cross));
fprintf('  Test 3 passed: cpsd cross-spectrum (median relErr=%.6f)\n', median(relErr_cross));

%% Test 4: Time-series pwelch PSD comparison
rng(200);
x_ts = randn(2000, 1) + 0.5 * sin(2 * pi * 0.1 * (1:2000)');
Lsub4 = 256; Psub4 = 128; nfft4 = 256;
win4 = hann(Lsub4);

[Pxx_mw, F_mw4] = pwelch(x_ts, win4, Psub4, nfft4, Fs);
Pxx_mw_noDC = Pxx_mw(2:end);

result_ts = sidFreqMap(x_ts, [], 'Algorithm', 'welch', 'SegmentLength', length(x_ts), ...
    'SubSegmentLength', Lsub4, 'SubOverlap', Psub4, 'NFFT', nfft4, ...
    'SampleTime', Ts);
Pxx_sid = result_ts.NoiseSpectrum(:, 1);

% MathWorks pwelch uses PSD convention (V^2/Hz = V^2*Ts)
% sid uses raw: Pxx_raw = (1/J/S1) * sum(|X|^2)
% To match: Pxx_sid * Ts ≈ Pxx_mw
Pxx_sid_scaled = Pxx_sid * Ts;
relErr_psd = abs(Pxx_sid_scaled - Pxx_mw_noDC) ./ max(abs(Pxx_mw_noDC), eps);
assert(median(relErr_psd) < 0.05, ...
    'Time-series pwelch median relErr too large: %.4f', median(relErr_psd));
fprintf('  Test 4 passed: pwelch PSD (median relErr=%.6f)\n', median(relErr_psd));

%% Test 5: Custom sub-segment parameters
rng(300);
N5 = 4000;
u5 = randn(N5, 1);
y5 = filter([1 0.5], [1 -0.8], u5) + 0.05 * randn(N5, 1);

Lsub5 = 128; Psub5 = 64; nfft5 = 256;
win5 = hann(Lsub5);

[Txy5_mw, ~] = tfestimate(u5, y5, win5, Psub5, nfft5, Fs);
Txy5_mw_noDC = Txy5_mw(2:end);

result5 = sidFreqMap(y5, u5, 'Algorithm', 'welch', 'SegmentLength', N5, ...
    'SubSegmentLength', Lsub5, 'SubOverlap', Psub5, 'NFFT', nfft5, ...
    'SampleTime', Ts);
G5_sid = result5.Response(:, 1);

relErr5 = abs(G5_sid - Txy5_mw_noDC) ./ max(abs(Txy5_mw_noDC), eps);
assert(median(relErr5) < 0.02, ...
    'Custom params tfestimate median relErr too large: %.4f', median(relErr5));
fprintf('  Test 5 passed: custom sub-segment params (median relErr=%.6f)\n', median(relErr5));

%% Test 6: MIMO (2-output, 1-input)
rng(400);
N6 = 4000;
u6 = randn(N6, 1);
y6_1 = filter([1], [1 -0.9], u6) + 0.1 * randn(N6, 1);
y6_2 = filter([0.5 0.3], [1 -0.7 0.2], u6) + 0.1 * randn(N6, 1);
y6 = [y6_1, y6_2];

Lsub6 = 256; Psub6 = 128; nfft6 = 512;
win6 = hann(Lsub6);

% MathWorks: separate tfestimate per output channel
[Txy6_1_mw, ~] = tfestimate(u6, y6_1, win6, Psub6, nfft6, Fs);
[Txy6_2_mw, ~] = tfestimate(u6, y6_2, win6, Psub6, nfft6, Fs);
Txy6_1_noDC = Txy6_1_mw(2:end);
Txy6_2_noDC = Txy6_2_mw(2:end);

result6 = sidFreqMap(y6, u6, 'Algorithm', 'welch', 'SegmentLength', N6, ...
    'SubSegmentLength', Lsub6, 'SubOverlap', Psub6, 'NFFT', nfft6, ...
    'SampleTime', Ts);

% result6.Response is (nf x 1 x 2 x 1) for single segment
G6_ch1 = squeeze(result6.Response(:, 1, 1, 1));
G6_ch2 = squeeze(result6.Response(:, 1, 2, 1));

relErr6_1 = abs(G6_ch1 - Txy6_1_noDC) ./ max(abs(Txy6_1_noDC), eps);
relErr6_2 = abs(G6_ch2 - Txy6_2_noDC) ./ max(abs(Txy6_2_noDC), eps);
assert(median(relErr6_1) < 0.05, ...
    'MIMO ch1 median relErr too large: %.4f', median(relErr6_1));
assert(median(relErr6_2) < 0.05, ...
    'MIMO ch2 median relErr too large: %.4f', median(relErr6_2));
fprintf('  Test 6 passed: MIMO tfestimate (ch1=%.6f, ch2=%.6f)\n', ...
    median(relErr6_1), median(relErr6_2));

%% Test 7: Non-unit sample time (Ts=0.001, Fs=1000)
rng(500);
Fs7 = 1000; Ts7 = 1 / Fs7;
N7 = 4000;
u7 = randn(N7, 1);
y7 = filter([1], [1 -0.99], u7) + 0.1 * randn(N7, 1);

Lsub7 = 256; Psub7 = 128; nfft7 = 512;
win7 = hann(Lsub7);

[Txy7_mw, F7_mw] = tfestimate(u7, y7, win7, Psub7, nfft7, Fs7);
Txy7_mw_noDC = Txy7_mw(2:end);
F7_mw_noDC = F7_mw(2:end);

result7 = sidFreqMap(y7, u7, 'Algorithm', 'welch', 'SegmentLength', N7, ...
    'SubSegmentLength', Lsub7, 'SubOverlap', Psub7, 'NFFT', nfft7, ...
    'SampleTime', Ts7);
G7_sid = result7.Response(:, 1);

% Frequency axis check: sid FrequencyHz should match MathWorks F (after DC skip)
freqErr = max(abs(result7.FrequencyHz - F7_mw_noDC));
assert(freqErr < 1e-6, 'Frequency axis mismatch: %.6f Hz', freqErr);

% Response comparison
relErr7 = abs(G7_sid - Txy7_mw_noDC) ./ max(abs(Txy7_mw_noDC), eps);
assert(median(relErr7) < 0.02, ...
    'Non-unit Ts median relErr too large: %.4f', median(relErr7));
fprintf('  Test 7 passed: non-unit Ts (median relErr=%.6f, freqErr=%.2e)\n', ...
    median(relErr7), freqErr);

fprintf('test_compareWelch: ALL TESTS PASSED\n');
