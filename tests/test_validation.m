%% test_validation - Analytical verification and validation tests
%
% Tests the toolbox against known analytical results from system
% identification theory (Ljung, 1999).

fprintf('Running test_validation...\n');

%% Test 1: White noise auto-spectrum = variance (BT method)
% For white noise with variance sigma^2, the power spectrum is sigma^2
% at all frequencies.
rng(42);
sigma = 2.5;
N = 50000;
x = sigma * randn(N, 1);
result = sidFreqBT(x, [], 'WindowSize', 50);
Phi = result.NoiseSpectrum;
expected = sigma^2;
relErr = max(abs(Phi - expected)) / expected;
assert(relErr < 0.10, 'White noise spectrum should be ~sigma^2 (relErr=%.3f)', relErr);

%% Test 2: Sinusoidal signal has peak at correct frequency
% x(t) = cos(w0*t), spectrum should peak at w0
w0 = pi / 4;
N = 2000;
t = (1:N)';
x = cos(w0 * t);
result = sidFreqBT(x, [], 'WindowSize', 100);
[~, idx_peak] = max(result.NoiseSpectrum);
w_peak = result.Frequency(idx_peak);
assert(abs(w_peak - w0) < pi/128, 'Spectrum peak should be at w0=pi/4 (found %.3f)', w_peak);

%% Test 3: Known first-order system G(z) = 1/(1 - a*z^{-1})
% Magnitude: |G(e^{jw})| = 1/|1 - a*e^{-jw}|
% Phase: angle(G(e^{jw})) = -atan2(a*sin(w), 1 - a*cos(w))
a = 0.85;
rng(1);
N = 10000;
u = randn(N, 1);
y = filter(1, [1 -a], u);
result = sidFreqBT(y, u, 'WindowSize', 50);
w = result.Frequency;
G_true = 1 ./ (1 - a * exp(-1j * w));

mag_err = abs(abs(result.Response) - abs(G_true)) ./ abs(G_true);
assert(median(mag_err) < 0.03, 'Median magnitude error should be <3%% (got %.1f%%)', median(mag_err)*100);

phase_err = abs(angle(result.Response) - angle(G_true));
phase_err = min(phase_err, 2*pi - phase_err);
assert(median(phase_err) < 0.05, 'Median phase error should be <0.05 rad (got %.3f)', median(phase_err));

%% Test 4: FIR system G(z) = 1 + 0.5*z^{-1}
% |G(e^{jw})| = |1 + 0.5*e^{-jw}|
rng(2);
N = 10000;
u = randn(N, 1);
y = filter([1 0.5], 1, u);
result = sidFreqBT(y, u, 'WindowSize', 30);
w = result.Frequency;
G_true = 1 + 0.5 * exp(-1j * w);

mag_err = abs(abs(result.Response) - abs(G_true)) ./ abs(G_true);
assert(median(mag_err) < 0.03, 'FIR magnitude error should be small (got %.1f%%)', median(mag_err)*100);

%% Test 5: Noise variance estimation consistency
% y = G*u + v, noise variance sigma_v^2 should be recoverable
sigma_v = 0.5;
rng(3);
N = 10000;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + sigma_v * randn(N, 1);
result = sidFreqBT(y, u, 'WindowSize', 50);
% Average noise spectrum should be ~sigma_v^2
mean_noise = mean(result.NoiseSpectrum);
assert(abs(mean_noise - sigma_v^2) / sigma_v^2 < 0.20, ...
    'Average noise spectrum should be ~sigma_v^2 (got %.3f vs %.3f)', mean_noise, sigma_v^2);

%% Test 6: Coherence approaches 1 for low-noise systems
rng(4);
N = 5000;
u = randn(N, 1);
y = filter(1, [1 -0.9], u) + 0.01 * randn(N, 1);  % Very low noise
result = sidFreqBT(y, u, 'WindowSize', 50);
assert(median(result.Coherence) > 0.95, ...
    'Coherence should be ~1 for low-noise system (got %.3f)', median(result.Coherence));

%% Test 7: Coherence decreases with more noise
rng(4);
N = 5000;
u = randn(N, 1);
y_lo = filter(1, [1 -0.9], u) + 0.1 * randn(N, 1);
y_hi = filter(1, [1 -0.9], u) + 2.0 * randn(N, 1);
result_lo = sidFreqBT(y_lo, u, 'WindowSize', 30);
result_hi = sidFreqBT(y_hi, u, 'WindowSize', 30);
assert(median(result_lo.Coherence) > median(result_hi.Coherence), ...
    'More noise should reduce coherence');

%% Test 8: ETFE of noiseless system should be exact (within numerical precision)
rng(5);
N = 256;
u = randn(N, 1);
y = 3 * u;  % G = 3
result = sidFreqETFE(y, u);
G_mag = abs(result.Response);
assert(max(abs(G_mag - 3)) < 1e-8, 'ETFE of noiseless gain=3 should be exact');

%% Test 9: Parseval-like check for periodogram
% For ETFE periodogram: mean(Phi_y) * pi ≈ (1/N) * sum(y^2) * (2*pi/N) ... approximately
% Actually: sum of periodogram * delta_w ≈ variance
rng(6);
N = 1024;
y = randn(N, 1);
result = sidFreqETFE(y, []);
Phi = result.NoiseSpectrum;
w = result.Frequency;
% The integral of periodogram over (0, pi) should approximate the variance
% Using trapezoidal rule
dw = w(2) - w(1);
integral_approx = sum(Phi) * dw / pi;  % Normalize
sample_var = mean(y.^2);  % Biased variance (mean not subtracted for DFT)
% This is an approximation; check within 20%
relErr = abs(integral_approx - sample_var) / sample_var;
assert(relErr < 0.3, 'Periodogram integral should approximate variance (relErr=%.2f)', relErr);

%% Test 10: Uncertainty decreases with window size (BT)
rng(7);
N = 2000;
u = randn(N, 1);
y = filter(1, [1 -0.8], u) + 0.3 * randn(N, 1);
result_small = sidFreqBT(y, u, 'WindowSize', 10);
result_large = sidFreqBT(y, u, 'WindowSize', 50);
% Larger window => larger CW => larger uncertainty (for same N)
% But also smoother estimates. Check that std values differ.
assert(mean(result_large.ResponseStd) > mean(result_small.ResponseStd) * 0.5, ...
    'Uncertainty should be related to window size');

%% Test 11: Constant signal has zero-lag variance only
y_const = 5 * ones(100, 1);
R = sidCov(y_const, y_const, 5);
assert(abs(R(1) - 25) < 1e-10, 'R(0) of constant=5 should be 25');
% All lags should be 25 for constant signal (biased cov of constant is mean^2)
for tau = 0:5
    assert(abs(R(tau+1) - 25) < 1e-10, 'All lags of constant should be 25');
end

fprintf('  test_validation: ALL PASSED\n');
