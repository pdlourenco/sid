%% test_sidWindowedDFT - Unit tests for windowed Fourier transform of covariances
%
% Tests sidWindowedDFT(R, W, freqs, useFFT) for both FFT and direct paths,
% consistency between them, and spectral properties.

fprintf('Running test_sidWindowedDFT...\n');

%% Test 1: FFT vs direct path consistency (scalar)
rng(42);
N = 500;
x = randn(N, 1);
M = 30;
R = sidCov(x, x, M);
W = sidHannWin(M);
freqs = (1:128)' * pi / 128;

Phi_fft = sidWindowedDFT(R, W, freqs, true);
Phi_direct = sidWindowedDFT(R, W, freqs, false);

relErr = max(abs(Phi_fft - Phi_direct)) / max(abs(Phi_direct));
assert(relErr < 1e-10, 'FFT and direct should agree for scalar (relErr=%.2e)', relErr);

%% Test 2: Output is real for auto-covariance (scalar)
% For real x, R_xx is real and symmetric, so Phi_xx should be real
assert(max(abs(imag(Phi_fft))) < 1e-10, 'Auto-spectrum should be real-valued');

%% Test 3: Auto-spectrum is non-negative
% Windowed periodogram of auto-covariance should be non-negative
assert(all(real(Phi_fft) > -1e-10), 'Auto-spectrum should be non-negative');

%% Test 4: Output size for scalar signal
assert(isequal(size(Phi_fft), [128, 1]), 'Scalar output should be (nf x 1)');

%% Test 5: Output size for matrix signal
rng(7);
N = 200;
x = randn(N, 2);
z = randn(N, 3);
M = 20;
R = sidCov(x, z, M);  % (M+1 x 2 x 3)
W = sidHannWin(M);
freqs = (1:64)' * pi / 64;
Phi = sidWindowedDFT(R, W, freqs, false);
assert(isequal(size(Phi), [64, 2, 3]), 'Matrix output should be (nf x p x q)');

%% Test 6: FFT vs direct for matrix signal
freqs128 = (1:128)' * pi / 128;
R2 = sidCov(x, x, M);  % (M+1 x 2 x 2)
Phi_fft2 = sidWindowedDFT(R2, W, freqs128, true);
Phi_dir2 = sidWindowedDFT(R2, W, freqs128, false);
relErr2 = max(abs(Phi_fft2(:) - Phi_dir2(:))) / max(abs(Phi_dir2(:)));
assert(relErr2 < 1e-10, 'FFT and direct should agree for matrix (relErr=%.2e)', relErr2);

%% Test 7: White noise spectrum should be approximately flat
rng(99);
N = 10000;
x = randn(N, 1);
M = 50;
R = sidCov(x, x, M);
W = sidHannWin(M);
freqs = (1:128)' * pi / 128;
Phi = real(sidWindowedDFT(R, W, freqs, true));
% For unit variance white noise, spectrum should be ~1 everywhere
assert(max(abs(Phi - 1)) < 0.15, 'White noise spectrum should be ~1 (flat)');

%% Test 8: Custom frequencies (non-default grid)
freqs_custom = [0.1; 0.5; 1.0; 2.0; 3.0];
rng(42);
x = randn(200, 1);
M = 20;
R = sidCov(x, x, M);
W = sidHannWin(M);
Phi = sidWindowedDFT(R, W, freqs_custom, false);
assert(isequal(size(Phi), [5, 1]), 'Custom freq output should be (5 x 1)');
assert(all(real(Phi) > -1e-10), 'Auto-spectrum at custom freqs should be non-negative');

fprintf('  test_sidWindowedDFT: ALL PASSED\n');
