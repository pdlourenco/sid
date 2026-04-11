%% test_sidWindowedDFT - Unit tests for windowed Fourier transform of covariances
%
% Tests sidWindowedDFT(R, W, freqs, useFFT) for both FFT and direct paths,
% consistency between them, and spectral properties.

fprintf('Running test_sidWindowedDFT...\n');
runner__nPassed = 0;

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
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: FFT vs direct path consistency (scalar).\n');

%% Test 2: Output is real for auto-covariance (scalar)
% For real x, R_xx is real and symmetric, so Phi_xx should be real
assert(max(abs(imag(Phi_fft))) < 1e-10, 'Auto-spectrum should be real-valued');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: output is real for auto-covariance (scalar).\n');

%% Test 3: Auto-spectrum is non-negative
% Windowed periodogram of auto-covariance should be non-negative
assert(all(real(Phi_fft) > -1e-10), 'Auto-spectrum should be non-negative');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: auto-spectrum is non-negative.\n');

%% Test 4: Output size for scalar signal
assert(isequal(size(Phi_fft), [128, 1]), 'Scalar output should be (nf x 1)');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: output size for scalar signal.\n');

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
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: output size for matrix signal.\n');

%% Test 6: FFT vs direct for matrix signal
freqs128 = (1:128)' * pi / 128;
R2 = sidCov(x, x, M);  % (M+1 x 2 x 2)
Phi_fft2 = sidWindowedDFT(R2, W, freqs128, true);
Phi_dir2 = sidWindowedDFT(R2, W, freqs128, false);
relErr2 = max(abs(Phi_fft2(:) - Phi_dir2(:))) / max(abs(Phi_dir2(:)));
assert(relErr2 < 1e-10, 'FFT and direct should agree for matrix (relErr=%.2e)', relErr2);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: FFT vs direct for matrix signal.\n');

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
assert(max(abs(Phi - 1)) < 0.25, 'White noise spectrum should be ~1 (flat)');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: white noise spectrum should be approximately flat.\n');

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
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: custom frequencies (non-default grid).\n');

%% Test 9: FFT vs direct path for large M, scalar auto-covariance
% Regression test for the L = 2*nf FFT-length hardcode (SPEC.md §2.5.1).
% The previous implementation chose L = 256 irrespective of M, so for
% M >= nf = 128 the positive lags s(2..M+1) and the wrapped negative
% lags s(L..L-M+1) collided in the buffer. For M >= 2*nf = 256, the
% MATLAB assignment auto-extended s beyond length L and the subsequent
% fft(s, L) truncated it back, silently corrupting the output. Covers
% the boundary M=127 (which was correct), the silent-overlap region
% M in [128, 255], and the truncation region M >= 256.
rng(42);
M_list = [127, 128, 150, 200, 255, 256, 300, 500];
freqs128 = (1:128)' * pi / 128;
for ii = 1:length(M_list)
    M_i = M_list(ii);
    N_i = 2 * M_i + 50;  % comfortably > 2M so freq_bt would not clamp
    x = randn(N_i, 1);
    R = sidCov(x, x, M_i);
    W = sidHannWin(M_i);
    Phi_fft = sidWindowedDFT(R, W, freqs128, true);
    Phi_dir = sidWindowedDFT(R, W, freqs128, false);
    relErr = max(abs(Phi_fft - Phi_dir)) / max(abs(Phi_dir));
    assert(relErr < 1e-10, ...
        'M=%d (scalar auto): FFT vs direct relErr=%.2e', M_i, relErr);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: FFT vs direct for large M (scalar auto-covariance).\n');

%% Test 10: FFT vs direct path for large M, scalar cross-covariance
% Same regression envelope with a cross-covariance input. The collision
% is more visible here than for auto-covariance at M=128 because
% R_yu(tau) != conj(R_yu(tau)) in general, so the overlap write is
% not harmless even when W(M+1) = 0.
rng(7);
M_list_cross = [128, 200, 500];
for ii = 1:length(M_list_cross)
    M_i = M_list_cross(ii);
    N_i = 2 * M_i + 50;
    y = randn(N_i, 1);
    u = randn(N_i, 1);
    Ryu = sidCov(y, u, M_i);
    Ruy = sidCov(u, y, M_i);
    W = sidHannWin(M_i);
    Phi_fft = sidWindowedDFT(Ryu, W, freqs128, true,  Ruy);
    Phi_dir = sidWindowedDFT(Ryu, W, freqs128, false, Ruy);
    relErr = max(abs(Phi_fft - Phi_dir)) / max(abs(Phi_dir));
    assert(relErr < 1e-10, ...
        'M=%d (cross): FFT vs direct relErr=%.2e', M_i, relErr);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: FFT vs direct for large M (cross-covariance).\n');

%% Test 11: FFT vs direct path for large M, matrix auto-covariance
rng(11);
M_list_mat = [128, 200, 500];
for ii = 1:length(M_list_mat)
    M_i = M_list_mat(ii);
    N_i = 2 * M_i + 50;
    x = randn(N_i, 2);
    R = sidCov(x, x, M_i);  % (M+1 x 2 x 2)
    W = sidHannWin(M_i);
    Phi_fft = sidWindowedDFT(R, W, freqs128, true);
    Phi_dir = sidWindowedDFT(R, W, freqs128, false);
    relErr = max(abs(Phi_fft(:) - Phi_dir(:))) / max(abs(Phi_dir(:)));
    assert(relErr < 1e-10, ...
        'M=%d (matrix): FFT vs direct relErr=%.2e', M_i, relErr);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: FFT vs direct for large M (matrix auto-covariance).\n');

fprintf('test_sidWindowedDFT: %d/%d passed\n', runner__nPassed, runner__nPassed);
