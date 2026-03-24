%% test_sidDFT - Unit tests for the DFT computation
%
% Tests sidDFT(x, freqs, useFFT) for correctness of both FFT and direct
% paths, consistency between paths, and known analytical DFTs.

fprintf('Running test_sidDFT...\n');

%% Test 1: Direct DFT of a known signal
% x = [1; 0; 0; 0] (impulse), DFT should be exp(-j*w*1) = exp(-jw) for all w
x = [1; 0; 0; 0];
freqs = [pi/4; pi/2; pi];
X = sidDFT(x, freqs, false);
expected = exp(-1j * freqs);
assert(max(abs(X - expected)) < 1e-12, 'DFT of impulse at t=1 should be exp(-jw)');

%% Test 2: DFT of constant signal
% x = ones(N, 1): X(w) = sum_{t=1}^{N} exp(-jwt) = exp(-jw) * (1-exp(-jwN))/(1-exp(-jw))
N = 16;
x = ones(N, 1);
w = pi/3;
X = sidDFT(x, w, false);
expected = sum(exp(-1j * w * (1:N)'));
assert(abs(X - expected) < 1e-10, 'DFT of constant should match direct sum');

%% Test 3: FFT vs direct DFT consistency on default 128-point grid
rng(42);
N = 200;
x = randn(N, 1);
freqs = (1:128)' * pi / 128;
X_fft = sidDFT(x, freqs, true);
X_direct = sidDFT(x, freqs, false);
% Tolerance needs to be moderate due to FFT binning/interpolation
relErr = max(abs(X_fft - X_direct)) / max(abs(X_direct));
assert(relErr < 0.05, 'FFT and direct DFT should agree on default grid (relErr=%.4f)', relErr);

%% Test 4: Multi-channel signal
rng(7);
N = 50;
x = randn(N, 3);
freqs = [0.5; 1.0; 2.0];
X = sidDFT(x, freqs, false);
assert(isequal(size(X), [3, 3]), 'Output should be (nf x p)');
% Verify first channel manually
X1 = zeros(3, 1);
t = (1:N)';
for k = 1:3
    X1(k) = sum(x(:,1) .* exp(-1j * freqs(k) * t));
end
assert(max(abs(X(:,1) - X1)) < 1e-10, 'Multi-channel DFT should match per-channel computation');

%% Test 5: Output dimensions
rng(1);
x = randn(100, 2);
freqs = (1:64)' * pi / 64;
X = sidDFT(x, freqs, false);
assert(isequal(size(X), [64, 2]), 'Output size should be (nf x p)');

%% Test 6: Single frequency
x = randn(50, 1);
freqs = pi/2;
X = sidDFT(x, freqs, false);
assert(isequal(size(X), [1, 1]), 'Single frequency should return scalar');
expected = sum(x .* exp(-1j * pi/2 * (1:50)'));
assert(abs(X - expected) < 1e-10, 'Single frequency DFT should match');

%% Test 7: Parseval-like energy check (direct path)
rng(88);
N = 64;
x = randn(N, 1);
% DFT at Fourier frequencies: w_k = 2*pi*k/N for k=0..N-1
% Parseval: (1/N) sum |X(w_k)|^2 = sum |x(t)|^2
freqs_full = (0:N-1)' * 2 * pi / N;
freqs_full = freqs_full(2:end);  % exclude DC since sidDFT expects (0,pi]
X = sidDFT(x, freqs_full, false);
% Add DC manually: X(0) = sum(x)
Xdc = sum(x);
energy_freq = (abs(Xdc)^2 + sum(abs(X).^2)) / N;
energy_time = sum(x.^2);
assert(abs(energy_freq - energy_time) < 1e-8, 'Parseval theorem should hold');

fprintf('  test_sidDFT: ALL PASSED\n');
