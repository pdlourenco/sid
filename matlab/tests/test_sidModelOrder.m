%% test_sidModelOrder - Unit tests for model order estimation
%
% Tests sidModelOrder for known systems of various orders, SISO and MIMO,
% gap and threshold detection, output struct fields, and edge cases.

fprintf('Running test_sidModelOrder...\n');
runner__nPassed = 0;

%% Test 1: Known 2nd-order SISO system
% Generate data from a 2nd-order discrete system and verify n = 2.
rng(100);
N = 2000;
u = randn(N, 1);

% 2nd-order system: y(k) = G(z) u(k) where G(z) = (0.5z) / (z^2 - 1.2z + 0.5)
% Poles at 0.6 +/- 0.2i => stable
a1 = -1.2; a2 = 0.5; b1 = 0.5;
y = zeros(N, 1);
for k = 3:N
    y(k) = -a1 * y(k-1) - a2 * y(k-2) + b1 * u(k-1) + 0.01 * randn;
end

G = sidFreqBT(y, u, 'WindowSize', 60);
[n, sv] = sidModelOrder(G);

assert(n == 2, 'Expected n = 2 for 2nd-order system, got n = %d', n);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: 2nd-order SISO system detected (n = %d).\n', n);

%% Test 2: Known 4th-order SISO system
% Two resonant modes at the same radius but well-separated frequencies.
rng(200);
N = 5000;
u = randn(N, 1);

% Section 1: poles at 0.85*exp(+/-j*0.4)  (low-frequency resonance)
% Section 2: poles at 0.85*exp(+/-j*2.0)  (high-frequency resonance)
r_pole = 0.85;
a_poly = conv([1, -2*r_pole*cos(0.4), r_pole^2], ...
              [1, -2*r_pole*cos(2.0), r_pole^2]);
y = filter(1, a_poly, u) + 0.005 * randn(N, 1);

G = sidFreqBT(y, u, 'WindowSize', 100);
[n, sv] = sidModelOrder(G);

assert(n == 4, 'Expected n = 4 for 4th-order system, got n = %d', n);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: 4th-order SISO system detected (n = %d).\n', n);

%% Test 3: Known 1st-order system (edge case)
rng(300);
N = 1000;
u = randn(N, 1);
y = filter([0.8], [1, -0.9], u) + 0.005 * randn(N, 1);

G = sidFreqBT(y, u, 'WindowSize', 40);
[n, ~] = sidModelOrder(G);

assert(n == 1, 'Expected n = 1 for 1st-order system, got n = %d', n);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: 1st-order system detected (n = %d).\n', n);

%% Test 4: MIMO system (2x2, 3rd-order)
rng(400);
N = 3000;
ny = 2; nu = 2; nx = 3;
A_true = [0.7 0.2 0; -0.2 0.8 0.1; 0 -0.1 0.6];
B_true = [1 0; 0 0.5; 0.3 0.2];
C_obs = [1 0 0.5; 0 1 0];
x = zeros(nx, 1);
y = zeros(N, ny);
u = randn(N, nu);
for k = 1:N
    y(k, :) = (C_obs * x)' + 0.01 * randn(1, ny);
    x = A_true * x + B_true * u(k, :)';
end

G = sidFreqBT(y, u, 'WindowSize', 60);
[n, sv] = sidModelOrder(G);

assert(n >= 2 && n <= 4, ...
    'Expected n close to 3 for 3rd-order MIMO system, got n = %d', n);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: MIMO system order detected (n = %d, true = 3).\n', n);

%% Test 5: Threshold method
rng(500);
N = 2000;
u = randn(N, 1);
y = zeros(N, 1);
for k = 3:N
    y(k) = -a1 * y(k-1) - a2 * y(k-2) + b1 * u(k-1) + 0.01 * randn;
end

G = sidFreqBT(y, u, 'WindowSize', 60);
[n_thresh, ~] = sidModelOrder(G, 'Threshold', 0.01);

% Loose sanity range only. On a noisy BT estimate the singular-value tail is
% dense around a relative threshold of 0.01, so the exact count is sensitive
% to the BLAS/SVD engine (MATLAB ~4, Octave ~7); the bound is wide enough to
% absorb that while still catching a broken method. The threshold method's
% overcounting on noisy tails is part of the deferred noise-floor work (#139).
assert(n_thresh >= 1 && n_thresh <= 12, ...
    'Threshold method returned unreasonable n = %d', n_thresh);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: threshold method returned n = %d.\n', n_thresh);

%% Test 6: Output struct fields
rng(600);
N = 1000;
u = randn(N, 1);
y = filter([1], [1, -0.9], u) + 0.01 * randn(N, 1);
G = sidFreqBT(y, u);
[n, sv] = sidModelOrder(G);

assert(isstruct(sv), 'Second output must be a struct.');
assert(isfield(sv, 'SingularValues'), 'Missing field: SingularValues');
assert(isfield(sv, 'Horizon'), 'Missing field: Horizon');
assert(isvector(sv.SingularValues), 'SingularValues must be a vector.');
assert(isscalar(sv.Horizon), 'Horizon must be a scalar.');
assert(all(sv.SingularValues >= 0), 'Singular values must be non-negative.');
assert(issorted(sv.SingularValues(end:-1:1)), ...
    'Singular values must be in non-increasing order.');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: output struct fields correct.\n');

%% Test 7: Custom horizon
rng(700);
N = 1000;
u = randn(N, 1);
y = filter([0.5], [1, -0.8], u) + 0.01 * randn(N, 1);
G = sidFreqBT(y, u);
[n, sv] = sidModelOrder(G, 'Horizon', 20);

assert(sv.Horizon == 20, 'Expected Horizon = 20, got %d', sv.Horizon);
assert(length(sv.SingularValues) > 0, 'SingularValues should not be empty.');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: custom horizon accepted (r = %d).\n', sv.Horizon);

%% Test 8: Plot option runs without error
rng(800);
N = 500;
u = randn(N, 1);
y = filter([1], [1, -0.85], u) + 0.01 * randn(N, 1);
G = sidFreqBT(y, u);
try
    sidModelOrder(G, 'Plot', true);
    close all;
catch e
    if isempty(strfind(e.message, 'figure')) && ...
       isempty(strfind(e.message, 'display')) && ...
       isempty(strfind(e.message, 'DISPLAY')) && ...
       isempty(strfind(e.message, 'gnuplot'))
        rethrow(e);
    end
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: Plot option runs without error.\n');

%% Test 9: Input validation - bad struct
passed = false;
try
    sidModelOrder(42);
catch e
    if ~isempty(strfind(e.identifier, 'sid:'))
        passed = true;
    end
end
assert(passed, 'Should error on non-struct input.');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: input validation rejects non-struct.\n');

%% Test 10: Noisy data - still detects correct order with enough data
rng(1000);
N = 5000;
u = randn(N, 1);
y = zeros(N, 1);
for k = 3:N
    y(k) = -a1 * y(k-1) - a2 * y(k-2) + b1 * u(k-1) + 0.1 * randn;
end

G = sidFreqBT(y, u, 'WindowSize', 100);
[n, ~] = sidModelOrder(G);

assert(n == 2, ...
    'Expected n = 2 for noisy 2nd-order system with N=5000, got n = %d', n);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: noisy data still detects n = %d.\n', n);

%% Test 11: Gap method — verify SingularValues and gap location (SPEC §8.12.12)
% n = argmax_k (sigma_k / sigma_{k+1})
% For a 2nd-order system, the gap should be at k=2.
rng(1100);
N = 3000;
u = randn(N, 1);
% 2nd-order: y(k) = -a1*y(k-1) - a2*y(k-2) + b1*u(k-1)
a1 = -1.5; a2 = 0.7; b1 = 1;
y = zeros(N, 1);
for k = 3:N
    y(k) = -a1 * y(k-1) - a2 * y(k-2) + b1 * u(k-1) + 0.01 * randn;
end

G = sidFreqBT(y, u, 'WindowSize', 80);
[n_est, info] = sidModelOrder(G);

assert(n_est == 2, 'Expected n=2, got n=%d', n_est);

% SingularValues should be returned and have correct properties
assert(isfield(info, 'SingularValues'), 'Should have SingularValues');
sv = info.SingularValues;
assert(all(sv >= 0), 'Singular values should be non-negative');
assert(issorted(flipud(sv(:))), 'Singular values should be non-increasing');

% Verify the gap is at position n_est: sigma_2/sigma_3 should be the max ratio
ratios = sv(1:end-1) ./ max(sv(2:end), eps);
[~, gap_idx] = max(ratios);
assert(gap_idx == n_est, ...
    'Gap should be at index %d, found at %d', n_est, gap_idx);

% The gap ratio should be large (clear separation)
assert(ratios(n_est) > 5, ...
    'Gap ratio at n=%d should be large, got %.2f', n_est, ratios(n_est));
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: gap at k=%d (ratio=%.1f).\n', ...
    gap_idx, ratios(n_est));

%% Test 12: Threshold method — verify count formula (SPEC §8.12.12)
% n = count of sigma_k / sigma_1 > threshold
[n_thr, info_thr] = sidModelOrder(G, 'Threshold', 0.01);
sv = info_thr.SingularValues;

% Manual: count how many sigma_k / sigma_1 exceed threshold
manual_n = sum(sv / sv(1) > 0.01);
assert(n_thr == manual_n, ...
    'Threshold method: n=%d should match manual count %d', ...
    n_thr, manual_n);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 12 passed: threshold method (n=%d).\n', n_thr);

%% Test 13: White noise — no system, should return small n
rng(1300);
N = 2000;
u = randn(N, 1);
y = randn(N, 1);  % pure white noise, no relation to u

G_wn = sidFreqBT(y, u, 'WindowSize', 60);
[n_wn, info_wn] = sidModelOrder(G_wn);

% With no system, singular values should decay slowly (no clear gap)
% n can vary but should be well below the Hankel dimension
sv_wn = info_wn.SingularValues;
nSV = length(sv_wn);
assert(n_wn < nSV, ...
    'White noise: n=%d should be less than Hankel dim %d', n_wn, nSV);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 13 passed: white noise gives n=%d (Hankel dim=%d).\n', ...
    n_wn, nSV);

%% Test 14: Threshold tuning — different thresholds give different n
rng(1400);
N = 3000;
u = randn(N, 1);
a1 = -1.5; a2 = 0.7; b1 = 1;
y = zeros(N, 1);
for k = 3:N
    y(k) = -a1 * y(k-1) - a2 * y(k-2) + b1 * u(k-1) + 0.01 * randn;
end
G_th = sidFreqBT(y, u, 'WindowSize', 80);

% Strict threshold → smaller n
[n_strict, ~] = sidModelOrder(G_th, 'Threshold', 0.5);
% Loose threshold → larger n
[n_loose, ~] = sidModelOrder(G_th, 'Threshold', 0.001);

assert(n_strict <= n_loose, ...
    'Strict threshold n=%d should be <= loose n=%d', n_strict, n_loose);
assert(n_strict >= 1, 'Strict threshold should give at least n=1');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 14 passed: threshold tuning (strict=%d, loose=%d).\n', ...
    n_strict, n_loose);

%% Test 15: Strictly-proper plant with a finite zero -> n = 2 (issue #139)
% Exact analytic G(z) = (z + 0.5) / (z^2 - 1.2728 z + 0.81) on the (0, pi]
% grid (true n = 2, zero at -0.5). The pre-fix code kept the lag-0 IFFT
% sample and built the Hankel of z^{-1}G(z), estimating n = 3.
nf = 128;
w = (1:nf)' * pi / nf;
zc = exp(1j * w);
Gz = (zc + 0.5) ./ (zc.^2 - 1.2728 * zc + 0.81);
resFZ = struct('Frequency', w, 'Response', reshape(Gz, nf, 1, 1), ...
    'Method', 'sidFreqBT');
% Order via the threshold method: deterministic and cross-language-consistent
% (the gap method's noise-floor handling is a separate follow-up; see the
% note in sidModelOrder). The lag-0 fix is verified directly on the third
% singular value, which pre-fix was structural (~0.4/6.1) and post-fix
% collapses to numerical noise (~1e-4).
[nFZ, svFZ] = sidModelOrder(resFZ, 'Threshold', 1e-2);
assert(nFZ == 2, ...
    'finite-zero strictly proper plant: expected n=2, got %d', nFZ);
sFZ = svFZ.SingularValues;
assert(sFZ(3) / sFZ(1) < 1e-3, ...
    'spurious structural singular value: %.2e', sFZ(3) / sFZ(1));
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 15 passed: strictly-proper finite-zero plant n=2 (#139).\n');

fprintf('test_sidModelOrder: %d/%d passed\n', runner__nPassed, runner__nPassed);
