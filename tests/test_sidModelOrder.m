%% test_sidModelOrder - Unit tests for model order estimation
%
% Tests sidModelOrder for known systems of various orders, SISO and MIMO,
% gap and threshold detection, output struct fields, and edge cases.

fprintf('Running test_sidModelOrder...\n');

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
fprintf('  Test 2 passed: 4th-order SISO system detected (n = %d).\n', n);

%% Test 3: Known 1st-order system (edge case)
rng(300);
N = 1000;
u = randn(N, 1);
y = filter([0.8], [1, -0.9], u) + 0.005 * randn(N, 1);

G = sidFreqBT(y, u, 'WindowSize', 40);
[n, ~] = sidModelOrder(G);

assert(n == 1, 'Expected n = 1 for 1st-order system, got n = %d', n);
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

assert(n_thresh >= 1 && n_thresh <= 6, ...
    'Threshold method returned unreasonable n = %d', n_thresh);
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
fprintf('  Test 10 passed: noisy data still detects n = %d.\n', n);

fprintf('test_sidModelOrder: all tests passed.\n');
