%% test_sidCov - Unit tests for biased cross-covariance estimation
%
% Tests sidCov(x, z, maxLag) for correctness against hand-computed values,
% white noise properties, and output dimensions.

fprintf('Running test_sidCov...\n');

%% Test 1: Known hand-computed values (auto-covariance)
% x = [1; 2; 3; 4], N = 4
% R(0) = (1*1 + 2*2 + 3*3 + 4*4) / 4 = 30/4 = 7.5
% R(1) = (2*1 + 3*2 + 4*3) / 4 = 20/4 = 5.0
% R(2) = (3*1 + 4*2) / 4 = 11/4 = 2.75
x = [1; 2; 3; 4];
R = sidCov(x, x, 2);
assert(abs(R(1) - 7.5) < 1e-12, 'R(0) should be 7.5');
assert(abs(R(2) - 5.0) < 1e-12, 'R(1) should be 5.0');
assert(abs(R(3) - 2.75) < 1e-12, 'R(2) should be 2.75');

%% Test 2: Known hand-computed values (cross-covariance)
% x = [1; 2; 3; 4], z = [4; 3; 2; 1], N = 4
% R_xz(0) = (1*4 + 2*3 + 3*2 + 4*1) / 4 = 20/4 = 5.0
% R_xz(1) = (2*4 + 3*3 + 4*2) / 4 = 25/4 = 6.25
% R_xz(2) = (3*4 + 4*3) / 4 = 24/4 = 6.0
z = [4; 3; 2; 1];
R = sidCov(x, z, 2);
assert(abs(R(1) - 5.0) < 1e-12, 'R_xz(0) should be 5.0');
assert(abs(R(2) - 6.25) < 1e-12, 'R_xz(1) should be 6.25');
assert(abs(R(3) - 6.0) < 1e-12, 'R_xz(2) should be 6.0');

%% Test 3: Output size for scalar signals is (maxLag+1 x 1)
R = sidCov(x, x, 3);
assert(isequal(size(R), [4, 1]), 'Scalar output should be (maxLag+1 x 1)');

%% Test 4: Output size for multi-channel signals
% x: (10 x 2), z: (10 x 3), maxLag = 4 => R: (5 x 2 x 3)
rng(42);
xm = randn(10, 2);
zm = randn(10, 3);
R = sidCov(xm, zm, 4);
assert(isequal(size(R), [5, 2, 3]), 'Multi-channel output should be (5 x 2 x 3)');

%% Test 5: Biased estimator divides by N (not N-|tau|)
% For x = [1; 0; 0; 0], auto-cov at lag 0:
% Biased: R(0) = 1/4 = 0.25
% Unbiased would be: 1/4 = 0.25 (same at lag 0)
% At lag 1: biased = 0/4 = 0, unbiased = 0/3 = 0
% At lag 3: x(4)*x(1) / N = 0*1/4 = 0 (biased) vs 0/1 (unbiased)
x_imp = [1; 0; 0; 0];
R = sidCov(x_imp, x_imp, 3);
assert(abs(R(1) - 0.25) < 1e-12, 'Biased R(0) for impulse should be 1/N');
assert(abs(R(2)) < 1e-12, 'R(1) should be 0 for impulse');

%% Test 6: White noise auto-covariance (statistical)
rng(123);
N = 100000;
x_wn = randn(N, 1);
R = sidCov(x_wn, x_wn, 10);
% R(0) should be approximately 1 (variance of standard normal)
assert(abs(R(1) - 1.0) < 0.02, 'R(0) of white noise should be ~1');
% R(tau>0) should be approximately 0
for tau = 1:10
    assert(abs(R(tau+1)) < 0.02, 'R(%d) of white noise should be ~0', tau);
end

%% Test 7: Lag 0 auto-covariance equals biased variance
rng(99);
x = randn(200, 1) * 3 + 5;  % mean=5, std=3
R = sidCov(x, x, 0);
biasedVar = sum(x.^2) / length(x);  % NOT mean-subtracted
assert(abs(R(1) - biasedVar) < 1e-10, 'R(0) should equal (1/N)*sum(x^2)');

%% Test 8: maxLag = 0 returns single value
R = sidCov([1; 2; 3], [1; 2; 3], 0);
assert(length(R) == 1, 'maxLag=0 should return single value');

fprintf('  test_sidCov: ALL PASSED\n');
