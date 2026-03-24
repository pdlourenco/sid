%% test_sidHannWin - Unit tests for the Hann lag window function
%
% Tests sidHannWin(M) for correctness of boundary values, shape,
% symmetry, and known analytical values.

fprintf('Running test_sidHannWin...\n');

%% Test 1: Boundary values W(0)=1, W(M)=0
for M = [2, 5, 10, 30, 100]
    W = sidHannWin(M);
    assert(abs(W(1) - 1) < 1e-15, 'W(0) should be 1 for M=%d', M);
    assert(abs(W(end)) < 1e-15, 'W(M) should be 0 for M=%d', M);
end

%% Test 2: Output size is (M+1 x 1) column vector
for M = [2, 7, 50]
    W = sidHannWin(M);
    assert(isequal(size(W), [M+1, 1]), 'Size should be [%d,1] for M=%d', M+1, M);
end

%% Test 3: All values in [0, 1]
W = sidHannWin(50);
assert(all(W >= 0) && all(W <= 1), 'All window values should be in [0,1]');

%% Test 4: Known values for M=2
W = sidHannWin(2);
expected = [1; 0.5*(1+cos(pi/2)); 0.5*(1+cos(pi))];  % [1, 0.5, 0]
assert(max(abs(W - expected)) < 1e-15, 'M=2: values should be [1, 0.5, 0]');

%% Test 5: Known values for M=4
W = sidHannWin(4);
expected = zeros(5, 1);
for tau = 0:4
    expected(tau+1) = 0.5 * (1 + cos(pi * tau / 4));
end
assert(max(abs(W - expected)) < 1e-15, 'M=4: values should match formula');

%% Test 6: Monotonically decreasing
W = sidHannWin(20);
assert(all(diff(W) <= 0), 'Hann window should be monotonically decreasing for lags 0..M');

%% Test 7: Mid-point value for even M
M = 10;
W = sidHannWin(M);
% At tau = M/2: W = 0.5*(1 + cos(pi/2)) = 0.5
assert(abs(W(M/2 + 1) - 0.5) < 1e-15, 'Mid-point should be 0.5 for even M');

fprintf('  test_sidHannWin: ALL PASSED\n');
