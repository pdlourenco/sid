%% test_util_msd - Unit tests for the LTI SMD plant helper.
%
% Verifies util_msd (living in matlab/examples/) produces the expected
% discrete state-space matrices across a range of chain sizes n=1..5.
%
% The n=3 parity test against the legacy private helper sidTestMSD
% that used to live here was removed when sidTestMSD was deleted; its
% regression role is now carried by the cross-language check in
% testdata/reference_test_msd.json, which this helper must still
% match bit-for-bit under the MATLAB port.

fprintf('Running test_util_msd...\n');
runner__nPassed = 0;

% Make util_msd* discoverable (they live in matlab/examples/, not on
% the test runner's search path by default when invoked standalone).
% runAllTests.m already addpaths matlab/examples/, so this branch is
% only exercised when the file is run directly.
if isempty(which('util_msd'))
    test__here = fileparts(mfilename('fullpath'));
    if isempty(test__here)
        test__here = pwd;
    end
    addpath(fullfile(fileparts(test__here), 'examples'));
end

%% Test 1: n=1 SDOF matches the analytic 2x2 ZOH result
m1 = 1.0; k1 = 100.0; c1 = 0.5;
[Ad1, Bd1] = util_msd(m1, k1, c1, 1.0, 0.01);
Ac1 = [0, 1; -k1/m1, -c1/m1];
Bc1 = [0; 1/m1];
Ad1_exp = expm(Ac1 * 0.01);
Bd1_exp = Ac1 \ (Ad1_exp - eye(2)) * Bc1;
assert(max(max(abs(Ad1 - Ad1_exp))) < 1e-14, 'n=1 Ad analytic mismatch');
assert(max(max(abs(Bd1 - Bd1_exp))) < 1e-14, 'n=1 Bd analytic mismatch');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: n=1 SDOF analytic match.\n');

%% Test 2: n=2 undamped eigenvalues match modal analysis
m2 = [1; 1]; k2 = [100; 80]; c2 = [0; 0]; F2 = eye(2);
Ts2 = 0.001;
[Ad2, ~] = util_msd(m2, k2, c2, F2, Ts2);
K_ref = [180 -80; -80 80];
omega_sq = sort(eig(K_ref));
omega = sqrt(omega_sq);
disc_args = sort(unique(abs(angle(eig(Ad2)))));
expected = omega * Ts2;
assert(max(abs(disc_args - expected)) < 1e-5, ...
    'n=2 eigenvalue mismatch: %.2e', max(abs(disc_args - expected)));
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: n=2 modal eigenvalues.\n');

%% Test 3: n=5 chain shapes and stability
n5 = 5;
m5 = ones(n5, 1);
k5 = linspace(100, 60, n5)';
c5 = 0.5 * ones(n5, 1);
F5 = [1; zeros(n5 - 1, 1)];
[Ad5, Bd5] = util_msd(m5, k5, c5, F5, 0.01);
assert(isequal(size(Ad5), [2*n5 2*n5]), 'n=5 Ad shape');
assert(isequal(size(Bd5), [2*n5 1]),    'n=5 Bd shape');
assert(all(abs(eig(Ad5)) < 1), 'n=5 chain not strictly stable');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: n=5 shapes and stability.\n');

%% Test 4: rejects mismatched shapes
threw = false;
try
    util_msd([1; 1], 100, [0.5; 0.5], eye(2), 0.01);  % k scalar instead of (2,)
catch
    threw = true;
end
assert(threw, 'util_msd should reject mismatched k length');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: rejects mismatched shapes.\n');

fprintf('test_util_msd: %d/%d passed\n', runner__nPassed, runner__nPassed);
