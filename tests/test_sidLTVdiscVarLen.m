%% test_sidLTVdiscVarLen - Unit tests for variable-length trajectory support
%
% Tests the Phase 8a extension of sidLTVdisc that accepts cell arrays
% of trajectories with different horizons.

fprintf('Running test_sidLTVdiscVarLen...\n');

%% Test 1: Cell array input runs without error
rng(2000);
p = 2; q = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
sigma = 0.02;

% 3 trajectories with different lengths
Ns = [20, 30, 40];
X_cell = cell(3, 1);
U_cell = cell(3, 1);
for l = 1:3
    Nl = Ns(l);
    U_cell{l} = randn(Nl, q);
    X_cell{l} = zeros(Nl+1, p);
    X_cell{l}(1, :) = randn(1, p);
    for k = 1:Nl
        X_cell{l}(k+1, :) = (A_true * X_cell{l}(k, :)' + B_true * U_cell{l}(k, :)')' + sigma * randn(1, p);
    end
end

result = sidLTVdisc(X_cell, U_cell, 'Lambda', 1e4);
assert(isstruct(result), 'Result should be a struct');
assert(strcmp(result.Method, 'sidLTVdisc'), 'Method should be sidLTVdisc');
assert(result.DataLength == 40, 'DataLength should be max horizon = 40');
assert(result.NumTrajectories == 3, 'NumTrajectories should be 3');
assert(isequal(size(result.A), [p, p, 40]), 'A dimensions should be (p x p x N)');
assert(isequal(size(result.B), [p, q, 40]), 'B dimensions should be (p x q x N)');
fprintf('  Test 1 passed: cell array input works.\n');

%% Test 2: Both paths recover the system well with uniform-length data
% Note: The 3D path normalizes by sqrt(N), the cell path by sqrt(L).
% With the same lambda, the effective regularization differs, so results
% are not bit-identical. We verify both recover the system independently.
rng(2100);
p = 2; q = 1; N = 30; L = 15;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
sigma = 0.02;

X_3d = zeros(N+1, p, L);
U_3d = randn(N, q, L);
X_cell = cell(L, 1);
U_cell = cell(L, 1);

for l = 1:L
    X_3d(1, :, l) = randn(1, p);
    for k = 1:N
        X_3d(k+1, :, l) = (A_true * X_3d(k, :, l)' + B_true * U_3d(k, :, l)')' + sigma * randn(1, p);
    end
    X_cell{l} = X_3d(:, :, l);
    U_cell{l} = U_3d(:, :, l);
end

res_3d   = sidLTVdisc(X_3d, U_3d, 'Lambda', 1e4);
res_cell = sidLTVdisc(X_cell, U_cell, 'Lambda', 1e4);

A_mean_3d = mean(res_3d.A, 3);
A_mean_cell = mean(res_cell.A, 3);
errA_3d   = norm(A_mean_3d - A_true, 'fro') / norm(A_true, 'fro');
errA_cell = norm(A_mean_cell - A_true, 'fro') / norm(A_true, 'fro');
assert(errA_3d < 0.15, '3D path should recover A well: %.4f', errA_3d);
assert(errA_cell < 0.15, 'Cell path should recover A well: %.4f', errA_cell);
fprintf('  Test 2 passed: both paths recover LTI system (3d=%.4f, cell=%.4f).\n', errA_3d, errA_cell);

%% Test 3: Mixed lengths - correct dimensions
rng(2200);
p = 1; q = 1;
Ns = [10, 20, 30, 15, 25];
L = length(Ns);
X_cell = cell(L, 1);
U_cell = cell(L, 1);
for l = 1:L
    X_cell{l} = randn(Ns(l)+1, p);
    U_cell{l} = randn(Ns(l), q);
end

result = sidLTVdisc(X_cell, U_cell, 'Lambda', 100);
assert(result.DataLength == 30, 'DataLength should be max(Ns) = 30');
assert(result.NumTrajectories == 5, 'NumTrajectories should be 5');
assert(isequal(size(result.A), [p, p, 30]), 'A size should be (1 x 1 x 30)');
assert(isequal(size(result.B), [p, q, 30]), 'B size should be (1 x 1 x 30)');
fprintf('  Test 3 passed: mixed-length dimensions correct.\n');

%% Test 4: Short + long trajectories - no crash, reasonable results
rng(2300);
p = 2; q = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
sigma = 0.02;

% 1 very short trajectory (N=3) + 4 long ones (N=40)
X_cell = cell(5, 1);
U_cell = cell(5, 1);

% Short trajectory
U_cell{1} = randn(3, q);
X_cell{1} = zeros(4, p);
X_cell{1}(1, :) = randn(1, p);
for k = 1:3
    X_cell{1}(k+1, :) = (A_true * X_cell{1}(k, :)' + B_true * U_cell{1}(k, :)')' + sigma * randn(1, p);
end

% Long trajectories
for l = 2:5
    U_cell{l} = randn(40, q);
    X_cell{l} = zeros(41, p);
    X_cell{l}(1, :) = randn(1, p);
    for k = 1:40
        X_cell{l}(k+1, :) = (A_true * X_cell{l}(k, :)' + B_true * U_cell{l}(k, :)')' + sigma * randn(1, p);
    end
end

result = sidLTVdisc(X_cell, U_cell, 'Lambda', 1e4);
assert(result.DataLength == 40, 'DataLength should be 40');
A_mean = mean(result.A, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
assert(errA < 0.2, 'Short+long should still recover A: %.4f', errA);
fprintf('  Test 4 passed: short+long trajectories work (errA=%.4f).\n', errA);

%% Test 5: Single trajectory in cell array
rng(2400);
p = 2; q = 1; N = 20;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
U_cell = {randn(N, q)};
X_cell = {zeros(N+1, p)};
X_cell{1}(1, :) = randn(1, p);
for k = 1:N
    X_cell{1}(k+1, :) = (A_true * X_cell{1}(k, :)' + B_true * U_cell{1}(k, :)')';
end

result = sidLTVdisc(X_cell, U_cell, 'Lambda', 1e6);
assert(result.NumTrajectories == 1, 'NumTrajectories should be 1');
assert(result.DataLength == N, 'DataLength should be N');
fprintf('  Test 5 passed: single trajectory in cell array works.\n');

%% Test 6: Input validation for cell arrays
% Mismatched p across trajectories
try
    sidLTVdisc({randn(10, 2), randn(10, 3)}, {randn(9, 1), randn(9, 1)}, 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for mismatched p');
catch e
    assert(strcmp(e.identifier, 'sid:dimMismatch'), ...
        'Expected sid:dimMismatch, got %s', e.identifier);
end

% Mismatched q across trajectories
try
    sidLTVdisc({randn(10, 2), randn(10, 2)}, {randn(9, 1), randn(9, 2)}, 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for mismatched q');
catch e
    assert(strcmp(e.identifier, 'sid:dimMismatch'), ...
        'Expected sid:dimMismatch, got %s', e.identifier);
end

% X is cell but U is not
try
    sidLTVdisc({randn(10, 2)}, randn(9, 1), 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for mismatched types');
catch e
    assert(strcmp(e.identifier, 'sid:badInput'), ...
        'Expected sid:badInput, got %s', e.identifier);
end

% Empty cell arrays
try
    sidLTVdisc({}, {}, 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for empty cell');
catch e
    assert(strcmp(e.identifier, 'sid:badInput'), ...
        'Expected sid:badInput, got %s', e.identifier);
end

% U/X row mismatch within a trajectory
try
    sidLTVdisc({randn(10, 2)}, {randn(8, 1)}, 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for row mismatch');
catch e
    assert(strcmp(e.identifier, 'sid:dimMismatch'), ...
        'Expected sid:dimMismatch, got %s', e.identifier);
end

% Different number of trajectories
try
    sidLTVdisc({randn(10, 2), randn(10, 2)}, {randn(9, 1)}, 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for different L');
catch e
    assert(strcmp(e.identifier, 'sid:dimMismatch'), ...
        'Expected sid:dimMismatch, got %s', e.identifier);
end

fprintf('  Test 6 passed: input validation for cell arrays correct.\n');

%% Test 7: L-curve auto lambda with cell array input
rng(2500);
p = 2; q = 1;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
sigma = 0.02;
Ns = [20, 30, 25];
X_cell = cell(3, 1);
U_cell = cell(3, 1);
for l = 1:3
    Nl = Ns(l);
    U_cell{l} = randn(Nl, q);
    X_cell{l} = zeros(Nl+1, p);
    X_cell{l}(1, :) = randn(1, p);
    for k = 1:Nl
        X_cell{l}(k+1, :) = (A_true * X_cell{l}(k, :)' + B_true * U_cell{l}(k, :)')' + sigma * randn(1, p);
    end
end

result = sidLTVdisc(X_cell, U_cell, 'Lambda', 'auto');
assert(all(result.Lambda > 0), 'Auto lambda should be positive');
fprintf('  Test 7 passed: L-curve auto lambda with cell arrays works.\n');

fprintf('test_sidLTVdiscVarLen: ALL TESTS PASSED\n');
