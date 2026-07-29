%% test_sidLTVcosmicSolve - Unit test for the COSMIC solver ill-conditioning warning
%
% Exercises the SPEC §8.3.4 diagnostic (#120): the forward pass warns
% (sid:singularLbd) when an information block Lbd(k-1) is near-singular, while
% still returning a finite result. This is a defensive path not reachable from
% normal sidLTVdisc inputs, so it is driven directly via the private solver
% (on the path here through runAllTests' private_test_shim). Port-symmetric with
% test_ltv_cosmic_solve.py.

fprintf('Running test_sidLTVcosmicSolve...\n');
runner__nPassed = 0;

%% Test 1: Near-singular forward-pass block warns (SPEC §8.3.4, #120)
p = 1; q = 1; N = 2; d = p + q;
S = zeros(d, d, N);
S(:, :, 1) = diag([1.0, 1e-17]);   % rcond ~1e-17 < eps -> warns
S(:, :, 2) = 5.0 * eye(d);
T = ones(d, p, N);
lambda = ones(N - 1, 1);

warning('error', 'sid:singularLbd');
threw = false;
try
    sidLTVcosmicSolve(S, T, lambda, N, p, q);
catch e
    threw = ~isempty(strfind(e.identifier, 'singularLbd'));
end
warning('on', 'sid:singularLbd');
assert(threw, 'Expected sid:singularLbd warning on a near-singular block');

% With the warning left enabled, the solver still returns a finite result.
[C, ~] = sidLTVcosmicSolve(S, T, lambda, N, p, q);
assert(all(isfinite(C(:))), 'Solver must still return a finite result');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: sid:singularLbd warns on near-singular block.\n');

fprintf('test_sidLTVcosmicSolve: %d/%d passed\n', runner__nPassed, runner__nPassed);
