%% test_sidFreqDomainSim - Unit tests for the frequency-domain simulation helper
%
% Tests sidFreqDomainSim (SPEC §15), the IFFT-based frequency-domain model
% simulation used by sidResidual / sidCompare. Port-symmetric with
% test_freq_domain_sim.py (#124). Also cross-validated via
% reference_freq_domain_sim.json; these are independent-oracle unit tests.
%
% sidFreqDomainSim is a private helper — runAllTests.m puts matlab/sid/private
% on the path (via the private_test_shim), so it is callable here.

fprintf('Running test_sidFreqDomainSim...\n');
runner__nPassed = 0;

%% Test 1: Constant gain scales a zero-mean input (SPEC §15, #124)
rng(124);
N = 256;
u = randn(N, 1); u = u - mean(u);
freqs = (1:N/2)' * pi / (N/2);       % dense grid over (0, pi]
G = 2.0 * ones(numel(freqs), 1);     % constant gain
Y = sidFreqDomainSim(G, freqs, u, N);
relErr = max(abs(Y(:) - 2 * u(:))) / max(abs(2 * u(:)));
assert(relErr < 1e-10, ...
    'Constant gain should scale input by 2 (relErr=%.2e)', relErr);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: constant gain scales input (relErr=%.2e).\n', relErr);

%% Test 2: Out-of-grid frequencies contribute nothing (SPEC §15, #124)
rng(125);
N = 128;
u = randn(N, 1);
freqs = linspace(0.05, pi / 3, 30)';  % low third of (0, pi] only
G = ones(numel(freqs), 1);
Y = sidFreqDomainSim(G, freqs, u, N);
assert(sum(Y(:).^2) < sum(u(:).^2), ...
    'Out-of-grid (high-freq) content should be removed');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: out-of-grid frequencies zeroed.\n');

fprintf('test_sidFreqDomainSim: %d/%d passed\n', runner__nPassed, runner__nPassed);
