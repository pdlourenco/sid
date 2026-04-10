%% test_util_msd_ltv - Unit tests for LTV SMD plant construction.

fprintf('Running test_util_msd_ltv...\n');
runner__nPassed = 0;

if isempty(which('util_msd_ltv'))
    test__here = fileparts(mfilename('fullpath'));
    if isempty(test__here)
        test__here = pwd;
    end
    addpath(fullfile(fileparts(test__here), 'examples'));
end

%% Test 1: LTV collapses to LTI when every input is constant
m = [1; 1]; k = [100; 80]; c = [0.5; 0.5];
F = [1; 0]; Ts = 0.01; N = 50;
[Ad, Bd] = util_msd_ltv(m, k, c, F, Ts, N);
[Ad_ref, Bd_ref] = util_msd(m, k, c, F, Ts);
assert(isequal(size(Ad), [4 4 N]), 'Ad shape');
assert(isequal(size(Bd), [4 1 N]), 'Bd shape');
for kk = 1:N
    assert(max(max(abs(Ad(:,:,kk) - Ad_ref))) < 1e-14, ...
        'LTI replication slice %d', kk);
    assert(max(max(abs(Bd(:,:,kk) - Bd_ref))) < 1e-14, ...
        'LTI replication slice %d (Bd)', kk);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: LTV collapses to LTI.\n');

%% Test 2: time-invariant inputs without N raises
threw = false;
try
    util_msd_ltv([1], [100], [0.5], 1, 0.01);
catch
    threw = true;
end
assert(threw, 'Expected error when N omitted for all-constant inputs');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: missing N rejected.\n');

%% Test 3: ramping stiffness tracked per step
N3 = 100;
m3 = [1; 1]; c3 = [0.5; 0.5];
k_tv = zeros(2, N3);
k_tv(1, :) = linspace(100, 400, N3);
k_tv(2, :) = 80;
F3 = [1; 0];
[Ad3, ~] = util_msd_ltv(m3, k_tv, c3, F3, 0.01);
[Ad3_0, ~] = util_msd(m3, k_tv(:, 1),    c3, F3, 0.01);
[Ad3_N, ~] = util_msd(m3, k_tv(:, end),  c3, F3, 0.01);
assert(max(max(abs(Ad3(:, :, 1)  - Ad3_0))) < 1e-14, 'first slice mismatch');
assert(max(max(abs(Ad3(:, :, N3) - Ad3_N))) < 1e-14, 'last slice mismatch');
assert(max(max(abs(Ad3(:, :, round(N3/2)) - Ad3_0))) > 1e-6, ...
    'middle slice should differ from first');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: ramping stiffness tracked.\n');

%% Test 4: step change in k produces piecewise-constant sequence
N4 = 40;
k_step = zeros(1, N4);
k_step(1, 1:N4/2)   = 100;
k_step(1, N4/2+1:N4) = 50;
[Ad4, ~] = util_msd_ltv(1, k_step, 0.5, 1, 0.01);
assert(max(max(abs(Ad4(:,:,1) - Ad4(:,:,N4/2)))) < 1e-14, ...
    'pre-step segment not constant');
assert(max(max(abs(Ad4(:,:,N4/2+1) - Ad4(:,:,N4)))) < 1e-14, ...
    'post-step segment not constant');
assert(max(max(abs(Ad4(:,:,N4/2) - Ad4(:,:,N4/2+1)))) > 1e-6, ...
    'step discontinuity missing');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: step-change discontinuity preserved.\n');

%% Test 5: time-varying F accepted
N5 = 20;
F_tv = zeros(2, 1, N5);
F_tv(:, 1, :) = repmat(linspace(0, 1, N5), 2, 1);
[Ad5, Bd5] = util_msd_ltv([1; 1], [100; 80], [0.5; 0.5], F_tv, 0.01);
assert(isequal(size(Ad5), [4 4 N5]), 'Ad5 shape');
assert(isequal(size(Bd5), [4 1 N5]), 'Bd5 shape');
% Dynamics are LTI so all Ad slices must match.
for kk = 2:N5
    assert(max(max(abs(Ad5(:,:,kk) - Ad5(:,:,1)))) < 1e-14, ...
        'Ad slice drift with time-varying F');
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: time-varying F accepted.\n');

fprintf('test_util_msd_ltv: %d/%d passed\n', runner__nPassed, runner__nPassed);
