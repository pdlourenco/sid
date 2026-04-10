%% test_util_msd_nl - Unit tests for the Duffing-capable RK4 simulator.

fprintf('Running test_util_msd_nl...\n');
runner__nPassed = 0;

if isempty(which('util_msd_nl'))
    test__here = fileparts(mfilename('fullpath'));
    if isempty(test__here)
        test__here = pwd;
    end
    addpath(fullfile(fileparts(test__here), 'examples'));
end

%% Test 1: linear case (k_cubic = 0) matches the exact ZOH propagation
m = 1; k = 100; kc = 0; c = 0.5; F = 1; Ts = 0.005;
rng(7);
N = 200;
u = randn(N, 1);
x_rk  = util_msd_nl(m, k, kc, c, F, Ts, u, [], 4);
[Ad, Bd] = util_msd(m, k, c, F, Ts);
x_zoh = zeros(N + 1, 2);
for kk = 1:N
    x_zoh(kk + 1, :) = (Ad * x_zoh(kk, :)' + Bd * u(kk))';
end
err_rms = sqrt(mean(mean((x_rk - x_zoh).^2)));
assert(err_rms < 1e-8, 'RK4 vs ZOH mismatch: %.2e', err_rms);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: linear RK4 matches ZOH (err=%.2e).\n', err_rms);

%% Test 2: Duffing deviates measurably from linear at large amplitude
k_cub = 1000;
rng(42);
N2 = 1000;
u2 = 3.0 * randn(N2, 1);
x_lin = util_msd_nl(m, k, 0,     c, F, 0.01, u2, [], 4);
x_duf = util_msd_nl(m, k, k_cub, c, F, 0.01, u2, [], 4);
assert(all(isfinite(x_lin(:))), 'linear blew up');
assert(all(isfinite(x_duf(:))), 'duffing blew up');
rms_diff = sqrt(mean((x_lin(:, 1) - x_duf(:, 1)).^2));
rms_lin  = sqrt(mean(x_lin(:, 1).^2));
assert(rms_diff / rms_lin > 0.05, ...
    'Duffing response too close to linear: %.3f%%', 100 * rms_diff / rms_lin);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: Duffing deviation %.1f%%.\n', ...
    100 * rms_diff / rms_lin);

%% Test 3: zero input, zero initial state -> stays at origin
m3 = [1; 1]; k3 = [100; 80]; kc3 = [500; 0]; c3 = [0.5; 0.5];
F3 = [1; 0];
N3 = 200;
x3 = util_msd_nl(m3, k3, kc3, c3, F3, 0.01, zeros(N3, 1), [], 1);
assert(max(max(abs(x3))) < 1e-14, 'zero-input trajectory not zero');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: zero-input trajectory is zero.\n');

%% Test 4: canonical notebook parameters remain stable
rng(0);
N4 = 2000;
u4 = 2.0 * randn(N4, 1);
x4 = util_msd_nl(1, 100, 1000, 0.5, 1, 0.01, u4, [], 4);
assert(all(isfinite(x4(:))), 'canonical run blew up');
assert(max(abs(x4(:, 1))) < 10, 'canonical run drifted too far');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: canonical run stable, peak |x|=%.3f.\n', ...
    max(abs(x4(:, 1))));

%% Test 5: mismatched input width rejected
threw = false;
try
    util_msd_nl(1, 100, 0, 0.5, 1, 0.01, zeros(10, 2));
catch
    threw = true;
end
assert(threw, 'util_msd_nl should reject mismatched u width');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: mismatched input width rejected.\n');

fprintf('test_util_msd_nl: %d/%d passed\n', runner__nPassed, runner__nPassed);
