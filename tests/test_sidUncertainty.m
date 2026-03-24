%% test_sidUncertainty - Unit tests for asymptotic standard deviations
%
% Tests sidUncertainty(G, PhiV, Coh, N, W) for correct formulas,
% edge cases, and dimensional consistency.

fprintf('Running test_sidUncertainty...\n');

%% Test 1: Noise spectrum std formula
% PhiVStd = sqrt(2*CW/N) * |PhiV|
% For W = [1; 0.5; 0]: CW = 1^2 + 2*(0.5^2 + 0^2) = 1 + 0.5 = 1.5
W = [1; 0.5; 0];
N = 100;
PhiV = [2.0; 3.0; 5.0];
CW = 1.5;
expected_PhiVStd = sqrt(2 * CW / N) * abs(PhiV);
[~, PhiVStd] = sidUncertainty([], PhiV, [], N, W);
assert(max(abs(PhiVStd - expected_PhiVStd)) < 1e-12, 'PhiVStd formula should match');

%% Test 2: Transfer function std formula (SISO)
G = [1+1j; 2-0.5j; 0.5+0.3j];
PhiV = [1; 1; 1];
Coh = [0.9; 0.5; 0.99];
N = 1000;
W = sidHannWin(10);
CW = W(1)^2 + 2 * sum(W(2:end).^2);

expected_GVar = (CW / N) .* abs(G).^2 .* (1 - Coh) ./ Coh;
expected_GStd = sqrt(expected_GVar);

[GStd, ~] = sidUncertainty(G, PhiV, Coh, N, W);
assert(max(abs(GStd - expected_GStd)) < 1e-12, 'GStd formula should match for SISO');

%% Test 3: Time series mode (G = [])
[GStd, PhiVStd] = sidUncertainty([], [1; 2], [], 100, sidHannWin(5));
assert(isempty(GStd), 'GStd should be empty when G is empty');
assert(length(PhiVStd) == 2, 'PhiVStd should still be computed');

%% Test 4: MIMO mode (Coh = []) returns NaN for GStd
G_mimo = randn(10, 2, 3) + 1j * randn(10, 2, 3);
PhiV_mimo = abs(randn(10, 2, 2));
[GStd, ~] = sidUncertainty(G_mimo, PhiV_mimo, [], 200, sidHannWin(10));
assert(isequal(size(GStd), size(G_mimo)), 'GStd should be same size as G');
assert(all(isnan(GStd(:))), 'MIMO GStd should be all NaN');

%% Test 5: Uncertainty decreases with more data
W = sidHannWin(20);
G = 1 + 0.5j;
PhiV = 1.0;
Coh = 0.8;
[GStd1, PhiVStd1] = sidUncertainty(G, PhiV, Coh, 100, W);
[GStd2, PhiVStd2] = sidUncertainty(G, PhiV, Coh, 10000, W);
assert(GStd2 < GStd1, 'More data should reduce G uncertainty');
assert(PhiVStd2 < PhiVStd1, 'More data should reduce PhiV uncertainty');

%% Test 6: Higher coherence reduces G uncertainty
[GStd_hi, ~] = sidUncertainty(G, PhiV, 0.99, 1000, W);
[GStd_lo, ~] = sidUncertainty(G, PhiV, 0.3, 1000, W);
assert(GStd_hi < GStd_lo, 'Higher coherence should reduce G uncertainty');

%% Test 7: Window norm CW computation
% For M=2, W = [1, 0.5, 0]: CW = 1 + 2*(0.25 + 0) = 1.5
W = sidHannWin(2);
[~, PhiVStd] = sidUncertainty([], [1], [], 100, W);
expected = sqrt(2 * 1.5 / 100) * 1;
assert(abs(PhiVStd - expected) < 1e-12, 'CW should be 1.5 for M=2 Hann window');

%% Test 8: Zero coherence handling (clamped to eps)
G = 1 + 0j;
PhiV = 1.0;
Coh = 0;
[GStd, ~] = sidUncertainty(G, PhiV, Coh, 1000, sidHannWin(10));
assert(isfinite(GStd), 'Zero coherence should not produce Inf (clamped to eps)');

fprintf('  test_sidUncertainty: ALL PASSED\n');
