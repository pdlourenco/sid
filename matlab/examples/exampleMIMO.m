%% exampleMIMO - Multi-Input Multi-Output frequency response estimation
%
% This example demonstrates MIMO system identification with sidFreqBT.
% Key differences from SISO: Response is a 3D array, NoiseSpectrum is a
% spectral matrix, and Coherence is not available.

runner__nCompleted = 0;

%% 2-output, 1-input system
% Two independent channels driven by the same input:
%   G_1(z) = 1 / (1 - 0.5 z^{-1})
%   G_2(z) = 0.3 / (1 - 0.7 z^{-1})

rng(10);
N = 3000;
u = randn(N, 1);
y1 = filter(1, [1 -0.5], u) + 0.1 * randn(N, 1);
y2 = filter(0.3, [1 -0.7], u) + 0.1 * randn(N, 1);
y = [y1, y2];   % (N x 2) output matrix

result = sidFreqBT(y, u, 'WindowSize', 30);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: 2-output, 1-input system.\n', runner__nCompleted);

%% Inspect MIMO result dimensions
% Response is (nf x ny x nu) = (128 x 2 x 1) for this system.
% Coherence is empty for MIMO.

fprintf('Response size:      [%s]\n', num2str(size(result.Response)));
fprintf('NoiseSpectrum size: [%s]\n', num2str(size(result.NoiseSpectrum)));
fprintf('Coherence is empty: %d\n', isempty(result.Coherence));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Inspect MIMO result dimensions.\n', runner__nCompleted);

%% Plot individual channels
% sidBodePlot only shows the first channel, so we plot both manually.

w = result.Frequency;
G1_true = 1 ./ (1 - 0.5 * exp(-1j * w));
G2_true = 0.3 ./ (1 - 0.7 * exp(-1j * w));

figure;
subplot(2,1,1);
semilogx(w, 20*log10(abs(result.Response(:, 1))), 'b', 'DisplayName', 'Estimated');
hold on;
semilogx(w, 20*log10(abs(G1_true)), 'k--', 'DisplayName', 'True');
ylabel('Magnitude (dB)');
title('Channel 1: G_1(z) = 1/(1 - 0.5z^{-1})');
legend('show');
grid on;
hold off;

subplot(2,1,2);
semilogx(w, 20*log10(abs(result.Response(:, 2))), 'r', 'DisplayName', 'Estimated');
hold on;
semilogx(w, 20*log10(abs(G2_true)), 'k--', 'DisplayName', 'True');
ylabel('Magnitude (dB)');
xlabel('Frequency (rad/sample)');
title('Channel 2: G_2(z) = 0.3/(1 - 0.7z^{-1})');
legend('show');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Plot individual channels.\n', runner__nCompleted);

%% Noise spectral matrix
% For a 2-output system, NoiseSpectrum is (nf x 2 x 2): a Hermitian
% positive semi-definite matrix at each frequency. The diagonal elements
% are the individual output noise spectra.

nf = length(w);
diag11 = real(result.NoiseSpectrum(:, 1, 1));
diag22 = real(result.NoiseSpectrum(:, 2, 2));

figure;
semilogx(w, 10*log10(diag11), 'b', 'DisplayName', '\Phi_{v,11}');
hold on;
semilogx(w, 10*log10(diag22), 'r', 'DisplayName', '\Phi_{v,22}');
xlabel('Frequency (rad/sample)');
ylabel('Noise Spectrum (dB)');
title('Diagonal Elements of Noise Spectral Matrix');
legend('show');
grid on;
hold off;

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: Noise spectral matrix.\n', runner__nCompleted);

%% 2-output, 2-input system
% Full 2x2 transfer matrix: each output depends on both inputs.
%   y1 = G_11 * u1 + G_12 * u2 + v1
%   y2 = G_21 * u1 + G_22 * u2 + v2

rng(20);
N = 3000;
u2in = randn(N, 2);
y2out_1 = filter(1, [1 -0.5], u2in(:,1)) + filter(0.3, [1 -0.3], u2in(:,2)) ...
    + 0.1 * randn(N, 1);
y2out_2 = filter(0.5, [1 -0.7], u2in(:,1)) + filter(1, [1 -0.4], u2in(:,2)) ...
    + 0.1 * randn(N, 1);
y2out = [y2out_1, y2out_2];

result_22 = sidFreqBT(y2out, u2in, 'WindowSize', 30);

% Response is now (nf x 2 x 2)
fprintf('\n2x2 MIMO Response size: [%s]\n', num2str(size(result_22.Response)));

% Plot the full transfer matrix
figure;
titles = {'G_{11}', 'G_{12}'; 'G_{21}', 'G_{22}'};
for iy = 1:2
    for iu = 1:2
        subplot(2, 2, (iy-1)*2 + iu);
        semilogx(w, 20*log10(abs(result_22.Response(:, iy, iu))), 'b');
        ylabel('Magnitude (dB)');
        if iy == 2; xlabel('Frequency (rad/sample)'); end
        title(titles{iy, iu});
        grid on;
    end
end

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: 2-output, 2-input system.\n', runner__nCompleted);

%% MIMO uncertainty
% In v1.0, ResponseStd is NaN for MIMO (no asymptotic formula implemented).

fprintf('\nMIMO ResponseStd contains NaN: %d\n', all(isnan(result_22.ResponseStd(:))));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: MIMO uncertainty.\n', runner__nCompleted);

fprintf('exampleMIMO: %d/%d sections completed\n', runner__nCompleted, runner__nCompleted);
