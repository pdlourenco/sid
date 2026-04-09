%% test_sidPlotting - Unit tests for sidBodePlot and sidSpectrumPlot
%
% Beyond smoke testing, these tests validate that the correct data is
% plotted on the axes (line data, axis scales, labels, confidence bands).
% Uses xvfb in CI for headless execution.

fprintf('Running test_sidPlotting...\n');
runner__nPassed = 0;

%% Setup: create test result structs
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1 0.5], [1 -0.8], u) + 0.1 * randn(N, 1);
result_siso = sidFreqBT(y, u);
result_ts = sidFreqBT(randn(300, 1), []);

%% Test 1: sidBodePlot returns correct handles
h = sidBodePlot(result_siso);
assert(isfield(h, 'fig'), 'Should have fig handle');
assert(isfield(h, 'axMag'), 'Should have axMag handle');
assert(isfield(h, 'axPhase'), 'Should have axPhase handle');
assert(isfield(h, 'lineMag'), 'Should have lineMag handle');
assert(isfield(h, 'linePhase'), 'Should have linePhase handle');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: sidBodePlot returns correct handles.\n');

%% Test 2: Magnitude line data matches 20*log10(|G|)
h = sidBodePlot(result_siso, 'Confidence', 0);
yPlotted = get(h.lineMag, 'YData');
G = result_siso.Response(:, 1);
expectedMagDB = 20 * log10(abs(G));
relErr = max(abs(yPlotted(:) - expectedMagDB(:)) ./ max(abs(expectedMagDB(:)), 1e-300));
assert(relErr < 1e-10, 'Magnitude data mismatch (relErr=%.2e)', relErr);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: magnitude data matches 20*log10(|G|).\n');

%% Test 3: Phase line data matches angle(G) in degrees
h = sidBodePlot(result_siso, 'Confidence', 0);
yPlotted = get(h.linePhase, 'YData');
expectedPhase = angle(result_siso.Response(:, 1)) * 180 / pi;
relErr = max(abs(yPlotted(:) - expectedPhase(:)) ./ max(abs(expectedPhase(:)), 1e-300));
assert(relErr < 1e-10, 'Phase data mismatch (relErr=%.2e)', relErr);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: phase data matches angle(G) in degrees.\n');

%% Test 4: Frequency axis uses rad/s by default
h = sidBodePlot(result_siso, 'Confidence', 0);
xPlotted = get(h.lineMag, 'XData');
expectedFreq = result_siso.Frequency / result_siso.SampleTime;
relErr = max(abs(xPlotted(:) - expectedFreq(:)) ./ max(abs(expectedFreq(:)), 1e-300));
assert(relErr < 1e-10, 'Frequency axis rad/s mismatch (relErr=%.2e)', relErr);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: frequency axis uses rad/s by default.\n');

%% Test 5: Frequency axis uses Hz when requested
h = sidBodePlot(result_siso, 'Confidence', 0, 'FrequencyUnit', 'Hz');
xPlotted = get(h.lineMag, 'XData');
expectedFreqHz = result_siso.FrequencyHz;
relErr = max(abs(xPlotted(:) - expectedFreqHz(:)) ./ max(abs(expectedFreqHz(:)), 1e-300));
assert(relErr < 1e-10, 'Frequency axis Hz mismatch (relErr=%.2e)', relErr);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: frequency axis uses Hz when requested.\n');

%% Test 6: Both axes use log x-scale
h = sidBodePlot(result_siso);
assert(strcmp(get(h.axMag, 'XScale'), 'log'), 'Magnitude axis should be log');
assert(strcmp(get(h.axPhase, 'XScale'), 'log'), 'Phase axis should be log');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: both axes use log x-scale.\n');

%% Test 7: Confidence bands present when Confidence > 0
h = sidBodePlot(result_siso, 'Confidence', 3);
kids = get(h.axMag, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
hasPatch = any(strcmp(types, 'patch'));
assert(hasPatch, 'Magnitude axis should have confidence band (patch)');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: confidence bands present when Confidence > 0.\n');

%% Test 8: No confidence bands when Confidence = 0
h = sidBodePlot(result_siso, 'Confidence', 0);
kids = get(h.axMag, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
hasPatch = any(strcmp(types, 'patch'));
assert(~hasPatch, 'No confidence band expected when Confidence=0');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: no confidence bands when Confidence = 0.\n');

%% Test 9: sidBodePlot errors on time-series result
threw = false;
try
    sidBodePlot(result_ts);
catch e
    if strcmp(e.identifier, 'sid:noResponse')
        threw = true;
    end
end
assert(threw, 'Should have thrown sid:noResponse');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: sidBodePlot errors on time-series result.\n');

%% Test 10: sidBodePlot with custom axes
fig = figure;
ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
h = sidBodePlot(result_siso, 'Axes', [ax1 ax2], 'Confidence', 0);
assert(h.axMag == ax1, 'Should use provided axMag');
assert(h.axPhase == ax2, 'Should use provided axPhase');
close(fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: sidBodePlot with custom axes.\n');

%% Test 11: sidSpectrumPlot data matches 10*log10(NoiseSpectrum)
h = sidSpectrumPlot(result_siso, 'Confidence', 0);
yPlotted = get(h.line, 'YData');
PhiV = result_siso.NoiseSpectrum(:, 1);
expectedDB = 10 * log10(max(PhiV, 1e-20));
relErr = max(abs(yPlotted(:) - expectedDB(:)) ./ max(abs(expectedDB(:)), 1e-300));
assert(relErr < 1e-6, 'Spectrum data mismatch (relErr=%.2e)', relErr);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: spectrum data matches 10*log10(NoiseSpectrum).\n');

%% Test 12: sidSpectrumPlot uses log x-scale
h = sidSpectrumPlot(result_siso);
assert(strcmp(get(h.ax, 'XScale'), 'log'), 'Spectrum x-axis should be log');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 12 passed: sidSpectrumPlot uses log x-scale.\n');

%% Test 13: sidSpectrumPlot works for time-series
h = sidSpectrumPlot(result_ts);
nLines = numel(findobj(h.ax, 'Type', 'line'));
assert(nLines > 0, 'Time-series spectrum should have a plotted line');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 13 passed: sidSpectrumPlot works for time-series.\n');

%% Test 14: sidSpectrumPlot confidence bands
h = sidSpectrumPlot(result_siso, 'Confidence', 3);
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
hasPatch = any(strcmp(types, 'patch'));
assert(hasPatch, 'Spectrum should have confidence band (patch)');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 14 passed: sidSpectrumPlot confidence bands present.\n');

%% Test 15: sidBodePlot with ETFE result
result_etfe = sidFreqETFE(y, u);
h = sidBodePlot(result_etfe, 'Confidence', 0);
yPlotted = get(h.lineMag, 'YData');
G_etfe = result_etfe.Response(:, 1);
expectedMagDB = 20 * log10(abs(G_etfe));
relErr = max(abs(yPlotted(:) - expectedMagDB(:)) ./ max(abs(expectedMagDB(:)), 1e-300));
assert(relErr < 1e-10, 'ETFE magnitude data mismatch (relErr=%.2e)', relErr);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 15 passed: sidBodePlot with ETFE result data validated.\n');

fprintf('test_sidPlotting: %d/%d passed\n', runner__nPassed, runner__nPassed);
