%% test_sidMapPlot - Unit tests for sidMapPlot
%
% Validates plotted data content, axis scales, and error handling.
% Uses xvfb in CI for headless execution.

fprintf('Running test_sidMapPlot...\n');
runner__nPassed = 0;

%% Setup: create test result structs
rng(42);
N = 2000;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);
result_siso = sidFreqMap(y, u, 'SegmentLength', 256);
result_ts = sidFreqMap(randn(1000, 1), [], 'SegmentLength', 128);

%% Test 1: Returns correct handles
h = sidMapPlot(result_siso);
assert(isfield(h, 'fig'), 'Should have fig handle');
assert(isfield(h, 'ax'), 'Should have ax handle');
assert(isfield(h, 'surf'), 'Should have surf handle');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: returns correct handles.\n');

%% Test 2: Frequency axis is log scale
h = sidMapPlot(result_siso);
assert(strcmp(get(h.ax, 'YScale'), 'log'), 'Frequency axis should be log');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: frequency axis is log scale.\n');

%% Test 3: Pcolormesh present for magnitude plot
h = sidMapPlot(result_siso, 'PlotType', 'magnitude');
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
hasSurface = any(strcmp(types, 'surface'));
assert(hasSurface, 'Should have pcolor surface object');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: pcolormesh present for magnitude plot.\n');

%% Test 4: Phase plot runs and has surface
h = sidMapPlot(result_siso, 'PlotType', 'phase');
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
assert(any(strcmp(types, 'surface')), 'Phase plot should have surface');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: phase plot has surface.\n');

%% Test 5: Coherence plot runs
h = sidMapPlot(result_siso, 'PlotType', 'coherence');
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
assert(any(strcmp(types, 'surface')), 'Coherence plot should have surface');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: coherence plot has surface.\n');

%% Test 6: Noise plot for time series
h = sidMapPlot(result_ts, 'PlotType', 'noise');
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
assert(any(strcmp(types, 'surface')), 'Noise plot should have surface');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: noise plot for time series.\n');

%% Test 7: Hz frequency unit
h = sidMapPlot(result_siso, 'FrequencyUnit', 'Hz');
yLabel = get(get(h.ax, 'YLabel'), 'String');
assert(~isempty(strfind(lower(yLabel), 'hz')), ...
    'Y-label should mention Hz, got ''%s''', yLabel);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: Hz frequency unit label.\n');

%% Test 8: Custom CLim applied
h = sidMapPlot(result_siso, 'CLim', [-40 10]);
clim = get(h.ax, 'CLim');
assert(abs(clim(1) - (-40)) < 1e-10 && abs(clim(2) - 10) < 1e-10, ...
    'CLim should be [-40 10], got [%.1f %.1f]', clim(1), clim(2));
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: custom CLim applied.\n');

%% Test 9: Error on invalid result struct
threw = false;
try
    sidMapPlot(struct('Method', 'wrong'));
catch e
    if strcmp(e.identifier, 'sid:invalidResult')
        threw = true;
    end
end
assert(threw, 'Should have thrown sid:invalidResult');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: error on invalid result struct.\n');

%% Test 10: Error on magnitude plot for time series
threw = false;
try
    sidMapPlot(result_ts, 'PlotType', 'magnitude');
catch e
    if strcmp(e.identifier, 'sid:noResponse')
        threw = true;
    end
end
assert(threw, 'Should have thrown sid:noResponse');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: error on magnitude for time series.\n');

%% Test 11: Error on invalid PlotType
threw = false;
try
    sidMapPlot(result_siso, 'PlotType', 'invalid');
catch e
    if strcmp(e.identifier, 'sid:invalidPlotType')
        threw = true;
    end
end
assert(threw, 'Should have thrown sid:invalidPlotType');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: error on invalid PlotType.\n');

fprintf('test_sidMapPlot: %d/%d passed\n', runner__nPassed, runner__nPassed);
