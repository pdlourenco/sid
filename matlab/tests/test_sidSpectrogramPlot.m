%% test_sidSpectrogramPlot - Unit tests for sidSpectrogramPlot
%
% Validates plotted data, axis scales, labels, and error handling.
% Uses xvfb in CI for headless execution.

fprintf('Running test_sidSpectrogramPlot...\n');
runner__nPassed = 0;

%% Setup: create test result struct
rng(42);
result = sidSpectrogram(randn(1000, 1), 'WindowLength', 64);
result_mc = sidSpectrogram(randn(500, 2), 'WindowLength', 64);

%% Test 1: Returns correct handles
h = sidSpectrogramPlot(result);
assert(isfield(h, 'fig'), 'Should have fig handle');
assert(isfield(h, 'ax'), 'Should have ax handle');
assert(isfield(h, 'surf'), 'Should have surf handle');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: returns correct handles.\n');

%% Test 2: Has pcolor surface object
h = sidSpectrogramPlot(result);
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
hasSurface = any(strcmp(types, 'surface'));
assert(hasSurface, 'Should have pcolor surface object');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: has pcolor surface object.\n');

%% Test 3: Default frequency scale is linear
h = sidSpectrogramPlot(result);
assert(strcmp(get(h.ax, 'YScale'), 'linear'), ...
    'Default frequency scale should be linear');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: default frequency scale is linear.\n');

%% Test 4: Log frequency scale when requested
h = sidSpectrogramPlot(result, 'FrequencyScale', 'log');
assert(strcmp(get(h.ax, 'YScale'), 'log'), ...
    'Frequency scale should be log when requested');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: log frequency scale when requested.\n');

%% Test 5: X-axis label contains 'time'
h = sidSpectrogramPlot(result);
xLabel = get(get(h.ax, 'XLabel'), 'String');
assert(~isempty(strfind(lower(xLabel), 'time')), ...
    'X-label should mention time, got ''%s''', xLabel);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: x-axis label contains ''time''.\n');

%% Test 6: Y-axis label contains 'freq'
h = sidSpectrogramPlot(result);
yLabel = get(get(h.ax, 'YLabel'), 'String');
assert(~isempty(strfind(lower(yLabel), 'freq')), ...
    'Y-label should mention frequency, got ''%s''', yLabel);
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: y-axis label contains ''freq''.\n');

%% Test 7: Custom CLim applied
h = sidSpectrogramPlot(result, 'CLim', [-60 0]);
clim = get(h.ax, 'CLim');
assert(abs(clim(1) - (-60)) < 1e-10 && abs(clim(2) - 0) < 1e-10, ...
    'CLim should be [-60 0], got [%.1f %.1f]', clim(1), clim(2));
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: custom CLim applied.\n');

%% Test 8: Multi-channel with channel selection
h = sidSpectrogramPlot(result_mc, 'Channel', 2);
kids = get(h.ax, 'Children');
types = arrayfun(@(c) get(c, 'Type'), kids, 'UniformOutput', false);
assert(any(strcmp(types, 'surface')), ...
    'Channel 2 plot should have surface');
close(h.fig);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: multi-channel with channel selection.\n');

%% Test 9: Error on invalid result struct
threw = false;
try
    sidSpectrogramPlot(struct('Method', 'wrong'));
catch e
    if strcmp(e.identifier, 'sid:invalidResult')
        threw = true;
    end
end
assert(threw, 'Should have thrown sid:invalidResult');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: error on invalid result struct.\n');

%% Test 10: Error on invalid channel
threw = false;
try
    sidSpectrogramPlot(result, 'Channel', 5);
catch e
    if strcmp(e.identifier, 'sid:invalidChannel')
        threw = true;
    end
end
assert(threw, 'Should have thrown sid:invalidChannel');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: error on invalid channel.\n');

fprintf('test_sidSpectrogramPlot: %d/%d passed\n', runner__nPassed, runner__nPassed);
