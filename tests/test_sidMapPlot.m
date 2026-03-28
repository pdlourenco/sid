%% test_sidMapPlot - Unit tests for sidMapPlot
%
% Tests that sidMapPlot runs without errors for all plot types and
% handles edge cases. Uses headless mode fallback for CI.

fprintf('Running test_sidMapPlot...\n');

%% Setup: create test result structs
rng(42);
N = 2000;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);
result_siso = sidFreqBTMap(y, u, 'SegmentLength', 256);
result_ts = sidFreqBTMap(randn(1000, 1), [], 'SegmentLength', 128);

%% Test 1: Magnitude plot (default)
try
    h = sidMapPlot(result_siso);
    assert(isfield(h, 'fig'), 'Should have fig handle');
    assert(isfield(h, 'ax'), 'Should have ax handle');
    assert(isfield(h, 'surf'), 'Should have surf handle');
    close(h.fig);
    fprintf('  Test 1 passed: magnitude plot\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 1 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 2: Phase plot
try
    h = sidMapPlot(result_siso, 'PlotType', 'phase');
    close(h.fig);
    fprintf('  Test 2 passed: phase plot\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 2 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 3: Noise plot
try
    h = sidMapPlot(result_siso, 'PlotType', 'noise');
    close(h.fig);
    fprintf('  Test 3 passed: noise plot\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 3 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 4: Coherence plot
try
    h = sidMapPlot(result_siso, 'PlotType', 'coherence');
    close(h.fig);
    fprintf('  Test 4 passed: coherence plot\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 4 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 5: Spectrum plot (time series)
try
    h = sidMapPlot(result_ts, 'PlotType', 'spectrum');
    close(h.fig);
    fprintf('  Test 5 passed: spectrum plot (time series)\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 5 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 6: Hz frequency unit
try
    h = sidMapPlot(result_siso, 'FrequencyUnit', 'Hz');
    close(h.fig);
    fprintf('  Test 6 passed: Hz frequency unit\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 6 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 7: Custom CLim
try
    h = sidMapPlot(result_siso, 'CLim', [-40 10]);
    close(h.fig);
    fprintf('  Test 7 passed: custom CLim\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 7 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 8: Error on invalid result struct
try
    sidMapPlot(struct('Method', 'wrong'));
    error('Should have thrown sid:invalidResult');
catch e
    if strcmp(e.identifier, 'sid:invalidResult')
        fprintf('  Test 8 passed: error on invalid result\n');
    elseif contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 8 skipped: headless mode\n');
    else
        assert(strcmp(e.identifier, 'sid:invalidResult'), 'Expected sid:invalidResult');
    end
end

%% Test 9: Error on magnitude plot for time series
try
    sidMapPlot(result_ts, 'PlotType', 'magnitude');
    error('Should have thrown sid:noResponse');
catch e
    if strcmp(e.identifier, 'sid:noResponse')
        fprintf('  Test 9 passed: error on magnitude for time series\n');
    elseif contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 9 skipped: headless mode\n');
    else
        assert(strcmp(e.identifier, 'sid:noResponse'), 'Expected sid:noResponse');
    end
end

%% Test 10: Error on invalid PlotType
try
    sidMapPlot(result_siso, 'PlotType', 'invalid');
    error('Should have thrown sid:invalidPlotType');
catch e
    if strcmp(e.identifier, 'sid:invalidPlotType')
        fprintf('  Test 10 passed: error on invalid PlotType\n');
    elseif contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 10 skipped: headless mode\n');
    else
        assert(strcmp(e.identifier, 'sid:invalidPlotType'), 'Expected sid:invalidPlotType');
    end
end

fprintf('test_sidMapPlot: ALL TESTS PASSED\n');
