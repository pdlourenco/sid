%% test_sidPlotting - Unit tests for sidBodePlot and sidSpectrumPlot
%
% Tests that plotting functions run without errors, return correct handles,
% and handle edge cases. Uses -nodisplay headless mode.

fprintf('Running test_sidPlotting...\n');

%% Setup: create test result structs
rng(42);
N = 500;
u = randn(N, 1);
y = filter([1], [1 -0.8], u) + 0.1 * randn(N, 1);
result_siso = sidFreqBT(y, u);
result_ts = sidFreqBT(randn(300, 1), []);

%% Test 1: sidBodePlot runs without error for SISO
try
    h = sidBodePlot(result_siso);
    assert(isfield(h, 'fig'), 'Should have fig handle');
    assert(isfield(h, 'axMag'), 'Should have axMag handle');
    assert(isfield(h, 'axPhase'), 'Should have axPhase handle');
    assert(isfield(h, 'lineMag'), 'Should have lineMag handle');
    assert(isfield(h, 'linePhase'), 'Should have linePhase handle');
    close(h.fig);
catch e
    % In headless environments, figure creation may fail
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 2: sidBodePlot errors on time-series result
try
    sidBodePlot(result_ts);
    error('Should have thrown sid:noResponse');
catch e
    if strcmp(e.identifier, 'sid:noResponse')
        % Expected
    else
        % Might be a display error in headless mode
        if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
                && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
            assert(strcmp(e.identifier, 'sid:noResponse'), 'Expected sid:noResponse error');
        end
    end
end

%% Test 3: sidSpectrumPlot runs without error
try
    h = sidSpectrumPlot(result_ts);
    assert(isfield(h, 'fig'), 'Should have fig handle');
    assert(isfield(h, 'ax'), 'Should have ax handle');
    assert(isfield(h, 'line'), 'Should have line handle');
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 4: sidSpectrumPlot with SISO result (noise spectrum)
try
    h = sidSpectrumPlot(result_siso);
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 5: sidBodePlot with Hz frequency unit
try
    h = sidBodePlot(result_siso, 'FrequencyUnit', 'Hz');
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 6: sidBodePlot with no confidence band
try
    h = sidBodePlot(result_siso, 'Confidence', 0);
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 7: sidBodePlot with custom color and line width
try
    h = sidBodePlot(result_siso, 'Color', [1 0 0], 'LineWidth', 2);
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 8: Plotting ETFE result
result_etfe = sidFreqETFE(y, u);
try
    h = sidBodePlot(result_etfe);
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

%% Test 9: Plotting BTFDR result
result_btfdr = sidFreqBTFDR(y, u);
try
    h = sidBodePlot(result_btfdr);
    close(h.fig);
catch e
    if ~contains(e.message, 'figure') && ~contains(e.message, 'display') ...
            && ~contains(e.message, 'DISPLAY') && ~contains(e.message, 'java')
        rethrow(e);
    end
    fprintf('    (Skipping display-dependent test in headless mode)\n');
end

fprintf('  test_sidPlotting: ALL PASSED\n');
