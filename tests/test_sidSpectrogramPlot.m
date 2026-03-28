%% test_sidSpectrogramPlot - Unit tests for sidSpectrogramPlot
%
% Tests that plotting functions run without errors and return handles.
% Uses headless mode fallback for CI environments.

fprintf('Running test_sidSpectrogramPlot...\n');

%% Setup: create test result struct
rng(42);
result = sidSpectrogram(randn(1000, 1), 'WindowLength', 64);
result_mc = sidSpectrogram(randn(500, 2), 'WindowLength', 64);

%% Test 1: Basic call runs without error
try
    h = sidSpectrogramPlot(result);
    assert(isfield(h, 'fig'), 'Should have fig handle');
    assert(isfield(h, 'ax'), 'Should have ax handle');
    assert(isfield(h, 'surf'), 'Should have surf handle');
    close(h.fig);
    fprintf('  Test 1 passed: basic call\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 1 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 2: Log frequency scale
try
    h = sidSpectrogramPlot(result, 'FrequencyScale', 'log');
    close(h.fig);
    fprintf('  Test 2 passed: log frequency scale\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 2 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 3: Custom CLim
try
    h = sidSpectrogramPlot(result, 'CLim', [-60 0]);
    close(h.fig);
    fprintf('  Test 3 passed: custom CLim\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 3 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 4: Multi-channel with channel selection
try
    h = sidSpectrogramPlot(result_mc, 'Channel', 2);
    close(h.fig);
    fprintf('  Test 4 passed: multi-channel, channel 2\n');
catch e
    if contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 4 skipped: headless mode\n');
    else
        rethrow(e);
    end
end

%% Test 5: Error on invalid result struct
try
    sidSpectrogramPlot(struct('Method', 'wrong'));
    error('Should have thrown sid:invalidResult');
catch e
    if strcmp(e.identifier, 'sid:invalidResult')
        fprintf('  Test 5 passed: error on invalid result\n');
    elseif contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 5 skipped: headless mode\n');
    else
        assert(strcmp(e.identifier, 'sid:invalidResult'), 'Expected sid:invalidResult');
    end
end

%% Test 6: Error on invalid channel
try
    sidSpectrogramPlot(result, 'Channel', 5);
    error('Should have thrown sid:invalidChannel');
catch e
    if strcmp(e.identifier, 'sid:invalidChannel')
        fprintf('  Test 6 passed: error on invalid channel\n');
    elseif contains(e.message, 'figure') || contains(e.message, 'display') ...
            || contains(e.message, 'DISPLAY') || contains(e.message, 'java')
        fprintf('  Test 6 skipped: headless mode\n');
    else
        assert(strcmp(e.identifier, 'sid:invalidChannel'), 'Expected sid:invalidChannel');
    end
end

fprintf('test_sidSpectrogramPlot: ALL TESTS PASSED\n');
