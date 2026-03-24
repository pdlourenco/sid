%% test_sidValidate - Unit tests for input parsing and validation
%
% Tests sidValidate and sidValidateData for correct error handling,
% default values, and name-value parsing.

fprintf('Running test_sidValidate...\n');

%% Test 1: Default window size
y = randn(1000, 1);
[~, ~, M, freqs, Ts, iTS] = sidValidate(y, []);
assert(M == 30, 'Default M should be min(floor(N/10), 30) = 30 for N=1000');
assert(length(freqs) == 128, 'Default should be 128 frequencies');
assert(Ts == 1.0, 'Default Ts should be 1.0');
assert(iTS == true, 'Empty u means time series');

%% Test 2: Default window size for short data
y = randn(100, 1);
[~, ~, M] = sidValidate(y, []);
assert(M == 10, 'Default M should be floor(100/10) = 10 for N=100');

%% Test 3: Positional syntax
y = randn(200, 1);
u = randn(200, 1);
[~, ~, M, freqs] = sidValidate(y, u, 15);
assert(M == 15, 'Positional M should be 15');

%% Test 4: Positional with frequencies
w = [0.1; 0.5; 1.0];
[~, ~, M, freqs] = sidValidate(y, u, 20, w);
assert(M == 20, 'Positional M should be 20');
assert(isequal(freqs, w), 'Positional freqs should match');

%% Test 5: Name-value syntax
[~, ~, M, freqs, Ts] = sidValidate(y, u, 'WindowSize', 25, 'SampleTime', 0.01);
assert(M == 25, 'Name-value M should be 25');
assert(Ts == 0.01, 'Name-value Ts should be 0.01');

%% Test 6: Row vector gets converted to column
y_row = [1 2 3 4 5 6 7 8 9 10];
[y_out, ~] = sidValidate(y_row, []);
assert(isequal(size(y_out), [10, 1]), 'Row vector should be converted to column');

%% Test 7: Error on too-short data
try
    sidValidate([1], []);
    error('Should have thrown sid:tooShort');
catch e
    assert(strcmp(e.identifier, 'sid:tooShort'), 'Expected sid:tooShort error');
end

%% Test 8: Error on complex data
try
    sidValidate([1+1j; 2; 3], []);
    error('Should have thrown sid:complexData');
catch e
    assert(strcmp(e.identifier, 'sid:complexData'), 'Expected sid:complexData error');
end

%% Test 9: Error on NaN data
try
    sidValidate([1; NaN; 3], []);
    error('Should have thrown sid:nonFinite');
catch e
    assert(strcmp(e.identifier, 'sid:nonFinite'), 'Expected sid:nonFinite error');
end

%% Test 10: Error on Inf data
try
    sidValidate([1; Inf; 3], []);
    error('Should have thrown sid:nonFinite');
catch e
    assert(strcmp(e.identifier, 'sid:nonFinite'), 'Expected sid:nonFinite error');
end

%% Test 11: Error on size mismatch
try
    sidValidate(randn(100,1), randn(50,1));
    error('Should have thrown sid:sizeMismatch');
catch e
    assert(strcmp(e.identifier, 'sid:sizeMismatch'), 'Expected sid:sizeMismatch error');
end

%% Test 12: Error on bad frequencies
try
    sidValidate(randn(100,1), [], 10, [-0.1; 1.0]);
    error('Should have thrown sid:badFreqs');
catch e
    assert(strcmp(e.identifier, 'sid:badFreqs'), 'Expected sid:badFreqs error');
end

%% Test 13: Error on bad window size
try
    sidValidate(randn(100,1), [], 1);
    error('Should have thrown sid:badWindowSize');
catch e
    assert(strcmp(e.identifier, 'sid:badWindowSize'), 'Expected sid:badWindowSize error');
end

%% Test 14: Warning on window size exceeding N/2
lastwarn('');
y = randn(20, 1);
[~, ~, M] = sidValidate(y, [], 15);
[warnMsg, warnId] = lastwarn;
assert(M == 10, 'M should be reduced to N/2 = 10');
assert(strcmp(warnId, 'sid:windowReduced'), 'Should warn about window reduction');

%% Test 15: Error on negative sample time
try
    sidValidate(randn(100,1), [], 'SampleTime', -1);
    error('Should have thrown sid:badTs');
catch e
    assert(strcmp(e.identifier, 'sid:badTs'), 'Expected sid:badTs error');
end

%% Test 16: Error on unknown option
try
    sidValidate(randn(100,1), [], 'BadOption', 5);
    error('Should have thrown sid:unknownOption');
catch e
    assert(strcmp(e.identifier, 'sid:unknownOption'), 'Expected sid:unknownOption error');
end

%% Test 17: sidValidateData basic functionality
[y_out, u_out, N, ny, nu, iTS] = sidValidateData(randn(50, 2), randn(50, 3));
assert(N == 50, 'N should be 50');
assert(ny == 2, 'ny should be 2');
assert(nu == 3, 'nu should be 3');
assert(iTS == false, 'Not a time series');

%% Test 18: sidValidateData time series
[~, ~, ~, ~, nu, iTS] = sidValidateData(randn(50, 1), []);
assert(nu == 0, 'nu should be 0 for time series');
assert(iTS == true, 'Should be time series');

fprintf('  test_sidValidate: ALL PASSED\n');
