%% test_template - Template for sid test files
%
% Copy this file to create new tests. Follow the structure below so the
% test runner (runAllTests.m) can track individual test case results.
%
% Convention:
%   - Initialize runner__nPassed = 0 at the top
%   - Increment runner__nPassed after each passing test case
%   - Print 'Test N passed: description' for every test case
%   - End with 'test_<name>: N/N passed' summary line
%   - Use assert() for test assertions (throws on failure, stops the file)

fprintf('Running test_template...\n');
runner__nPassed = 0;

%% Test 1: Example - basic assertion
x = 1 + 1;
assert(x == 2, 'Expected 2, got %d', x);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: basic assertion.\n');

%% Test 2: Example - tolerance check
actual = 0.1 + 0.2;
expected = 0.3;
assert(abs(actual - expected) < 1e-10, ...
    'Expected %.15g, got %.15g', expected, actual);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: tolerance check.\n');

fprintf('test_template: %d/%d passed\n', runner__nPassed, runner__nPassed);
