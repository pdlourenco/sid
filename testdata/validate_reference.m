function validate_reference()
%VALIDATE_REFERENCE Validate sid outputs against reference data.
%
%   validate_reference
%
%   Reads each reference_*.json file, calls the corresponding sid function
%   with the stored input data, and verifies outputs match within tolerance.
%   Designed to run under any engine (Octave, Python, Julia) to validate
%   cross-language numerical equivalence.
%
%   Usage:
%     run('testdata/validate_reference.m')

fprintf('=== Cross-language reference validation ===\n\n');

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
sidDir = fullfile(rootDir, 'matlab', 'sid');
addpath(sidDir);
% MATLAB ignores addpath on directories named 'private'. Copy to a
% temporary non-private-named directory so private helpers are accessible.
privateDir = fullfile(sidDir, 'private');
shimDir = fullfile(thisDir, 'private_shim');
if exist(shimDir, 'dir')
    rmdir(shimDir, 's');
end
mkdir(shimDir);
copyfile(fullfile(privateDir, '*.m'), shimDir);
addpath(shimDir);
cleanupObj = onCleanup(@() cleanupShim(shimDir));

files = dir(fullfile(thisDir, 'reference_*.json'));
if isempty(files)
    fprintf('No reference JSON files found in %s\n', thisDir);
    fprintf('Run generate_reference.m under MATLAB first.\n');
    fprintf('\n=== SKIP (no reference data) ===\n');
    return;
end

nPass = 0;
nFail = 0;
failures = {};

for i = 1:numel(files)
    name = files(i).name;
    filepath = fullfile(thisDir, name);
    fprintf('  %s\n', name);

    ref = jsondecode(fileread(filepath));

    % Call the sid function with stored inputs and params
    result = callSidFunction(ref.function_name, ref.input, ref.params);

    % Compare each output field within tolerance
    [ok, msgs] = compareOutputs(result, ref.output, ref.tolerance);

    if ok
        nPass = nPass + 1;
        fprintf('    PASS\n');
    else
        nFail = nFail + 1;
        fprintf('    FAIL\n');
        for m = 1:numel(msgs)
            fprintf('    %s\n', msgs{m});
        end
        failures{end+1} = name;
    end
end

fprintf('\n=== %d passed, %d failed ===\n', nPass, nFail);
if nFail > 0
    error('Validation failed for: %s', strjoin(failures, ', '));
end

end


% -----------------------------------------------------------------------
%  Local functions
% -----------------------------------------------------------------------

function cleanupShim(shimDir)
%CLEANUPSHIM Remove the temporary shim directory from the path and disk.
    if exist(shimDir, 'dir')
        rmpath(shimDir);
        rmdir(shimDir, 's');
    end
end


function result = callSidFunction(funcName, input, params)
%CALLSIDFUNCTION Dispatch to the right sid function with stored args.
    args = structToNameValue(params);

    switch funcName
        case {'sidFreqBT', 'sidFreqETFE', 'sidFreqBTFDR'}
            result = feval(funcName, input.y, input.u, args{:});
        case 'sidSpectrogram'
            result = sidSpectrogram(input.x, args{:});
        case 'sidFreqMap'
            result = sidFreqMap(input.y, input.u, args{:});
        case 'sidLTVdisc'
            result = feval(funcName, input.X, input.U, args{:});
        case 'internals'
            % Internal helpers validated separately — skip dispatch
            result = struct();
        otherwise
            error('Unknown function: %s', funcName);
    end
end


function args = structToNameValue(s)
%STRUCTTONAMEVALUE Convert struct fields to {'Name', value, ...} cell array.
    fields = fieldnames(s);
    args = cell(1, 2 * numel(fields));
    for i = 1:numel(fields)
        args{2*i - 1} = fields{i};
        args{2*i}     = s.(fields{i});
    end
end


function [ok, messages] = compareOutputs(result, expected, tolerance)
%COMPAREOUTPUTS Check each expected output field against actual result.
    ok = true;
    messages = {};

    fields = fieldnames(expected);
    for i = 1:numel(fields)
        name = fields{i};
        expVal = expected.(name);

        % Map JSON output fields to result struct fields.
        % Response is stored as separate _real / _imag components.
        if strcmp(name, 'Response_real')
            actVal = real(result.Response);
            tolField = 'Response_rel';
        elseif strcmp(name, 'Response_imag')
            actVal = imag(result.Response);
            tolField = 'Response_rel';
        else
            if ~isfield(result, name)
                ok = false;
                messages{end+1} = sprintf('  %s: field missing from result', name);
                continue;
            end
            actVal = result.(name);
            tolField = [name '_rel'];
        end

        % Flatten to vectors — handles plain arrays and cell arrays
        expVec = flatten(expVal);
        actVec = flatten(actVal);

        if numel(expVec) ~= numel(actVec)
            ok = false;
            messages{end+1} = sprintf('  %s: size mismatch (%d vs %d)', ...
                name, numel(expVec), numel(actVec));
            continue;
        end

        % Relative error (guard against division by zero)
        denom = max(abs(expVec), 1e-300);
        relErr = max(abs(actVec - expVec) ./ denom);

        % Look up tolerance (default 1e-6 for cross-engine comparison)
        if isfield(tolerance, tolField)
            tol = tolerance.(tolField);
        else
            tol = 1e-6;
        end

        if relErr > tol
            ok = false;
            messages{end+1} = sprintf( ...
                '  %s: max relative error %.2e exceeds tolerance %.2e', ...
                name, relErr, tol);
        else
            fprintf('    %s: max relative error %.2e (tol %.2e)\n', ...
                name, relErr, tol);
        end
    end
end


function v = flatten(val)
%FLATTEN Convert any numeric value (array, cell array) to a column vector.
    if iscell(val)
        v = cell2mat(cellfun(@(x) x(:), val(:), 'UniformOutput', false));
    else
        v = val(:);
    end
end
