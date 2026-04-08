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
    if isfield(ref, 'params')
        params = ref.params;
    else
        params = struct();
    end
    result = callSidFunction(ref.function_name, ref.input, params);

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
            if isfield(input, 'u')
                u = input.u;
            else
                u = [];
            end
            result = feval(funcName, input.y, u, args{:});
        case 'sidLTVdisc'
            result = feval(funcName, input.X, input.U, args{:});
        case 'internals'
            x = input.x; z = input.z; M = input.M;
            R_xx = sidCov(x, x, M);
            R_xz = sidCov(x, z, M);
            W    = sidHannWin(M);
            freqs = (1:128)' * pi / 128;
            Phi_xx = sidWindowedDFT(R_xx, W, freqs, true, R_xx);
            R_zx   = sidCov(z, x, M);
            Phi_xz = sidWindowedDFT(R_xz, W, freqs, true, R_zx);
            result = struct( ...
                'R_xx', R_xx, 'R_xz', R_xz, 'W', W, ...
                'Phi_xx_real', real(Phi_xx), ...
                'Phi_xx_imag', imag(Phi_xx), ...
                'Phi_xz_real', real(Phi_xz), ...
                'Phi_xz_imag', imag(Phi_xz));
        otherwise
            error('Unknown function: %s', funcName);
    end
end


function args = structToNameValue(s)
%STRUCTTONAMEVALUE Convert struct fields to {'Name', value, ...} cell array.
    fields = fieldnames(s);
    % Filter out metadata-only fields not accepted by sid functions
    meta = {'mode'};
    fields = fields(~ismember(fields, meta));
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

        % Look up tolerances (default rtol 1e-6, atol 0 for cross-engine)
        if isfield(tolerance, tolField)
            rtol = tolerance.(tolField);
        else
            rtol = 1e-6;
        end

        atolField = strrep(tolField, '_rel', '_atol');
        if isfield(tolerance, atolField)
            atol = tolerance.(atolField);
        else
            atol = 0;
        end

        % allclose check: |actual - expected| <= atol + rtol * |expected|
        absDiff = abs(actVec - expVec);
        thresh  = atol + rtol * abs(expVec);
        worstIdx = find(absDiff - thresh == max(absDiff - thresh), 1);

        % Report the effective relative error for logging (using the
        % allclose denominator so small entries don't dominate).
        denom = max(abs(expVec), max(atol, 1e-300));
        relErr = max(absDiff ./ denom);

        if any(absDiff > thresh)
            ok = false;
            messages{end+1} = sprintf( ...
                '  %s: element %d |diff|=%.2e exceeds atol(%.0e)+rtol(%.0e)*|exp|(%.2e) = %.2e', ...
                name, worstIdx, absDiff(worstIdx), atol, rtol, ...
                abs(expVec(worstIdx)), thresh(worstIdx));
        else
            fprintf('    %s: max relative error %.2e (rtol %.0e, atol %.0e)\n', ...
                name, relErr, rtol, atol);
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
