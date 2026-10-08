function bad = refAtolViolations(output, tolerance)
%REFATOLVIOLATIONS Output fields that need an absolute tolerance floor but lack one.
%
%   bad = refAtolViolations(output, tolerance)
%
%   Standing rule 1 (ADR-0002, carried forward by ADR-0007): a relative-only
%   comparison is vacuous on near-zero entries, so an output field that holds
%   an exact zero, or whose smallest finite magnitude is below 1e-12 of its
%   largest, must carry <key>_atol. <key> resolves as the validators resolve
%   it: Response_real and Response_imag share Response_*, every other field
%   uses its own name. Returns the offending field names (cell array).
%
%   Used by generate_reference.m (refuses to write) and validate_reference.m
%   (fails the vector); python/tests/test_cross_validation.py mirrors it.

    bad = {};
    fields = fieldnames(output);
    for i = 1:numel(fields)
        name = fields{i};
        if any(strcmp(name, {'Response_real', 'Response_imag'}))
            key = 'Response';
        else
            key = name;
        end
        v = output.(name);
        if iscell(v) && all(cellfun(@(x) isnumeric(x) || islogical(x), v(:)))
            v = cell2mat(cellfun(@(x) double(x(:)), v(:), 'UniformOutput', false));
        end
        if ~(isnumeric(v) || islogical(v))
            error('refAtolViolations:unsupported', ...
                'Output field %s is not numeric; the floor rule cannot check it.', name);
        end
        v = abs(double(v(:)));
        v = v(isfinite(v));
        if isempty(v)
            continue;
        end
        nearZero = min(v) == 0 || min(v) < 1e-12 * max(v);
        if nearZero && ~isfield(tolerance, [key '_atol'])
            bad{end + 1} = name; %#ok<AGROW>
        end
    end
end
