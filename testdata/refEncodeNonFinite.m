function json = refEncodeNonFinite(json)
%REFENCODENONFINITE Store NaN/Inf/-Inf in reference JSON as strings (ADR-0007).
%
%   json = refEncodeNonFinite(json)
%
%   Takes the output of jsonencode(data, 'ConvertInfAndNaN', false), where
%   non-finite numbers appear as the bare tokens NaN, Infinity and -Infinity
%   (not valid JSON), and rewrites each token in value position as the JSON
%   string "NaN", "Inf" or "-Inf". Errors if any bare token or a null is left,
%   so a vector can never be written with an ambiguous non-finite value
%   (jsonencode's default writes NaN, Inf and -Inf all as null).
%
%   The inverse is refDecodeNonFinite. Both are used by generate_reference.m
%   and validate_reference.m; the validator's self-test round-trips them.

    json = regexprep(json, '([\[,:])(-?)Infinity(?=[,\]\}])', '$1"$2Inf"');
    json = regexprep(json, '([\[,:])NaN(?=[,\]\}])', '$1"NaN"');
    if ~isempty(regexp(json, '[\[,:](null|NaN|-?Infinity)[,\]\}]', 'once'))
        error('refEncodeNonFinite:unencoded', ...
            'Reference JSON still holds a null or a bare NaN/Infinity token.');
    end
end
