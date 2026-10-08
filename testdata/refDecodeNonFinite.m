function json = refDecodeNonFinite(json)
%REFDECODENONFINITE Map "NaN"/"Inf"/"-Inf" strings in reference JSON to numbers.
%
%   json = refDecodeNonFinite(json)
%
%   Inverse of refEncodeNonFinite (ADR-0007): rewrites each JSON string
%   "NaN", "Inf" or "-Inf" in value position as the bare token NaN, Infinity
%   or -Infinity, which jsondecode (MATLAB and Octave) parses as a number in
%   place, so arrays keep their shape. Decoding the strings directly would
%   turn every array holding one into a cell array.

    json = regexprep(json, '([\[,:]\s*)"(-?)Inf"(?=\s*[,\]\}])', '$1$2Infinity');
    json = regexprep(json, '([\[,:]\s*)"NaN"(?=\s*[,\]\}])', '$1NaN');
end
