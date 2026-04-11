function tf = sidIsDefaultFreqs(freqs, nf)
% SIDISDEFAULTFREQS Check if frequency vector matches the default grid.
%
%   tf = sidIsDefaultFreqs(freqs, nf)
%
%   Returns true if freqs matches the default 128-point linearly spaced
%   grid in (0, pi]. Used to decide whether the FFT fast path can be
%   applied.
%
%   INPUTS:
%     freqs - (nf x 1) frequency vector in rad/sample.
%     nf    - Length of the frequency vector.
%
%   OUTPUTS:
%     tf - Logical scalar, true if freqs matches the default grid.
%
%   EXAMPLES:
%     freqs = (1:128)' * pi / 128;
%     sidIsDefaultFreqs(freqs, 128)  % returns true
%     sidIsDefaultFreqs([0.1; 0.5], 2)  % returns false
%
%   SPECIFICATION:
%     SPEC.md §2.2 — Default Frequency Grid
%
%   See also: sidFreqBT, sidFreqETFE
%
%   Changelog:
%   2026-04-04: Extracted from sidFreqBT and sidFreqETFE.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid
%  -----------------------------------------------------------------------

    if nf ~= 128
        tf = false;
        return;
    end
    defaultFreqs = (1:128)' * pi / 128;
    tf = max(abs(freqs(:) - defaultFreqs)) < 1e-12;
end
