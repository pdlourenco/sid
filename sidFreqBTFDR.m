function result = sidFreqBTFDR(y, u, varargin)
%SIDFREQBTFDR Blackman-Tukey spectral analysis with frequency-dependent resolution.
%
%   result = sidFreqBTFDR(y, u)
%   result = sidFreqBTFDR(y, u, 'Resolution', R, 'Frequencies', w)
%
%   Like sidFreqBT, but the window size varies across frequencies.
%   The user specifies a resolution parameter (in rad/sample) instead
%   of a fixed window size.
%
%   This is an open-source replacement for the System Identification
%   Toolbox function 'spafdr'.
%
%   See SPEC.md section 5 for algorithm details.
%
%   NOT YET IMPLEMENTED — placeholder for Phase 7.
%
%   See also: sidFreqBT, sidFreqETFE

    error('sid:notImplemented', 'sidFreqBTFDR is not yet implemented.');
end
