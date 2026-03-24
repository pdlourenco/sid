function result = sidFreqETFE(y, u, varargin)
%SIDFREQETFE Empirical transfer function estimate.
%
%   result = sidFreqETFE(y, u)
%   result = sidFreqETFE(y, u, 'Smoothing', S)
%
%   Estimates the frequency response as the ratio of output and input
%   discrete Fourier transforms. Provides maximum frequency resolution
%   but high variance. Optional smoothing reduces variance.
%
%   This is an open-source replacement for the System Identification
%   Toolbox function 'etfe'.
%
%   See SPEC.md section 4 for algorithm details.
%
%   NOT YET IMPLEMENTED — placeholder for Phase 7.
%
%   See also: sidFreqBT, sidFreqBTFDR

    error('sid:notImplemented', 'sidFreqETFE is not yet implemented.');
end
