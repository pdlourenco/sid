function sidResultTypes()
% SIDRESULTTYPES Central reference for all sid result struct types.
%
%   sidResultTypes
%   help sidResultTypes
%
%   This file is a documentation-only reference. It lists every result
%   struct returned by the sid public API, together with its fields,
%   dimensions, and the functions that produce or consume it.
%
%   Use this file as a single lookup point when writing code that reads
%   fields from sid result structs.
%
%   See the Python equivalent in python/sid/_results.py (frozen dataclasses
%   with per-field type annotations and docstrings).
%
%   INPUTS:
%     (none — this function takes no arguments)
%
%   OUTPUTS:
%     (none — prints a help reminder to the console)
%
%   EXAMPLES:
%     help sidResultTypes     % view all result struct definitions
%     sidResultTypes          % prints a reminder to use help
%
% =========================================================================
%  1. FreqResult  (sidFreqBT, sidFreqBTFDR, sidFreqETFE)
% =========================================================================
%
%   Produced by:  sidFreqBT, sidFreqBTFDR, sidFreqETFE
%   Consumed by:  sidBodePlot, sidSpectrumPlot, sidCompare, sidResidual,
%                 sidLTVdiscTune (frequency method)
%
%   Field               Dimensions          Description
%   .................................................................
%   .Frequency          (nf x 1)            Frequency vector, rad/sample
%   .FrequencyHz        (nf x 1)            Frequency vector, Hz
%   .Response           (nf x ny x nu)      Complex frequency response.
%                                            [] in time-series mode.
%   .ResponseStd        (nf x ny x nu)      Standard deviation of Response.
%                                            [] in time-series mode.
%   .NoiseSpectrum      (nf x ny x ny)      Noise spectrum (or output
%                                            spectrum in time-series mode).
%   .NoiseSpectrumStd   (nf x ny x ny)      Standard deviation of
%                                            NoiseSpectrum.
%   .Coherence          (nf x 1)            Squared coherence (SISO only).
%                                            [] for MIMO or time-series.
%   .SampleTime         scalar               Sample time in seconds.
%   .WindowSize         scalar or (nf x 1)  Lag window size M. Scalar for
%                                            BT/ETFE; per-freq vector for
%                                            BTFDR.
%   .DataLength         scalar               Number of samples N.
%   .NumTrajectories    scalar               Number of trajectories L.
%   .Method             char                 'sidFreqBT', 'sidFreqBTFDR',
%                                            or 'sidFreqETFE'.
%
% =========================================================================
%  2. FreqMapResult  (sidFreqMap)
% =========================================================================
%
%   Produced by:  sidFreqMap
%   Consumed by:  sidMapPlot
%
%   Field               Dimensions              Description
%   .................................................................
%   .Time               (K x 1)                 Center time of each segment
%                                                (seconds).
%   .Frequency          (nf x 1)                Frequency vector, rad/sample.
%   .FrequencyHz        (nf x 1)                Frequency vector, Hz.
%   .Response           (nf x K [x ny x nu])    Time-varying complex response.
%                                                [] in time-series mode.
%   .ResponseStd        (nf x K [x ny x nu])    Std dev of Response.
%                                                [] in time-series mode.
%   .NoiseSpectrum      (nf x K [x ny x ny])    Noise spectrum.
%   .NoiseSpectrumStd   (nf x K [x ny x ny])    Std dev of NoiseSpectrum.
%   .Coherence          (nf x K)                Squared coherence (SISO).
%                                                [] for MIMO or time-series.
%   .SampleTime         scalar                   Sample time in seconds.
%   .SegmentLength      scalar                   Segment length L.
%   .Overlap            scalar                   Overlap P.
%   .WindowSize         scalar                   BT lag window M, or []
%                                                for Welch.
%   .Algorithm          char                     'bt' or 'welch'.
%   .NumTrajectories    scalar                   Number of trajectories.
%   .Method             char                     'sidFreqMap'.
%
% =========================================================================
%  3. SpectrogramResult  (sidSpectrogram)
% =========================================================================
%
%   Produced by:  sidSpectrogram
%   Consumed by:  sidSpectrogramPlot
%
%   Field               Dimensions              Description
%   .................................................................
%   .Time               (K x 1)                 Center time of each segment
%                                                (seconds).
%   .Frequency          (n_bins x 1)            Frequency vector, Hz.
%   .FrequencyRad       (n_bins x 1)            Frequency vector, rad/s.
%   .Power              (n_bins x K x n_ch)     Power spectral density.
%   .PowerDB            (n_bins x K x n_ch)     10*log10(Power).
%   .Complex            (n_bins x K x n_ch)     Complex STFT coefficients.
%   .SampleTime         scalar                   Sample time in seconds.
%   .WindowLength       scalar                   Segment length L.
%   .Overlap            scalar                   Overlap P.
%   .NFFT               scalar                   FFT length.
%   .NumTrajectories    scalar                   Number of trajectories.
%   .Method             char                     'sidSpectrogram'.
%
% =========================================================================
%  4. LTVResult  (sidLTVdisc)
% =========================================================================
%
%   Produced by:  sidLTVdisc, sidLTVdiscTune (bestResult output)
%   Consumed by:  sidLTVdiscFrozen, sidCompare, sidResidual
%
%   Field               Dimensions          Description
%   .................................................................
%   .A                  (p x p x N)         Time-varying dynamics matrices.
%   .B                  (p x q x N)         Time-varying input matrices.
%   .Lambda             (N-1 x 1)           Regularization values used.
%   .Cost               (1 x 3)             [total, data_fidelity,
%                                            regularization].
%   .DataLength         scalar               Number of time steps N.
%   .StateDim           scalar               State dimension p.
%   .InputDim           scalar               Input dimension q.
%   .NumTrajectories    scalar               Number of trajectories L.
%   .Algorithm          char                 'cosmic'.
%   .Preconditioned     logical              Preconditioning flag.
%   .Method             char                 'sidLTVdisc'.
%
%   When 'Uncertainty' is true, the following fields are added:
%
%   .AStd               (p x p x N)         Std dev of A(k) entries.
%   .BStd               (p x q x N)         Std dev of B(k) entries.
%   .P                  (d x d x N)         Row-wise posterior covariance,
%                                            d = p + q.
%   .NoiseCov           (p x p)             Noise covariance (provided or
%                                            estimated).
%   .NoiseCovEstimated  logical              true if estimated from
%                                            residuals.
%   .NoiseVariance      scalar               trace(NoiseCov) / p.
%   .DegreesOfFreedom   scalar               Effective d.o.f. (NaN if
%                                            NoiseCov was provided).
%
% =========================================================================
%  5. LTVIOResult  (sidLTVdiscIO)
% =========================================================================
%
%   Produced by:  sidLTVdiscIO
%   Consumed by:  sidCompare, sidResidual
%
%   Field               Dimensions          Description
%   .................................................................
%   .A                  (n x n x N)         Estimated dynamics matrices.
%   .B                  (n x q x N)         Estimated input matrices.
%   .X                  (N+1 x n x L) or    Estimated state trajectories.
%                       cell {L x 1}
%   .H                  (py x n)            Observation matrix (copy).
%   .R                  (py x py)           Noise covariance used.
%   .Cost               (n_iter x 1)        Cost J at each iteration.
%   .Iterations         scalar               Number of alternating iters.
%   .Lambda             (N-1 x 1)           Regularisation used.
%   .DataLength         scalar               Number of time steps N.
%   .StateDim           scalar               State dimension n.
%   .OutputDim          scalar               Output dimension py.
%   .InputDim           scalar               Input dimension q.
%   .NumTrajectories    scalar               Number of trajectories L.
%   .Algorithm          char                 'cosmic'.
%   .Method             char                 'sidLTVdiscIO'.
%
% =========================================================================
%  6. FrozenResult  (sidLTVdiscFrozen)
% =========================================================================
%
%   Produced by:  sidLTVdiscFrozen
%   Consumed by:  sidBodePlot, sidMapPlot (via manual extraction)
%
%   Field               Dimensions          Description
%   .................................................................
%   .Frequency          (nf x 1)            Frequency vector, rad/sample.
%   .FrequencyHz        (nf x 1)            Frequency vector, Hz.
%   .TimeSteps          (nk x 1)            Selected time step indices
%                                            (1-based).
%   .Response           (nf x p x q x nk)   Complex frozen transfer
%                                            function G(w, k).
%   .ResponseStd        (nf x p x q x nk)   Std dev of Response.
%                                            [] if no uncertainty.
%   .SampleTime         scalar               Sample time in seconds.
%   .Method             char                 'sidLTVdiscFrozen'.
%
% =========================================================================
%  7. CompareResult  (sidCompare)
% =========================================================================
%
%   Produced by:  sidCompare
%
%   Field               Dimensions          Description
%   .................................................................
%   .Predicted          (N x ny)            Model-predicted output.
%   .Measured           (N x ny)            Measured output (copy).
%   .Fit                (1 x ny)            NRMSE fit percentage per
%                                            channel (100% = perfect).
%   .Residual           (N x ny)            Measured - Predicted.
%   .Method             char                 Method of the source model.
%
% =========================================================================
%  8. ResidualResult  (sidResidual)
% =========================================================================
%
%   Produced by:  sidResidual
%
%   Field               Dimensions          Description
%   .................................................................
%   .Residual           (N x ny)            Residual time series e(t).
%   .AutoCorr           (M+1 x 1)           Normalised autocorrelation
%                                            r_ee(tau).
%   .CrossCorr          (2M+1 x 1)          Normalised cross-correlation
%                                            r_eu(tau). [] for time-series.
%   .ConfidenceBound    scalar               99% bound 2.58 / sqrt(N).
%   .WhitenessPass      logical              true if autocorrelation
%                                            test passes.
%   .IndependencePass   logical              true if cross-correlation
%                                            test passes.
%   .DataLength         scalar               Number of samples N.
%
%   See also: sidFreqBT, sidFreqBTFDR, sidFreqETFE, sidFreqMap,
%             sidSpectrogram, sidLTVdisc, sidLTVdiscIO, sidLTVdiscFrozen,
%             sidCompare, sidResidual
%
%   Changelog:
%   2026-04-09: First version by Pedro Lourenco.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenco, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    fprintf('Run "help sidResultTypes" to view all result struct definitions.\n');
end
