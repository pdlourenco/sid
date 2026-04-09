# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""sid -- System Identification Toolbox (Python implementation).

MATLAB to Python function mapping
----------------------------------
sidFreqBT       -> sid.freq_bt
sidFreqETFE     -> sid.freq_etfe
sidFreqBTFDR    -> sid.freq_btfdr
sidFreqMap      -> sid.freq_map
sidSpectrogram  -> sid.spectrogram
sidLTVdisc      -> sid.ltv_disc
sidLTVdiscIO    -> sid.ltv_disc_io
sidLTVdiscTune  -> sid.ltv_disc_tune
sidLTVdiscFrozen -> sid.ltv_disc_frozen
sidLTIfreqIO    -> sid.lti_freq_io
sidLTVStateEst  -> sid.ltv_state_est
sidModelOrder   -> sid.model_order
sidDetrend      -> sid.detrend
sidResidual     -> sid.residual
sidCompare      -> sid.compare
sidBodePlot     -> sid.bode_plot
sidSpectrumPlot -> sid.spectrum_plot
sidMapPlot      -> sid.map_plot
sidSpectrogramPlot -> sid.spectrogram_plot
"""

from __future__ import annotations

__version__ = "1.0.0"

from sid._exceptions import SidError
from sid._results import FreqMapResult, FreqResult, FrozenResult, LTVResult, SpectrogramResult
from sid.detrend import detrend
from sid.freq_bt import freq_bt
from sid.freq_btfdr import freq_btfdr
from sid.freq_etfe import freq_etfe
from sid.freq_map import freq_map
from sid.ltv_disc import ltv_disc
from sid.ltv_disc_frozen import ltv_disc_frozen
from sid.ltv_disc_tune import ltv_disc_tune
from sid.spectrogram import spectrogram

__all__ = [
    "__version__",
    "FreqMapResult",
    "FreqResult",
    "FrozenResult",
    "LTVResult",
    "SidError",
    "SpectrogramResult",
    "detrend",
    "freq_bt",
    "freq_btfdr",
    "freq_etfe",
    "freq_map",
    "ltv_disc",
    "ltv_disc_frozen",
    "ltv_disc_tune",
    "spectrogram",
]
