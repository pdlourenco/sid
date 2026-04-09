# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for plotting functions.

All tests use matplotlib's Agg backend for headless execution.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # Must be before pyplot import
import matplotlib.pyplot as plt

import numpy as np
import pytest
from scipy.signal import lfilter

from sid.freq_bt import freq_bt
from sid.bode_plot import bode_plot
from sid.spectrum_plot import spectrum_plot
from sid.freq_map import freq_map
from sid.spectrogram import spectrogram as sid_spectrogram  # avoid name clash
from sid.map_plot import map_plot
from sid.spectrogram_plot import spectrogram_plot


# ---------------------------------------------------------------------------
# Bode plot tests (4)
# ---------------------------------------------------------------------------


class TestBodePlot:
    """Tests for bode_plot."""

    @pytest.fixture()
    def siso_result(self):
        """SISO frequency response result."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return freq_bt(y, u)

    @pytest.fixture()
    def timeseries_result(self):
        """Time-series result (no input)."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(300)
        return freq_bt(y, None)

    def test_bode_plot_siso(self, siso_result) -> None:
        """SISO result returns dict with 'fig', 'ax_mag', 'ax_phase'."""
        h = bode_plot(siso_result)
        assert "fig" in h
        assert "ax_mag" in h
        assert "ax_phase" in h
        assert isinstance(h["fig"], plt.Figure)
        assert isinstance(h["ax_mag"], plt.Axes)
        assert isinstance(h["ax_phase"], plt.Axes)
        plt.close(h["fig"])

    def test_bode_plot_confidence(self, siso_result) -> None:
        """confidence=0 and confidence=3 both run without error."""
        h0 = bode_plot(siso_result, confidence=0)
        assert "fig" in h0
        plt.close(h0["fig"])

        h3 = bode_plot(siso_result, confidence=3)
        assert "fig" in h3
        plt.close(h3["fig"])

    def test_bode_plot_hz(self, siso_result) -> None:
        """frequency_unit='Hz' runs without error."""
        h = bode_plot(siso_result, frequency_unit="Hz")
        assert "fig" in h
        plt.close(h["fig"])

    def test_bode_plot_timeseries_error(self, timeseries_result) -> None:
        """Time-series result (response=None) raises an error."""
        with pytest.raises(Exception):
            bode_plot(timeseries_result)


# ---------------------------------------------------------------------------
# Spectrum plot tests (3)
# ---------------------------------------------------------------------------


class TestSpectrumPlot:
    """Tests for spectrum_plot."""

    @pytest.fixture()
    def siso_result(self):
        """SISO frequency response result."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return freq_bt(y, u)

    @pytest.fixture()
    def timeseries_result(self):
        """Time-series result (no input)."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(300)
        return freq_bt(y, None)

    def test_spectrum_plot_siso(self, siso_result) -> None:
        """SISO result returns dict with 'fig', 'ax'."""
        h = spectrum_plot(siso_result)
        assert "fig" in h
        assert "ax" in h
        assert isinstance(h["fig"], plt.Figure)
        assert isinstance(h["ax"], plt.Axes)
        plt.close(h["fig"])

    def test_spectrum_plot_timeseries(self, timeseries_result) -> None:
        """Time-series result (u=None) works for output spectrum."""
        h = spectrum_plot(timeseries_result)
        assert "fig" in h
        assert "ax" in h
        plt.close(h["fig"])

    def test_spectrum_plot_confidence(self, siso_result) -> None:
        """confidence=0 and confidence=3 both work."""
        h0 = spectrum_plot(siso_result, confidence=0)
        assert "fig" in h0
        plt.close(h0["fig"])

        h3 = spectrum_plot(siso_result, confidence=3)
        assert "fig" in h3
        plt.close(h3["fig"])


# ---------------------------------------------------------------------------
# Map plot tests (3)
# ---------------------------------------------------------------------------


class TestMapPlot:
    """Tests for map_plot."""

    @pytest.fixture()
    def freq_map_result(self):
        """SISO freq_map result for plotting."""
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return freq_map(y, u, segment_length=512)

    @pytest.fixture()
    def freq_map_timeseries(self):
        """Time-series freq_map result (no input)."""
        rng = np.random.default_rng(42)
        N = 2000
        y = rng.standard_normal(N)
        return freq_map(y, None, segment_length=512)

    def test_map_plot_magnitude(self, freq_map_result) -> None:
        """Default magnitude plot returns dict with 'fig', 'ax'."""
        h = map_plot(freq_map_result)
        assert "fig" in h
        assert "ax" in h
        assert isinstance(h["fig"], plt.Figure)
        assert isinstance(h["ax"], plt.Axes)
        plt.close(h["fig"])

    def test_map_plot_coherence(self, freq_map_result) -> None:
        """plot_type='coherence' works for SISO result."""
        h = map_plot(freq_map_result, plot_type="coherence")
        assert "fig" in h
        assert "ax" in h
        plt.close(h["fig"])

    def test_map_plot_timeseries_noise(self, freq_map_timeseries) -> None:
        """Time-series freq_map with plot_type='noise' works."""
        h = map_plot(freq_map_timeseries, plot_type="noise")
        assert "fig" in h
        assert "ax" in h
        plt.close(h["fig"])


# ---------------------------------------------------------------------------
# Spectrogram plot tests (2)
# ---------------------------------------------------------------------------


class TestSpectrogramPlot:
    """Tests for spectrogram_plot."""

    @pytest.fixture()
    def spectrogram_result(self):
        """Spectrogram result for plotting."""
        rng = np.random.default_rng(42)
        N = 2000
        x = rng.standard_normal(N)
        return sid_spectrogram(x, window_length=256)

    def test_spectrogram_plot(self, spectrogram_result) -> None:
        """Spectrogram plot returns dict with 'fig', 'ax'."""
        h = spectrogram_plot(spectrogram_result)
        assert "fig" in h
        assert "ax" in h
        assert isinstance(h["fig"], plt.Figure)
        assert isinstance(h["ax"], plt.Axes)
        plt.close(h["fig"])

    def test_spectrogram_plot_log_scale(self, spectrogram_result) -> None:
        """frequency_scale='log' works."""
        h = spectrogram_plot(spectrogram_result, frequency_scale="log")
        assert "fig" in h
        assert "ax" in h
        plt.close(h["fig"])
