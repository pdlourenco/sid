# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for plotting functions.

Beyond smoke testing, these tests validate that the correct data is plotted
on the axes (line data, axis scales, labels, confidence bands, colormap data).
All tests use matplotlib's Agg backend for headless execution.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pytest
from scipy.signal import lfilter

from sid.bode_plot import bode_plot
from sid.freq_bt import freq_bt
from sid.freq_map import freq_map
from sid.map_plot import map_plot
from sid.spectrogram import spectrogram as sid_spectrogram
from sid.spectrogram_plot import spectrogram_plot
from sid.spectrum_plot import spectrum_plot


# ======================================================================
#  Bode plot
# ======================================================================


class TestBodePlot:
    """Tests for bode_plot — validates plotted data, not just handles."""

    @pytest.fixture()
    def siso_result(self):
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return freq_bt(y, u)

    @pytest.fixture()
    def ts_result(self):
        rng = np.random.default_rng(42)
        return freq_bt(rng.standard_normal(300), None)

    def test_returns_handles(self, siso_result) -> None:
        h = bode_plot(siso_result)
        assert isinstance(h["fig"], plt.Figure)
        assert isinstance(h["ax_mag"], plt.Axes)
        assert isinstance(h["ax_phase"], plt.Axes)
        plt.close(h["fig"])

    def test_magnitude_data_matches_response(self, siso_result) -> None:
        """The plotted y-data on the magnitude axis matches 20*log10(|G|)."""
        h = bode_plot(siso_result, confidence=0)
        line = h["ax_mag"].get_lines()[0]
        y_plotted = line.get_ydata()
        G = siso_result.response
        expected_mag_db = 20 * np.log10(np.abs(G))
        np.testing.assert_allclose(y_plotted, expected_mag_db, rtol=1e-10)
        plt.close(h["fig"])

    def test_phase_data_matches_response(self, siso_result) -> None:
        """The plotted y-data on the phase axis matches angle(G) in degrees."""
        h = bode_plot(siso_result, confidence=0)
        line = h["ax_phase"].get_lines()[0]
        y_plotted = line.get_ydata()
        expected_phase = np.angle(siso_result.response) * 180 / np.pi
        np.testing.assert_allclose(y_plotted, expected_phase, rtol=1e-10)
        plt.close(h["fig"])

    def test_frequency_axis_rad_s(self, siso_result) -> None:
        """Default frequency unit is rad/s with correct values."""
        h = bode_plot(siso_result, confidence=0)
        x_plotted = h["ax_mag"].get_lines()[0].get_xdata()
        expected_freq = siso_result.frequency / siso_result.sample_time
        np.testing.assert_allclose(x_plotted, expected_freq, rtol=1e-10)
        plt.close(h["fig"])

    def test_frequency_axis_hz(self, siso_result) -> None:
        """frequency_unit='Hz' uses Hz values."""
        h = bode_plot(siso_result, confidence=0, frequency_unit="Hz")
        x_plotted = h["ax_mag"].get_lines()[0].get_xdata()
        np.testing.assert_allclose(x_plotted, siso_result.frequency_hz, rtol=1e-10)
        plt.close(h["fig"])

    def test_log_x_scale(self, siso_result) -> None:
        """Both axes use log scale on x-axis."""
        h = bode_plot(siso_result)
        assert h["ax_mag"].get_xscale() == "log"
        assert h["ax_phase"].get_xscale() == "log"
        plt.close(h["fig"])

    def test_confidence_bands_present(self, siso_result) -> None:
        """With confidence > 0 and valid std, shaded bands are added."""
        h = bode_plot(siso_result, confidence=3)
        # fill_between creates PolyCollection objects in ax.collections
        collections_mag = h["ax_mag"].collections
        assert len(collections_mag) > 0, "Magnitude axis should have confidence band"
        plt.close(h["fig"])

    def test_no_confidence_bands_when_zero(self, siso_result) -> None:
        """With confidence=0, no shaded bands are added."""
        h = bode_plot(siso_result, confidence=0)
        collections_mag = h["ax_mag"].collections
        assert len(collections_mag) == 0, "No confidence band when confidence=0"
        plt.close(h["fig"])

    def test_timeseries_raises(self, ts_result) -> None:
        with pytest.raises(Exception):
            bode_plot(ts_result)

    def test_custom_axes(self, siso_result) -> None:
        """Plotting into provided axes works."""
        fig, (ax1, ax2) = plt.subplots(2, 1)
        h = bode_plot(siso_result, ax=(ax1, ax2), confidence=0)
        assert h["ax_mag"] is ax1
        assert h["ax_phase"] is ax2
        plt.close(fig)


# ======================================================================
#  Spectrum plot
# ======================================================================


class TestSpectrumPlot:
    """Tests for spectrum_plot — validates plotted spectrum data."""

    @pytest.fixture()
    def siso_result(self):
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return freq_bt(y, u)

    @pytest.fixture()
    def ts_result(self):
        rng = np.random.default_rng(42)
        return freq_bt(rng.standard_normal(300), None)

    def test_spectrum_data_matches(self, siso_result) -> None:
        """Plotted y-data matches 10*log10(noise_spectrum)."""
        h = spectrum_plot(siso_result, confidence=0)
        line = h["ax"].get_lines()[0]
        y_plotted = line.get_ydata()
        ns = siso_result.noise_spectrum
        expected_db = 10 * np.log10(np.maximum(ns, np.finfo(float).tiny))
        np.testing.assert_allclose(y_plotted, expected_db, rtol=1e-6)
        plt.close(h["fig"])

    def test_log_x_scale(self, siso_result) -> None:
        h = spectrum_plot(siso_result)
        assert h["ax"].get_xscale() == "log"
        plt.close(h["fig"])

    def test_timeseries_works(self, ts_result) -> None:
        """Time-series result plots output spectrum without error."""
        h = spectrum_plot(ts_result)
        assert len(h["ax"].get_lines()) > 0
        plt.close(h["fig"])

    def test_confidence_bands(self, siso_result) -> None:
        h = spectrum_plot(siso_result, confidence=3)
        assert len(h["ax"].collections) > 0, "Should have confidence band"
        plt.close(h["fig"])


# ======================================================================
#  Map plot
# ======================================================================


class TestMapPlot:
    """Tests for map_plot — validates colormap data and axes."""

    @pytest.fixture()
    def fmap_result(self):
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return freq_map(y, u, segment_length=512)

    @pytest.fixture()
    def fmap_ts(self):
        rng = np.random.default_rng(42)
        return freq_map(rng.standard_normal(2000), None, segment_length=512)

    def test_magnitude_plot(self, fmap_result) -> None:
        h = map_plot(fmap_result)
        assert h["ax"].get_yscale() == "log", "Frequency axis should be log"
        assert len(h["ax"].collections) > 0, "Should have pcolormesh"
        plt.close(h["fig"])

    def test_coherence_plot(self, fmap_result) -> None:
        h = map_plot(fmap_result, plot_type="coherence")
        assert len(h["ax"].collections) > 0
        plt.close(h["fig"])

    def test_noise_timeseries(self, fmap_ts) -> None:
        h = map_plot(fmap_ts, plot_type="noise")
        assert len(h["ax"].collections) > 0
        plt.close(h["fig"])

    def test_invalid_plot_type(self, fmap_result) -> None:
        with pytest.raises(Exception):
            map_plot(fmap_result, plot_type="invalid_type")


# ======================================================================
#  Spectrogram plot
# ======================================================================


class TestSpectrogramPlot:
    """Tests for spectrogram_plot — validates axes and colormap."""

    @pytest.fixture()
    def spec_result(self):
        rng = np.random.default_rng(42)
        return sid_spectrogram(rng.standard_normal(2000), window_length=256)

    def test_returns_handles(self, spec_result) -> None:
        h = spectrogram_plot(spec_result)
        assert isinstance(h["fig"], plt.Figure)
        assert isinstance(h["ax"], plt.Axes)
        assert len(h["ax"].collections) > 0, "Should have pcolormesh"
        plt.close(h["fig"])

    def test_linear_frequency_scale(self, spec_result) -> None:
        h = spectrogram_plot(spec_result, frequency_scale="linear")
        assert h["ax"].get_yscale() == "linear"
        plt.close(h["fig"])

    def test_log_frequency_scale(self, spec_result) -> None:
        h = spectrogram_plot(spec_result, frequency_scale="log")
        assert h["ax"].get_yscale() == "log"
        plt.close(h["fig"])

    def test_xlabel_is_time(self, spec_result) -> None:
        h = spectrogram_plot(spec_result)
        xlabel = h["ax"].get_xlabel().lower()
        assert "time" in xlabel, f"x-axis should mention time, got '{xlabel}'"
        plt.close(h["fig"])

    def test_ylabel_is_frequency(self, spec_result) -> None:
        h = spectrogram_plot(spec_result)
        ylabel = h["ax"].get_ylabel().lower()
        assert "freq" in ylabel, f"y-axis should mention frequency, got '{ylabel}'"
        plt.close(h["fig"])
