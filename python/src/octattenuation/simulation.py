"""Synthetic linear-intensity generators for examples and verification."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .core import _as_vector, _uniform_spacing, left_cumulative_integral


FloatArray = NDArray[np.float64]


def simulate_beer_lambert_intensity(
    z: ArrayLike,
    attenuation: ArrayLike,
    amplitude: float = 1.0,
) -> FloatArray:
    r"""Generate :math:`I=A\mu(z)\exp[-2\int_0^z\mu(u)du]`.

    This is the ideal single-scattering model used for the numerical A-line
    demonstration. The returned signal is linear intensity.
    """

    depth = _as_vector(z, "z")
    mu = _as_vector(attenuation, "attenuation")
    if depth.shape != mu.shape:
        raise ValueError("z and attenuation must have the same shape")
    if np.any(np.diff(depth) <= 0):
        raise ValueError("z must be strictly increasing")
    if np.any(mu < 0):
        raise ValueError("attenuation cannot be negative")
    if not np.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("amplitude must be positive and finite")
    spacing = _uniform_spacing(depth)
    optical_depth = left_cumulative_integral(mu, spacing)
    return amplitude * mu * np.exp(-2.0 * optical_depth)


def simulate_first_order_intensity(
    z: ArrayLike,
    mu0: float,
    mu1: ArrayLike,
    amplitude: float = 1.0,
) -> FloatArray:
    """Generate intensity from the manuscript's local first-order model."""

    depth = _as_vector(z, "z")
    correction = _as_vector(mu1, "mu1")
    if depth.shape != correction.shape:
        raise ValueError("z and mu1 must have the same shape")
    if np.any(np.diff(depth) <= 0):
        raise ValueError("z must be strictly increasing")
    if not np.isfinite(mu0) or mu0 <= 0:
        raise ValueError("mu0 must be positive and finite")
    if not np.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("amplitude must be positive and finite")
    spacing = _uniform_spacing(depth)
    xi = depth - depth[0]
    integral = left_cumulative_integral(correction, spacing)
    signal = amplitude * (mu0 + correction - 2.0 * mu0 * integral) * np.exp(-2.0 * mu0 * xi)
    if np.any(signal < 0):
        raise ValueError("the selected correction makes the first-order intensity negative")
    return signal


def finite_depth_profiles(sample_count: int = 1000, spacing: float = 0.002) -> Tuple[FloatArray, FloatArray, FloatArray]:
    """Return deterministic step-like and smooth profiles for the public demo."""

    if sample_count < 10 or spacing <= 0:
        raise ValueError("sample_count and spacing must be positive")
    z = np.arange(sample_count, dtype=float) * spacing
    if sample_count == 1000 and np.isclose(spacing, 0.002):
        step = np.concatenate(
            (
                np.full(150, 0.6),
                np.full(10, 0.2),
                np.full(140, 0.6),
                np.full(300, 0.4),
                np.full(100, 0.8),
                np.full(10, 1.0),
                np.full(290, 0.8),
            )
        )
        smooth = np.concatenate((np.full(300, 1.0), np.full(300, 1.6), np.full(400, 1.2)))
        smooth = smooth + 0.05 * np.sin(0.2 * np.arange(1, 1001))
    else:
        step = np.where(z < 0.3 * z[-1], 0.6, np.where(z < 0.7 * z[-1], 0.4, 0.8))
        smooth = 0.9 + 0.35 * np.sin(2.0 * np.pi * z / max(z[-1], spacing))
    return z, step, smooth
