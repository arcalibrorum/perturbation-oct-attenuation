"""Core perturbation-based OCT attenuation reconstruction routines.

All public functions operate on linear OCT intensity. Depth is expressed in
millimetres and attenuation in inverse millimetres.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def _bounded_golden_section(objective, lower: float, upper: float, tolerance: float = 1e-11) -> float:
    """Minimize a scalar objective on a closed interval without extra dependencies."""

    ratio = (np.sqrt(5.0) - 1.0) / 2.0
    left, right = float(lower), float(upper)
    x1 = right - ratio * (right - left)
    x2 = left + ratio * (right - left)
    f1, f2 = objective(x1), objective(x2)
    for _ in range(1000):
        if right - left <= tolerance * max(1.0, abs(x1), abs(x2)):
            break
        if f1 <= f2:
            right, x2, f2 = x2, x1, f1
            x1 = right - ratio * (right - left)
            f1 = objective(x1)
        else:
            left, x1, f1 = x1, x2, f2
            x2 = left + ratio * (right - left)
            f2 = objective(x2)
    return float((left + right) / 2.0)


@dataclass(frozen=True)
class ZerothOrderFit:
    """Parameters and diagnostics for the local zeroth-order fit."""

    mu0: float
    amplitude: float
    rmse: float
    r_squared: float
    sample_count: int


@dataclass(frozen=True)
class ReconstructionResult:
    """Zeroth-order, first-order, and combined attenuation estimates."""

    mu0: FloatArray
    mu1: FloatArray
    total: FloatArray
    fits: Tuple[ZerothOrderFit, ...] = ()


def _as_vector(values: ArrayLike, name: str) -> FloatArray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size < 2:
        raise ValueError(f"{name} must be a one-dimensional array with at least two samples")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _validate_depth_and_intensity(z: ArrayLike, intensity: ArrayLike) -> Tuple[FloatArray, FloatArray]:
    depth = _as_vector(z, "z")
    signal = _as_vector(intensity, "intensity")
    if depth.shape != signal.shape:
        raise ValueError("z and intensity must have the same shape")
    if np.any(np.diff(depth) <= 0):
        raise ValueError("z must be strictly increasing")
    if np.any(signal < 0):
        raise ValueError("linear OCT intensity cannot be negative")
    return depth, signal


def _uniform_spacing(z: FloatArray) -> float:
    increments = np.diff(z)
    spacing = float(np.median(increments))
    if not np.allclose(increments, spacing, rtol=1e-7, atol=max(1e-12, abs(spacing) * 1e-10)):
        raise ValueError("the current implementation requires a uniformly sampled depth grid")
    return spacing


def left_cumulative_integral(values: ArrayLike, spacing: float) -> FloatArray:
    """Return the left-rectangle cumulative integral with integral[0] = 0.

    This is the cumulative-sum convention used for the local integral in the
    manuscript. ``spacing`` must use the same depth unit as the attenuation.
    """

    vector = _as_vector(values, "values")
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be positive and finite")
    result = np.zeros_like(vector)
    result[1:] = spacing * np.cumsum(vector[:-1])
    return result


def fit_zeroth_order(
    z: ArrayLike,
    intensity: ArrayLike,
    mu_bounds: Tuple[float, float] = (1e-5, 50.0),
    fit_window: Optional[Tuple[float, float]] = None,
) -> ZerothOrderFit:
    r"""Fit :math:`I=A\mu_0\exp(-2\mu_0\xi)` by least squares.

    The amplitude is eliminated analytically for each candidate ``mu0`` and
    the resulting one-dimensional bounded objective is minimized. This is
    mathematically equivalent to the two-parameter least-squares problem in
    the manuscript while avoiding a toolbox-specific nonlinear fitter.
    """

    depth, signal = _validate_depth_and_intensity(z, intensity)
    local_depth = depth - depth[0]
    lower, upper = map(float, mu_bounds)
    if not (0 < lower < upper and np.isfinite(upper)):
        raise ValueError("mu_bounds must satisfy 0 < lower < upper")

    mask = np.ones(depth.size, dtype=bool)
    if fit_window is not None:
        start, stop = map(float, fit_window)
        if not start < stop:
            raise ValueError("fit_window must satisfy start < stop")
        mask = (local_depth >= start) & (local_depth <= stop)
    if np.count_nonzero(mask) < 3:
        raise ValueError("the fitting window must contain at least three samples")

    x_fit = local_depth[mask]
    y_fit = signal[mask]

    def amplitude_for(mu0: float) -> float:
        basis = mu0 * np.exp(-2.0 * mu0 * x_fit)
        denominator = float(np.dot(basis, basis))
        if denominator <= np.finfo(float).tiny:
            return 0.0
        return max(0.0, float(np.dot(basis, y_fit) / denominator))

    def objective(mu0: float) -> float:
        amplitude = amplitude_for(mu0)
        residual = y_fit - amplitude * mu0 * np.exp(-2.0 * mu0 * x_fit)
        return float(np.dot(residual, residual))

    mu0 = _bounded_golden_section(objective, lower, upper)
    amplitude = amplitude_for(mu0)
    prediction = amplitude * mu0 * np.exp(-2.0 * mu0 * x_fit)
    residual = y_fit - prediction
    sse = float(np.dot(residual, residual))
    rmse = float(np.sqrt(sse / y_fit.size))
    centered = y_fit - float(np.mean(y_fit))
    sst = float(np.dot(centered, centered))
    r_squared = float(1.0 - sse / sst) if sst > 0 else float("nan")
    return ZerothOrderFit(mu0, amplitude, rmse, r_squared, int(y_fit.size))


def first_order_correction(
    z: ArrayLike,
    intensity: ArrayLike,
    mu0: float,
    amplitude: float,
) -> FloatArray:
    r"""Evaluate the local first-order correction in the manuscript.

    The local coordinate and cumulative integral are reset at the first depth
    sample. ``amplitude`` is :math:`A=2 C L_0 R`.
    """

    depth, signal = _validate_depth_and_intensity(z, intensity)
    spacing = _uniform_spacing(depth)
    mu0 = float(mu0)
    amplitude = float(amplitude)
    if not np.isfinite(mu0) or mu0 <= 0:
        raise ValueError("mu0 must be positive and finite")
    if not np.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("amplitude must be positive and finite")

    xi = depth - depth[0]
    cumulative_intensity = left_cumulative_integral(signal, spacing)
    bracket = signal / amplitude + (2.0 * mu0 / amplitude) * cumulative_intensity - mu0
    return bracket * np.exp(2.0 * mu0 * xi)


def reconstruct_interval(
    z: ArrayLike,
    intensity: ArrayLike,
    fit: Optional[ZerothOrderFit] = None,
    mu_bounds: Tuple[float, float] = (1e-5, 50.0),
    fit_window: Optional[Tuple[float, float]] = None,
) -> ReconstructionResult:
    """Reconstruct one local depth interval using fitted or fixed parameters."""

    depth, signal = _validate_depth_and_intensity(z, intensity)
    local_fit = fit or fit_zeroth_order(depth, signal, mu_bounds, fit_window)
    correction = first_order_correction(depth, signal, local_fit.mu0, local_fit.amplitude)
    baseline = np.full_like(depth, local_fit.mu0)
    return ReconstructionResult(baseline, correction, baseline + correction, (local_fit,))


FixedParameter = Union[ZerothOrderFit, Tuple[float, float]]


def _coerce_fit(value: FixedParameter, sample_count: int) -> ZerothOrderFit:
    if isinstance(value, ZerothOrderFit):
        return value
    mu0, amplitude = value
    return ZerothOrderFit(float(mu0), float(amplitude), float("nan"), float("nan"), sample_count)


def reconstruct_piecewise(
    z: ArrayLike,
    intensity: ArrayLike,
    interval_edges: ArrayLike,
    fixed_parameters: Optional[Sequence[FixedParameter]] = None,
    fit_windows: Optional[Sequence[Optional[Tuple[float, float]]]] = None,
    mu_bounds: Tuple[float, float] = (1e-5, 50.0),
) -> ReconstructionResult:
    """Apply the local reconstruction independently in consecutive intervals.

    ``interval_edges`` are physical depths. A sample exactly on an internal
    edge belongs to the deeper interval. ``fit_windows`` use local coordinates
    relative to the corresponding interval's first sampled depth.
    """

    depth, signal = _validate_depth_and_intensity(z, intensity)
    edges = _as_vector(interval_edges, "interval_edges")
    if np.any(np.diff(edges) <= 0):
        raise ValueError("interval_edges must be strictly increasing")
    interval_count = edges.size - 1
    if fixed_parameters is not None and len(fixed_parameters) != interval_count:
        raise ValueError("fixed_parameters must contain one pair per interval")
    if fit_windows is not None and len(fit_windows) != interval_count:
        raise ValueError("fit_windows must contain one entry per interval")

    baseline = np.full_like(depth, np.nan)
    correction = np.full_like(depth, np.nan)
    fits = []
    for index in range(interval_count):
        start = int(np.searchsorted(depth, edges[index], side="left"))
        side = "right" if index == interval_count - 1 else "left"
        stop = int(np.searchsorted(depth, edges[index + 1], side=side))
        if stop - start < 3:
            raise ValueError(f"interval {index} contains fewer than three samples")
        local_z = depth[start:stop]
        local_i = signal[start:stop]
        if fixed_parameters is None:
            window = None if fit_windows is None else fit_windows[index]
            fit = fit_zeroth_order(local_z, local_i, mu_bounds, window)
        else:
            fit = _coerce_fit(fixed_parameters[index], local_z.size)
        local = reconstruct_interval(local_z, local_i, fit)
        baseline[start:stop] = local.mu0
        correction[start:stop] = local.mu1
        fits.append(fit)

    return ReconstructionResult(baseline, correction, baseline + correction, tuple(fits))


def reconstruct_bscan(
    z: ArrayLike,
    intensity: ArrayLike,
    layer_boundaries: ArrayLike,
    fixed_parameters: Sequence[FixedParameter],
) -> ReconstructionResult:
    """Reconstruct a B-scan with fixed layer-wise zeroth-order parameters.

    ``intensity`` has shape ``(depth, lateral)``. ``layer_boundaries`` may be a
    vector shared by all columns or an array of shape ``(layers + 1, lateral)``.
    The vessel mask, if any, is not needed by this routine.
    """

    depth = _as_vector(z, "z")
    image = np.asarray(intensity, dtype=float)
    if image.ndim != 2 or image.shape[0] != depth.size:
        raise ValueError("intensity must have shape (depth, lateral)")
    if not np.all(np.isfinite(image)) or np.any(image < 0):
        raise ValueError("intensity must be finite and nonnegative")
    if np.any(np.diff(depth) <= 0):
        raise ValueError("z must be strictly increasing")
    _uniform_spacing(depth)

    boundaries = np.asarray(layer_boundaries, dtype=float)
    if boundaries.ndim == 1:
        boundaries = np.repeat(boundaries[:, None], image.shape[1], axis=1)
    if boundaries.shape != (len(fixed_parameters) + 1, image.shape[1]):
        raise ValueError("layer_boundaries must have shape (layers + 1, lateral)")

    baseline = np.full_like(image, np.nan)
    correction = np.full_like(image, np.nan)
    for column in range(image.shape[1]):
        result = reconstruct_piecewise(
            depth,
            image[:, column],
            boundaries[:, column],
            fixed_parameters=fixed_parameters,
        )
        baseline[:, column] = result.mu0
        correction[:, column] = result.mu1
    return ReconstructionResult(baseline, correction, baseline + correction, tuple())


def fit_control_baselines(
    z: ArrayLike,
    intensity: ArrayLike,
    layer_boundaries: ArrayLike,
    control_columns: ArrayLike,
    fit_windows: Optional[Sequence[Optional[Tuple[float, float]]]] = None,
    mu_bounds: Tuple[float, float] = (1e-5, 50.0),
) -> Tuple[ZerothOrderFit, ...]:
    """Estimate fixed layer parameters from remote control columns.

    Each selected control A-line is fitted independently inside its layer
    boundaries. The returned ``mu0`` and ``amplitude`` are the medians of the
    valid control-column fits, matching the no-vessel-prior workflow. The
    ``sample_count`` field records the number of accepted column fits.
    """

    depth = _as_vector(z, "z")
    image = np.asarray(intensity, dtype=float)
    if image.ndim != 2 or image.shape[0] != depth.size:
        raise ValueError("intensity must have shape (depth, lateral)")
    if not np.all(np.isfinite(image)) or np.any(image < 0):
        raise ValueError("intensity must be finite and nonnegative")
    boundaries = np.asarray(layer_boundaries, dtype=float)
    if boundaries.ndim == 1:
        layer_count = boundaries.size - 1
        boundaries = np.repeat(boundaries[:, None], image.shape[1], axis=1)
    elif boundaries.ndim == 2:
        layer_count = boundaries.shape[0] - 1
    else:
        raise ValueError("layer_boundaries must be one- or two-dimensional")
    if boundaries.shape != (layer_count + 1, image.shape[1]):
        raise ValueError("layer_boundaries must have shape (layers + 1, lateral)")
    if fit_windows is not None and len(fit_windows) != layer_count:
        raise ValueError("fit_windows must contain one entry per layer")

    selected = np.asarray(control_columns)
    if selected.dtype == bool:
        if selected.ndim != 1 or selected.size != image.shape[1]:
            raise ValueError("a Boolean control_columns mask must match the lateral dimension")
        columns = np.flatnonzero(selected)
    else:
        columns = np.asarray(selected, dtype=int).ravel()
    if columns.size == 0 or np.any(columns < 0) or np.any(columns >= image.shape[1]):
        raise ValueError("control_columns must identify at least one valid column")

    aggregate = []
    for layer in range(layer_count):
        fits = []
        for column in columns:
            start = int(np.searchsorted(depth, boundaries[layer, column], side="left"))
            side = "right" if layer == layer_count - 1 else "left"
            stop = int(np.searchsorted(depth, boundaries[layer + 1, column], side=side))
            if stop - start < 3:
                continue
            window = None if fit_windows is None else fit_windows[layer]
            try:
                fit = fit_zeroth_order(depth[start:stop], image[start:stop, column], mu_bounds, window)
            except (ValueError, RuntimeError, FloatingPointError):
                continue
            if fit.mu0 > 0 and fit.amplitude > 0 and np.isfinite(fit.mu0 + fit.amplitude):
                fits.append(fit)
        if not fits:
            raise RuntimeError(f"no valid control-column fit was obtained for layer {layer}")
        aggregate.append(
            ZerothOrderFit(
                float(np.median([fit.mu0 for fit in fits])),
                float(np.median([fit.amplitude for fit in fits])),
                float(np.median([fit.rmse for fit in fits])),
                float(np.median([fit.r_squared for fit in fits])),
                len(fits),
            )
        )
    return tuple(aggregate)


def conventional_depth_resolved(intensity: ArrayLike, spacing: float, axis: int = 0) -> FloatArray:
    """Evaluate the conventional distal-tail estimator from the manuscript."""

    signal = np.asarray(intensity, dtype=float)
    if signal.ndim not in (1, 2):
        raise ValueError("intensity must be one- or two-dimensional")
    if not np.all(np.isfinite(signal)) or np.any(signal < 0):
        raise ValueError("intensity must be finite and nonnegative")
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be positive and finite")
    if axis < 0:
        axis += signal.ndim
    if axis < 0 or axis >= signal.ndim:
        raise ValueError("axis is out of range")

    reversed_signal = np.flip(signal, axis=axis)
    tail = np.flip(np.cumsum(reversed_signal, axis=axis), axis=axis)
    denominator = 2.0 * spacing * tail
    result = np.full_like(signal, np.nan, dtype=float)
    np.divide(signal, denominator, out=result, where=denominator > 0)
    return result
