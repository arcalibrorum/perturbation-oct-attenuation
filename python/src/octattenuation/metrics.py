"""Evaluation metrics used by the no-vessel-prior component ablation."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from numpy.typing import ArrayLike


def vessel_contrast(image: ArrayLike, vessel_mask: ArrayLike, background_mask: ArrayLike) -> float:
    """Median vessel attenuation minus median surrounding-retina attenuation."""

    values = np.asarray(image, dtype=float)
    vessel = np.asarray(vessel_mask, dtype=bool)
    background = np.asarray(background_mask, dtype=bool)
    if values.shape != vessel.shape or values.shape != background.shape:
        raise ValueError("image and masks must have identical shapes")
    if not np.any(vessel) or not np.any(background):
        raise ValueError("both masks must contain at least one sample")
    return float(np.nanmedian(values[vessel]) - np.nanmedian(values[background]))


def vessel_roi_mae(image: ArrayLike, reference: ArrayLike, vessel_mask: ArrayLike) -> float:
    """Mean absolute error within a post hoc vessel ROI."""

    values = np.asarray(image, dtype=float)
    target = np.asarray(reference, dtype=float)
    vessel = np.asarray(vessel_mask, dtype=bool)
    if values.shape != target.shape or values.shape != vessel.shape:
        raise ValueError("image, reference, and vessel_mask must have identical shapes")
    if not np.any(vessel):
        raise ValueError("vessel_mask must contain at least one sample")
    return float(np.nanmean(np.abs(values[vessel] - target[vessel])))


def ablation_metrics(
    image: ArrayLike,
    reference: ArrayLike,
    vessel_mask: ArrayLike,
    background_mask: ArrayLike,
) -> Tuple[float, float]:
    """Return vessel contrast and vessel-ROI MAE as defined in the manuscript."""

    return (
        vessel_contrast(image, vessel_mask, background_mask),
        vessel_roi_mae(image, reference, vessel_mask),
    )
