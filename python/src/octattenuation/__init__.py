"""Perturbation-based reconstruction of depth-resolved OCT attenuation."""

from .core import (
    ReconstructionResult,
    ZerothOrderFit,
    conventional_depth_resolved,
    fit_control_baselines,
    first_order_correction,
    fit_zeroth_order,
    left_cumulative_integral,
    reconstruct_bscan,
    reconstruct_interval,
    reconstruct_piecewise,
)
from .metrics import ablation_metrics, vessel_contrast, vessel_roi_mae
from .simulation import (
    finite_depth_profiles,
    simulate_beer_lambert_intensity,
    simulate_first_order_intensity,
)

__all__ = [
    "ReconstructionResult",
    "ZerothOrderFit",
    "ablation_metrics",
    "conventional_depth_resolved",
    "finite_depth_profiles",
    "fit_control_baselines",
    "first_order_correction",
    "fit_zeroth_order",
    "left_cumulative_integral",
    "reconstruct_bscan",
    "reconstruct_interval",
    "reconstruct_piecewise",
    "simulate_beer_lambert_intensity",
    "simulate_first_order_intensity",
    "vessel_contrast",
    "vessel_roi_mae",
]

__version__ = "0.1.0"
