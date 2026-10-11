"""Synthetic fixed-baseline B-scan and no-prior component demonstration."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import octattenuation as oa


def build_case() -> tuple:
    rng = np.random.default_rng(7)
    z = np.arange(360, dtype=float) * 0.002
    x = np.linspace(-1.0, 1.0, 121)
    upper = 0.08 + 0.018 * (1.0 - x**2)
    middle = 0.34 + 0.025 * (1.0 - x**2)
    lower = np.full_like(x, z[-1])
    boundaries = np.vstack((upper, middle, lower))
    parameters = ((0.35, 2.0), (0.75, 1.2))
    intensity = np.full((z.size, x.size), 1e-5)
    reference = np.full_like(intensity, np.nan)
    vessel_mask = np.zeros_like(intensity, dtype=bool)

    for column, lateral in enumerate(x):
        for layer in range(2):
            start = np.searchsorted(z, boundaries[layer, column], side="left")
            stop = np.searchsorted(z, boundaries[layer + 1, column], side="right" if layer == 1 else "left")
            local_z = z[start:stop]
            mu0, amplitude = parameters[layer]
            mu1 = np.zeros(local_z.size)
            if layer == 0:
                vessel = np.exp(-((lateral / 0.18) ** 2)) * np.exp(-((local_z - 0.20) / 0.055) ** 2)
                mu1 = 0.20 * vessel
                vessel_mask[start:stop, column] = vessel > 0.55
            intensity[start:stop, column] = oa.simulate_first_order_intensity(local_z, mu0, mu1, amplitude)
            reference[start:stop, column] = mu0 + mu1

    intensity *= np.maximum(0.0, 1.0 + 0.01 * rng.standard_normal(intensity.shape))
    control_columns = np.flatnonzero(np.abs(x) > 0.65)
    fitted_parameters = oa.fit_control_baselines(z, intensity, boundaries, control_columns, mu_bounds=(0.05, 2.0))
    reconstruction = oa.reconstruct_bscan(z, intensity, boundaries, fitted_parameters)
    surrounding = np.zeros_like(vessel_mask)
    surrounding[:, (np.abs(x) > 0.30) & (np.abs(x) < 0.55)] = vessel_mask.any(axis=1)[:, None]
    return z, x, intensity, reference, vessel_mask, surrounding, reconstruction, fitted_parameters


def run(output: Path, make_plot: bool = True) -> None:
    output.mkdir(parents=True, exist_ok=True)
    z, x, intensity, reference, vessel, surrounding, result, fitted_parameters = build_case()
    contrast, mae = oa.ablation_metrics(result.total, reference, vessel, surrounding)
    np.savez_compressed(
        output / "fixed_baseline_case.npz",
        z_mm=z,
        x=x,
        intensity=intensity,
        reference=reference,
        mu0=result.mu0,
        mu1=result.mu1,
        total=result.total,
        vessel_mask=vessel,
        surrounding_mask=surrounding,
        vessel_contrast=contrast,
        vessel_roi_mae=mae,
        fitted_mu0=np.array([fit.mu0 for fit in fitted_parameters]),
        fitted_amplitude=np.array([fit.amplitude for fit in fitted_parameters]),
    )
    print(f"vessel contrast: {contrast:.6f} mm^-1")
    print(f"vessel ROI MAE: {mae:.6f} mm^-1")

    if not make_plot:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Install the optional plotting dependency with: pip install -e '.[plot]'") from exc

    figure, axes = plt.subplots(1, 4, figsize=(12, 3), constrained_layout=True)
    extent = (x[0], x[-1], z[-1], z[0])
    panels = (
        (intensity, "Linear intensity"),
        (result.mu0, r"$\hat{\mu}_0$"),
        (np.abs(result.mu1), r"$|\hat{\mu}_1|$"),
        (result.total, r"$\hat{\mu}_0+\hat{\mu}_1$"),
    )
    for axis, (image, title) in zip(axes, panels):
        handle = axis.imshow(image, aspect="auto", extent=extent, cmap="viridis")
        axis.set_title(title)
        axis.set_xlabel("Lateral coordinate")
        figure.colorbar(handle, ax=axis, shrink=0.75)
    axes[0].set_ylabel("Depth (mm)")
    figure.savefig(output / "fixed_baseline_bscan_demo.png", dpi=200)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("output/fixed_baseline_bscan"))
    parser.add_argument("--no-plot", action="store_true")
    arguments = parser.parse_args()
    run(arguments.output, make_plot=not arguments.no_plot)
