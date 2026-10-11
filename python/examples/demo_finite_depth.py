"""Finite-depth A-line demonstration corresponding to the manuscript workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import octattenuation as oa


def run(output: Path, make_plot: bool = True) -> None:
    output.mkdir(parents=True, exist_ok=True)
    z, step_mu, smooth_mu = oa.finite_depth_profiles()
    edges = np.array([0.0, 0.6, 1.2, z[-1]])
    fit_windows = {
        "step": ((0.0, 0.298), (0.0, 0.598), (0.0, 0.198)),
        "smooth": ((0.0, 0.598), (0.0, 0.598), (0.0, 0.798)),
    }
    cases = {
        "step": oa.simulate_beer_lambert_intensity(z, step_mu, amplitude=1.67),
        "smooth": oa.simulate_beer_lambert_intensity(z, smooth_mu, amplitude=1.0),
    }
    references = {"step": step_mu, "smooth": smooth_mu}

    rows = []
    results = {}
    for name, intensity in cases.items():
        proposed = oa.reconstruct_piecewise(z, intensity, edges, fit_windows=fit_windows[name])
        conventional = oa.conventional_depth_resolved(intensity, z[1] - z[0])
        results[name] = (intensity, proposed, conventional)
        for index in range(z.size):
            rows.append(
                (
                    name,
                    z[index],
                    intensity[index],
                    references[name][index],
                    proposed.mu0[index],
                    proposed.mu1[index],
                    proposed.total[index],
                    conventional[index],
                )
            )

    header = "case,z_mm,intensity,reference_mu,mu0,mu1,total,conventional"
    np.savetxt(
        output / "finite_depth_results.csv",
        np.asarray(rows, dtype=object),
        delimiter=",",
        fmt=["%s"] + ["%.12g"] * 7,
        header=header,
        comments="",
    )

    if not make_plot:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Install the optional plotting dependency with: pip install -e '.[plot]'") from exc

    figure, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    for column, name in enumerate(("step", "smooth")):
        intensity, proposed, conventional = results[name]
        axes[0, column].semilogy(z, intensity, color="0.25")
        axes[0, column].set_title(f"{name.capitalize()} linear OCT intensity")
        axes[0, column].set_xlabel("Depth (mm)")
        axes[0, column].set_ylabel("Intensity (a.u.)")
        axes[1, column].plot(z, references[name], "k--", label="Reference")
        axes[1, column].plot(z, proposed.total, label="Perturbation-based")
        axes[1, column].plot(z, conventional, label="Conventional", alpha=0.75)
        axes[1, column].set_ylim(0, 3)
        axes[1, column].set_xlabel("Depth (mm)")
        axes[1, column].set_ylabel("Attenuation (mm$^{-1}$)")
    axes[1, 0].legend(frameon=False)
    figure.savefig(output / "finite_depth_demo.png", dpi=200)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("output/finite_depth"))
    parser.add_argument("--no-plot", action="store_true")
    arguments = parser.parse_args()
    run(arguments.output, make_plot=not arguments.no_plot)
