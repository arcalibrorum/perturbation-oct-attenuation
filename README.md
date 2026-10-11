# Perturbation-Based OCT Attenuation Reconstruction

Reference MATLAB and Python implementations accompanying the manuscript
**“Perturbation-Based Reconstruction of Depth-Resolved OCT Attenuation without Requiring a Fully Decayed Signal Tail.”**

The package is the method-level reference implementation of the proposed reconstruction. It reconstructs effective depth-resolved OCT attenuation from **linear OCT intensity** and exposes the same zeroth-order fit, local first-order correction, piecewise reset, and conventional comparison described in the manuscript.

The examples are deliberately compact demonstrations of how to use the method. They are not intended to regenerate every manuscript figure or redistribute the Monte Carlo and experimental datasets.

## Central method claims represented in code

- attenuation at the current depth is reconstructed without integrating the complete distal signal tail;
- a fitted zeroth-order exponential provides the local baseline;
- the analytical first-order term restores local deviations from that baseline;
- the coordinate, boundary condition, and cumulative integral are reset in every depth interval;
- fixed layer-wise parameters can be reused across a B-scan while the first-order term is evaluated independently in every A-line;
- the no-vessel-prior workflow estimates baseline parameters from remote control columns and does not accept a vessel mask during reconstruction;
- experimental outputs are interpreted as effective OCT attenuation under the model assumptions.

## What is included

- matched MATLAB and Python public APIs;
- local zeroth-order exponential fitting;
- analytical first-order correction;
- piecewise reconstruction with a reset at every interval boundary;
- B-scan reconstruction with fixed layer-wise parameters;
- control-column parameter estimation without a vessel mask;
- the conventional depth-resolved estimator;
- synthetic examples and shared cross-language reference values;
- dependency-free numerical tests in both languages.

The repository contains synthetic data only. The Monte Carlo and experimental OCT datasets used in the manuscript are not redistributed and remain available from the authors on reasonable request.

## Numerical conventions

| Quantity | Convention |
|---|---|
| OCT signal | Linear intensity; never logarithmic or decibel values |
| Depth `z` | mm |
| Attenuation `mu0`, `mu1`, `total` | mm⁻¹ |
| Local amplitude `A` | Effective scale in `I = A mu0 exp(-2 mu0 xi)` |
| Local coordinate `xi` | Reset to zero at each interval upper boundary |
| Local integral | Left cumulative sum on a uniform axial grid; first value is zero |
| Interval edge | An internal-edge sample belongs to the deeper interval |
| B-scan layout | `(depth, lateral)` |

These choices match the reconstruction description in the manuscript. Display compression, confocal focusing, sensitivity roll-off, and lateral smoothing of displayed baseline maps are not applied to reconstruction inputs by this package.

Dataset-specific depths, attenuation values, control columns, layer boundaries, and random seeds appear only in the example scripts. They are not embedded in the reconstruction functions. The default positive fitting range, `1e-5` to `50 mm⁻¹`, is an exposed numerical safeguard and can be replaced through the `mu_bounds` argument.

## Python quick start

Python 3.9 or newer is required.

```bash
python -m pip install -e .
python -m unittest discover -s python/tests -v
python python/examples/demo_finite_depth.py --output output/finite_depth
python python/examples/demo_fixed_baseline_bscan.py --output output/fixed_baseline
```

Plot generation is optional:

```bash
python -m pip install -e '.[plot]'
```

Minimal reconstruction:

```python
import numpy as np
import octattenuation as oa

z, reference_mu, _ = oa.finite_depth_profiles()
intensity = oa.simulate_beer_lambert_intensity(z, reference_mu, amplitude=1.67)
edges = [0.0, 0.6, 1.2, z[-1]]
fit_windows = [(0.0, 0.298), (0.0, 0.598), (0.0, 0.198)]

result = oa.reconstruct_piecewise(z, intensity, edges, fit_windows=fit_windows)
mu0 = result.mu0
mu1 = result.mu1
mu = result.total
mu_dr = oa.conventional_depth_resolved(intensity, spacing=0.002)
```

For the no-vessel-prior workflow:

```python
fits = oa.fit_control_baselines(z, bscan, layer_boundaries, control_columns)
result = oa.reconstruct_bscan(z, bscan, layer_boundaries, fits)
```

`control_columns` are chosen remotely from the target structure. A vessel mask is not accepted by either function. Vessel and surrounding-tissue masks are supplied only afterward to `ablation_metrics`.

## MATLAB quick start

MATLAB R2021b or newer is recommended. The reconstruction functions use base MATLAB only.

```matlab
addpath('matlab', 'matlab/tests', 'matlab/examples');
run_all_tests
demo_finite_depth
demo_fixed_baseline_bscan
```

Minimal reconstruction:

```matlab
z = (0:999).' .* 0.002;  % mm
reference_mu = 0.6 .* ones(size(z));
intensity = octattenuation.simulate_beer_lambert_intensity(z, reference_mu, 1.67);
edges = [0, 0.6, 1.2, z(end)];

result = octattenuation.reconstruct_piecewise(z, intensity, edges);
mu_dr = octattenuation.conventional_depth_resolved(intensity, 0.002);
```

For the no-vessel-prior workflow:

```matlab
fits = octattenuation.fit_control_baselines(z, bscan, layer_boundaries, control_columns);
result = octattenuation.reconstruct_bscan(z, bscan, layer_boundaries, fits);
```

## Manuscript correspondence

The main mappings are:

| Manuscript operation | Python | MATLAB |
|---|---|---|
| Zeroth-order least-squares fit | `fit_zeroth_order` | `octattenuation.fit_zeroth_order` |
| Local first-order solution | `first_order_correction` | `octattenuation.first_order_correction` |
| Final local estimate | `reconstruct_interval` | `octattenuation.reconstruct_interval` |
| Piecewise coordinate/integral reset | `reconstruct_piecewise` | `octattenuation.reconstruct_piecewise` |
| Fixed layer-wise B-scan reconstruction | `reconstruct_bscan` | `octattenuation.reconstruct_bscan` |
| Remote-control median fitting | `fit_control_baselines` | `octattenuation.fit_control_baselines` |
| Conventional distal-tail estimator | `conventional_depth_resolved` | `octattenuation.conventional_depth_resolved` |
| Post hoc vessel contrast and MAE | `ablation_metrics` | `octattenuation.ablation_metrics` |

See [docs/paper_mapping.md](docs/paper_mapping.md) for the equation-level and method-workflow mapping.

The correspondence document is an implementation audit, not a promise of pixel-for-pixel figure reproduction. The public scientific contract is the mathematical model, numerical conventions, and reconstruction workflow.

## Scope and limitations

- The output is an **effective OCT attenuation estimate** under the model assumptions; it is not automatically an intrinsic microscopic optical coefficient.
- The first-order approximation is most reliable when local deviations from the fitted zeroth-order exponential are moderate.
- The method reduces dependence on unreliable distal signal below the current depth; it does not eliminate multiple scattering.
- Layer boundaries and fitting intervals remain application-specific inputs.
- Fit quality should be inspected. Weak or contaminated fitting intervals must not be interpreted as independent optical-property truth.
- The public synthetic demonstrations validate the implementation and workflow. They are not substitutes for the manuscript's Monte Carlo or experimental datasets.

## Validation status

- Python package build and isolated import: passed.
- Python numerical suite: 6 tests passed.
- MATLAB R2021b numerical suite: all tests passed.
- Both synthetic workflows completed in Python and MATLAB without reconstruction errors.
- Cross-language finite-depth outputs agreed to a maximum absolute difference below `1.4e-8` across `mu0`, `mu1`, and the combined estimate.
- The shared first-order reference case agreed within `2e-10` in both test suites.

The repository also includes public continuous-integration workflows for Python and MATLAB.

## Repository layout

```text
python/src/octattenuation/  Python package
python/examples/            Python demonstrations
python/tests/               Python numerical tests
matlab/+octattenuation/     MATLAB package
matlab/examples/            MATLAB demonstrations
matlab/tests/               MATLAB numerical tests
shared/test_data/           Cross-language golden case
docs/                       Equation and data-contract documentation
```

## Citation

Please cite the accompanying manuscript if this implementation contributes to published work. Citation metadata are provided in `CITATION.cff` and can be updated with the final DOI after publication.

After the public repository and archival release are created, the manuscript can link to this implementation in its code-availability statement. The archived software DOI can then be added to `CITATION.cff` so that use of the implementation can be credited independently while the article remains the primary scientific citation.

## License

BSD 3-Clause License. See [LICENSE](LICENSE).
