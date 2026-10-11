# Manuscript-to-Code Mapping

This document audits the method-level correspondence between the public implementation and the mathematical and experimental operations described in the accompanying manuscript. It does not define a figure-reproduction package. LaTeX equation labels are used because printed equation numbers may change during production.

## Core equations

| Manuscript item | Implemented operation | Python | MATLAB | Verification |
|---|---|---|---|---|
| Beer–Lambert model, `eq:beer_lambert`, and OCT signal relation, `eq:oct_signal_relation` | Generate ideal single-scattering linear intensity from a prescribed attenuation profile | `simulate_beer_lambert_intensity` | `octattenuation.simulate_beer_lambert_intensity` | Finite-depth example |
| Local amplitude definition, `A_ab = 2 C_ab L_0,ab R_ab` | Represent multiplicative system, incident-irradiance, and backscatter scaling by one fitted positive parameter | `ZerothOrderFit.amplitude` | `fit_result.amplitude` | Exact-exponential fit test |
| Zeroth-order fit, `eq:zeroth_order_fit` | Minimize `sum[I_k - A mu0 exp(-2 mu0 xi_k)]^2` | `fit_zeroth_order` | `octattenuation.fit_zeroth_order` | Exact-exponential fit test |
| Local first-order solution following `eq:mu1_ode` | Evaluate `mu1 = [I/A + (2 mu0/A) integral(I) - mu0] exp(2 mu0 xi)` | `first_order_correction` | `octattenuation.first_order_correction` | Forward/inverse first-order test and shared CSV |
| Local boundary condition | Enforce `mu1(0) = I(0)/A - mu0` because the cumulative integral is exactly zero at the first sample | `first_order_correction` | `octattenuation.first_order_correction` | First sample in shared CSV |
| Final estimate, `eq:reconstructed_mu` | Return `mu0`, `mu1`, and `mu0 + mu1` separately | `reconstruct_interval` | `octattenuation.reconstruct_interval` | Unit tests |
| Piecewise implementation paragraph following `eq:zeroth_order_fit` | Reset local coordinate, boundary condition, and cumulative integral at every interval | `reconstruct_piecewise` | `octattenuation.reconstruct_piecewise` | Boundary-reset test |
| Conventional comparison, `eq:conventional_dr` | Evaluate `I_i / [2 dz sum_(j=i)^N I_j]` | `conventional_depth_resolved` | `octattenuation.conventional_depth_resolved` | Closed-form array test and shared CSV |

The local integral uses a left cumulative sum on a uniform depth grid:

```text
integral[0] = 0
integral[i] = dz * sum(intensity[0:i])
```

This convention is implemented by `left_cumulative_integral` in both languages.
The implementation evaluates the truncated first-order approximation with the bookkeeping parameter set to `epsilon = 1`; it does not include second- or higher-order terms.

## Parameter fitting correspondence

The manuscript writes the zeroth-order estimate as a two-parameter nonlinear least-squares problem. The code solves the same objective by variable projection:

1. choose a bounded candidate `mu0`;
2. compute the least-squares-optimal positive amplitude `A` analytically;
3. minimize the residual sum of squares over `mu0` by the same golden-section algorithm in both languages;
4. report `mu0`, `A`, RMSE, R², and the fit sample count.

This formulation changes the numerical solver, not the fitted model or objective, and removes the need for MATLAB Optimization Toolbox.

## Connection to manuscript use cases

| Manuscript use case | Public method path | Scope of the public example |
|---|---|---|
| Finite-depth numerical A-lines | `demo_finite_depth` | Demonstrates the central distal-tail failure mode and the complete method call sequence using deterministic step/sinusoidal profiles and interface-avoiding fitting windows. |
| Monte Carlo evaluation | Core functions accept linear-intensity profiles and physical layer boundaries | Photon histories and manuscript input profiles are not part of the method-level release. Display-only focusing, roll-off, and compression are not reconstruction operations. |
| No-vessel-prior component assessment | `fit_control_baselines` → `reconstruct_bscan` → `ablation_metrics`; see `demo_fixed_baseline_bscan` | Demonstrates that control columns determine fixed layer-wise parameters while evaluation masks are introduced only after reconstruction. Its synthetic metrics are illustrative and are not the manuscript's Monte Carlo metrics. |
| Layered phantom or tissue B-scan | Fit representative layer intervals and pass the fixed parameters to `reconstruct_bscan` | The same method path supports the phantom and bovine analyses without embedding dataset-specific values in the core functions. |
| Experimental first-order assessment | Reuse fixed baseline parameters with `reconstruct_bscan` and analyze `abs(result.mu1)` afterward | Surface alignment, profile smoothing, control subtraction, and percentile visualization are downstream assessment choices rather than reconstruction inputs. |

## Deliberately excluded operations

- logarithmic or decibel conversion before reconstruction;
- confocal focusing and sensitivity roll-off used only to render selected simulation panels;
- undocumented offsets, overrides, or empirical attenuation compensation;
- independent per-column parameters in the final experimental maps;
- vessel masks or vessel boundaries in the no-prior reconstruction;
- interpretation of experimental fitting benchmarks as independently measured intrinsic optical coefficients.

## Shared numerical case

`shared/test_data/cross_language_case.csv` contains one deterministic local first-order case. Both test suites reconstruct it using `mu0 = 0.7 mm^-1` and `A = 1.8` and compare the complete correction vector with the same stored values.
