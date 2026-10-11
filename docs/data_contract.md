# Input and Output Contract

## Required preprocessing

1. Supply linear OCT intensity after background subtraction and any documented system-response correction.
2. Do not supply logarithmic intensity, decibel display values, or gamma-corrected images.
3. Use a strictly increasing, uniformly sampled physical depth vector in millimetres.
4. Keep the array layout as depth-by-lateral for B-scans.
5. Define layer or interval boundaries in physical depth, not row indices.

## Required application-specific information

- interval or layer boundaries;
- fitting windows when the full interval contains reflective interfaces or unreliable samples;
- representative A-lines or remote control columns for fixed layer-wise fitting;
- optional post-reconstruction ROIs for evaluation.

## Outputs

- `mu0`: fitted or fixed layer-wise zeroth-order baseline, in mm⁻¹;
- `mu1`: signed local first-order correction, in mm⁻¹;
- `total`: `mu0 + mu1`, in mm⁻¹;
- `fits`: fitted attenuation, effective amplitude, RMSE, R², and sample count where applicable;
- `conventional_depth_resolved`: distal-tail comparison estimate, in mm⁻¹.

Samples outside supplied interval boundaries are returned as `NaN`, making the reconstruction support explicit.

## Fit-quality checks

Before interpreting a result, inspect the fitting window, parameter positivity, RMSE, R², and sensitivity to small changes in the window. A low R² can indicate interface contamination, insufficient decay, noise-floor dominance, or model mismatch. The software reports the fit; it does not convert a weak fit into an independent optical-property measurement.

## Reproducibility metadata to retain

For each result, retain:

- input dataset identifier and preprocessing version;
- axial spacing and physical units;
- interval boundaries and fitting windows;
- control or representative column indices;
- fitted `mu0`, amplitude, RMSE, R², and sample count;
- software version and random seed for synthetic or stochastic processing.
