"""Numerical tests for the Python implementation."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path

import numpy as np

import octattenuation as oa


class CoreTests(unittest.TestCase):
    def test_zeroth_order_fit_recovers_exponential(self) -> None:
        z = np.arange(501, dtype=float) * 0.001
        expected_mu0 = 0.8
        expected_amplitude = 2.4
        intensity = expected_amplitude * expected_mu0 * np.exp(-2 * expected_mu0 * z)
        fit = oa.fit_zeroth_order(z, intensity, mu_bounds=(0.01, 3.0))
        self.assertAlmostEqual(fit.mu0, expected_mu0, places=8)
        self.assertAlmostEqual(fit.amplitude, expected_amplitude, places=8)
        self.assertGreater(fit.r_squared, 0.999999999)

    def test_first_order_inverse_matches_forward_model(self) -> None:
        z = np.arange(501, dtype=float) * 0.001
        mu0, amplitude = 0.8, 2.4
        expected = 0.05 * np.sin(2 * np.pi * z / z[-1])
        intensity = oa.simulate_first_order_intensity(z, mu0, expected, amplitude)
        recovered = oa.first_order_correction(z, intensity, mu0, amplitude)
        np.testing.assert_allclose(recovered, expected, atol=8.1e-4, rtol=0)

    def test_piecewise_resets_local_coordinate(self) -> None:
        z = np.arange(601, dtype=float) * 0.001
        edges = np.array([0.0, 0.3, 0.6])
        intensity = np.empty_like(z)
        first = z < 0.3
        intensity[first] = 2.0 * 0.4 * np.exp(-2 * 0.4 * z[first])
        local = z[~first] - z[~first][0]
        intensity[~first] = 1.2 * 0.9 * np.exp(-2 * 0.9 * local)
        result = oa.reconstruct_piecewise(z, intensity, edges, fixed_parameters=((0.4, 2.0), (0.9, 1.2)))
        self.assertAlmostEqual(result.mu1[np.searchsorted(z, 0.3)], 0.0, places=12)
        self.assertTrue(np.all(np.isfinite(result.total)))

    def test_conventional_estimator(self) -> None:
        intensity = np.array([4.0, 3.0, 2.0, 1.0])
        expected = intensity / (0.2 * np.array([10.0, 6.0, 3.0, 1.0]))
        np.testing.assert_allclose(oa.conventional_depth_resolved(intensity, 0.1), expected)

    def test_control_baselines_use_only_selected_columns(self) -> None:
        z = np.arange(101, dtype=float) * 0.005
        boundaries = np.array([0.0, 0.25, 0.5])
        image = np.empty((z.size, 5))
        for column in range(image.shape[1]):
            image[z < 0.25, column] = 1.6 * 0.4 * np.exp(-2 * 0.4 * z[z < 0.25])
            local = z[z >= 0.25] - z[z >= 0.25][0]
            image[z >= 0.25, column] = 1.1 * 0.9 * np.exp(-2 * 0.9 * local)
        image[:, 2] *= 1.4
        fits = oa.fit_control_baselines(z, image, boundaries, [0, 1, 3, 4], mu_bounds=(0.05, 2.0))
        self.assertAlmostEqual(fits[0].mu0, 0.4, places=7)
        self.assertAlmostEqual(fits[1].mu0, 0.9, places=7)
        self.assertEqual(fits[0].sample_count, 4)

    def test_shared_cross_language_case(self) -> None:
        path = Path(__file__).resolve().parents[2] / "shared" / "test_data" / "cross_language_case.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        z = np.array([float(row["z_mm"]) for row in rows])
        intensity = np.array([float(row["intensity"]) for row in rows])
        expected = np.array([float(row["mu1_expected"]) for row in rows])
        recovered = oa.first_order_correction(z, intensity, 0.7, 1.8)
        np.testing.assert_allclose(recovered, expected, rtol=2e-10, atol=2e-12)


if __name__ == "__main__":
    unittest.main()
