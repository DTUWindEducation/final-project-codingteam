import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add src directory to path for importing bem
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))

from bem.solver import BEMSolver, plot_cp_ct, plot_cp_ct_contours, plot_bem_vs_measured


# === FIXTURES === #

@pytest.fixture
def synthetic_geometry():
    return pd.DataFrame({
        "r": [1, 2],
        "chord": [2, 2],
        "twist_deg": [5, 5],
        "airfoil_id": [1, 1]
    })


@pytest.fixture
def synthetic_polars():
    return {
        0: (
            np.array([-10, 0, 10, 20]),
            np.array([0.5, 0.8, 1.2, 1.5]),  # More realistic Cl
            np.array([0.01, 0.02, 0.03, 0.04])
        )
    }


@pytest.fixture
def solver(synthetic_geometry, synthetic_polars):
    return BEMSolver(synthetic_geometry, synthetic_polars)


# === UNIT TESTS === #

def test_interpolate_polar(solver):
    cl, cd = solver.interpolate_polar(0, 10)
    assert isinstance(cl, float)
    assert isinstance(cd, float)
    assert cl > 0
    assert cd > 0

def test_interpolate_polar_out_of_bounds(solver):
    cl, cd = solver.interpolate_polar(0, 100)  # Extrapolation
    assert isinstance(cl, float)
    assert isinstance(cd, float)

def test_induction_factors(solver):
    a, a_prime = solver.compute_induction_factors(0.1, np.radians(30), 1.0, 0.5)
    assert 0 <= a <= 1
    assert -0.5 <= a_prime <= 1

def test_compute_performance(solver):
    T, M, P = solver.compute_performance(8, 1.0, 2.0)
    assert isinstance(T, float) and np.isfinite(T)
    assert isinstance(M, float) and np.isfinite(M)
    assert isinstance(P, float) and np.isfinite(P)

def test_r_equals_zero_handling():
    geometry = pd.DataFrame({
        "r": [0, 2],
        "chord": [2, 2],
        "twist_deg": [5, 5],
        "airfoil_id": [1, 1]
    })
    polars = {
        0: (
            np.array([-10, 0, 10, 20]),
            np.array([0.5, 0.8, 1.2, 1.5]),
            np.array([0.01, 0.02, 0.03, 0.04])
        )
    }
    solver = BEMSolver(geometry, polars)
    T, M, P = solver.compute_performance(8.0, 1.0, 2.0)
    assert isinstance(T, float) and np.isfinite(T)
    assert isinstance(M, float) and np.isfinite(M)
    assert isinstance(P, float) and np.isfinite(P)


# === SURFACE TEST === #

def test_cp_ct_surface(solver):
    pitch_range = (-2, 2)
    lambda_range = (6, 8)
    V0 = 8.0
    R = 120.0
    results = solver.compute_cp_ct_surface(pitch_range, lambda_range, V0, R, n_points=5)
    pitch_grid, lambda_grid, cp_grid, ct_grid, max_cp_point, max_ct_point = results
    assert cp_grid.shape == ct_grid.shape
    assert pitch_grid.shape == lambda_grid.shape
    assert np.isfinite(cp_grid).all()
    assert np.isfinite(ct_grid).all()


# === PLOTTING TESTS (Smoke tests) === #

def test_plot_cp_ct(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    V0 = np.array([6, 8, 10])
    P = np.array([300, 600, 900])  # kW
    T = np.array([100, 200, 300])  # kN
    plot_cp_ct(V0, P, T, R=120)

def test_plot_cp_ct_contours(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    pitch_vals = np.linspace(-2, 2, 5)
    lambda_vals = np.linspace(6, 8, 5)
    pitch_grid, lambda_grid = np.meshgrid(pitch_vals, lambda_vals)
    cp_grid = np.random.rand(*pitch_grid.shape)
    ct_grid = np.random.rand(*lambda_grid.shape)
    plot_cp_ct_contours(pitch_grid, lambda_grid, cp_grid, ct_grid)

def test_plot_bem_vs_measured(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    df = pd.DataFrame({
        "V0": [6, 8, 10],
        "power_kw": [100, 300, 600],
        "thrust_kN": [50, 150, 250]
    })
    plot_bem_vs_measured(df, [110, 310, 610], [55, 155, 255])
