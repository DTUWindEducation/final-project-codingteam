import sys
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Add the src directory to the Python path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))

from bem.solver import BEMSolver, plot_cp_ct, plot_cp_ct_contours
from bem.geometry import load_geometry
from bem.airfoil_data import read_polar


# === FIXTURES === #

@pytest.fixture
def bem_fixture():
    geom_path = Path("./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat")
    polar_path = Path("./inputs/IEA-15-240-RWT/Airfoils/polar").glob("*.dat")
    geometry = load_geometry(geom_path)
    polars = read_polar(polar_path)
    return BEMSolver(geometry, polars)

@pytest.fixture
def synthetic_bem():
    geometry = pd.DataFrame({
        "r": [0, 20],
        "chord": [2, 2],
        "twist_deg": [5, 5],
        "airfoil_id": [1, 1]
    })
    polars = {
        0: (
            np.array([0, 10, 20]),
            np.array([0.5, 1.0, 1.2]),
            np.array([0.01, 0.02, 0.03])
        )
    }
    return BEMSolver(geometry, polars)


# === UNIT TESTS === #

def test_interpolate_polar(bem_fixture):
    cl, cd = bem_fixture.interpolate_polar(0, 5.0)
    assert isinstance(cl, (float, np.floating))
    assert isinstance(cd, (float, np.floating))
    assert cl > 0

def test_interpolate_out_of_range(bem_fixture):
    cl, cd = bem_fixture.interpolate_polar(0, 100.0)  # extrapolate
    assert isinstance(cl, float) and isinstance(cd, float)

def test_compute_induction_factors():
    sigma = 0.05
    phi = np.radians(30)
    Cn, Ct = 1.0, 0.8
    bem = BEMSolver(None, None)
    a, a_prime = bem.compute_induction_factors(sigma, phi, Cn, Ct)
    assert 0 <= a <= 1
    assert -0.5 <= a_prime <= 1

def test_compute_performance(bem_fixture):
    V0 = 8.0
    omega = 1.0
    pitch = 2.0
    T, M, P = bem_fixture.compute_performance(V0, omega, pitch)
    assert isinstance(T, (float, np.floating))
    assert isinstance(M, (float, np.floating))
    assert isinstance(P, (float, np.floating))
    assert T > 0
    assert P > 0

def test_compute_performance_with_synthetic_data(synthetic_bem):
    T, M, P = synthetic_bem.compute_performance(8.0, 1.0, 2.0)
    assert T > 0 and M > 0 and P > 0

def test_r_equals_zero_handling(synthetic_bem):
    # Verifies that r=0 doesn't break the computation
    T, M, P = synthetic_bem.compute_performance(8.0, 1.0, 2.0)
    assert T > 0 and M > 0 and P > 0


# === PLOTTING === #

def test_plot_cp_ct(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    V0 = np.array([6, 8, 10])
    P = np.array([300, 600, 900])  # kW
    T = np.array([100, 200, 300])  # kN
    plot_cp_ct(V0, P, T, R=120)  # smoke test

def test_plot_bem_vs_measured(monkeypatch):
    from bem.solver import plot_bem_vs_measured
    monkeypatch.setattr(plt, "show", lambda: None)

    df = pd.DataFrame({
        "V0": [6, 8, 10],
        "power_kw": [100, 300, 600],
        "thrust_kN": [50, 150, 250]
    })
    plot_bem_vs_measured(df, [110, 310, 610], [55, 155, 255])  # smoke test

def test_plot_cp_ct_surface_runs(monkeypatch, synthetic_bem):
    monkeypatch.setattr(plt, "show", lambda: None)  # prevent GUI popup

    pitch_range = np.linspace(-5, 5, 5)
    tsr_range = np.linspace(7, 11, 5)

    pitch_grid, lambda_grid, cp_grid, ct_grid, max_cp_point, max_ct_point = synthetic_bem.compute_cp_ct_surface(
        pitch_range=(pitch_range[0], pitch_range[-1]),
        lambda_range=(tsr_range[0], tsr_range[-1]),
        V0=10,
        R=synthetic_bem.geometry["r"].values[-1]
    )

    plot_cp_ct_contours(
        pitch_grid, lambda_grid, cp_grid, ct_grid,
        design_point=(0, 9),
        max_cp_point=max_cp_point,
        max_ct_point=max_ct_point
    )