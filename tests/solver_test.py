import sys
from pathlib import Path

# Add the src directory to the Python path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))

import pytest
import numpy as np
from bem.solver import BEMSolver, plot_cp_ct
from bem.geometry import blad_geometry
from bem.airfoil_data import read_polar

@pytest.fixture
def bem_fixture():
    geom_path = Path("./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat")
    polar_path = Path("./inputs/IEA-15-240-RWT/Airfoils/polar").glob("*.dat")
    geometry = blad_geometry(geom_path)
    polars = read_polar(polar_path)
    return BEMSolver(geometry, polars)


def test_interpolate_polar(bem_fixture):
    cl, cd = bem_fixture.interpolate_polar(0, 5.0)
    assert isinstance(cl, (float, np.floating))
    assert isinstance(cd, (float, np.floating))
    assert cl > 0


def test_compute_induction_factors():
    sigma = 0.05
    phi = np.radians(30)
    Cn, Ct = 1.0, 0.8
    bem = BEMSolver(None, None)
    a, a_prime = bem.compute_induction_factors(sigma, phi, Cn, Ct)
    assert 0 <= a <= 1
    assert -0.5 <= a_prime <= 1

def test_compute_performance(bem_fixture):
    V0 = 8.0  # m/s
    omega = 1.0  # rad/s
    pitch = 2.0  # deg
    T, M, P = bem_fixture.compute_performance(V0, omega, pitch)
    assert isinstance(T, (float, np.floating))
    assert isinstance(M, (float, np.floating))
    assert isinstance(P, (float, np.floating))
    assert T > 0
    assert P > 0

# def test_plot_cp_ct(monkeypatch):
#     monkeypatch.setattr("matplotlib.pyplot.show", lambda: None)
#     V0 = np.array([6, 8, 10])
#     P = np.array([300, 600, 900])  # kW
#     T = np.array([100, 200, 300])  # kN
#     plot_cp_ct(V0, P, T, R=120)  # Just ensure it runs without error
def test_plot_cp_ct(monkeypatch):
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, "show", lambda: None)  # Prevent plot from opening
    V0 = np.array([6, 8, 10])
    P = np.array([300, 600, 900])  # kW
    T = np.array([100, 200, 300])  # kN
    plot_cp_ct(V0, P, T, R=120)  # Just check that it runs

