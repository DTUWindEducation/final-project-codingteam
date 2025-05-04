"""Main script to run and visualize the full BEM workflow."""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Add module path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))

from bem.solver import BEMSolver, plot_cp_ct, plot_bem_vs_measured, plot_cp_ct_contours
from bem.geometry import load_geometry, plot_geometry
from bem.airfoil_data import read_polar, plot_airfoils, read_coords
from bem.operation import blad_operation, plot_operational_subplots

# Minimal input paths
base_dir = Path("./inputs/IEA-15-240-RWT")
geometry_path = base_dir / "IEA-15-240-RWT_AeroDyn15_blade.dat"
polars_path = (base_dir / "Airfoils" / "polar").glob("*.dat")
coords_path = (base_dir / "Airfoils" / "coord").glob("*.txt")
operation_path = base_dir / "IEA_15MW_RWT_Onshore.opt"


def run_demo():
    """Run the BEM demo with IEA-15-240-RWT data."""
    
    # Load data
    polars = read_polar(polars_path)
    geometry = load_geometry(geometry_path)
    operation_df = blad_operation(operation_path)
    coords = read_coords(coords_path)

    # Plot input data
    plot_geometry(geometry)                 
    plot_operational_subplots(operation_df)
    plot_airfoils(coords)

    bem = BEMSolver(geometry, polars)   # Class excutes BEM solver

    power_output, thrust_output = [], []

    # Compute performance for each row in the operation DataFrame
    for _, row in operation_df.iterrows():
        v_0 = row["V0"]
        omega_rps = row["omega_rpm"] * 2 * np.pi / 60
        pitch_deg = row["pitch_deg"]
        thrust, _, power = bem.compute_performance(v_0, omega_rps, pitch_deg)
        thrust_output.append(thrust / 1000)
        power_output.append(power / 1000)
    # Fast grid: lower resolution
    pitch_range = (-5, 5)
    lambda_range = (6, 11)
    V0 = 10.59  # rated wind speed
    R = 240 / 2 # rotor radius

    pitch_grid, lambda_grid, cp_grid, ct_grid, max_cp_point, max_ct_point = bem.compute_cp_ct_surface(
    pitch_range, lambda_range, V0, R, n_points=20
    )
    # Plot cp and ct contours
    plot_cp_ct_contours(pitch_grid, lambda_grid, cp_grid, ct_grid,
                        design_point=(0.0, 9.0),
                        max_cp_point=max_cp_point,
                        max_ct_point=max_ct_point)

    # Plot power and thrust output
    plot_bem_vs_measured(operation_df, power_output, thrust_output)
    # Plot cp and ct for each operation point as a function of V0
    plot_cp_ct(operation_df["V0"].values, power_output, thrust_output, R=240 / 2)


if __name__ == "__main__":
    run_demo()
