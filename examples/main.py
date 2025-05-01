"""Main script to run and visualize the full BEM workflow."""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Add module path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))
from bem.solver import BEMSolver, plot_cp_ct, plot_bem_vs_measured
from bem.geometry import load_geometry, plot_geometry
from bem.airfoil_data import read_polar, plot_airfoils, read_coords
from bem.operation import blad_operation, plot_operational_subplots

def main():
    base_dir = Path("./inputs/IEA-15-240-RWT")
    polars_path = (base_dir / "Airfoils" / "polar").glob("*.dat")
    coords_path = (base_dir / "Airfoils" / "coord").glob("*.txt")
    geometry_path = base_dir / "IEA-15-240-RWT_AeroDyn15_blade.dat"
    operation_path = base_dir / "IEA_15MW_RWT_Onshore.opt"

    # Load data
    polars = read_polar(polars_path)
    geometry = load_geometry(geometry_path)
    operation_df = blad_operation(operation_path)
    coords = read_coords(coords_path)

    # Plot data
    plot_geometry(geometry)
    plot_operational_subplots(operation_df)
    plot_airfoils(coords)

    # Run BEM model
    bem = BEMSolver(geometry, polars)
    power_output, thrust_output = [], []

    for _, row in operation_df.iterrows():
        v_0 = row["V0"]
        omega_rps = row["omega_rpm"] * 2 * np.pi / 60
        pitch_deg = row["pitch_deg"]
        thrust, _, power = bem.compute_performance(v_0, omega_rps, pitch_deg)
        thrust_output.append(thrust / 1000)
        power_output.append(power / 1000)

    # Plots
    plot_bem_vs_measured(operation_df, power_output, thrust_output)
    plot_cp_ct(operation_df["V0"].values, power_output, thrust_output, R=240 / 2)

# if __name__ == "__main__":
main()
