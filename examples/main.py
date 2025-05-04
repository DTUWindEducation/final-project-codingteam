# """Main script to run and visualize the full BEM workflow."""

# import sys
# from pathlib import Path
# import numpy as np
# import matplotlib.pyplot as plt

# # Add module path
# sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))
# from bem.solver import BEMSolver, plot_cp_ct, plot_bem_vs_measured, plot_cp_ct_contours
# from bem.geometry import load_geometry, plot_geometry
# from bem.airfoil_data import read_polar, plot_airfoils, read_coords
# from bem.operation import blad_operation, plot_operational_subplots

# def main():
#     base_dir = Path("./inputs/IEA-15-240-RWT")
#     polars_path = (base_dir / "Airfoils" / "polar").glob("*.dat")
#     coords_path = (base_dir / "Airfoils" / "coord").glob("*.txt")
#     geometry_path = base_dir / "IEA-15-240-RWT_AeroDyn15_blade.dat"
#     operation_path = base_dir / "IEA_15MW_RWT_Onshore.opt"

#     # Load data
#     polars = read_polar(polars_path)
#     geometry = load_geometry(geometry_path)
#     operation_df = blad_operation(operation_path)
#     coords = read_coords(coords_path)

#     # Plot data
#     plot_geometry(geometry)
#     plot_operational_subplots(operation_df)
#     plot_airfoils(coords)

#     # Run BEM model

#     bem = BEMSolver(geometry, polars)
#     power_output, thrust_output = [], []

#     for _, row in operation_df.iterrows():
#         v_0 = row["V0"]
#         omega_rps = row["omega_rpm"] * 2 * np.pi / 60
#         pitch_deg = row["pitch_deg"]
#         thrust, _, power = bem.compute_performance(v_0, omega_rps, pitch_deg)
#         thrust_output.append(thrust / 1000)
#         power_output.append(power / 1000)

#     # CP-CT surface analysis
#     pitch_range = (-5, 5)
#     lambda_range = (7, 12)
#     V0 = 10.59  # m/s # rated wind speed
#     R = 240 / 2

#     pitch_grid, lambda_grid, cp_grid, ct_grid, max_cp_point, max_ct_point = bem.compute_cp_ct_surface(
#     pitch_range, lambda_range, V0, R, n_points=60
#     )
#     plot_cp_ct_contours(pitch_grid, lambda_grid, cp_grid, ct_grid,
#                         design_point=(0.0, 9.0),
#                         max_cp_point=max_cp_point,
#                         max_ct_point=max_ct_point)

#     # Plots
#     plot_bem_vs_measured(operation_df, power_output, thrust_output)
#     plot_cp_ct(operation_df["V0"].values, power_output, thrust_output, R=240 / 2)

# # if __name__ == "__main__":
# main()

def run_demo():
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
    
    # Load data
    polars = read_polar(polars_path)
    geometry = load_geometry(geometry_path)
    operation_df = blad_operation(operation_path)
    coords = read_coords(coords_path)

    # Plot data
    plot_geometry(geometry)
    plot_operational_subplots(operation_df)
    plot_airfoils(coords)
    #Load subset
    
    bem = BEMSolver(geometry, polars)
    power_output, thrust_output = [], []
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
    R = 240 / 2

    pitch_grid, lambda_grid, cp_grid, ct_grid, max_cp_point, max_ct_point = bem.compute_cp_ct_surface(
    pitch_range, lambda_range, V0, R, n_points=20
    )
    plot_cp_ct_contours(pitch_grid, lambda_grid, cp_grid, ct_grid,
                        design_point=(0.0, 9.0),
                        max_cp_point=max_cp_point,
                        max_ct_point=max_ct_point)

    # Plots
    plot_bem_vs_measured(operation_df, power_output, thrust_output)
    plot_cp_ct(operation_df["V0"].values, power_output, thrust_output, R=240 / 2)


if __name__ == "__main__":
    run_demo()
