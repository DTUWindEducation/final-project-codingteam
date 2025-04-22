import sys
from pathlib import Path

# Add the src directory to Python path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))
from solver import BEM
from operation import blad_operation
from geometry import blad_geometry
from airfoil_data import read_polar
from solver import plot_cp_ct

import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
# Load data
geometry = blad_geometry(Path("./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat"))
polars = read_polar(Path("./inputs/IEA-15-240-RWT/Airfoils/polar").glob("*.dat"))
operation = blad_operation(Path("./inputs/IEA-15-240-RWT/IEA_15MW_RWT_Onshore.opt"))

# Initialize BEM model
bem = BEM(geometry, polars)

# Preallocate arrays
P_out, T_out = [], []
for idx, row in operation.iterrows():
    V0 = row["V0"]
    omega_rps = row["omega_rpm"] * 2 * np.pi / 60  # convert to rad/s
    pitch = row["pitch_deg"]
    T, M, P = bem.compute_performance(V0, omega_rps, pitch)
    T_out.append(T / 1000)  # N -> kN
    P_out.append(P / 1000)  # W -> kW

# Plot results
plt.figure()
plt.plot(operation["V0"], P_out, label="Predicted Power [kW]")
plt.plot(operation["V0"], T_out, label="Predicted Thrust [kN]")
plt.xlabel("Wind Speed [m/s]")
plt.ylabel("Performance")
plt.title("BEM Power and Thrust Curves")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot results with reference (measured) data
plt.figure()
plt.plot(operation["V0"], P_out, label="Predicted Power [kW]", markersize=2)
plt.plot(operation["V0"], T_out, label="Predicted Thrust [kN]", markersize=2)

# Plot reference (measured) values
plt.plot(operation["V0"], operation["power_kw"], 'o--', label="Measured Power [kW]", markersize=4)
plt.plot(operation["V0"], operation["thrust_kN"], 'o--', label="Measured Thrust [kN]", markersize=4)

plt.xlabel("Wind Speed [m/s]")
plt.ylabel("Performance")
plt.title("BEM vs Measured Power and Thrust Curves")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plot_cp_ct(operation["V0"].values, P_out, T_out, R=240/2)

