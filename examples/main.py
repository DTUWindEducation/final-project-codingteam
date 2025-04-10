import pandas as pd
from pathlib import Path
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
data_folder = Path("./inputs/IEA-15-240-RWT/Airfoils")


class Airfoil:
    def __init__(self, polar_files, coord_files, name=None):
        self.name = name
        self.polar_files = polar_files
        self.coord_files = coord_files
        self.polar_data = {}  # Store polar data here

    def read_all_polar(self):
        """
        Reads all polar files and stores in self.polar_data
        """
        for polar_file in self.polar_files:
            data = []
            try:
                with open(polar_file, "r") as f:
                    for i, line in enumerate(f):
                        if i < 20:
                            continue
                        if line.strip().startswith("!") or not line.strip():
                            continue
                        try:
                            alpha, cl, cd, cm = map(float, line.strip().split()[:4])
                            data.append((alpha, cl, cd, cm))
                        except ValueError:
                            continue

                df = pd.DataFrame(data, columns=["Alpha", "Cl", "Cd", "Cm"])
                self.polar_data[polar_file.name] = df

            except FileNotFoundError:
                print(f"File not found: {polar_file}")

    def interpolate_polar(self, alpha):
        """
        Interpolates Cl and Cd for the given alpha from all polar files.
        Returns a dictionary: { filename: (Cl, Cd) }
        """
        results = {}
        for name, df in self.polar_data.items():
            try:
                cl_interp = interp1d(df["Alpha"], df["Cl"], kind="linear", fill_value ="extrapolate")
                cd_interp = interp1d(df["Alpha"], df["Cd"], kind="linear", fill_value ="extrapolate")
                cl_val = float(cl_interp(alpha))
                cd_val = float(cd_interp(alpha))
                results[name] = (cl_val, cd_val)
            except Exception as e:
                print(f"Interpolation failed for {name}: {e}")
                results[name] = (None, None)
        return results
    def read_all_coords(self):
        """
        Reads all coordinate files and returns a dictionary of DataFrames.
        """
        coords = {}
        for coord_file in self.coord_files:
            try:
                df = pd.read_csv(coord_file, delim_whitespace=True, header=None, skiprows=8)
                coords[coord_file.name] = df
            except FileNotFoundError:
                print(f"File not found: {coord_file}")
        return coords
    
    
# Get all files
polar_files = list(data_folder.glob("IEA-15-240-RWT_AeroDyn15_Polar_*.dat"))
coord_files = list(data_folder.glob("IEA-15-240-RWT_AF*_Coords*.txt"))

# Safety check
if not polar_files or not coord_files:
    print("No files found.")
    exit()

# Create airfoil object
airfoil = Airfoil(polar_files=polar_files, coord_files=coord_files)

# Load data once
airfoil.read_all_polar()
airfoil.read_all_coords()
print(airfoil.read_all_coords())
# Interpolate Cl and Cd at a specific Alpha (e.g. 10 degrees)
alpha_target = 10
interpolated_values = airfoil.interpolate_polar(alpha=alpha_target)

# Print results
for fname, (cl, cd) in interpolated_values.items():
    print(f"{fname}: Alpha={alpha_target}° → Cl={cl:.3f}, Cd={cd:.3f}")

figure, ax = plt.subplots()
for fname, df in airfoil.polar_data.items():
    ax.plot(df["Alpha"], df["Cl"], label=f"{fname} Cl")
    ax.plot(df["Alpha"], df["Cd"], label=f"{fname} Cd") 
ax.axvline(x=alpha_target, color='r', linestyle='--', label=f"Alpha={alpha_target}°")
ax.set_xlabel("Alpha (degrees)")
ax.set_ylabel("Coefficient")
ax.set_title("Polar Data")
#ax.legend()
plt.show()
