import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

# Define paths
data_folder = Path("./inputs/IEA-15-240-RWT/Airfoils")
polar_paths = list((data_folder/"polar").glob("*.dat"))  # Polar = .dat
coord_paths = list((data_folder/"coord").glob("*.txt"))  # Shape = .txt
def read_coords(coord_paths):
    """
    Reads all coordinate files and returns a dictionary of DataFrames.
    """
    coords = {}
    for path in coord_paths:
        try:
            df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=8)
            coords[path.name] = df
        except Exception as e:
            print(f" Failed to read {path.name}: {e}")
    return coords


def read_polar(polar_paths):
    """
    Reads all polar files and returns a dictionary:
    {index: (alpha_array, Cl_array, Cd_array)}
    """
    polars = {}
    for path in polar_paths:
        try:
            df = pd.read_csv(
                path,
                sep=r"\s+",
                header=None,
                skiprows=20,
                usecols=[0, 1, 2],
                names=["alpha", "Cl", "Cd"]
            )
            # Extract airfoil index from filename
            idx = int(path.stem.split("_")[-1])  # 'Polar_00' → 0
            polars[idx] = (df["alpha"].values, df["Cl"].values, df["Cd"].values)
        except Exception as e:
            print(f" Failed to read {path.name}: {e}")
    return polars
polars_data = read_polar(polar_paths)
alpha, Cl, Cd = polars_data[0]  # AF00
print(Cl[:5])



# Call functions
coords_data = read_coords(coord_paths)
polars_data = read_polar(polar_paths)

# Quick test print
print(f"Loaded {len(coords_data)} coord files and {len(polars_data)} polar files.")  


def plot_airfoils(coords_data):
    """
    Plot all airfoil shapes from parsed coordinate DataFrames.
    """
    plt.figure(figsize=(10, 6))
    for filename, df in coords_data.items():
        label = filename.split('_')[2]  # e.g., 'AF00'
        plt.plot(df.iloc[:, 0], df.iloc[:, 1], label=label)
    plt.xlabel('x/c')
    plt.ylabel('y/c')
    plt.title('Airfoil Shapes')
    plt.axis('equal')
    # plt.legend(ncol=2, fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_airfoils(coords_data) 