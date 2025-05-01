"""Functions for loading and plotting airfoil data."""

import re

import matplotlib.pyplot as plt
import pandas as pd


def read_coords(coord_paths):
    """
    Reads all coordinate files and returns a dictionary of DataFrames.
    """
    coords_dict = {}
    for path in coord_paths:
        try:
            coord_df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=8)
            coords_dict[path.name] = coord_df
        except (OSError, ValueError) as error:
            print(f" Failed to read {path.name}: {error}")
    return coords_dict


def read_polar(polar_file_paths):
    """
    Reads polar files and returns a dictionary with alpha, Cl, Cd arrays.
    """
    polars_dict = {}
    for file_path in polar_file_paths:
        try:
            with open(file_path, 'r', encoding='utf-8') as polar_file:
                lines = polar_file.readlines()

            num_data_lines = None
            for i, line in enumerate(lines):
                if 'NumAlf' in line:
                    match = re.search(r'\d+', line)
                    if match:
                        num_data_lines = int(match.group())
                        data_start = i + 2
                        break

            if num_data_lines is None:
                raise ValueError("Couldn't find 'NumAlf' line in file")

            data_rows = []
            for line in lines[data_start:data_start + num_data_lines]:
                if line.strip().startswith("!"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        alpha = float(parts[0])
                        lift = float(parts[1])
                        drag = float(parts[2])
                        data_rows.append((alpha, lift, drag))
                    except ValueError:
                        continue

            polar_df = pd.DataFrame(data_rows, columns=["alpha", "cl", "cd"])
            airfoil_idx = int(file_path.stem.split("_")[-1])
            polars_dict[airfoil_idx] = (
                polar_df["alpha"].values,
                polar_df["cl"].values,
                polar_df["cd"].values
            )
        except (OSError, ValueError) as error:
            print(f" Failed to read {file_path.name}: {error}")
    return polars_dict


def plot_airfoils(coords_dict):
    """
    Plot all airfoil shapes from parsed coordinate DataFrames.
    """
    plt.figure(figsize=(10, 6))
    for filename, coord_df in coords_dict.items():
        label = filename.split('_')[2]  # e.g., 'AF00'
        plt.plot(coord_df.iloc[:, 0], coord_df.iloc[:, 1], label=label)
    plt.xlabel('x/c')
    plt.ylabel('y/c')
    plt.title('Airfoil Shapes')
    plt.axis('equal')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# if __name__ == "__main__":
#     base_folder = Path("inputs/IEA-15-240-RWT/Airfoils")
#     polar_files = list((base_folder / "polar").glob("*.dat"))
#     coord_files = list((base_folder / "coord").glob("*.txt"))

#     coords = read_coords(coord_files)
#     polars = read_polar(polar_files)
