"""Functions for loading and plotting blade geometry data."""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def load_geometry(path_geometry):
    """
    Load blade geometry and return a DataFrame with:
    - r: spanwise position [m]
    - chord: chord length [m]
    - twist_deg: twist angle [deg]
    - airfoil_id: airfoil index (1–50)
    """
    path_geometry = Path(path_geometry)
    geometry_df = pd.read_csv(path_geometry, sep=r'\s+', header=None, skiprows=6)

    geometry_df.columns = [
        "r", "BlCrvAC", "BlSwpAC", "BlCrvAng", "twist_deg",
        "chord", "airfoil_id", "BlCb", "BlCenBn", "BlCenBt"
    ]

    geometry_df = geometry_df[["r", "chord", "twist_deg", "airfoil_id"]]
    geometry_df["airfoil_id"] = geometry_df["airfoil_id"].astype(int)

    return geometry_df


def plot_geometry(geometry_df):
    """
    Plot blade geometry: chord and twist angle vs. spanwise position.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(geometry_df["r"], geometry_df["chord"], label="Chord [m]")
    plt.plot(geometry_df["r"], geometry_df["twist_deg"], label="Twist [deg]")
    plt.xlabel("Spanwise position r [m]")
    plt.ylabel("Value")
    plt.title("Blade Geometry")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# # Uncomment to run
# if __name__ == "__main__":
#     geometry_path = Path("inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat")
#     blade_geometry = load_geometry(geometry_path)
#     plot_geometry(blade_geometry)
