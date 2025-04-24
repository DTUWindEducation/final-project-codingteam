import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

geometry=Path("./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat")  # Geometry = .dat

def blad_geometry(geometry):
    """
    Load blade geometry and return a DataFrame with:
    - r: spanwise position [m]
    - chord: chord length [m]
    - twist_deg: twist angle [deg]
    - airfoil_id: airfoil index (1–50)
    """
    df = pd.read_csv(geometry, sep=r'\s+', header=None, skiprows=6)
    
    df.columns = [
        "r", "BlCrvAC", "BlSwpAC", "BlCrvAng", "twist_deg",
        "chord", "airfoil_id", "BlCb", "BlCenBn", "BlCenBt"
    ]

    # Clean up DataFrame
    df = df[["r", "chord", "twist_deg", "airfoil_id"]]
    df["airfoil_id"] = df["airfoil_id"].astype(int)

    return df

bladegeometry = blad_geometry(geometry)
print(bladegeometry)


def plot_geometry(bladegeometry):
    """
    Plot blade geometry: chord and twist angle vs. spanwise position.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(bladegeometry["r"], bladegeometry["chord"], label="Chord [m]")
    plt.plot(bladegeometry["r"], bladegeometry["twist_deg"], label="Twist [deg]")
    plt.xlabel("Spanwise position r [m]")
    plt.ylabel("Value")
    plt.title("Blade Geometry")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
   
plot_geometry(bladegeometry)   