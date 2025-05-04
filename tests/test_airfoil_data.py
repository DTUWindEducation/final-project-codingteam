import sys
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add src/bem to Python path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))

from airfoil_data import read_coords, read_polar, plot_airfoils

# === TEST read_coords === #

def test_read_coords_with_mock_file(tmp_path):
    coord_content = "\n" * 8 + "0.0 0.0\n1.0 1.0\n"
    coord_file = tmp_path / "IEA-15-240_AF00_coords.txt"
    coord_file.write_text(coord_content)

    result = read_coords([coord_file])
    assert "IEA-15-240_AF00_coords.txt" in result
    df = result["IEA-15-240_AF00_coords.txt"]
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (2, 2)


# === TEST read_polar === #

def test_read_polar_with_mock_file(tmp_path):
    polar_content = (
        "HEADER\n"
        "NumAlf 2\n\n"
        "0.0 1.0 0.01\n"
        "5.0 1.2 0.02\n"
    )
    polar_file = tmp_path / "IEA-15-240-RWT_AeroDyn15_Polar_00.dat"
    polar_file.write_text(polar_content)

    result = read_polar([polar_file])
    assert 0 in result
    alpha, cl, cd = result[0]
    assert isinstance(alpha, np.ndarray)
    assert cl[0] == 1.0
    assert cd[1] == 0.02


# === TEST plot_airfoils === #

def test_plot_airfoils(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)

    mock_df = pd.DataFrame({0: [0.0, 1.0], 1: [0.0, 1.0]})
    coords_dict = {"IEA-15-240_AF00_coords.txt": mock_df}
    plot_airfoils(coords_dict)  # Smoke test
