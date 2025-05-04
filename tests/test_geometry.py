import sys
from pathlib import Path
import pytest
import pandas as pd
import matplotlib.pyplot as plt

# Add src/bem to the Python path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src" / "bem"))

from geometry import load_geometry, plot_geometry

# === TEST load_geometry === #

def test_load_geometry_with_mock_file(tmp_path):
    # Simulate file with 6 skipped header lines and one valid data line
    content = "\n" * 6 + "1 0 0 0 5 2.5 10 0 0 0\n"
    geom_file = tmp_path / "mock_blade.dat"
    geom_file.write_text(content)

    df = load_geometry(geom_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["r", "chord", "twist_deg", "airfoil_id"]
    assert df.iloc[0]["r"] == 1
    assert df.iloc[0]["chord"] == 2.5
    assert df.iloc[0]["twist_deg"] == 5
    assert df.iloc[0]["airfoil_id"] == 10


# === TEST plot_geometry === #

def test_plot_geometry(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)

    mock_df = pd.DataFrame({
        "r": [1, 2],
        "chord": [2.0, 2.5],
        "twist_deg": [5.0, 6.0],
        "airfoil_id": [1, 1]
    })
    plot_geometry(mock_df)  # Smoke test: just make sure it runs
