import pytest
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from bem.geometry import load_geometry, plot_geometry


@pytest.fixture
def sample_geometry_df():
    return pd.DataFrame({
        "r": [1.0, 2.0, 3.0],
        "chord": [2.5, 2.0, 1.5],
        "twist_deg": [10.0, 7.5, 5.0],
        "airfoil_id": [1, 2, 3]
    })


def test_load_geometry(tmp_path):
    # Create mock blade geometry file
    fake_data = (
        "\n" * 6 +
        "1.0 0 0 0 10.0 2.5 1 0 0 0\n"
        "2.0 0 0 0 7.5 2.0 2 0 0 0\n"
        "3.0 0 0 0 5.0 1.5 3 0 0 0\n"
    )
    test_file = tmp_path / "test_blade.dat"
    test_file.write_text(fake_data)

    df = load_geometry(test_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["r", "chord", "twist_deg", "airfoil_id"]
    assert df.shape == (3, 4)
    assert df["airfoil_id"].dtype == int


def test_plot_geometry(monkeypatch, sample_geometry_df):
    # Prevent plt.show() from blocking test
    monkeypatch.setattr(plt, "show", lambda: None)
    plot_geometry(sample_geometry_df)  # Smoke test: should not raise error
