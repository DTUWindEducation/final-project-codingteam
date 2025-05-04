import pytest
import pandas as pd
from pathlib import Path
from bem.operation import blad_operation, plot_operational_subplots
import matplotlib.pyplot as plt

def test_blad_operation_with_mock_file(tmp_path):
    mock_data = "V0 pitch_deg omega_rpm power_kw thrust_kN\n10 2 12 1000 200"
    test_file = tmp_path / "mock.opt"
    test_file.write_text(mock_data)

    df = blad_operation(test_file)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert list(df.columns) == ["V0", "pitch_deg", "omega_rpm", "power_kw", "thrust_kN"]
    assert df.iloc[0]["V0"] == 10

def test_plot_operational_subplots(monkeypatch, tmp_path):
    monkeypatch.setattr(plt, "show", lambda: None)
    mock_data = "V0 pitch_deg omega_rpm power_kw thrust_kN\n10 2 12 1000 200"
    test_file = tmp_path / "mock.opt"
    test_file.write_text(mock_data)
    df = blad_operation(test_file)
    plot_operational_subplots(df)
