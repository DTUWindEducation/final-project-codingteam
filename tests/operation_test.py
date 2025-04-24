import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd
from operation import blad_operation, plot_operational_subplots

def test_blad_operation():
    file_path = Path("./inputs/IEA-15-240-RWT/IEA_15MW_RWT_Onshore.opt")
    df = blad_operation(file_path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert list(df.columns) == ["V0", "pitch_deg", "omega_rpm", "power_kw", "thrust_kN"]

def test_plot_operational_subplots(monkeypatch):
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, "show", lambda: None)
    file_path = Path("./inputs/IEA-15-240-RWT/IEA_15MW_RWT_Onshore.opt")
    df = blad_operation(file_path)
    plot_operational_subplots(df)
