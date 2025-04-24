import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
operation = Path("./inputs/IEA-15-240-RWT/IEA_15MW_RWT_Onshore.opt")  # Operation = .txt

def blad_operation(operation):
    # Read the file, assuming whitespace-separated and no header
    df = pd.read_csv(operation, sep=r'\s+', header=None, skiprows=1)

    # Assign expected column names (based on your image)
    df.columns = [
        "V0",           # wind speed [m/s]
        "pitch_deg",    # pitch [deg]
        "omega_rpm",    # rotational speed [rpm]
        "power_kw",     # aero power [kW]
        "thrust_kN"     # aero thrust [kN]
    ]

    return df

oper_df = blad_operation(operation)
#print(oper_df)


def plot_operational_subplots(df):
    """
    Plot 2x2 subplots for turbine operation data:
    - Power curve
    - Pitch schedule
    - Thrust curve
    - Rotor speed schedule
    """
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Power curve
    axs[0, 0].plot(df["V0"], df["power_kw"], marker='o', color='blue')
    axs[0, 0].set_title("Power Curve")
    axs[0, 0].set_xlabel("Wind Speed [m/s]")
    axs[0, 0].set_ylabel("Power [kW]")
    axs[0, 0].grid(True)

    # Pitch angle schedule
    axs[0, 1].plot(df["V0"], df["pitch_deg"], marker='s', color='orange')
    axs[0, 1].set_title("Pitch Schedule")
    axs[0, 1].set_xlabel("Wind Speed [m/s]")
    axs[0, 1].set_ylabel("Pitch [deg]")
    axs[0, 1].grid(True)

    # Thrust curve
    axs[1, 0].plot(df["V0"], df["thrust_kN"], marker='^', color='green')
    axs[1, 0].set_title("Thrust Curve")
    axs[1, 0].set_xlabel("Wind Speed [m/s]")
    axs[1, 0].set_ylabel("Thrust [kN]")
    axs[1, 0].grid(True)

    # Rotor speed
    axs[1, 1].plot(df["V0"], df["omega_rpm"], marker='d', color='purple')
    axs[1, 1].set_title("Rotor Speed Schedule")
    axs[1, 1].set_xlabel("Wind Speed [m/s]")
    axs[1, 1].set_ylabel("Rotor Speed [rpm]")
    axs[1, 1].grid(True)

    plt.tight_layout()
    plt.show()


plot_operational_subplots(oper_df)
