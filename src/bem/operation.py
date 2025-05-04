"""Functions for loading and plotting wind turbine operational data."""

import pandas as pd
import matplotlib.pyplot as plt


def blad_operation(filepath):
    """
    Load wind turbine operation data from a whitespace-separated file.

    Args:
        filepath (str): Path to the data file.

    Returns:
        pd.DataFrame: DataFrame with columns for wind speed, pitch, 
                      rotational speed, power, and thrust.
    """
    df_operation = pd.read_csv(filepath, sep=r'\s+', header=None, skiprows=1)
    df_operation.columns = [
        "V0",           # wind speed [m/s]
        "pitch_deg",    # pitch [deg]
        "omega_rpm",    # rotational speed [rpm]
        "power_kw",     # aero power [kW]
        "thrust_kN"     # aero thrust [kN]
    ]
    return df_operation


def plot_operational_subplots(df_operation):
    """
    Plot 2x2 subplots for turbine operation data.
    """
    _, axs = plt.subplots(2, 2, figsize=(12, 8))

    axs[0, 0].plot(df_operation["V0"], df_operation["power_kw"], marker='o', color='blue')
    axs[0, 0].set_title("Power Curve")
    axs[0, 0].set_xlabel("Wind Speed [m/s]")
    axs[0, 0].set_ylabel("Power [kW]")
    axs[0, 0].grid(True)

    axs[0, 1].plot(df_operation["V0"], df_operation["pitch_deg"], marker='s', color='orange')
    axs[0, 1].set_title("Pitch Schedule")
    axs[0, 1].set_xlabel("Wind Speed [m/s]")
    axs[0, 1].set_ylabel("Pitch [deg]")
    axs[0, 1].grid(True)

    axs[1, 0].plot(df_operation["V0"], df_operation["thrust_kN"], marker='^', color='green')
    axs[1, 0].set_title("Thrust Curve")
    axs[1, 0].set_xlabel("Wind Speed [m/s]")
    axs[1, 0].set_ylabel("Thrust [kN]")
    axs[1, 0].grid(True)

    axs[1, 1].plot(df_operation["V0"], df_operation["omega_rpm"], marker='d', color='purple')
    axs[1, 1].set_title("Rotor Speed Schedule")
    axs[1, 1].set_xlabel("Wind Speed [m/s]")
    axs[1, 1].set_ylabel("Rotor Speed [rpm]")
    axs[1, 1].grid(True)

    plt.tight_layout()
    plt.show()