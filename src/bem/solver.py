"""
Solver module for Blade Element Momentum (BEM) method for wind turbine analysis.
"""
# pylint: disable=C0103, R0914
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class BEMSolver:
    """
    Blade Element Momentum (BEM) method for wind turbine analysis.

    Attributes:
        R (float): Radius of the rotor (m).
        omega (float): Angular velocity of the rotor (rad/s).
        theta_p (float): Pitch angle of the blade (rad).
        beta (float): Angle of attack (rad).
        rho (float): Air density (kg/m^3).
        v0 (float): Wind speed (m/s).
        B (int): Number of blades.
        c (float): Chord length of the blade (m).
    """

    def __init__(self, blade_geometry, airfoil_polars, n_blades=3, rho=1.225):
        self.geometry = blade_geometry
        self.polars = airfoil_polars
        self.B = n_blades
        self.rho = rho


    def interpolate_polar(self, airfoil_id, alpha_deg):
        '''Interpolate the lift and drag coefficients for a given angle of attack.
        Args:
            airfoil_id (int): The ID of the airfoil.
            alpha_deg (float): Angle of attack in degrees.'''
        alpha, cl, cd = self.polars[airfoil_id]
        cl_interp = interp1d(alpha, cl, bounds_error=False, fill_value="extrapolate")
        cd_interp = interp1d(alpha, cd, bounds_error=False, fill_value="extrapolate")
        # return cl_interp(alpha_deg), cd_interp(alpha_deg)
        return cl_interp(alpha_deg).item(), cd_interp(alpha_deg).item()


    def compute_induction_factors(self, sigma, phi, Cn, Ct):
        '''Compute the induction factors a and a' based on the BEM theory.
        Args:
            sigma (float): Solidity of the rotor.
            phi (float): Angle between the wind and the rotor plane.
            Cn (float): Normal force coefficient.
            Ct (float): Tangential force coefficient.'''
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)

        a = 1 / (4 * sin_phi**2 / (sigma * Cn) + 1)
        a_prime = 1 / (4 * sin_phi * cos_phi / (sigma * Ct) - 1)
        return a, a_prime


    def compute_performance(self, V0, omega, pitch_deg):
        '''Compute the thrust, torque, and power of the rotor. '''
        dr = np.gradient(self.geometry["r"])
        r_vals = self.geometry["r"].values
        chord = self.geometry["chord"].values
        twist = self.geometry["twist_deg"].values
        af_id = self.geometry["airfoil_id"].astype(int).values

        a = np.zeros_like(r_vals)
        a_prime = np.zeros_like(r_vals)
        max_iter = 100
        tol = 1e-5

        dT = np.zeros_like(r_vals)
        dM = np.zeros_like(r_vals)

        for i, r in enumerate(r_vals):
            if r == 0:
                continue 

            local_chord = chord[i]
            local_twist = np.radians(twist[i])
            airfoil = af_id[i] - 1  # convert to 0-based
            sigma = self.B * local_chord / (2 * np.pi * r)
            # Init
            a_i, a_p_i = 0.0, 0.0
            for _ in range(max_iter):
                phi = np.arctan2((1 - a_i) * V0, (1 + a_p_i) * omega * r)
                alpha = np.degrees(phi - (np.radians(pitch_deg) + local_twist))
                cl, cd = self.interpolate_polar(airfoil, alpha)
                cn = cl * np.cos(phi) + cd * np.sin(phi)
                ct = cl * np.sin(phi) - cd * np.cos(phi)
                a_new, a_p_new = self.compute_induction_factors(sigma, phi, cn, ct)
                if np.abs(a_new - a_i) < tol and np.abs(a_p_new - a_p_i) < tol:
                    break
                a_i, a_p_i = a_new, a_p_new
            a[i] = a_i
            a_prime[i] = a_p_i

            # Local forces
            dT[i] = 4 * np.pi * r * self.rho * V0**2 * a_i * (1 - a_i) #* dr[i]
            dM[i] = 4 * np.pi * r**3 * self.rho * V0 * omega * a_p_i * (1 - a_i) #* dr[i]
        T = np.trapezoid(dT, r_vals)
        M = np.trapezoid(dM, r_vals)

        P = M * omega
        return T, M, P


    def compute_cp_ct_surface(self, pitch_range, lambda_range, V0, R, n_points=20):
        """
        Compute CP and CT over a grid of pitch angles and tip speed ratios (lambda).
    
        Args:
        pitch_range (tuple): (min_pitch_deg, max_pitch_deg)
        lambda_range (tuple): (min_lambda, max_lambda)
        V0 (float): Wind speed [m/s]
        R (float): Rotor radius [m]
        n_points (int): Number of points in each dimension

        Returns:
        Tuple of (pitch_grid, lambda_grid, cp_grid, ct_grid)
        """
        pitch_vals = np.linspace(*pitch_range, n_points)
        lambda_vals = np.linspace(*lambda_range, n_points)
        pitch_grid, lambda_grid = np.meshgrid(pitch_vals, lambda_vals)

        cp_grid = np.zeros_like(pitch_grid)
        ct_grid = np.zeros_like(lambda_grid)
        A = np.pi * R**2

        for i in range(n_points):
            for j in range(n_points):
                pitch = pitch_grid[i, j]
                lam = lambda_grid[i, j]
                omega = lam * V0 / R
                T, _, P = self.compute_performance(V0, omega, pitch)
                cp_grid[i, j] = P / (0.5 * self.rho * A * V0**3)
                ct_grid[i, j] = T  / (0.5 * self.rho * A * V0**2)
        # Max CP
        max_cp_idx = np.unravel_index(np.argmax(cp_grid), cp_grid.shape)
        max_cp = cp_grid[max_cp_idx]
        max_cp_pitch = pitch_grid[max_cp_idx]
        max_cp_lambda = lambda_grid[max_cp_idx]

        # Max CT
        max_ct_idx = np.unravel_index(np.argmax(ct_grid), ct_grid.shape)
        max_ct = ct_grid[max_ct_idx]
        max_ct_pitch = pitch_grid[max_ct_idx]
        max_ct_lambda = lambda_grid[max_ct_idx]
        max_cp_point = (max_cp_pitch, max_cp_lambda)
        max_ct_point = (max_ct_pitch, max_ct_lambda)

        print(f"Max CP = {max_cp:.3f} at pitch = {max_cp_pitch:.2f} deg, TSR = {max_cp_lambda:.2f}")
        print(f"Max CT = {max_ct:.3f} at pitch = {max_ct_pitch:.2f} deg, TSR = {max_ct_lambda:.2f}")

        return pitch_grid, lambda_grid, cp_grid, ct_grid, max_cp_point, max_ct_point


def plot_cp_ct_contours(pitch_grid, lambda_grid, cp_grid, ct_grid, 
                        design_point=None, max_cp_point=None, max_ct_point=None):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    cp_levels = np.linspace(np.min(cp_grid), np.max(cp_grid), 10)
    ct_levels = np.linspace(np.min(ct_grid), np.max(ct_grid), 10)

    cp_contour = axs[0].contour(pitch_grid, lambda_grid, cp_grid, levels=cp_levels, cmap='viridis')
    axs[0].clabel(cp_contour, inline=True, fontsize=8)
    axs[0].set_title("Power Coefficient")
    axs[0].set_xlabel("Blade Pitch Angle, [°]")
    axs[0].set_ylabel("Tip Speed Ratio, [λ]")

    ct_contour = axs[1].contour(pitch_grid, lambda_grid, ct_grid, levels=ct_levels, cmap='viridis')
    axs[1].clabel(ct_contour, inline=True, fontsize=8)
    axs[1].set_title("Thrust Coefficient")
    axs[1].set_xlabel("Blade Pitch Angle, [°]")
    axs[1].set_ylabel("Tip Speed Ratio, [λ]")

    if design_point:
        for ax in axs:
            ax.plot(design_point[0], design_point[1], 'ks', label="Design Point")

    if max_cp_point:
        axs[0].plot(max_cp_point[0], max_cp_point[1], 'ro', label="Max $C_P$")
        axs[0].legend()

    if max_ct_point:
        axs[1].plot(max_ct_point[0], max_ct_point[1], 'bo', label="Max $C_T$")
        axs[1].legend()

    plt.tight_layout()
    plt.show()


def plot_bem_vs_measured(operation_df, power_output, thrust_output):
    """
    Create subplots comparing predicted vs measured power and thrust curves.

    Args:
        operation_df (pd.DataFrame): DataFrame with measured data.
        power_output (list or array): Predicted power [kW].
        thrust_output (list or array): Predicted thrust [kN].
    """
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    axs[0].plot(operation_df["V0"], power_output, label="Predicted Power [kW]")
    axs[0].plot(operation_df["V0"], operation_df["power_kw"], label="Measured Power [kW]")
    axs[0].set_ylabel("Power [kW]")
    axs[0].set_title("BEM vs Measured: Power")
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(operation_df["V0"], thrust_output, label="Predicted Thrust [kN]")
    axs[1].plot(operation_df["V0"], operation_df["thrust_kN"], label="Measured Thrust [kN]")
    axs[1].set_xlabel("Wind Speed [m/s]")
    axs[1].set_ylabel("Thrust [kN]")
    axs[1].set_title("BEM vs Measured: Thrust")
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()

def plot_cp_ct(V0_array, P_array, T_array, R, rho=1.225):
    """
    Compute and plot power and thrust coefficients.     
    """
    A = np.pi * R**2
    P_array = np.array(P_array) * 1000  # kW to W
    T_array = np.array(T_array) * 1000  # kN to N

    cp = P_array / (0.5 * rho * A * V0_array**3)
    ct = T_array / (0.5 * rho * A * V0_array**2)
    plt.figure()
    plt.plot(V0_array, cp, label="$C_P$ (Power Coefficient)")
    plt.plot(V0_array, ct, label="$C_T$ (Thrust Coefficient)")
    plt.xlabel("Wind Speed [m/s]")
    plt.ylabel("Coefficient Value")
    plt.title("Power and Thrust Coefficients")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
