import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class BEM:
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

