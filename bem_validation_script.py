

import numpy as np

def calculate_propeller_thrust(theta_deg_start,theta_deg_end, R, R_root, c0,c1, RPM, V_c, N_b, a_w, cla0, alpha_l0_deg, N_elements=20):
    """
    Calculates the thrust coefficient (C_T) for a propeller using BEM theory.
    """
    # 1. Convert units and define globals
    Omega = RPM * (2 * np.pi / 60)       # Rotational speed [rad/s]       # Twist angle [rad]
    alpha_l0 = np.radians(alpha_l0_deg)  # Zero-lift angle of attack [rad]
    #cla = np.deg2rad(cla0)
    cla = cla0 * (180 / np.pi)


    # Climb velocity coefficient
    lambda_c = V_c / (Omega * R)

    # 2. Discretize the blade into N elements
    r_root_nd = R_root / R               # the wind start by %
    dr = (1.0 - r_root_nd) / N_elements  # each element %

    #  array each element
    r_array = np.linspace(r_root_nd + dr/2, 1.0 - dr/2, N_elements)

    c_array= np.linspace(c0,c1,N_elements) # un chord

    Teta_arr= np.linspace(np.deg2rad(theta_deg_start),np.deg2rad(theta_deg_end)  ,N_elements) #un angle


    # 3. Calculate solidity (sigma)
    # Since chord is constant in this validation case, sigma is constant


    total_CT = 0.0 # Initialize total thrust coefficient

    # 4. The Core Loop: Solve for each element
    for k in range(N_elements):
        r= r_array[k]
        c= c_array[k]
        sigma = (N_b * c) / (np.pi * R)

        theta= Teta_arr[k]



        A = (2 * r) / a_w
        B = (0.5 * sigma * cla * r) - (2 * r * lambda_c)
        C_const = -0.5 * sigma * cla * (r**2) * (theta + alpha_l0)

        # x = (-B + sqrt(B^2 - 4AC)) / 2A
        discriminant = B**2 - 4 * A * C_const

        if discriminant >= 0:
            lambda_val = (-B + np.sqrt(discriminant)) / (2 * A)
        else:
            lambda_val = 0

        # Calculate elemental thrust coefficient dC_T
        dCT = 0.5 * sigma * cla * (theta * r**2 - lambda_val * r + alpha_l0 * r**2) * dr


        total_CT += dCT

    return total_CT

if __name__ == "__main__":
    # --- TABLE 1: Global Parameters ---
    R = 1.829           # Radius [m]
    R_root = 0.178      # Root cut-off [m]
    c = 0.1524          # Chord [m]
    RPM = 600           # Rotational speed [RPM]
    V_c = 0.0           # Climb velocity [m/s]
    N_b = 3           # Number of blades
    a_w = 0.61            # Wake contraction ratio


    cla = 0.1096         #in Dgree  # Lift-curve slope [1/rad]
    alpha_l0 = 0.0     # Zero-lift angle of attack [deg]



    for j in [1,8,15]:
        print("Ct:",calculate_propeller_thrust(j,j, R, R_root, c,c, RPM, V_c, N_b, a_w, cla, alpha_l0))

