# -*- coding: utf-8 -*-
"""Propeller Design using BEM + NeuralFoil (NACA 9112)
Spitfire Mk 24 geometry: R=1.587 m, 5 blades, 1240 RPM
Optimise: piecewise-linear N-point chord and twist distributions for max eta at 250 kts
"""

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import aerosandbox as asb
import neuralfoil as nf

# ============================================================
# 1. ISA Conditions at 10,000 ft
# ============================================================
alt_m = 10_000 * 0.3048
T = 288.15 - 0.0065 * alt_m
rho = 1.225 * (T / 288.15) ** 4.256
a_sound = np.sqrt(1.4 * 287 * T)
mu = 1.458e-6 * T**1.5 / (T + 110.4)
nu = mu / rho

# ============================================================
# 2. Propeller Parameters — Spitfire Mk 24
# ============================================================
R = 1.587   # radius [m]
R_root = 0.150   # hub radius [m]
N_b = 5   # blades
RPM = 1240   # propeller RPM
Omega = RPM * 2 * np.pi / 60
N_EL = 100   # spanwise elements
a_w = 0.61   # wake contraction ratio

V_DESIGN = 250 * 0.51444   # design speed: 250 kts [m/s]
N_CTRL = 7
GEOM_CTRL_STATIONS = np.linspace(0.0, 1.0, N_CTRL)  # active blade span: root -> tip

CHORD_GUESS_ROOT = 0.25
CHORD_GUESS_TIP = 0.12

CHORD_MIN_ROOT = 0.05
CHORD_MIN_TIP = 0.02
CHORD_MAX_ROOT = 0.30
CHORD_MAX_TIP = 0.10

# Twist is referenced to the design-point inflow angle instead of fixed
# absolute bounds. This avoids forcing the optimizer into strongly negative AoA
# at high advance ratio.
TWIST_ALPHA_GUESS_DEG = 4.0
TWIST_ALPHA_MIN_DEG = -2.0
TWIST_ALPHA_MAX_DEG = 12.0
TWIST_ABS_MIN_DEG = 5.0
TWIST_ABS_MAX_DEG = 95.0

ETA_PENALTY = 1e6

# ============================================================
# 3. NeuralFoil Polars for NACA 9112
# Pre-computed at 75% span (most aerodynamically loaded region)
# ============================================================
airfoil = asb.Airfoil("naca9112")
ALPHA_POLAR = np.linspace(-10, 25, 71)

lambda_c_design = V_DESIGN / (Omega * R)
r75 = 0.75 * R
V_loc_75 = Omega * R * np.sqrt((r75 / R) ** 2 + lambda_c_design**2)
c_ref = 0.20   # reference chord estimate for Re [m]
Re_ref = V_loc_75 * c_ref / nu
M_ref = V_loc_75 / a_sound

print("Reference polar conditions (75% span, 250 kts):")
print(f"  V_local = {V_loc_75:.1f} m/s  M = {M_ref:.3f}  Re = {Re_ref/1e6:.2f}M")
print("  Running NeuralFoil for NACA 9112...")

aero = nf.get_aero_from_airfoil(airfoil=airfoil, alpha=ALPHA_POLAR, Re=Re_ref, model_size="large")
CL_raw = np.array(aero["CL"])
CD_raw = np.array(aero["CD"])

CL_func = interp1d(ALPHA_POLAR, CL_raw, bounds_error=False, fill_value=(CL_raw[0], CL_raw[-1]))
CD_func = interp1d(ALPHA_POLAR, CD_raw, bounds_error=False, fill_value=(CD_raw[0], CD_raw[-1]))

print(f"  CL_max = {CL_raw.max():.3f}  at alpha = {ALPHA_POLAR[CL_raw.argmax()]:.1f} deg")
print(f"  Best L/D = {(CL_raw / CD_raw).max():.1f}  at alpha = {ALPHA_POLAR[(CL_raw / CD_raw).argmax()]:.1f} deg")
phi_ctrl_deg = np.degrees(np.arctan2(lambda_c_design, (R_root / R) + GEOM_CTRL_STATIONS * (1.0 - R_root / R)))
print(
    "  Design inflow angle range (no induced flow): "
    f"{phi_ctrl_deg[0]:.1f} deg at root to {phi_ctrl_deg[-1]:.1f} deg at tip"
)


# ============================================================
# 4. Geometry + BEM Solver (Full — with Drag + Prandtl-Glauert per element)
# ============================================================
def build_control_labels():
    labels = []
    for station in GEOM_CTRL_STATIONS:
        if np.isclose(station, 0.0):
            labels.append("Root")
        elif np.isclose(station, 1.0):
            labels.append("Tip")
        else:
            labels.append(f"{int(round(station * 100))}% span")
    return labels


GEOM_CTRL_LABELS = build_control_labels()


def unpack_design_variables(x):
    n_ctrl = len(GEOM_CTRL_STATIONS)
    if len(x) != 2 * n_ctrl:
        raise ValueError(f"Expected {2 * n_ctrl} design variables, got {len(x)}")
    chord_ctrl = np.asarray(x[:n_ctrl], dtype=float)
    twist_ctrl_deg = np.asarray(x[n_ctrl:], dtype=float)
    return chord_ctrl, twist_ctrl_deg


def build_blade_geometry(chord_ctrl, twist_ctrl_deg, N):
    r_nd = R_root / R
    dr = (1.0 - r_nd) / N
    r_arr = np.linspace(r_nd + dr / 2, 1.0 - dr / 2, N)

    # Piecewise-linear interpolation between the optimized control points.
    span_frac = np.linspace(0.5 / N, 1.0 - 0.5 / N, N)
    c_arr = np.interp(span_frac, GEOM_CTRL_STATIONS, chord_ctrl)
    theta_arr_deg = np.interp(span_frac, GEOM_CTRL_STATIONS, twist_ctrl_deg)
    theta_arr_rad = np.deg2rad(theta_arr_deg)
    return r_arr, dr, c_arr, theta_arr_deg, theta_arr_rad


def solve_inflow(r, sigma, theta_rad, lambda_c):
    def residual(lambda_i):
        lam = lambda_c + lambda_i
        phi = np.arctan(lam / r)
        alpha_deg = np.degrees(theta_rad - phi)

        # Prandtl-Glauert correction for local Mach (CL only)
        V_loc = Omega * R * np.sqrt(r**2 + lam**2)
        M_loc = V_loc / a_sound
        pg = 1.0 / np.sqrt(max(1.0 - M_loc**2, 0.05))

        CL = float(CL_func(alpha_deg)) * pg
        CD = float(CD_func(alpha_deg))

        CT_mom = 4.0 * lambda_i * lam * r
        CT_blade = 0.5 * sigma * (r**2 + lam**2) * (CL * np.cos(phi) - CD * np.sin(phi))
        return CT_mom - CT_blade

    try:
        lambda_i = brentq(residual, 0.0, 5.0, maxiter=100, xtol=1e-8)
        return lambda_c + lambda_i
    except Exception:
        return lambda_c


def design_inflow_angles_deg():
    r_root_nd = R_root / R
    r_ctrl_nd = r_root_nd + GEOM_CTRL_STATIONS * (1.0 - r_root_nd)
    return np.degrees(np.arctan2(lambda_c_design, r_ctrl_nd))


def build_initial_guess():
    chord_guess = np.linspace(CHORD_GUESS_ROOT, CHORD_GUESS_TIP, N_CTRL)
    twist_guess = np.clip(
        design_inflow_angles_deg() + TWIST_ALPHA_GUESS_DEG,
        TWIST_ABS_MIN_DEG,
        TWIST_ABS_MAX_DEG,
    )
    return np.concatenate([chord_guess, twist_guess])


def build_bounds():
    chord_lower = np.linspace(CHORD_MIN_ROOT, CHORD_MIN_TIP, N_CTRL)
    chord_upper = np.linspace(CHORD_MAX_ROOT, CHORD_MAX_TIP, N_CTRL)
    phi_ctrl_deg = design_inflow_angles_deg()
    twist_lower = np.clip(
        phi_ctrl_deg + TWIST_ALPHA_MIN_DEG,
        TWIST_ABS_MIN_DEG,
        TWIST_ABS_MAX_DEG,
    )
    twist_upper = np.clip(
        phi_ctrl_deg + TWIST_ALPHA_MAX_DEG,
        TWIST_ABS_MIN_DEG,
        TWIST_ABS_MAX_DEG,
    )

    chord_bounds = list(zip(chord_lower, chord_upper))
    twist_bounds = list(zip(twist_lower, twist_upper))
    return chord_bounds + twist_bounds


def build_monotonic_constraints():
    constraints = []
    n_ctrl = len(GEOM_CTRL_STATIONS)

    for i in range(n_ctrl - 1):
        constraints.append({"type": "ineq", "fun": lambda x, i=i: x[i] - x[i + 1]})

    offset = n_ctrl
    for i in range(n_ctrl - 1):
        constraints.append(
            {"type": "ineq", "fun": lambda x, i=i, offset=offset: x[offset + i] - x[offset + i + 1]}
        )

    return constraints


def bem_solver(design_vars, V_c, N=N_EL):
    """
    Full BEM with drag.
    Inflow ratio lambda solved iteratively via brentq at each element.
    Returns CT, CQ, propulsive efficiency eta.
    """
    lambda_c = V_c / (Omega * R)
    chord_ctrl, twist_ctrl_deg = unpack_design_variables(design_vars)
    r_arr, dr, c_arr, _, theta_arr = build_blade_geometry(chord_ctrl, twist_ctrl_deg, N)

    total_CT = 0.0
    total_CQ = 0.0

    for k in range(N):
        r = r_arr[k]
        c = c_arr[k]
        theta = theta_arr[k]
        sigma = N_b * c / (np.pi * R)

        lam = solve_inflow(r, sigma, theta, lambda_c)

        phi = np.arctan(lam / r)
        alpha_deg = np.degrees(theta - phi)
        V_loc = Omega * R * np.sqrt(r**2 + lam**2)
        M_loc = V_loc / a_sound
        pg = 1.0 / np.sqrt(max(1.0 - M_loc**2, 0.05))

        CL = float(CL_func(alpha_deg)) * pg
        CD = float(CD_func(alpha_deg))

        dCT = 0.5 * sigma * (r**2 + lam**2) * (CL * np.cos(phi) - CD * np.sin(phi)) * dr
        dCQ = 0.5 * sigma * (r**2 + lam**2) * (CL * np.sin(phi) + CD * np.cos(phi)) * r * dr

        total_CT += dCT
        total_CQ += dCQ

    eta = (total_CT * lambda_c / total_CQ) if (total_CT > 0.0 and total_CQ > 1e-9) else np.nan
    return total_CT, total_CQ, eta


# ============================================================
# 5. Optimisation — Max eta at 250 kts
# ============================================================
def objective(x):
    CT, CQ, eta = bem_solver(x, V_DESIGN)
    if not np.isfinite(eta):
        penalty = ETA_PENALTY
        if np.isfinite(CT) and CT < 0.0:
            penalty += 1e4 * abs(CT)
        if np.isfinite(CQ) and CQ <= 0.0:
            penalty += 1e4 * abs(CQ)
        return penalty
    return -eta


print(
    f"\nOptimising piecewise-linear propeller geometry for eta with {N_CTRL} control points "
    f"per distribution (design point: 250 kts)..."
)

result = minimize(
    objective,
    x0=build_initial_guess(),
    method="SLSQP",
    bounds=build_bounds(),
    constraints=build_monotonic_constraints(),
)

if not result.success:
    print(f"WARNING: Optimizer did not fully converge: {result.message}")

opt_design = result.x
opt_chord_ctrl, opt_twist_ctrl = unpack_design_variables(opt_design)
opt_CT, opt_CQ, opt_eta = bem_solver(opt_design, V_DESIGN)

print(f"\n{'='*52}")
print("OPTIMAL PROPELLER GEOMETRY  (NACA 9112, MAX ETA)")
print(f"{'='*52}")
print(f"  Radius      : {R:.3f} m")
print(f"  Blades      : {N_b}")
print(f"  RPM         : {RPM}")
print(f"  Ctrl pts    : {N_CTRL}")
print("  Variation   : piecewise-linear between adjacent control points")
print(f"  Design CT   : {opt_CT:.5f}")
print(f"  Design CQ   : {opt_CQ:.5f}")
print(f"  Design eta  : {opt_eta:.4f}")
for label, chord, twist in zip(GEOM_CTRL_LABELS, opt_chord_ctrl, opt_twist_ctrl):
    print(f"  {label:8}: chord = {chord*100:.1f} cm   twist = {twist:.1f} deg")

# ============================================================
# 6. Performance Across Speed Range
# ============================================================
print(f"\n{'='*65}")
print("PERFORMANCE ACROSS SPEED RANGE")
print(f"{'='*65}")
print(f"{'Speed':>8} {'CT':>9} {'CQ':>9} {'eta':>7} {'Thrust [N]':>12} {'Power [kW]':>11}")
print("-" * 65)

for kts in [200, 250, 300]:
    V = kts * 0.51444
    CT, CQ, eta = bem_solver(opt_design, V)
    T_force = CT * rho * np.pi * R**2 * (Omega * R)**2
    P = CQ * rho * np.pi * R**2 * (Omega * R)**2 * R * Omega / 1000
    print(f"{kts:>7} kts  {CT:>9.5f}  {CQ:>9.5f}  {eta:>7.3f}  {T_force:>12.1f}  {P:>11.1f}")

# ============================================================
# 7. Spanwise Distributions for Plotting
# ============================================================
N_plot = 100
lambda_c = V_DESIGN / (Omega * R)
r_arr, _, c_arr, th_arr, th_arr_rad = build_blade_geometry(opt_chord_ctrl, opt_twist_ctrl, N_plot)
ctrl_r_phys = R_root + GEOM_CTRL_STATIONS * (R - R_root)

alpha_dist, phi_dist = [], []
for k in range(N_plot):
    r = r_arr[k]
    sigma = N_b * c_arr[k] / (np.pi * R)
    theta = th_arr_rad[k]
    lam = solve_inflow(r, sigma, theta, lambda_c)
    phi = np.arctan(lam / r)
    alpha = np.degrees(theta - phi)
    alpha_dist.append(alpha)
    phi_dist.append(np.degrees(phi))

r_phys = r_arr * R   # physical radius [m]

# ============================================================
# 8. Plots
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Chord distribution
ax = axes[0, 0]
ax.plot(r_phys, np.array(c_arr) * 100, color="steelblue", linewidth=2)
ax.scatter(ctrl_r_phys, opt_chord_ctrl * 100, color="navy", s=30, zorder=3)
ax.fill_between(r_phys, np.array(c_arr) * 100, alpha=0.2, color="steelblue")
ax.set_xlabel("Radius [m]")
ax.set_ylabel("Chord [cm]")
ax.set_title("Chord Distribution")
ax.grid(True, alpha=0.3)

# Twist distribution
ax = axes[0, 1]
ax.plot(r_phys, th_arr, color="tomato", linewidth=2)
ax.scatter(ctrl_r_phys, opt_twist_ctrl, color="darkred", s=30, zorder=3)
ax.set_xlabel("Radius [m]")
ax.set_ylabel("Twist angle [deg]")
ax.set_title("Twist Distribution")
ax.grid(True, alpha=0.3)

# Inflow angle
ax = axes[1, 0]
ax.plot(r_phys, phi_dist, color="seagreen", linewidth=2, label="Inflow angle phi")
ax.plot(r_phys, th_arr, color="tomato", linewidth=1.5, linestyle="--", label="Blade pitch theta")
ax.set_xlabel("Radius [m]")
ax.set_ylabel("Angle [deg]")
ax.set_title("Inflow vs Pitch Angle")
ax.legend()
ax.grid(True, alpha=0.3)

# Effective AoA
ax = axes[1, 1]
ax.plot(r_phys, alpha_dist, color="purple", linewidth=2)
ax.axhline(
    ALPHA_POLAR[CL_raw.argmax()],
    color="tomato",
    linestyle="--",
    linewidth=1.5,
    label=f"Stall alpha = {ALPHA_POLAR[CL_raw.argmax()]:.1f} deg",
)
ax.set_xlabel("Radius [m]")
ax.set_ylabel("Effective AoA [deg]")
ax.set_title("Effective Angle of Attack along Blade")
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle(
    f"NACA 9112 Propeller — R={R}m, {N_b} blades, {RPM} RPM\n"
    f"{N_CTRL}-point piecewise-linear chord/twist control (max eta)   "
    f"Root->Tip chord {opt_chord_ctrl[0]*100:.1f}->{opt_chord_ctrl[-1]*100:.1f} cm   "
    f"Twist {opt_twist_ctrl[0]:.1f}->{opt_twist_ctrl[-1]:.1f} deg",
    fontsize=11,
    fontweight="bold",
)
plt.tight_layout()
plt.savefig("propeller_design_eta.png", dpi=150, bbox_inches="tight")
print("\nPlot saved: propeller_design_eta.png")
