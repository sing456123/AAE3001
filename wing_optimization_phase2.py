# -*- coding: utf-8 -*-
"""Phase 2: 3D Wing Optimisation using Prandtl Lifting Line Theory
Fixed  : NACA 9112, span=0.50 m, root chord=0.12 m, no sweep
Optimise: taper ratio (0.3–1.0), tip twist (−5° to 0° washout)
Design point: 250 kts, 10,000 ft
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import aerosandbox as asb
import neuralfoil as nf

# ============================================================
# 1. ISA Conditions at 10,000 ft
# ============================================================
alt_m   = 10_000 * 0.3048
T       = 288.15 - 0.0065 * alt_m
a_sound = np.sqrt(1.4 * 287 * T)
mu      = 1.458e-6 * T**1.5 / (T + 110.4)
rho     = 1.225 * (T / 288.15) ** 4.256
nu      = mu / rho

# ============================================================
# 2. Fixed Wing Parameters
# ============================================================
NACA    = "naca9112"
SPAN    = 0.50    # m
C_ROOT  = 0.12    # m
N_LLT   = 40      # number of spanwise stations

# ============================================================
# 3. Extract 2D Airfoil Data from NeuralFoil + Prandtl-Glauert
# ============================================================
airfoil = asb.Airfoil(NACA)
ALPHA_2D = np.linspace(-15, 20, 71)

def get_2d_data(speed_kts):
    V    = speed_kts * 0.51444
    mach = V / a_sound
    Re   = V * C_ROOT / nu
    pg   = 1.0 / np.sqrt(1.0 - mach**2)

    aero  = nf.get_aero_from_airfoil(airfoil=airfoil, alpha=ALPHA_2D,
                                      Re=Re, model_size="large")
    cl_pg = np.array(aero["CL"]) * pg

    # Fit linear region for lift slope and zero-lift angle
    mask     = (ALPHA_2D > -8) & (ALPHA_2D < 6)
    coeffs   = np.polyfit(ALPHA_2D[mask], cl_pg[mask], 1)
    a_2d     = float(coeffs[0]) * (180.0 / np.pi)   # [1/rad]
    alpha_L0 = np.deg2rad(-coeffs[1] / coeffs[0])   # [rad]
    cl_max   = float(np.max(cl_pg))

    return {"kts": speed_kts, "mach": mach, "Re": Re,
            "a_2d": a_2d, "alpha_L0": alpha_L0, "cl_max_2d": cl_max}

print("Extracting NACA 9112 2D data from NeuralFoil...")
SPEEDS_KTS = [200, 250, 300]
cond = {kts: get_2d_data(kts) for kts in SPEEDS_KTS}

for kts, d in cond.items():
    print(f"  {kts} kts  M={d['mach']:.3f}  Re={d['Re']/1e6:.2f}M  "
          f"a_2d={d['a_2d']:.3f} /rad  "
          f"α_L0={np.degrees(d['alpha_L0']):.2f}°  "
          f"CL_max_2D={d['cl_max_2d']:.4f}")

# ============================================================
# 4. Prandtl Lifting Line Theory Solver
# ============================================================
def llt_cl_max(taper, tip_twist_deg, a_2d, alpha_L0, cl_max_2d,
               b=SPAN, c_root=C_ROOT, N=N_LLT):
    """
    Solves Prandtl LLT for 3D wing CL_max.
    Returns the wing CL when the first spanwise station reaches the 2D CL_max.

    Lift distribution: Γ(θ) = 2bV Σ Aₙ sin(nθ),  y = -b/2 cos(θ)
    Wing CL = π AR A₁
    """
    AR      = 2*b / (c_root * (1 + taper))
    theta   = np.linspace(np.pi/(N+1), N*np.pi/(N+1), N)
    n_terms = np.arange(1, N+1)

    # Local chord — linear taper; |2y/b| = |cos θ|
    c_loc = c_root * (1 - (1 - taper) * np.abs(np.cos(theta)))

    # Twist distribution — linear from 0 at root to tip_twist at tip
    twist = np.deg2rad(tip_twist_deg) * np.abs(np.cos(theta))

    # LLT influence matrix
    # Σ Aₙ sin(nθᵢ)[4b/(cᵢ a₀) + n/sin θᵢ] = α_eff(θᵢ)   [radians]
    M = np.array([
        [np.sin(n * theta[i]) * (4*b / (c_loc[i] * a_2d) + n / np.sin(theta[i]))
         for n in n_terms]
        for i in range(N)
    ])

    # Solve separately for unit AoA and twist contributions
    A_a = np.linalg.solve(M, np.ones(N))    # per 1 rad of (α_root − α_L0)
    A_t = np.linalg.solve(M, twist)         # twist contribution

    # Local section CL = 4b/c Σ Aₙ sin(nθ)
    local_cl_a = np.array([4*b/c_loc[i] * np.dot(A_a, np.sin(n_terms*theta[i])) for i in range(N)])
    local_cl_t = np.array([4*b/c_loc[i] * np.dot(A_t, np.sin(n_terms*theta[i])) for i in range(N)])

    # AoA at which each station stalls: local_cl_a·α + local_cl_t = cl_max_2d
    with np.errstate(divide='ignore', invalid='ignore'):
        alpha_stall = np.where(local_cl_a > 1e-6,
                               (cl_max_2d - local_cl_t) / local_cl_a,
                               np.inf)

    valid = alpha_stall[alpha_stall > 0]
    if len(valid) == 0:
        return 0.0, 0.0, None, None

    alpha_root_stall = float(np.min(valid))
    A1_stall  = alpha_root_stall * A_a[0] + A_t[0]
    CL_wing   = float(np.pi * AR * A1_stall)

    local_cl  = alpha_root_stall * local_cl_a + local_cl_t  # for plotting
    y_span    = -b/2 * np.cos(theta)

    return max(0.0, CL_wing), AR, local_cl, y_span

# ============================================================
# 5. Optimisation at 250 kts Design Point
# ============================================================
d250 = cond[250]

def objective(x):
    taper, tip_twist = x
    cl, *_ = llt_cl_max(taper, tip_twist,
                         d250["a_2d"], d250["alpha_L0"], d250["cl_max_2d"])
    return -cl

print("\nOptimising wing geometry (design point: 250 kts)...")
result = minimize(objective, x0=[0.6, -2.0], method="SLSQP",
                  bounds=[(0.3, 1.0), (-5.0, 0.0)])

opt_taper, opt_twist = result.x
opt_cl, opt_AR, local_cl, y_span = llt_cl_max(
    opt_taper, opt_twist, d250["a_2d"], d250["alpha_L0"], d250["cl_max_2d"])
opt_c_tip = C_ROOT * opt_taper

print(f"\n{'='*50}")
print("OPTIMAL 3D WING GEOMETRY  (NACA 9112)")
print(f"{'='*50}")
print(f"  Span          : {SPAN*100:.0f} cm")
print(f"  Root chord    : {C_ROOT*100:.0f} cm")
print(f"  Tip chord     : {opt_c_tip*100:.1f} cm")
print(f"  Taper ratio   : {opt_taper:.3f}")
print(f"  Tip twist     : {opt_twist:.2f}° (washout)")
print(f"  Aspect ratio  : {opt_AR:.2f}")

# ============================================================
# 6. Evaluate CL_wing Across All Speeds
# ============================================================
print(f"\n{'='*50}")
print("3D WING CL_max ACROSS SPEED RANGE")
print(f"{'='*50}")
print(f"{'Speed':>8} {'Mach':>6} {'Re':>8} {'2D CL_max':>10} {'Wing CL_max':>12}")
print("-"*50)
for kts in SPEEDS_KTS:
    d  = cond[kts]
    cl_w, ar_w, *_ = llt_cl_max(opt_taper, opt_twist,
                                  d["a_2d"], d["alpha_L0"], d["cl_max_2d"])
    print(f"{kts:>7} kts  {d['mach']:>5.3f}  {d['Re']/1e6:>6.2f}M  "
          f"{d['cl_max_2d']:>10.4f}  {cl_w:>12.4f}")

# ============================================================
# 7. Plots
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# --- Left: Spanwise CL distribution ---
ax = axes[0]
y_nd = y_span / (SPAN/2)   # normalised −1 to +1
ax.plot(y_nd, local_cl, color="steelblue", linewidth=2.5)
ax.axhline(d250["cl_max_2d"], color="tomato", linewidth=1.5, linestyle="--",
           label=f"2D CL_max = {d250['cl_max_2d']:.3f}")
ax.fill_between(y_nd, local_cl, alpha=0.2, color="steelblue")
ax.set_xlabel("Spanwise position  2y/b")
ax.set_ylabel("Local section CL")
ax.set_title("Spanwise CL Distribution at Stall\n(optimal geometry, 250 kts)")
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(-1, 1)

# --- Right: Wing planform ---
ax2 = axes[1]
half_b = SPAN / 2
# Chord runs in x, span in y; leading edge at x=0 (no sweep)
planform_y  = [-half_b, 0,      half_b]
planform_TE = [opt_c_tip, C_ROOT, opt_c_tip]
ax2.fill_betweenx(planform_y, 0, planform_TE,
                   alpha=0.25, color="steelblue", label="Wing planform")
ax2.plot([0, 0], [-half_b, half_b], color="steelblue", linewidth=2, label="Leading edge")
ax2.plot([opt_c_tip, C_ROOT, opt_c_tip], [-half_b, 0, half_b],
          color="steelblue", linewidth=2, linestyle="--", label="Trailing edge")

ax2.annotate(f"Root = {C_ROOT*100:.0f} cm", xy=(C_ROOT/2, 0.01),
             ha="center", fontsize=9, color="tomato")
ax2.annotate(f"Tip = {opt_c_tip*100:.1f} cm", xy=(opt_c_tip/2, half_b + 0.01),
             ha="center", fontsize=9, color="steelblue")
ax2.annotate(f"λ = {opt_taper:.3f}", xy=(C_ROOT*0.8, -half_b*0.6),
             fontsize=10, color="gray")
ax2.annotate(f"AR = {opt_AR:.2f}", xy=(C_ROOT*0.8, -half_b*0.8),
             fontsize=10, color="gray")

ax2.set_xlabel("Chord [m]")
ax2.set_ylabel("Semi-span [m]")
ax2.set_title(f"Optimised Wing Planform\n(Tip twist = {opt_twist:.2f}°  washout)")
ax2.set_aspect("equal")
ax2.legend(fontsize=8, loc="lower right")
ax2.grid(True, alpha=0.3)
ax2.invert_xaxis()

plt.suptitle(f"NACA 9112 · Span {SPAN*100:.0f} cm · Root chord {C_ROOT*100:.0f} cm · "
             f"Wing CL_max = {opt_cl:.4f}  @250 kts",
             fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig("wing_optimization_phase2.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nPlot saved: wing_optimization_phase2.png")
