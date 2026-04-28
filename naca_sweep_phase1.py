# -*- coding: utf-8 -*-
"""Phase 1: NACA 4-digit airfoil sweep for maximum CL
Conditions: 10,000 ft, 200-300 kts
NeuralFoil (incompressible) + Prandtl-Glauert compressibility correction
"""

import numpy as np
import pandas as pd
import aerosandbox as asb
import neuralfoil as nf

# ============================================================
# 1. ISA Conditions at 10,000 ft
# ============================================================
alt_m   = 10_000 * 0.3048
T       = 288.15 - 0.0065 * alt_m
rho     = 1.225 * (T / 288.15) ** 4.256
a_sound = np.sqrt(1.4 * 287 * T)
mu      = 1.458e-6 * T**1.5 / (T + 110.4)
nu      = mu / rho

# ============================================================
# 2. Flight Conditions
# ============================================================
CHORD      = 2.0   # Wing chord [m] — adjust to your wing design
SPEEDS_KTS = [200, 250, 300]

conditions = []
for spd in SPEEDS_KTS:
    V    = spd * 0.51444
    mach = V / a_sound
    Re   = V * CHORD / nu
    pg   = 1.0 / np.sqrt(1.0 - mach**2)   # Prandtl-Glauert factor
    conditions.append({"kts": spd, "V": V, "mach": mach, "Re": Re, "pg": pg})

print("Flight conditions at 10,000 ft:")
for c in conditions:
    print(f"  {c['kts']} kts  M={c['mach']:.3f}  Re={c['Re']/1e6:.2f}M  PG={c['pg']:.4f}")
print()

# ============================================================
# 3. NACA 4-digit Parameter Space
# ============================================================
THICKNESS_OPTIONS = [6, 8, 10, 12, 15, 18, 21, 24]
ALPHA = np.linspace(-5, 20, 51)

valid_profiles = [
    (m, p, tt)
    for m in range(10)
    for p in range(10)
    for tt in THICKNESS_OPTIONS
    if not (m == 0 and p != 0)   # NACA 00xx: symmetric, P must be 0
    and not (m > 0 and p == 0)   # cambered needs a camber position
]

print(f"Sweeping {len(valid_profiles)} NACA profiles across {len(SPEEDS_KTS)} speed conditions...")

# ============================================================
# 4. Main Sweep
# ============================================================
results = []

for i, (m, p, tt) in enumerate(valid_profiles):
    name = f"naca{m}{p}{tt:02d}"

    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(valid_profiles)} — {name.upper()}")

    try:
        airfoil = asb.Airfoil(name)
    except Exception:
        continue

    row = {"NACA": name.upper()}
    cl_max_all = []

    for cond in conditions:
        try:
            aero       = nf.get_aero_from_airfoil(airfoil=airfoil, alpha=ALPHA,
                                                   Re=cond["Re"], model_size="large")
            confidence = float(np.mean(aero["analysis_confidence"]))
            cl_pg      = np.array(aero["CL"]) * cond["pg"]   # compressibility correction
            idx        = int(np.argmax(cl_pg))
            cl_max     = float(cl_pg[idx])
            alpha_stall = float(ALPHA[idx])
        except Exception:
            cl_max = confidence = np.nan
            alpha_stall = np.nan

        row[f"CL_max_{cond['kts']}kts"]       = cl_max
        row[f"alpha_stall_{cond['kts']}kts"]  = alpha_stall
        cl_max_all.append(cl_max)

    row["CL_max_mean"] = float(np.nanmean(cl_max_all))
    row["CL_max_min"]  = float(np.nanmin(cl_max_all))   # worst case across speeds
    results.append(row)

# ============================================================
# 5. Results
# ============================================================
df = (pd.DataFrame(results)
        .dropna(subset=["CL_max_mean"])
        .sort_values("CL_max_mean", ascending=False)
        .reset_index(drop=True))

df.to_csv("naca_sweep_results.csv", index=False)

W = 82
print("\n" + "="*W)
print("TOP 10 NACA AIRFOILS FOR MAXIMUM CL  (10,000 ft, 200-300 kts, chord=2m)")
pg_str = "  ".join(f"{c['kts']}kts: PG={c['pg']:.3f}" for c in conditions)
print(f"Prandtl-Glauert corrected  |  {pg_str}")
print("="*W)
print(f"{'Rank':<5} {'NACA':<10} {'CL@200kts':>10} {'CL@250kts':>10} {'CL@300kts':>10} {'Mean CL':>9} {'Min CL':>8}")
print("-"*W)
for rank, (_, row) in enumerate(df.head(10).iterrows(), 1):
    print(f"{rank:<5} {row['NACA']:<10} "
          f"{row['CL_max_200kts']:>10.4f} "
          f"{row['CL_max_250kts']:>10.4f} "
          f"{row['CL_max_300kts']:>10.4f} "
          f"{row['CL_max_mean']:>9.4f} "
          f"{row['CL_max_min']:>8.4f}")

best = df.iloc[0]
print(f"\nBest airfoil : {best['NACA']}  (mean CL_max = {best['CL_max_mean']:.4f})")
print(f"Full results : naca_sweep_results.csv")
