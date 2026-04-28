# AAE3001 — Propeller and Wing Aerodynamic Design

A student project covering two aerodynamic design problems:
1. Propeller thrust analysis using Blade Element Momentum Theory (BEMT)
2. Wing airfoil selection for maximum lift coefficient using NeuralFoil

## Overview

| Script | Purpose |
|---|---|
| `bem_validation_script.py` | BEM propeller validation + C_T optimisation |
| `naca_sweep_phase1.py` | NACA 4-digit sweep for max CL (2D airfoil, Phase 1) |

## Requirements

```
numpy scipy neuralfoil aerosandbox pandas
```

Install with:

```bash
pip install numpy scipy neuralfoil aerosandbox pandas
```

## Usage

**Propeller BEM:**
```bash
python3 bem_validation_script.py
```
Prints C_T validation results for twist angles 1°, 8°, 15°, then runs the SLSQP optimiser for the drone propeller geometry.

**Wing airfoil sweep:**
```bash
python3 naca_sweep_phase1.py
```
Sweeps all 656 valid NACA 4-digit profiles and ranks them by maximum CL. Saves full results to `naca_sweep_results.csv`.

## Validation

The BEM implementation is validated against wind-tunnel experimental data from:

> *Blade Element Momentum Theory* (course reference, AAE3001)

A standard 3-blade rectangular propeller (R = 1.829 m, c = 0.1524 m, 600 RPM, NACA 0012) is used. The table below compares the computed C_T against experimental values for three uniform twist angles.

| Twist angle θ | C_T (BEM, computed) | C_T (Experiment) | Error  |
|:---:|:---:|:---:|:---:|
| 1°  | 0.00019 | 0.00017 | +9.7%  |
| 8°  | 0.00498 | 0.00442 | +12.7% |
| 15° | 0.01162 | 0.01078 | +7.7%  |

The BEM code consistently over-predicts C_T by ~8–13%. This is expected: the model assumes zero drag (`dD sin φ ≈ 0`) and includes no tip-loss correction, both of which reduce real-world thrust.

## Wing Airfoil Selection — Phase 1

### Flight Conditions

| Parameter | Value |
|---|---|
| Altitude | 10,000 ft (3,048 m) |
| Speed range | 200 – 300 kts |
| Mach range | 0.313 – 0.470 |
| Reynolds number (chord = 2 m) | 11.0M – 16.5M |
| Air density | 0.905 kg/m³ |

Compressibility is accounted for using the **Prandtl-Glauert correction** applied on top of NeuralFoil's incompressible predictions:

$$C_L^{\text{corrected}} = \frac{C_L^{\text{incompressible}}}{\sqrt{1 - M^2}}$$

| Speed | Mach | PG factor |
|:---:|:---:|:---:|
| 200 kts | 0.313 | 1.053 |
| 250 kts | 0.392 | 1.087 |
| 300 kts | 0.470 | 1.133 |

### Results — Top 10 NACA Airfoils by Maximum CL

656 valid NACA 4-digit profiles were evaluated across all three speed conditions.

| Rank | NACA | CL_max @ 200 kts | CL_max @ 250 kts | CL_max @ 300 kts | Mean CL |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | **NACA 9112** | 2.557 | 2.676 | 2.818 | **2.684** |
| 2 | NACA 9110 | 2.565 | 2.675 | 2.808 | 2.683 |
| 3 | NACA 8112 | 2.552 | 2.667 | 2.807 | 2.675 |
| 4 | NACA 7112 | 2.530 | 2.644 | 2.781 | 2.652 |
| 5 | NACA 9812 | 2.532 | 2.643 | 2.779 | 2.651 |
| 6 | NACA 8110 | 2.542 | 2.642 | 2.764 | 2.649 |
| 7 | NACA 9815 | 2.523 | 2.633 | 2.768 | 2.641 |
| 8 | NACA 9810 | 2.512 | 2.621 | 2.757 | 2.630 |
| 9 | NACA 8812 | 2.509 | 2.620 | 2.756 | 2.629 |
| 10 | NACA 9821 | 2.513 | 2.616 | 2.743 | 2.624 |

**Best airfoil: NACA 9112** — 9% camber, camber at 10% chord, 12% thickness.

> Note: Top-ranked profiles all use high camber (7–9%), which maximises CL but increases pitching moment. A camber limit of 4–5% is typical for practical fixed-wing aircraft.
