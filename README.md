# AAE3001 — Propeller Design using Blade Element Momentum Theory

A student project implementing Blade Element Momentum Theory (BEMT) to analyse and optimise propeller aerodynamic performance.

## Overview

This project is split into two parts:

1. **Validation** — computes the thrust coefficient C_T for a standard rectangular propeller and compares against wind-tunnel experimental data.
2. **Optimisation** — uses `scipy` to find the chord and twist distributions that maximise C_T for a small drone propeller under given geometric constraints.

## Requirements

```
numpy
scipy
```

Install with:

```bash
pip install numpy scipy
```

## Usage

```bash
python3 bem_validation_script.py
```

The script first prints the C_T validation results for twist angles of 1°, 8°, and 15°, then runs the SLSQP optimiser and prints the optimal propeller geometry.

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
