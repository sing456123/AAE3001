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
