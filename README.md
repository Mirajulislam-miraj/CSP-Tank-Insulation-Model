# CSP-Tank-Insulation-Model
# Thermal Loss Analysis of a Multi-Layer Insulated Molten Salt Storage Tank

A transient numerical model (MATLAB/Octave) that predicts heat loss, tank temperature, stored energy, and thermal efficiency of a large-scale hot molten-salt thermal energy storage tank used in Concentrated Solar Power (CSP) plants.

The tank's cylindrical side wall and flat top/bottom heads are insulated with a **carbon-steel shell + insulating firebrick (IFB) + mineral wool** stack. Conduction through the insulation is coupled to external natural convection and radiation using a nonlinear surface energy balance, solved at every time step with the **Newton–Raphson method**.
## Key Result

| Metric | Value |
|---|---|
| Tank diameter × height | 32 m × 15 m |
| Salt mass | 2.41 × 10⁷ kg (60% NaNO₃ / 40% KNO₃) |
| Total thermal capacity | 2762 MWh |
| 24-hour standing heat loss | 8.45 MWh |
| Thermal efficiency after 24 h | **99.69 %** |
| Dominant loss mechanism | Convection (≈1.85× radiation) |

## Repository Structure

.
miraj.m                # Main MATLAB/Octave simulation script

## Model Overview

1. **Geometry** — vertical cylindrical tank, side wall + flat top/bottom heads.
2. **Insulation stack** (series conduction resistances):
   | Carbon steel shell | 35 mm | 45 |
   | Insulating firebrick (IFB) | 400 mm | 0.35 |
   | Mineral wool | 150 mm | 0.04 |

3. **Surface energy balance** (solved via Newton–Raphson for outer surface temperature `Ts`):

   (T - Ts)/R = h·A·(Ts - Tamb) + ε·σ·A·(Ts⁴ - Tamb⁴)

4. **Transient update** — stored energy and salt temperature are updated every hour (`dt = 3600 s`) over a 24-hour horizon (`Nt = 24`).
5. **Output** — temperature, stored energy, convective/radiative loss breakdown, and thermal efficiency vs. time.

## Running the Simulation
Requires MATLAB:
This regenerates `results/results.csv` and `results/results_plots.png`.
## Author
**Mirajul Islam**
B.Sc. in Electrical and Electronic Engineering, Port City International University, Chittagong, Bangladesh
LinkedIn: https://www.linkedin.com/in/mirajul-islam-b7b700296
