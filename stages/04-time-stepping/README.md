# Stage 04: Time stepping (in progress)

**Question:** Can I reproduce Table 1 (second-order block) of Oosterlee, Leentvaar & Huang, "Accurate American Option Pricing by Grid Stretching and High Order Finite Differences"?

**Period:** ~28 Jul 2026 (file modified date at the time; reconstructed)

## Setup
European call, K = 15, sigma = 0.3, r = 0.05, rd = 0.03, T = 0.5, Smax = 3K. Second-order central differences in space and Crank-Nicolson in time on an equidistant grid, with as many time steps as space steps. Errors in u, Delta and Gamma at t = 0, infinity norm over the interior grid points, against the closed-form solution.

## Files
| File | Role |
|---|---|
| `Table1_CrankNicolson.m` | CN solver and error table |

## Results: paper vs mine
| Grid | u (paper / mine) | Delta (paper / mine) | Gamma (paper / mine) |
|---|---|---|---|
| 10x10 | 1.3e-1 / 9.0e-2 | 9.7e-2 / 5.3e-2 | 2.1e-2 / 1.1e-2 |
| 20x20 | 3.3e-2 / 1.5e-2 | 9.6e-3 / 2.3e-2 | 6.2e-3 / 6.2e-3 |
| 40x40 | 6.4e-3 / 4.9e-3 | 1.7e-3 / 6.4e-3 | 1.9e-3 / 1.8e-3 |

## Convergence on finer grids
| N | u error | ratio | Delta error | ratio | Gamma error | ratio |
|---|---|---|---|---|---|---|
| 10 | 9.0e-2 | - | 5.3e-2 | - | 1.1e-2 | - |
| 20 | 1.5e-2 | 5.8 | 2.3e-2 | 2.3 | 6.2e-3 | 1.8 |
| 40 | 4.9e-3 | 3.2 | 6.4e-3 | 3.6 | 1.8e-3 | 3.5 |
| 80 | 1.1e-3 | 4.4 | 1.7e-3 | 3.8 | 4.6e-4 | 3.8 |
| 160 | 2.9e-4 | 3.8 | 4.1e-4 | 4.1 | 1.1e-4 | 4.1 |
| 320 | 7.1e-5 | 4.1 | 1.0e-4 | 4.0 | 2.9e-5 | 3.9 |
| 640 | 1.8e-5 | 3.9 | 2.6e-5 | 4.0 | 7.1e-6 | 4.0 |

## Status
The solver is correct: error ratios per doubling approach 4 (second order) for u, Delta and Gamma. The paper's Table 1 numbers are NOT reproduced: only Gamma at 20x20 and 40x40 agree. Varying Smax (25-90) and the grid-size convention did not close the gap, and changing the number of time steps has no effect on the errors, so the discrepancy is spatial and probably due to a setup detail the paper does not state. Not resolved.

## Open items
- [ ] Add the BDF4 files (from the cn, bdf4 folder) to this stage
- [x] BDF4 files added (Hu & Gan Tables 4-5, above)
