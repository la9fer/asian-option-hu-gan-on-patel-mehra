# Stage 04: Time stepping (in progress)

**Question:** Can I reproduce Table 1 (second-order block) of Oosterlee, Leentvaar & Huang, "Accurate American Option Pricing by Grid Stretching and High Order Finite Differences"?

**Period:** ~28 Jul 2026 (file modified date; reconstructed)

## Setup
European call, K = 15, sigma = 0.3, r = 0.05, rd = 0.03, T = 0.5, Smax = 3K. Second-order central differences in space and Crank-Nicolson in time on an equidistant grid. Errors in u, Delta and Gamma at t = 0, infinity norm over the interior grid points, against the closed-form solution.

## Files
| File | Role |
|---|---|
| `Table1_CrankNicolson.m` | CN solver and error table |

## Results: paper vs mine
| N | u error | ratio | Delta error | ratio | Gamma error | ratio |
|---|---|---|---|---|---|---|
| 10 | 9.0e-2 | - | 5.3e-2 | - | 1.1e-2 | - |
| 20 | 1.5e-2 | 5.8 | 2.3e-2 | 2.3 | 6.2e-3 | 1.8 |
| 40 | 4.9e-3 | 3.2 | 6.4e-3 | 3.6 | 1.8e-3 | 3.5 |
| 80 | 1.1e-3 | 4.4 | 1.7e-3 | 3.8 | 4.6e-4 | 3.8 |
| 160 | 2.9e-4 | 3.8 | 4.1e-4 | 4.1 | 1.1e-4 | 4.1 |
| 320 | 7.1e-5 | 4.1 | 1.0e-4 | 4.0 | 2.9e-5 | 3.9 |
| 640 | 1.8e-5 | 3.9 | 2.6e-5 | 4.0 | 7.1e-6 | 4.0 |

## Convergence on finer grids (my MATLAB run)
[paste the extended table with N = 80...640 and the error ratios here]

## Status
The solver is correct: observed convergence is second order (error ratios approach 4). The paper's Table 1 numbers are NOT reproduced; only Gamma at 20x20 and 40x40 agree. Varying Smax (25-90) and the grid-size convention did not close the gap, and changing the number of time steps has no effect, so the discrepancy is spatial and likely due to a setup detail the paper does not state. I have not resolved it.

## Open items
- [ ] BDF4 files (from the `cn, bdf4` folder) to be added to this stage
- [ ] Reproduce the paper's fourth-order block of Table 1
