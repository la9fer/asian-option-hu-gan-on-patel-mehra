# Stage 04: Time stepping (in progress)

# Part 1: HONC Table 1 baseline

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

## Part 2: Hu & Gan (2018), Tables 4-5 (CN vs BDF4, unstretched grids)

**Paper:** Hu & Gan, "High order method for Black-Scholes PDE" (see references.md).
**Period:** ~11 Aug 2026 (file modified dates; reconstructed)
**Setup:** European put, Smax = 45, Smin = 0, sigma = 0.3, X = 15, r = 0.02, T = 0.5. Equidistant grid, same number of space and time steps. Space: fourth-order central FD (Table 4) and fourth-order compact FD (Table 5). Time: Crank-Nicolson vs BDF4. Error against the closed-form price.

| File | Role |
|---|---|
| `table4.m` | Table 4: fourth-order central FD in space |
| `table5.m` | Table 5: fourth-order compact FD in space |

### Table 4
| Scheme | Grid | Paper error | My error | Paper order | My order |
|---|---|---|---|---|---|
| CN | 10x10 | 0.0952 | 0.0952 | - | - |
| CN | 20x20 | 0.0158 | 0.0158 | 2.4547 | 2.5911 |
| CN | 40x40 | 0.0042 | 0.0042 | 1.9396 | 1.9268 |
| BDF4 | 10x10 | 0.0951 | 0.0951 | - | - |
| BDF4 | 20x20 | 0.0157 | 0.0157 | 2.4547 | 2.5977 |
| BDF4 | 40x40 | 0.0041 | 0.0041 | 1.9569 | 1.9269 |

### Table 5
| Scheme | Grid | Paper error | My error | Paper order | My order |
|---|---|---|---|---|---|
| CN | 20x20 | 0.0173 | 0.0173 | - | - |
| CN | 40x40 | 0.0043 | 0.0043 | 2.0058 | 2.0206 |
| CN | 80x80 | 0.0011 | 0.0011 | 1.9771 | 1.9404 |
| BDF4 | 20x20 | 0.0172 | 0.0172 | - | - |
| BDF4 | 40x40 | 0.0042 | 0.0042 | 2.0237 | 2.0213 |
| BDF4 | 80x80 | 0.0011 | 0.0011 | 1.9540 | 1.9412 |

Notes: the paper's first-step order in Table 4 (2.4547) is not consistent with its own printed errors (0.0952 to 0.0158 gives about 2.59); mine are computed from my errors. The paper's Table 5 labels the BDF4 grids 10/20/40, but its errors equal the CN values at 20/40/80, so I label them 20/40/80.

**Status:** reproduced. Errors match the paper's printed values. On unstretched grids the observed order stays near 2 for both CN and BDF4, so BDF4 gives no gain here because of the non-smooth payoff. This motivates stage 05.

## Status
The solver is correct: error ratios per doubling approach 4 (second order) for u, Delta and Gamma. The paper's Table 1 numbers are NOT reproduced: only Gamma at 20x20 and 40x40 agree. Varying Smax (25-90) and the grid-size convention did not close the gap, and changing the number of time steps has no effect on the errors, so the discrepancy is spatial and probably due to a setup detail the paper does not state. Not resolved.

## Open items
- [x] BDF4 files added (Hu & Gan Tables 4-5, above)
