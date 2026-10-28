Applying the Hu & Gan (2018) methodology (fourth-order compact finite differences, BDF4 time stepping, sinh grid stretching) to the Patel & Mehra formulation of the Asian option PDE (Rogers & Shi reduction). Written as a research log: each stage adds one piece, from finite-difference basics to the final solver.

## Headline results

| Result | Evidence |
|---|---|
| Matches semi-analytic reference prices | Max error 2.3e-6 on the 27 cases of Patel & Mehra Table 1; similar agreement across all 63 cases in Tables 1-3 (reference values: Zhang) |
| BDF4 needs about 20x fewer time steps than a Crank-Nicolson compact scheme at matched accuracy | N = 80, error about 4e-4: 320 steps (BDF4) vs 6,400 steps (CN with dt = h^2) |
| Stretching at the payoff kink gives clean fourth-order convergence | Native z-coordinate solver, sigma = 0.2, r = 0.09, K = 110, N = 320: error 1.62e-6 (stretched) vs 1.69e-5 (uniform); observed order about 4.0 (stretched) vs erratic (uniform) |
| Boundary error near x = 1 is a corner compatibility defect, not a scheme bug | Interior data force u_t(1,0) = 0 but the boundary data have nonzero slope at t = 0; confirmed with compatible vs incompatible boundary ramps |

## Scope and caveats

- The step-count comparison is the implementation-independent claim. Wall-clock comparisons here use dense-matrix implementations and are not claimed as general speedups.
- Results are for the test cases listed above; domain truncation (z_max) must be rechecked for other (sigma, T).
- Apparent errata found in benchmark sources (a likely digit transposition in one reference value, and an internal inconsistency between two of Patel & Mehra's tables) are recorded in LOG.md.

## Repository layout

- `stages/`: one folder per stage, each with code, a README and results
- `LOG.md`: dated timeline of the work
- `references.md`: papers used

| Stage | Topic |
|---|---|
| 01 | Background notes |
| 02 | Replication of Cen & Le (2011) |
| 03 | Compact schemes |
| 04 | Time stepping (Hu & Gan Tables 4-5, BDF3/BDF4) |
| 05 | Grid stretching (Hu & Gan Tables 6-7) |
| 06 | Asian PDE and Rogers-Shi reduction |
| 07 | Final solver and new results |

## Running the code

MATLAB (R20xx or later) and Python 3 with NumPy and SciPy. Each stage folder's README gives the script to run.

## References

See `references.md`. Key papers: Hu & Gan (2018); Patel & Mehra (2017-2019); Tangman, Gopaul & Bhuruth (2008); Rogers & Shi (1995); Zhang (2001).
