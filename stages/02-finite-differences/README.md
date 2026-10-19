# Stage 02: Finite differences on Black-Scholes

**Question:** Can I reproduce the published results of a finite difference scheme for a generalized Black-Scholes equation?

**Paper replicated:** Cen & Le (2011), "A robust and accurate finite difference method for a generalized Black-Scholes equation", J. Comput. Appl. Math. 235, 3728-3733.

**Period:** ~17-19 Jun 2026 (reconstructed from file modified dates)

## Method (from the paper)
- Central differences in space on a piecewise uniform mesh, refined near the strike K
- Implicit Euler in time
- Payoff max(x-K, 0) smoothed near x = K by a ninth-degree polynomial
- Truncated domain (0, Smax) with boundary conditions at both ends
- The paper proves second-order convergence in space

## Files
| File | Role |
|---|---|
| `pi_epsilon.m` | smoothed payoff |
| `space_mesh.m` | piecewise uniform spatial mesh |
| `solver.m` | implicit finite difference solver |
| `table1.m` | Test 1: sigma = 0.2(1 + t e^-x), K = 25, Smax = 100, T = 1, r = 0.06 |
| `table2.m` | Test 2: sigma = 0.4(2 + sin x), same other parameters |

## Results: paper vs mine
M = 1024, max-norm error against the N = 2048, M = 1024 reference solution (as in the paper).

| Test | N | Paper error | My error | Paper rate | My rate |
|---|---|---|---|---|---|
| 1 | 128 | 6.8027e-2 | 6.6890e-2 | 1.122 | 1.110 |
| 1 | 256 | 3.1258e-2 | 3.0990e-2 | 1.312 | 1.308 |
| 1 | 512 | 1.2588e-2 | 1.2518e-2 | 1.880 | 1.878 |
| 1 | 1024 | 3.4189e-3 | 3.4050e-3 | - | - |
| 2 | 128 | 5.7124e-2 | 5.6050e-2 | 1.551 | 1.540 |
| 2 | 256 | 1.9498e-2 | 1.9276e-2 | 1.919 | 1.912 |
| 2 | 512 | 5.1550e-3 | 5.1231e-3 | 2.284 | 2.280 |
| 2 | 1024 | 1.0586e-3 | 1.0550e-3 | - | - |

## Status
Reproduced. My errors are within about 2% of the paper's (largest gap at N = 128, shrinking as N grows) and my observed rates agree to about 0.01. Exact agreement was not obtained; I have not yet investigated the source of the small gap.

## Open items
- [ ] Find the source of the ~2% gap (mesh constants, smoothing parameter, how the reference solution is interpolated)
