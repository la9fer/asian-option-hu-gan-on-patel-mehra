# Research log

*Stages 01-06 were done before this repo existed. Dates are reconstructed from file timestamps and are approximate. Entries from 2026-10-20 onward are written as the work happens.*

## Reconstructed timeline
| Stage | Approx. period | Evidence | Status |
|---|---|---|---|
| 02 Finite differences (Cen & Le 2011) | ~17-19 Jun 2026 | .m files modified 17 and 19 Jun 2026 | Reproduced (errors within ~2% of the paper) |
| 04 Time stepping: Crank-Nicolson baseline (HONC Table 1) | ~28 Jul 2026 | Table1_CrankNicolson.m modified 28 Jul 2026 (re-saved later) | In progress: solver verified, paper's numbers not reproduced |

## Entries (newest first)
- 2026-10-20: added stage 02 (Cen & Le replication) with code, results table and README. Next: add the remaining stage folders one at a time.
- 2026-10-20: added stage 04 (Crank-Nicolson baseline for HONC Table 1). Solver converges at second order but does not match the paper's table; documented in the stage README. Next: `cn, bdf4` folder.
