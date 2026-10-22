# Research log

*Stages 01-06 were done before this repo existed. Dates are reconstructed from file timestamps and are approximate. Entries from 2026-10-20 onward are written as the work happens.*

## Reconstructed timeline
| Stage | Approx. period | Evidence | Status |
|---|---|---|---|
| 02 Finite differences (Cen & Le 2011) | ~17-19 Jun 2026 | .m files modified 17 and 19 Jun 2026 | Reproduced (errors within ~2% of the paper) |
| 04 Time stepping: HONC Table 1 baseline (CN) and Hu & Gan Tables 4-5 | ~28 Jul - 11 Aug 2026 | Table1_CrankNicolson.m modified 28 Jul 2026; table4.m, table5.m modified 11 Aug 2026 | HONC baseline: solver verified, paper's numbers not reproduced. Hu & Gan Tables 4-5: reproduced |
| 05 Grid stretching: Hu & Gan Tables 6-7 | ~9 Aug - 8 Sep 2026 | table6.m modified 9 Aug, table7.m 8 Sep 2026 | Partly: both tables match at 20 and 40, differ at 80 |

## Entries (newest first)
- 2026-10-23: added Hu & Gan (2018) Tables 4-7 replication code and results (stages 04 and 05). Next: check the 80x80 values in the paper, then the `novel` folder.
- 2026-10-20: added stage 02 (Cen & Le replication) with code, results table and README. Next: add the remaining stage folders one at a time.
- 2026-10-20: added stage 04 (Crank-Nicolson baseline for HONC Table 1). Solver converges at second order but does not match the paper's table; documented in the stage README. Next: `cn, bdf4` folder.
