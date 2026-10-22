# Stage 05: Grid stretching (in progress)

**Question:** Does the sinh grid stretching of Tavella & Randall restore fourth-order convergence with BDF4, as reported in Hu & Gan (2018), Tables 6-7?

**Period:** ~9 Aug - 8 Sep 2026 (file modified dates: table6.m 9 Aug, table7.m 8 Sep; reconstructed)

## Setup
European put, Smin = 0, sigma = 0.3, X = 15, r = 0.02, T = 0.5. Stretching around the strike, equidistant grid in the new variable, BDF4 in time, fourth-order central vs compact FD in space.
- Table 6 (paper Example 1): Smax = 45, xi = 12
- Table 7 (paper Example 2): Smax = 75, xi = 12 (the 20x20 and 40x40 rows reproduce the paper exactly)

| File | Role |
|---|---|
| `table6.m` | Table 6 replication |
| `table7.m` | Table 7 replication (extra grids 60, 160, 320) |

## Table 6
| Scheme | Grid | Paper error | My error | Paper order | My order | My time |
|---|---|---|---|---|---|---|
| Central | 20x20 | 0.0146 | 0.0146 | - | - | 0.0025 s |
| Central | 40x40 | 8.9287e-04 | 8.9204e-04 | 4.0435 | 4.0288 | 0.0053 s |
| Central | 80x80 | 6.0106e-05 | 5.5684e-05 | 3.8542 | 4.0018 | 0.0061 s |
| Compact | 20x20 | 0.0359 | 0.0359 | - | - | 0.0009 s |
| Compact | 40x40 | 0.0024 | 0.0024 | 3.8676 | 3.8889 | 0.0025 s |
| Compact | 80x80 | 1.5299e-04 | 2.1299e-04 | 3.9607 | 3.5082 | 0.0051 s |

## Table 7
Order = log(e_prev/e_cur) / log(N_cur/N_prev), relative to the previous row of my table. The paper's order in its 80x80 row is from 40 to 80; mine is from 60 to 80.
| Scheme | Grid | Paper error | My error | Paper order | My order | My time |
|---|---|---|---|---|---|---|
| Central | 20x20 | 0.0185 | 0.0185 | - | - | 0.1170 s |
| Central | 40x40 | 0.0011 | 0.0011 | 4.1010 | 4.0677 | 0.0063 s |
| Central | 60x60 | - | 2.1766e-04 | - | 4.00 | 0.0044 s |
| Central | 80x80 | 9.3713e-05 | 6.9059e-05 | 3.4261 | 3.99 | 0.0046 s |
| Central | 160x160 | - | 4.5301e-06 | - | 3.9302 | 0.0243 s |
| Central | 320x320 | - | 2.4110e-07 | - | 4.2319 | 0.5260 s |
| Compact | 20x20 | 0.0127 | 0.0127 | - | - | 0.0013 s |
| Compact | 40x40 | 6.4683e-04 | 6.4683e-04 | 4.4311 | 4.2980 | 0.0027 s |
| Compact | 60x60 | - | 4.0735e-05 | - | 6.82 | 0.0030 s |
| Compact | 80x80 | 3.4865e-05 | 1.4865e-05 | 4.3073 | 3.50 | 0.0043 s |
| Compact | 160x160 | - | 7.7943e-07 | - | 4.2533 | 0.0239 s |
| Compact | 320x320 | - | 5.9521e-08 | - | 3.7109 | 0.5044 s |

Times are single runs on my machine, in seconds (the first includes MATLAB start-up), and are not used for any comparison.

## Status
Tables 6 and 7 match the paper at 20x20 and 40x40 for both schemes. At 80x80 they differ: Table 6 central 5.57e-5 vs 6.01e-5 and compact 2.13e-4 vs 1.53e-4; Table 7 central 6.91e-5 vs 9.37e-5 and compact 1.49e-5 vs 3.49e-5. The two compact gaps are exactly 6.0e-5 and 2.0e-5 with otherwise identical digits, so I suspect a one-digit slip somewhere, but I have not determined whether it is in the paper or in my work. With stretching, the central scheme converges at about fourth order on every grid I ran (orders 3.9-4.2); the compact scheme is more erratic between 40 and 80 (3.5-6.8) and otherwise about 3.7-4.3.

## Open items
- [ ] Check the printed 80x80 values in the paper's Tables 6 and 7 against my transcription
- [ ] Find why the 80x80 rows differ
