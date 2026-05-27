# Results Snapshot

What the checked-in outputs actually show, stated plainly.

The project is not trying to prove Goldbach; it is measuring exact Goldbach counts against a heuristic model.

## One-line summary

After fixing the calibration of `h(N)`, the normalized Goldbach residual `z_h`
shows **no privileged residue class mod 30**. The Goldbach pair counts behave as
the Hardy-Littlewood heuristic predicts. This is a confirmed negative result.

## What was tested

For every even `N` up to 100,000:

- `r` = exact unordered Goldbach pair count (FFT convolution, verified against brute force).
- `h(N)` = Hardy-Littlewood expected count.
- `z_h` = normalized residual `(r - alpha*h) / sqrt(alpha*h)`.

Question: does `z_h` depend on `rho30 = N mod 30`?

## What was found

A single global calibration factor `alpha` leaves an apparent `rho30` signal:
between-chamber spread of `z_h` = 0.176, permutation-test `p < 0.001`.

That signal is **a calibration artifact, not Goldbach structure**. It is fully
removed by calibrating `h(N)` per residue class:

| Calibration | between-chamber spread of `z_h` |
|-------------|---------------------------------|
| single global `alpha` | 0.176 |
| per-`rho30` `alpha` | 0.000 |

Cause: `h_goldbach` omits the standard factor of 2 and uses a crude divisibility
boost instead of the full local singular series. A single `alpha` cannot absorb
a per-residue effect. Per-residue calibration, or using the full singular
series in `h`, makes the residual flat.

## Known open issues in the current pipeline

- `z_h` still drifts with `N` (slope about `-0.39` vs `log N`). Separate bug, same root cause: `h(N)` is missing terms. Not yet fixed.
- The `z_bucket` summary uses 0.5-wide bins, which is too coarse to display a 0.2 effect clearly. Use finer labels or report `mean_z` directly.

## Decimal mirror hits

[outputs/csv/goldbach_decimal_mirror_hits.csv](outputs/csv/goldbach_decimal_mirror_hits.csv): 13 hits for `r_G(N) = floor(N/10)`.

Base-10 artifact / negative control. Not a result.

## Scope

- Goldbach counts are exact on the tested ranges.
- The conjecture itself, every even `N > 2` is a sum of two primes, is **not** tested or proven here; it is only observed on the finite range, as already known far beyond this repo.
- This repo measures *how many* representations exist relative to the heuristic.
- Nothing here is a proof.
- Most valuable result: a clean, correctly explained null.
