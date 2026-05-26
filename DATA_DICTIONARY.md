# Data Dictionary

This repo is currently Goldbach-centric. The main CSVs come from `scripts/shattering_mirrors.py` and `scripts/shattering_compressed.py`.

## Core count columns

- `N`: even target being tested.
- `n`: shorthand for `N / 2`.
- `r`: exact unordered Goldbach count `r_G(N)`.
- `pair_fiber_size`: same count as `r`, reused when the compressed script expands pair-level views.

## Heuristic and residual columns

- `h`: raw Hardy-Littlewood-style heuristic count used by the scripts.
- `h_scale`: one global multiplicative calibration factor fitted over the run so the normalized residual has mean 0.
- `h_cal`: calibrated heuristic `h_scale * h`.
- `h_floor`: `floor(h)`. This is retained because older diagnostics used it.
- `eps_h`: raw integer residual `r - floor(h)`.
- `eps_bucket`: legacy column name for a clipped `eps_h` label. This is now a diagnostic only, not the main clustering axis.
- `z_h_raw`: raw normalized residual `(r - h) / sqrt(h)` before global calibration.
- `z_h`: calibrated normalized residual `(r - h_cal) / sqrt(h_cal)`.
- `z_bucket`: current column name for the discretized `z_h` label used in summaries.

## Structural columns

- `rho30`: `N mod 30`.
- `delta10`: `N - 10*r`. This is a decimal artifact diagnostic, not a theorem-bearing quantity.
- `boost_G`: local odd-prime singular-series multiplier used inside `h`.
- `native_cluster`: combined label `z=<z_bucket>;rho30=<residue>`.

## Mirror / string diagnostics

- `fake_mirror`: `r|delta10|r`.
- `h_mirror`: `r|eps_h|r`.
- `root_decimal`: `eps_h|delta10|eps_h`.
- `root_native`: `eps_h|rho30|eps_h`.
- `native_core`: `r|(rho30,eps_h)|r`.

These string columns are descriptive views over the computed quantities. They are not proof objects.

## Compressed pair-fiber columns

- `p`, `q`: a Goldbach prime pair with `p + q = N`.
- `a_q`: `q mod 30`.
- `b_p`: `p mod 30`.
- `pair_chamber`: `a_q|N mod 30|b_p`.
- `pair_i`: compressed residue signature written as `a_q + i*b_p [N mod 30]`.
- `z_N`: number-level label used by the compressed script for plotting and CSV summaries.
