# Data Dictionary

This repo is currently Goldbach-centric. The main CSVs come from
`scripts/shattering_mirrors.py` and `scripts/shattering_compressed.py`.

## Core count columns

- `N`: even target being tested.
- `n`: shorthand for `N / 2`.
- `r`: exact unordered Goldbach count `r_G(N)`.
- `pair_fiber_size`: same count as `r`, reused when the compressed script expands pair-level views.

## Heuristic and residual columns

- `h`: raw heuristic count built from a reciprocal-log density convolution and the Goldbach singular-series boost.
- `h_alpha`, `h_beta`: fitted two-parameter correction terms used in `h_cal = h * exp(alpha + beta / log N)`. These are repeated on every row for run provenance.
- `h_cal`: calibrated heuristic after the positive two-parameter correction.
- `var_scale`: fitted `c` in the working variance model `Var(r_G(N) | N) ~= c * h_cal(N)`.
- `z_h_raw`: raw normalized residual `(r - h) / sqrt(h)` before calibration.
- `z_h`: calibrated normalized residual `(r - h_cal) / sqrt(c * h_cal)`.
- `z_bucket`: discretized `z_h` label used in summaries.
- `h_floor`: `floor(h)`. Retained only for legacy diagnostics.
- `eps_h`: raw integer residual `r - floor(h)`. Legacy diagnostic only.
- `eps_bucket`: clipped `eps_h` label. Legacy diagnostic only, not the clustering axis.

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

These string columns are descriptive legacy views over the computed quantities.
They are not proof objects, and the main clustering logic no longer depends on
them.

## Compressed pair-fiber columns

- `p`, `q`: a Goldbach prime pair with `p + q = N`.
- `a_q`: `q mod 30`.
- `b_p`: `p mod 30`.
- `pair_chamber`: `a_q|N mod 30|b_p`.
- `pair_i`: compressed residue signature written as `a_q + i*b_p [N mod 30]`.
- `z_N`: number-level label used by the compressed script for plotting and CSV summaries.
