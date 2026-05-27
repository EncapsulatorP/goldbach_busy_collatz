# Results Snapshot

What the checked-in outputs actually show, stated plainly.

The project is not trying to prove Goldbach; it is measuring exact Goldbach
counts against a heuristic model and iterating on that model when the residuals
show obvious misspecification.

## One-line summary

The exact Goldbach counts are solid; the heuristic is better than before, but it
is still not finished. The current normalized residual is usable as a
diagnostic, not as a final null-result engine.

## What was tested

For every even `N` up to 100,000:

- `r` = exact unordered Goldbach pair count (FFT convolution, verified against brute force).
- `h(N)` = reciprocal-log density convolution times the singular-series boost.
- `h_cal(N)` = fitted positive correction `h * exp(alpha + beta / log N)`.
- `z_h` = normalized residual `(r - h_cal) / sqrt(c * h_cal)`.

Question: does `z_h` depend on `rho30 = N mod 30`?

## Current Diagnostics

On the current 100k run:

- `h_alpha ≈ 0.1422`
- `h_beta ≈ -1.6652`
- `var_scale ≈ 0.2177`
- `mean(z_h) ≈ -0.038`
- `sd(z_h) ≈ 1.01`
- drift slope of `z_h` against `log N` is still about `-0.129`
- the spread of `mean(z_h)` across `rho30` classes is still about `1.12`

So the improved normalization fixes the worst scale problem, but it does **not**
support the stronger claim that the residue signal is gone. The remaining signal
should be treated as heuristic error until a better model or a stronger
validation package says otherwise.

## Known open issues in the current pipeline

- `z_h` still drifts with `N`, so the two-parameter correction is an improvement, not an endpoint.
- `rho30`-conditioned means are still separated enough that residue effects cannot yet be called a null.
- The `z_bucket` summary uses 0.5-wide bins, which is fine for dashboards but too coarse for final inference.
- There is still no standalone replication note or benchmark-style statistical appendix.

## Clarifications

- The exact count `r` in this repo is an **unordered** Goldbach count, so the leading constant convention is different from the ordered-count formula often written with `2*C2`.
- Rewriting the boost over divisors of `N` rather than `N/2` is conceptually cleaner, but for even `N` it does not change the odd-prime support.
- The critique that the heuristic is incomplete is valid; the stronger claim that the current code was using the wrong odd-prime support is overstated.

## Decimal mirror hits

[outputs/csv/goldbach_decimal_mirror_hits.csv](outputs/csv/goldbach_decimal_mirror_hits.csv): 13 hits for `r_G(N) = floor(N/10)`.

Base-10 artifact / negative control. Not a result.

## Scope

- Goldbach counts are exact on the tested ranges.
- The conjecture itself, every even `N > 2` is a sum of two primes, is **not** tested or proven here; it is only observed on the finite range, as already known far beyond this repo.
- This repo measures *how many* representations exist relative to a heuristic model that is still under repair.
- Nothing here is a proof.
- Most valuable current result: the exact-count pipeline is trustworthy, while the heuristic layer still needs work.
