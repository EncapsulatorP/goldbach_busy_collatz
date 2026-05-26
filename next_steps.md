# Proposed Next Steps

This roadmap is ordered by value, not by novelty theater.

## 1. Keep the residual analysis honest

1. Treat `z_h` and `z_bucket` as the main residual diagnostics.
2. Keep `eps_h` and `eps_bucket` only as recorded raw diagnostics.
3. Add one or two explicit null summaries to `RESULTS.md`, for example cluster-share stability across `N` bins.

## 2. Tighten reproducibility

1. Keep `requirements.txt` current.
2. Expand the setup section in `README.md` only when a new dependency is added.
3. Add a small environment-report script if the dependency surface grows.

## 3. Grow the regression suite

1. Keep exact low-range Goldbach checks.
2. Keep the mod-30 pair-chamber invariants.
3. Add a regression check that the calibrated residual mean stays centered.

## 4. Separate source from artifacts

1. Decide which generated outputs belong in version control and remove the rest.
2. If large artifact check-in remains intentional, document why.
3. Prefer scripts and summaries over static bulk outputs.

## 5. Clarify interpretation boundaries

1. Keep README and RESULTS explicit about diagnostics versus claims.
2. Write down negative results when a pattern collapses under better normalization.
3. Treat visually striking plots as prompts for checks, not as findings.

## 6. Standardize script interfaces further

1. Continue moving defaults into `outputs/csv`, `outputs/plots`, and `outputs/html`.
2. Give every script separate CSV and plot output paths where that reduces confusion.
3. Add one short “reference run” command per script.

## 7. Benchmark scale deliberately

1. Measure runtime and memory at 10k, 100k, and 1M.
2. Cache reusable intermediate artifacts only after the measurements justify it.
3. Avoid premature chunking logic unless the benchmark data says it is needed.

## 8. Tau Collatz generalization

If this is added, it should be framed as a pedagogical re-implementation of the standard `qx+1` family with no novelty expected.

1. Do it only after the Goldbach diagnostics and repo hygiene are in good shape.
2. Label it as computational exploration from the start.
3. Do not let it displace higher-value maintenance and negative-result writeups.
