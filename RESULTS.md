# Results Snapshot

This document records the current interpretation of the checked-in outputs.

## Output inventory

- CSV files in [outputs/csv](outputs/csv)
- HTML files in [outputs/html](outputs/html)
- Plot PNG files in [outputs/plots](outputs/plots)

## Current headline

The main residual analysis is no longer based on `eps_h` buckets. The checked-in `shattering_mirrors` pipeline now applies one global calibration factor `alpha` to `h(N)` and clusters on normalized `z_h` buckets instead.

That change matters because the previous `eps_h` bucket story was dominated by scale drift. A clipped tail bucket recurring across the whole range is not a chamber discovery; it is a sign that the variable is the wrong axis.

## Key empirical counts

From the generated artifacts:

- [outputs/csv/goldbach_h_clusters.csv](outputs/csv/goldbach_h_clusters.csv): 49,999 rows
- [outputs/csv/goldbach_h_cluster_summary.csv](outputs/csv/goldbach_h_cluster_summary.csv): 30 rows
- [outputs/csv/goldbach_decimal_mirror_hits.csv](outputs/csv/goldbach_decimal_mirror_hits.csv): 13 rows

## Normalized residual summary

On the 100k run, the calibrated scale is approximately `alpha = 1.225600`.

Top rows in [outputs/csv/goldbach_h_cluster_summary.csv](outputs/csv/goldbach_h_cluster_summary.csv) now fall in narrow near-zero `z_bucket` bands rather than in a drifting `eps_h` tail. Representative rows are:

- `z_bucket=+0.0, rho30=12, count=1283, min_N=72, max_N=99912`
- `z_bucket=+0.0, rho30=14, count=1279, min_N=14, max_N=99884`
- `z_bucket=+0.0, rho30=26, count=1267, min_N=56, max_N=99956`

Interpretation: the residual mass lives close to the calibrated mean. The cluster table is useful as a descriptive occupancy map, but it does not isolate a privileged chamber.

## Decimal mirror hits

[outputs/csv/goldbach_decimal_mirror_hits.csv](outputs/csv/goldbach_decimal_mirror_hits.csv) contains 13 strict hits for `r_G(N) = floor(N/10)`.

This remains a base-10 artifact / negative control, not a theorem claim.

## Compressed views

The compressed script still visualizes pair-fiber residue structure, but it now inherits the normalized `native_cluster` labels from `shattering_mirrors.py`. That keeps the pair-level plots aligned with the corrected residual axis.

## Scope and caveat

- Goldbach counts here are exact on the tested finite ranges.
- The heuristic comparisons are empirical and calibration-assisted.
- A cleaner residual plot is still not a proof.
- The most valuable result so far is negative: after removing the raw-residual drift, the chamber story becomes much less dramatic.
