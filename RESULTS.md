# Results Snapshot

This document records factual outputs currently present in this repository.

## Output Inventory

- CSV files: 25 files in [outputs/csv](outputs/csv)
- HTML files: 4 files in [outputs/html](outputs/html)
- Plot PNG files: available in [outputs/plots](outputs/plots)

## Key Empirical Counts

From existing generated artifacts:

- [outputs/csv/goldbach_h_clusters.csv](outputs/csv/goldbach_h_clusters.csv): 49,999 rows
- [outputs/csv/goldbach_h_cluster_summary.csv](outputs/csv/goldbach_h_cluster_summary.csv): 30 rows
- [outputs/csv/goldbach_decimal_mirror_hits.csv](outputs/csv/goldbach_decimal_mirror_hits.csv): 13 rows
- [outputs/csv/shattering_compressed_100k_numbers.csv](outputs/csv/shattering_compressed_100k_numbers.csv): 26,517 rows
- [outputs/csv/shattering_compressed_100k_cluster_summary.csv](outputs/csv/shattering_compressed_100k_cluster_summary.csv): 8 rows

## Cluster/Family Highlights

### Native residual cluster summary (100k range)

Top rows in [outputs/csv/goldbach_h_cluster_summary.csv](outputs/csv/goldbach_h_cluster_summary.csv) include:

- `eps_bucket=>5, rho30=6, count=3321, min_N=246, max_N=99996`
- `eps_bucket=>5, rho30=0, count=3320, min_N=390, max_N=99990`
- `eps_bucket=>5, rho30=24, count=3320, min_N=324, max_N=99984`

These are descriptive chamber occupancies from the selected residual model.

### Decimal mirror hits

[outputs/csv/goldbach_decimal_mirror_hits.csv](outputs/csv/goldbach_decimal_mirror_hits.csv) contains 13 strict hits for the criterion used in code (`r_G(N) = floor(N/10)`).

This is treated as a base-10 artifact/diagnostic and not as a theorem claim.

### Family summaries

From [outputs/csv/goldbach_volume_5k_mod30_family_summary.csv](outputs/csv/goldbach_volume_5k_mod30_family_summary.csv) and [outputs/csv/goldbach_volume_5k_v2_family_summary.csv](outputs/csv/goldbach_volume_5k_v2_family_summary.csv):

- `shell_other`: 1321
- `cone_3`: 667
- `ring_5`: 334
- `mobius_3x5`: 166
- `spine_none`: 11

From [outputs/csv/goldbach_volume_10k_mod210_family_summary.csv](outputs/csv/goldbach_volume_10k_mod210_family_summary.csv):

- `cone_3`: 369
- `mobius_3x5`: 138
- `ring_5`: 46

## Visual Artifacts

Interactive pages:

- [outputs/html/goldbach_volume_5k.html](outputs/html/goldbach_volume_5k.html)
- [outputs/html/goldbach_volume_5k_mod30.html](outputs/html/goldbach_volume_5k_mod30.html)
- [outputs/html/goldbach_volume_5k_v2.html](outputs/html/goldbach_volume_5k_v2.html)
- [outputs/html/goldbach_volume_10k_mod210.html](outputs/html/goldbach_volume_10k_mod210.html)

Static dashboards and galleries:

- [outputs/plots/shattering_mirrors_100k_dashboard.png](outputs/plots/shattering_mirrors_100k_dashboard.png)
- [outputs/plots/shattering_mirrors_100k_gallery.png](outputs/plots/shattering_mirrors_100k_gallery.png)
- [outputs/plots/shattering_compressed_100k_dashboard.png](outputs/plots/shattering_compressed_100k_dashboard.png)
- [outputs/plots/shattering_compressed_100k_gallery.png](outputs/plots/shattering_compressed_100k_gallery.png)

## Scope and Caveat

- Goldbach computations here are empirical and heuristic-assisted.
- No proof of the Goldbach conjecture is claimed.
- Collatz `3m+1` and `5m+1` are acknowledged as context for integer-dynamics exploration, but this result set is Goldbach-centric.
