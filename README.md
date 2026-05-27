# Goldbach Diagnostics Workspace

Repository name aside, the active codebase is a Goldbach diagnostics workspace.
Collatz is historical context only, and the Busy Beaver script is a separate
toy encoding experiment rather than part of the main inference path.

## Mathematical status

- The Goldbach conjecture is still open.
- Nothing here is a proof.
- The scripts compute exact Goldbach counts on finite ranges and compare them to heuristic expectations.
- The checked-in Collatz references are framing only. The current code and outputs are overwhelmingly Goldbach-centric.

## Model Status

The main pipeline treats the heuristic as a model to be repaired, not as a
finished reference curve.

- Raw `h(N)` is built from a reciprocal-log density convolution times the Hardy-Littlewood singular-series boost.
- Calibrated `h_cal(N)` is a positive two-parameter correction of the form `h * exp(alpha + beta / log N)`.
- `z_h` is variance-scaled as `(r - h_cal) / sqrt(c * h_cal)` using an empirical fit for `c`.
- Legacy `eps_h`, `h_floor`, and mirror strings are retained only as diagnostics and backward-compatible CSV columns.

That is a stricter setup than the older `alpha * h` normalization, but it does
not by itself remove all model misspecification.

## Repository layout

```text
.
├── scripts/
│   ├── busy_beaver_waring_goldbach.py
│   ├── goldbach_automata.py
│   ├── goldbach_native_filter.py
│   ├── goldbach_volume.py
│   ├── shattering_compressed.py
│   └── shattering_mirrors.py
├── tests/
├── outputs/
│   ├── csv/
│   ├── html/
│   └── plots/
├── .github/workflows/ci.yml
├── BOUNDARY_RELATION.md
├── DATA_DICTIONARY.md
├── README.md
├── RESULTS.md
├── next_steps.md
└── requirements.txt
```

## Environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

The main scripts now default to the organized output folders.

```bash
python scripts/shattering_mirrors.py --max-n 100000 --plot --plot-prefix outputs/plots/shattering_mirrors_100k
python scripts/shattering_compressed.py \
  --max-n 100000 \
  --plot \
  --plot-prefix outputs/plots/shattering_compressed_100k \
  --numbers-out outputs/csv/shattering_compressed_100k_numbers.csv \
  --summary-out outputs/csv/shattering_compressed_100k_cluster_summary.csv \
  --signatures-out outputs/csv/shattering_compressed_100k_pair_signatures.csv

# Busy Beaver space-time tableau with a Waring-Goldbach-style encoding
python scripts/busy_beaver_waring_goldbach.py --machine bb3 --plot
python scripts/busy_beaver_waring_goldbach.py --machine bb4 --max-steps 200 --plot --plot-shells

# Goldbach z/rho automata search and dashboard
python scripts/goldbach_automata.py --max-n 100000

# Legacy 3D volume splash, rebuilt around the same raw heuristic
python scripts/goldbach_volume.py --max-n 5000 --plot --html
```

Run the regression checks with:

```bash
python -m unittest discover -s tests -v
```

## HTML splash previews

Click any screenshot to open the interactive HTML artifact:

[![Goldbach Volume 5k HTML preview](outputs/plots/goldbach_volume_5k_html_screenshot.png)](outputs/html/goldbach_volume_5k.html)

[![Goldbach Volume 5k mod30 HTML preview](outputs/plots/goldbach_volume_5k_mod30_html_screenshot.png)](outputs/html/goldbach_volume_5k_mod30.html)

[![Goldbach Volume 5k v2 HTML preview](outputs/plots/goldbach_volume_5k_v2_html_screenshot.png)](outputs/html/goldbach_volume_5k_v2.html)

[![Goldbach Volume 10k mod210 HTML preview](outputs/plots/goldbach_volume_10k_mod210_html_screenshot.png)](outputs/html/goldbach_volume_10k_mod210.html)

## Busy Beaver previews

The BB4 run needs a larger step cap than the default because the built-in 4-state
machine halts after 107 steps.

- BB4 dashboard: [outputs/plots/busy_beaver_waring_goldbach_bb4_dashboard.png](outputs/plots/busy_beaver_waring_goldbach_bb4_dashboard.png)
- BB4 shell diagram: [outputs/plots/busy_beaver_waring_goldbach_bb4_shells.png](outputs/plots/busy_beaver_waring_goldbach_bb4_shells.png)
- BB4 cell tableau: [outputs/csv/busy_beaver_waring_goldbach_bb4_cells.csv](outputs/csv/busy_beaver_waring_goldbach_bb4_cells.csv)

## Data model

The main dataset contains:

- exact counts `r`
- raw heuristic `h`
- fitted heuristic coefficients `h_alpha` and `h_beta`
- fitted variance scale `var_scale`
- calibrated heuristic `h_cal`
- calibrated normalized residual `z_h`
- legacy diagnostic residual `eps_h`
- label columns `z_bucket` and `native_cluster`

See [DATA_DICTIONARY.md](DATA_DICTIONARY.md) for column definitions.

## Scope boundaries

- Goldbach exact counting is the core result-producing path.
- The residue-family plots are descriptive diagnostics.
- The decimal mirror hits are negative controls and base-10 artifacts.
- The compressed pair-fiber plots are visualization aids, not theorem-bearing objects.
- The Busy Beaver script is a separate encoding toy, not evidence about Goldbach.
- Collatz is not part of the current measured pipeline; the repo name is historical.
- If one `z_bucket` label spans most of the range, that is a warning about the labeling choice, not a result by itself.

## References in repo

- [BOUNDARY_RELATION.md](BOUNDARY_RELATION.md)
- [RESULTS.md](RESULTS.md)
- [next_steps.md](next_steps.md)
