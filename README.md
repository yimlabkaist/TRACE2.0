# TRACE2.0 - Lineage Tracing

Code and materials from paper **"Directed evolution of CRISPR-Cas spacer acquisition
machinery enables noninvasive biological recording and lineage tracing in the
mammalian gut"** Yim and Jang et al., *Nature Chemical Biology* [In revision].

<p> The full paper and supplementary information can be accessed [here](URL_TO_PAPER). </p>
<p> Raw sequencing data can be found at NCBI SRA under [XXXXXXXXXXX](URL_TO_SRA). </p>

---

This repository provides the exact analysis code for extracting unique CRISPR arrays
from per-sample spacer-mapping files and building count / frequency / absolute-abundance
matrices across a time series, together with a runnable demo dataset and its outputs
(reviewer comment R1.20).

## Dependencies

- Python ≥ 3.8
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)

```bash
pip install pandas numpy
```

No other third-party packages are required.

## Contents

```
TRACE_crispr_array_demo/
├── code/
│   └── crispr_array_count_matrix.py     # the analysis script (pandas + numpy only)
├── demo_data/
│   ├── Ec257-beforegavage_uniq.txt      # day 0 (pre-gavage inoculum)
│   ├── Ec257-day1-rep4_uniq.txt
│   ├── Ec257-day10-rep4_uniq.txt
│   ├── Ec257-day20-rep4_uniq.txt
│   ├── Ec257-day30-rep4_uniq.txt
│   ├── Ec257-day40-rep4_uniq.txt
│   ├── Ec257-day50-rep4_uniq.txt
│   ├── Ec257-day60-rep4_uniq.txt
│   ├── Ec257-day80-rep4_uniq.txt
│   ├── Ec257-day100-rep4_uniq.txt
│   └── prop_array_expanded.csv          # per-sample array-expansion fractions
└── results/                             # produced by the command below
    ├── unique_arrays_counts.csv
    ├── unique_arrays_normfreq.csv
    ├── unique_arrays_absfreq.csv
    ├── unique_arrays_long.csv
    └── per_sample/<sample>_arrays.csv
```

## Input format

Each `*_uniq.txt` file is one sample (one mouse × one time point), comma-separated,
no header:

```
a1l1p1,ec257-v3r,4260750,4260782
a8l2p1,ec257-v3r,501296,501264
a8l2p2,ec257-v3r,3168840,3168808
```

| column | meaning |
|--------|---------|
| 1 | spacer ID `a{ARRAY}l{LEN}p{POS}`; the part before `p` (e.g. `a8l2`) identifies one observed array (one read); `p1,p2,…` are its ordered spacers (p1 = oldest/distal) |
| 2 | reference name (ignored) |
| 3 | spacer mapping **start** coordinate (defines spacer identity) |
| 4 | spacer mapping **end** coordinate (start > end ⇒ reverse strand) |

`prop_array_expanded.csv` has a `file` column (sample name) and a `prop. expanded`
column = the fraction of the total cell population that is array-bearing in that sample.

`Ec257-beforegavage` is treated as **day 0** (see `--day0-label`).

## How to run

```bash
cd TRACE_crispr_array_demo
python3 code/crispr_array_count_matrix.py \
    --input demo_data \
    --outdir results \
    --prop-expanded demo_data/prop_array_expanded.csv \
    --add-ancestor
```

Requirements: Python ≥ 3.8, `pandas`, `numpy`.

## Method

1. **Reconstruct arrays** — group spacers by array key (ID before `p`), order by
   position index, giving each observed array as an ordered tuple of spacer
   identities. A spacer's identity is its mapping start coordinate
   (`--spacer-identity start`, default; `start_end` for strand-aware).
2. **Unique arrays per sample** — collapse identical array tuples and count them
   → `per_sample/<sample>_arrays.csv` (`p1…pN`, `counts`, `norm_counts`).
3. **Count matrix** — union of all unique arrays across samples (the reference set)
   → `unique_arrays_counts.csv` (raw reads) and `unique_arrays_normfreq.csv`
   (within array-bearing population; each sample column sums to 1).
4. **Absolute abundance** — multiply each sample's normalized frequency by its
   `prop. expanded` → `unique_arrays_absfreq.csv` (fraction of the **total** cell
   population). With `--add-ancestor`, an ancestor row (identity 0, empty array)
   holds the unexpanded fraction `1 − prop. expanded`, so each column sums to 1
   over the whole population — ready for genealogy / Muller-plot reconstruction.
5. `unique_arrays_long.csv` — tidy long table with `array_id, p1…pN, sample, count,
   frequency, abs_frequency, day`.

## Options

| flag | default | description |
|------|---------|-------------|
| `--input` | (required) | directory or glob of `*_uniq.txt` files |
| `--outdir` | `crispr_array_results` | output directory |
| `--max-spacers` | `5` | number of spacer columns `p1…pN` kept |
| `--spacer-identity` | `start` | `start` or strand-aware `start_end` |
| `--prop-expanded` | none | CSV (`file`, `prop. expanded`) to compute absolute abundance |
| `--prop-col` | `prop. expanded` | column name in the prop-expanded CSV |
| `--add-ancestor` | off | add ancestor row (unexpanded fraction) |
| `--day0-label` | `beforegavage` | substring marking day 0 |

## Demo summary (this dataset)

10 samples (day 0, 1, 10, 20, 30, 40, 50, 60, 80, 100), replicate 4;
**82,350** unique CRISPR arrays across the series. The unexpanded (ancestor)
fraction declines from 0.95 at day 0 to ~0.01 by day 100 as recorded lineages expand.
