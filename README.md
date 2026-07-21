# TRACE2.0 - Lineage Tracing

Code and materials from paper **"Directed evolution of CRISPR-Cas spacer acquisition
machinery enables noninvasive biological recording and lineage tracing in the
mammalian gut"** Yim and Jang et al., *Nature Chemical Biology* [In revision].

<p> The full paper and supplementary information can be accessed [here](URL_TO_PAPER). </p>
<p> Raw sequencing data can be found at NCBI SRA under [PRJNA1498595](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1498595). </p>

---

This repository provides the analysis code for extracting unique CRISPR arrays
from per-sample spacer-mapping files and building count / frequency / absolute-abundance
matrices across a time series, together with a runnable demo dataset and its outputs.

## Dependencies

- Python ≥ 3.8
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/) and [seaborn](https://seaborn.pydata.org/) *(Muller plot only)*

```bash
pip install pandas numpy matplotlib seaborn
```

`crispr_array_count_matrix.py` needs only pandas + numpy; `plot_muller.py` additionally
uses matplotlib + seaborn. No other third-party packages are required.

## Contents

```
lineage_tracing_demo/
├── code/
│   ├── crispr_array_count_matrix.py     # unique-array extraction + count matrices (pandas + numpy)
│   └── plot_muller.py                   # Muller plot of recorded lineages (matplotlib + seaborn)
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
    ├── muller_plot.png                  # Muller plot from plot_muller.py
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
cd lineage_tracing_demo
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

## Muller plot (`plot_muller.py`)

Visualizes the recorded lineage dynamics over time from `unique_arrays_absfreq.csv`.

```bash
python3 code/plot_muller.py \
    --absfreq results/unique_arrays_absfreq.csv \
    --out results/muller_plot.png \
    --threshold 0.001 \
    --palette cubehelix_r
```

**Genealogy.** Spacers accumulate in reverse column order: `p1` is the **newest**
spacer (added at the left) and the right-most spacer is the **oldest**. Daughter
arrays therefore share their right-hand (suffix) spacers with the parent and differ
by a newly added left-hand spacer, so the parent of `(p1, p2, …, pL)` is `(p2, …, pL)`
(drop the newest spacer); the empty array is the common ancestor. For each array
passing `--threshold`, all of its suffixes are added as nodes so the tree stays fully
nested even when an intermediate array was not directly observed (inferred nodes have
zero population). Bands are stacked with each genotype bracketing its descendants.

Options: `--threshold` (abundance cutoff, default 0.1%), `--palette` (any seaborn
palette), `--smoothing` (gaussian std along the time axis; default 1, increase for
smoother band boundaries e.g. 2, 0 = none), `--no-normalize` (show absolute
abundance instead of normalizing each day to 1), `--day0-label`.

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

## Note

This demo uses a reduced set of time points (10 samples from a single replicate)
and therefore has a **lower temporal resolution than the analysis reported in the
paper**, which was based on a much denser time series. The demo is intended to make
the pipeline runnable end-to-end and to reproduce the analysis workflow; the
resulting lineage trajectories are correspondingly coarser than the published figures.

## Paper data (`paper_data/`)

Processed data underlying the lineage-tracing analysis in the paper. Raw sequencing
reads are deposited at NCBI SRA ([XXXXXXXXXXX](URL_TO_SRA)); this folder provides the
**unique CRISPR-array count tables** (per replicate) and the **Muller plots** derived
from them.

```
paper_data/
├── counts/
│   ├── rep1/ … rep4/
│   │   ├── unique_arrays_counts.csv.gz     array x sample raw read counts
│   │   ├── unique_arrays_normfreq.csv.gz   frequency within the array-bearing population
│   │   └── unique_arrays_absfreq.csv.gz    absolute abundance in the total cell population
│   │                                       (row array_id 0 = unexpanded ancestor)
└── figures/
    └── rep<N>_muller.png                    Muller plot of recorded lineages per replicate
```

Tables were generated with `lineage_tracing_demo/code/crispr_array_count_matrix.py`
(`--add-ancestor --no-long --gzip`) and Muller plots with
`lineage_tracing_demo/code/plot_muller.py` (`--threshold 0.001 --palette cubehelix_r
--smoothing 2`). Matrix rows are unique arrays defined by ordered spacer columns
`p1..pN` (`p1` = newest spacer, right-most = oldest); sample columns run day 0 → 100.
Unlike the low-resolution `lineage_tracing_demo`, these tables use the full,
densely-sampled time series analyzed in the paper.
