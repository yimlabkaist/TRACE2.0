#!/usr/bin/env python3
"""
crispr_array_count_matrix.py

Extract unique CRISPR arrays from per-sample spacer-mapping files and build a
count matrix across samples (for TRACE gut lineage-tracing analysis).

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
One text file per sample, comma-separated, no header, e.g. `*_uniq.txt`:

    a1l1p1,ec257-v3r,4260750,4260782
    a8l2p1,ec257-v3r,501296,501264
    a8l2p2,ec257-v3r,3168840,3168808
    ...

Column 1 : spacer ID = a{ARRAY}l{LEN}p{POS}
           - the substring before 'p'  (e.g. "a8l2") identifies one observed
             array (one read); multiple 'p' entries are the ordered spacers
             of that array (p1 = oldest / distal, p2 = next, ...).
Column 2 : reference name (ignored for counting).
Column 3 : spacer mapping START coordinate (defines spacer identity).
Column 4 : spacer mapping END coordinate (start>end => reverse strand).

--------------------------------------------------------------------------------
METHOD
--------------------------------------------------------------------------------
1. For each sample, group spacers by their array key (ID before 'p') and order
   them by position index p1,p2,... to reconstruct each observed array as an
   ordered tuple of spacer identities.
2. A spacer's identity is its mapping START coordinate by default
   (--spacer-identity start), matching the original TRACE analysis; use
   'start_end' for strand-aware identity based on the (start,end) pair.
3. Collapse identical array tuples within the sample and count them
   -> per-sample unique-array count table (columns p1..pN, counts, norm_counts).
4. Take the union of unique arrays across all samples (the reference array set),
   assign each a stable integer identity, and assemble two matrices
   (rows = unique arrays, columns = samples):
       - raw read counts
       - within-sample normalized frequencies
   Optionally prepend an empty "ancestor" array (row of NaNs) for downstream
   genealogy / Muller-plot reconstruction (--add-ancestor).

--------------------------------------------------------------------------------
OUTPUT (into --outdir)
--------------------------------------------------------------------------------
  per_sample/<sample>_arrays.csv        unique arrays + counts for each sample
  unique_arrays_counts.csv              array x sample raw count matrix
  unique_arrays_normfreq.csv            array x sample normalized-frequency matrix
                                        (within the array-bearing population)
  unique_arrays_absfreq.csv             array x sample absolute abundance in the
                                        TOTAL cell population (only with --prop-expanded):
                                        normfreq x (prop. expanded) per sample
  unique_arrays_long.csv                tidy long table (array_id, sample, day, count,
                                        frequency, abs_frequency)

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python3 crispr_array_count_matrix.py --input "Ec257-day*-rep4_uniq.txt" --outdir results
  python3 crispr_array_count_matrix.py --input ./data --outdir results --max-spacers 5
  python3 crispr_array_count_matrix.py --input "*.txt" --outdir results --spacer-identity start_end --add-ancestor

Requires: Python 3.8+, pandas, numpy  (no other dependencies)
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------- 
# parsing
# ----------------------------------------------------------------------------- 
def sample_name_from_path(path):
    """Derive a sample name from the filename (portion before '_uniq')."""
    base = os.path.basename(path)
    return base.split("_uniq")[0].rsplit(".", 1)[0] if "_uniq" in base else base.rsplit(".", 1)[0]


def day_of(sample_name, day0_label="beforegavage"):
    """Return the integer day for a sample (day0_label -> 0), else NaN."""
    if day0_label and day0_label in sample_name:
        return 0
    m = re.search(r"day(\d+)", sample_name)
    return int(m.group(1)) if m else float("nan")


def day_sort_key(sample_name, day0_label="beforegavage"):
    """Natural sort key: day0_label first (day 0), then by embedded day number."""
    d = day_of(sample_name, day0_label)
    return (0, d) if d == d else (1, sample_name)   # d==d is False only for NaN


def parse_sample_file(path, spacer_identity="start"):
    """
    Parse one sample file into a list of ordered array signatures.

    Returns a list of tuples, one per observed array (read); each tuple is the
    ordered sequence of spacer identities (int start, or 'start:end' strings).
    """
    arrays = defaultdict(list)   # array_key -> list of (pos_index, start, end)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 4:
                continue
            spacer_id, _ref, start, end = parts[0], parts[1], parts[2], parts[3]
            try:
                start, end = int(start), int(end)
            except ValueError:
                continue
            array_key = spacer_id.split("p")[0]           # e.g. "a8l2"
            m = re.search(r"p(\d+)\s*$", spacer_id)        # spacer position index
            pos_index = int(m.group(1)) if m else 1
            arrays[array_key].append((pos_index, start, end))

    signatures = []
    for _key, spacers in arrays.items():
        spacers.sort(key=lambda t: t[0])                   # order by p1,p2,...
        if spacer_identity == "start_end":
            sig = tuple(f"{s}:{e}" for _p, s, e in spacers)
        else:                                              # 'start' (default)
            sig = tuple(s for _p, s, _e in spacers)
        signatures.append(sig)
    return signatures


# ----------------------------------------------------------------------------- 
# per-sample unique arrays
# ----------------------------------------------------------------------------- 
def unique_array_table(signatures, max_spacers):
    """Collapse identical array signatures and count them -> DataFrame."""
    counts = defaultdict(int)
    for sig in signatures:
        counts[sig] += 1

    cols = [f"p{i}" for i in range(1, max_spacers + 1)]
    rows = []
    for sig, c in counts.items():
        padded = list(sig[:max_spacers]) + [np.nan] * (max_spacers - len(sig[:max_spacers]))
        rows.append(padded + [c])
    df = pd.DataFrame(rows, columns=cols + ["counts"])
    df = df.sort_values("counts", ascending=False).reset_index(drop=True)
    df["norm_counts"] = df["counts"] / df["counts"].sum()
    return df


# ----------------------------------------------------------------------------- 
# combined matrix
# ----------------------------------------------------------------------------- 
def build_matrix(per_sample, max_spacers, add_ancestor):
    """
    per_sample: dict sample_name -> per-sample unique-array DataFrame.
    Returns (counts_matrix, freq_matrix) indexed by a stable array identity,
    with the array definition columns (p1..pN) retained.
    """
    key_cols = [f"p{i}" for i in range(1, max_spacers + 1)]

    # union of all unique arrays (the reference set)
    ref = pd.concat([df[key_cols] for df in per_sample.values()], ignore_index=True)
    ref = ref.drop_duplicates().reset_index(drop=True)

    if add_ancestor:  # empty (all-NaN) array = common ancestor, identity 0
        ref = pd.concat([pd.DataFrame([[np.nan] * max_spacers], columns=key_cols), ref],
                        ignore_index=True)

    ref.index.name = "array_id"

    counts_mat = ref.copy()
    freq_mat = ref.copy()
    for sample, df in per_sample.items():
        merged = ref.merge(df, on=key_cols, how="left")
        counts_mat[sample] = merged["counts"].to_numpy()
        freq_mat[sample] = merged["norm_counts"].to_numpy()

    sample_cols = list(per_sample.keys())
    counts_mat[sample_cols] = counts_mat[sample_cols].fillna(0).astype(int)
    freq_mat[sample_cols] = freq_mat[sample_cols].fillna(0.0)
    return counts_mat, freq_mat


def apply_prop_expanded(freq_mat, sample_cols, max_spacers, prop_path, prop_col, add_ancestor):
    """
    Scale within-population normalized frequencies by each sample's
    'prop. expanded' (fraction of the total cell population that is array-bearing),
    giving the absolute abundance of each array in the TOTAL cell population.

        abs_freq[array, sample] = norm_freq[array, sample] * prop_expanded[sample]

    If add_ancestor, the ancestor row (array_id 0, empty array) receives the
    unexpanded fraction (1 - prop_expanded), so each sample column sums to ~1
    over the whole population (arrays + ancestor).
    """
    key_cols = [f"p{i}" for i in range(1, max_spacers + 1)]
    prop = pd.read_csv(prop_path)
    if "file" not in prop.columns or prop_col not in prop.columns:
        sys.exit(f"'{prop_path}' must contain a 'file' column and a '{prop_col}' column.")
    prop_map = dict(zip(prop["file"], prop[prop_col]))

    missing = [s for s in sample_cols if s not in prop_map]
    if missing:
        print(f"  [warn] no 'prop. expanded' value for: {', '.join(missing)} "
              f"(their absolute abundance will be 0)")

    abs_mat = freq_mat[key_cols].copy()
    for s in sample_cols:
        pe = prop_map.get(s, 0.0)
        abs_mat[s] = freq_mat[s].to_numpy() * pe

    if add_ancestor and 0 in abs_mat.index:
        # ancestor row = unexpanded fraction of the population
        for s in sample_cols:
            abs_mat.loc[0, s] = 1.0 - prop_map.get(s, 1.0)
    return abs_mat


def to_long(counts_mat, freq_mat, abs_mat, sample_cols, max_spacers, day0_label="beforegavage"):
    """Tidy long table: one row per (array, sample)."""
    key_cols = [f"p{i}" for i in range(1, max_spacers + 1)]
    c = counts_mat.reset_index().melt(id_vars=["array_id"] + key_cols,
                                      value_vars=sample_cols,
                                      var_name="sample", value_name="count")
    f = freq_mat.reset_index().melt(id_vars=["array_id"],
                                    value_vars=sample_cols,
                                    var_name="sample", value_name="frequency")
    long = c.merge(f, on=["array_id", "sample"])
    if abs_mat is not None:
        a = abs_mat.reset_index().melt(id_vars=["array_id"],
                                       value_vars=sample_cols,
                                       var_name="sample", value_name="abs_frequency")
        long = long.merge(a, on=["array_id", "sample"])
    long["day"] = long["sample"].map(lambda s: day_of(s, day0_label))
    return long


# ----------------------------------------------------------------------------- 
# main
# ----------------------------------------------------------------------------- 
def resolve_inputs(pattern):
    if os.path.isdir(pattern):
        files = glob.glob(os.path.join(pattern, "*_uniq.txt")) or glob.glob(os.path.join(pattern, "*.txt"))
    else:
        files = glob.glob(pattern)
    return sorted(files)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Extract unique CRISPR arrays and build a count matrix.")
    ap.add_argument("--input", required=True,
                    help="Directory, or glob pattern, of per-sample *_uniq.txt files.")
    ap.add_argument("--outdir", default="crispr_array_results", help="Output directory.")
    ap.add_argument("--max-spacers", type=int, default=5,
                    help="Number of spacer columns p1..pN to keep (default 5).")
    ap.add_argument("--spacer-identity", choices=["start", "start_end"], default="start",
                    help="Spacer identity: mapping start only (default) or strand-aware start:end.")
    ap.add_argument("--add-ancestor", action="store_true",
                    help="Prepend an empty ancestor array (identity 0) for genealogy reconstruction.")
    ap.add_argument("--prop-expanded", default=None,
                    help="Optional CSV with columns 'file' and 'prop. expanded' to scale "
                         "normalized frequencies into absolute abundance in the total cell population.")
    ap.add_argument("--prop-col", default="prop. expanded",
                    help="Column name in --prop-expanded holding the expansion fraction.")
    ap.add_argument("--day0-label", default="beforegavage",
                    help="Substring in a sample name that marks day 0 (default 'beforegavage').")
    args = ap.parse_args(argv)

    files = resolve_inputs(args.input)
    if not files:
        sys.exit(f"No input files matched: {args.input}")
    files.sort(key=lambda p: day_sort_key(sample_name_from_path(p), args.day0_label))

    os.makedirs(os.path.join(args.outdir, "per_sample"), exist_ok=True)

    per_sample = {}
    print(f"Parsing {len(files)} sample file(s):")
    for path in files:
        name = sample_name_from_path(path)
        sigs = parse_sample_file(path, spacer_identity=args.spacer_identity)
        tbl = unique_array_table(sigs, args.max_spacers)
        per_sample[name] = tbl
        tbl.to_csv(os.path.join(args.outdir, "per_sample", f"{name}_arrays.csv"), index=False)
        print(f"  {name:<28} reads={len(sigs):>8}  unique_arrays={len(tbl):>7}")

    counts_mat, freq_mat = build_matrix(per_sample, args.max_spacers, args.add_ancestor)
    sample_cols = list(per_sample.keys())

    abs_mat = None
    if args.prop_expanded:
        abs_mat = apply_prop_expanded(freq_mat, sample_cols, args.max_spacers,
                                      args.prop_expanded, args.prop_col, args.add_ancestor)

    counts_path = os.path.join(args.outdir, "unique_arrays_counts.csv")
    freq_path = os.path.join(args.outdir, "unique_arrays_normfreq.csv")
    long_path = os.path.join(args.outdir, "unique_arrays_long.csv")
    counts_mat.to_csv(counts_path)
    freq_mat.to_csv(freq_path)
    written = [counts_path, freq_path]
    if abs_mat is not None:
        abs_path = os.path.join(args.outdir, "unique_arrays_absfreq.csv")
        abs_mat.to_csv(abs_path)
        written.append(abs_path)
    to_long(counts_mat, freq_mat, abs_mat, sample_cols, args.max_spacers,
            args.day0_label).to_csv(long_path, index=False)
    written.append(long_path)

    print(f"\nTotal unique arrays across samples: {len(counts_mat)}")
    print("Wrote:")
    for p in written:
        print(f"  {p}")
    print(f"  {os.path.join(args.outdir, 'per_sample')}/<sample>_arrays.csv")


if __name__ == "__main__":
    main()
