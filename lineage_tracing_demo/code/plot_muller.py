#!/usr/bin/env python3
"""
plot_muller.py

Draw a Muller plot of CRISPR-recorded lineages from the absolute-abundance matrix
produced by `crispr_array_count_matrix.py` (unique_arrays_absfreq.csv).

--------------------------------------------------------------------------------
GENEALOGY
--------------------------------------------------------------------------------
Spacers accumulate in reverse column order: p1 is the NEWEST spacer (added at the
leader-proximal, left end) and the right-most spacer is the OLDEST. Daughter arrays
therefore share the right-hand (suffix) spacers with their parent and differ by a
newly added left-hand spacer. Accordingly, the parent of an array with ordered
spacer tuple (p1, p2, ..., pL) is the array (p2, ..., pL) — i.e. drop the newest
(left-most) spacer. The empty array () is the common ancestor (unexpanded cells).

For every array reaching the abundance threshold, all of its suffixes are added as
nodes so the genealogy is fully nested even when an intermediate array was never
directly observed (such inferred nodes get zero population).

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
unique_arrays_absfreq.csv  (from crispr_array_count_matrix.py --prop-expanded --add-ancestor)
    index  : array_id
    columns: p1..pN  (spacer identities; NaN-padded)  +  one column per sample
             (sample columns hold absolute abundance in the total cell population;
              each column sums to ~1 including the ancestor/unexpanded row)

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  python3 plot_muller.py --absfreq results/unique_arrays_absfreq.csv --out muller.png
  python3 plot_muller.py --absfreq results/unique_arrays_absfreq.csv \
                         --threshold 0.001 --palette cubehelix_r --no-normalize

Requires: Python >= 3.8, pandas, numpy, matplotlib, seaborn
"""

import argparse
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams["pdf.fonttype"] = 42
warnings.filterwarnings("ignore")


# ----------------------------------------------------------------------------- 
# genealogy construction (suffix-based)
# ----------------------------------------------------------------------------- 
def day_of(sample_name, day0_label="beforegavage"):
    if day0_label and day0_label in sample_name:
        return 0
    m = re.search(r"day(\d+)", sample_name)
    return int(m.group(1)) if m else float("nan")


def build_genealogy(absfreq_path, threshold, day0_label):
    """Return (populations_df, adjacency_df, ancestor_id, n_nodes)."""
    a = pd.read_csv(absfreq_path, index_col=0)
    pcols = [c for c in a.columns if c.startswith("p")]
    scols = [c for c in a.columns if not c.startswith("p")]
    days = [day_of(c, day0_label) for c in scols]

    # each observed array -> its ordered spacer tuple and abundance vector
    a["tuple"] = a.apply(lambda r: tuple(int(r[c]) for c in pcols if pd.notna(r[c])), axis=1)
    obs = {t: np.nan_to_num(a.loc[i, scols].to_numpy(float)) for i, t in zip(a.index, a["tuple"])}

    # seeds = arrays reaching the abundance threshold at any timepoint
    seeds = [t for t in obs if np.max(obs[t]) >= threshold]

    # nodes = seeds + ALL their suffixes (inferred intermediates), + ancestor ()
    nodes = {()}
    for t in seeds:
        for k in range(len(t) + 1):
            nodes.add(t[k:])
    nodes = sorted(nodes, key=lambda t: (len(t), t))
    node_id = {t: i for i, t in enumerate(nodes)}          # ancestor () -> id 0

    # parent = drop newest (left-most) spacer -> suffix t[1:]
    adjacency_df = pd.DataFrame(
        [{"Parent": node_id[t[1:]], "Identity": node_id[t]} for t in nodes if len(t) > 0]
    )

    zero = np.zeros(len(scols))
    rows = [(node_id[t], d, val) for t in nodes for d, val in zip(days, obs.get(t, zero))]
    populations_df = pd.DataFrame(rows, columns=["Identity", "Generation", "Population"]).astype(
        {"Identity": int, "Generation": int, "Population": float}
    )
    return populations_df, adjacency_df, 0, len(nodes)


# ----------------------------------------------------------------------------- 
# Muller stacking (self-contained; same nesting convention as ggmuller/pymuller)
# ----------------------------------------------------------------------------- 
def _strain_order(adjacency_df):
    """Depth-first order in which each genotype brackets its descendants."""
    children = adjacency_df.groupby("Parent")["Identity"].apply(lambda x: sorted(x))

    def inner(idt):
        kids = children.get(idt, [])
        if not kids:
            return [idt, idt]
        return [idt] + sum((inner(c) for c in kids), []) + [idt]

    order, seen = [], set()
    identities = set(adjacency_df["Identity"]) | set(adjacency_df["Parent"])
    for s in sorted(identities):
        if s not in seen:
            sub = inner(s)
            order += sub
            seen.update(sub)
    return np.array(order)


def _y_values(populations_df, adjacency_df, smoothing_std):
    """Stacked half-heights per genotype; children nested inside parents."""
    ordering = _strain_order(adjacency_df)
    pop_max = populations_df.groupby("Generation")["Population"].sum().max()
    span = populations_df["Generation"].max() - populations_df["Generation"].min()

    pivot = populations_df.pivot(index="Generation", columns="Identity", values="Population").sort_index()
    if smoothing_std and smoothing_std > 0:
        pivot = pivot.rolling(max(int(span), 1), 1, True, "gaussian").mean(std=smoothing_std)
    pivot = pivot.clip(0, pop_max)

    Y = pivot[ordering] / 2.0
    keep = [0]
    for i, c in enumerate(Y.columns[1:], 1):
        if c == Y.columns[i - 1]:      # a genotype's two brackets are adjacent (leaf) -> merge
            Y.iloc[:, i] *= 2
            keep.pop()
        keep.append(i)
    return Y.iloc[:, keep]


# ----------------------------------------------------------------------------- 
# plot
# ----------------------------------------------------------------------------- 
def plot_muller(populations_df, adjacency_df, ancestor_id, normalize,
                palette, smoothing_std, ax=None, title=None):
    pop = populations_df.copy()
    if normalize:
        pop["Population"] = pop.groupby("Generation")["Population"].transform(
            lambda x: x / x.sum() if x.sum() else x
        )

    Y = _y_values(pop, adjacency_df, smoothing_std)
    x = Y.index.values
    final_order = Y.columns.values

    color_labels = list(dict.fromkeys(final_order))            # unique arrays, stack order
    colors_list = sns.color_palette(palette, len(color_labels))
    cmap_id = dict(zip(color_labels, colors_list))
    colors = [(1, 1, 1) if i == ancestor_id else cmap_id[i] for i in final_order]

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 6))
    ax.stackplot(x, Y.to_numpy().T, colors=colors, linewidth=0)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1 if normalize else None)
    ax.set_xlabel("Day", fontsize=12)
    ax.set_ylabel("Relative abundance" + (" (normalized)" if normalize else ""), fontsize=12)
    ax.tick_params(direction="out", length=5, width=1, colors="black")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if title:
        ax.set_title(title, fontsize=10.5)
    return ax


def main(argv=None):
    ap = argparse.ArgumentParser(description="Muller plot of CRISPR-recorded lineages.")
    ap.add_argument("--absfreq", required=True, help="unique_arrays_absfreq.csv from the count-matrix step.")
    ap.add_argument("--out", default="muller.png", help="output figure (.png/.pdf/.svg).")
    ap.add_argument("--threshold", type=float, default=0.001,
                    help="keep arrays reaching this abundance at any timepoint (default 0.001 = 0.1%%).")
    ap.add_argument("--palette", default="cubehelix_r", help="seaborn color palette name.")
    ap.add_argument("--smoothing", type=float, default=1.0,
                    help="gaussian smoothing std along the time axis (default 1; "
                         "increase for smoother band boundaries, e.g. 2; set 0 to disable).")
    ap.add_argument("--no-normalize", action="store_true", help="show absolute abundance instead of normalizing to 1.")
    ap.add_argument("--day0-label", default="beforegavage", help="substring marking day 0.")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)

    pops, adj, anc, n = build_genealogy(args.absfreq, args.threshold, args.day0_label)
    title = (f"Muller plot of CRISPR-recorded lineages\n"
             f"{n} array nodes (\u2265{args.threshold:.3%} + inferred suffixes) | palette: {args.palette}")
    ax = plot_muller(pops, adj, anc, not args.no_normalize, args.palette, args.smoothing, title=title)
    plt.tight_layout()
    plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"saved -> {args.out}  ({n} nodes, threshold {args.threshold})")


if __name__ == "__main__":
    main()
