# plot_helpers.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2

def qqplot(pvalues, figsize=(6, 6)):
    """
    QQ plot for GWAS P-values.
    Also reports genomic inflation factor lambda_GC.
    """

    p = np.asarray(pvalues, dtype=float)

    p = p[
        np.isfinite(p)
        & (p > 0)
        & (p <= 1)
    ]

    p = np.sort(p)

    n = len(p)

    expected = -np.log10(
        (np.arange(1, n + 1) - 0.5) / n
    )

    observed = -np.log10(p)

    # Genomic inflation factor
    chi2_obs = chi2.isf(p, df=1)
    lambda_gc = np.median(chi2_obs) / chi2.ppf(0.5, df=1)

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(expected, observed, s=12, alpha=0.7)

    lim = max(expected.max(), observed.max())
    ax.plot([0, lim], [0, lim], linestyle="--")

    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"QQ plot   λGC = {lambda_gc:.3f}")

    fig.tight_layout()
    return fig, ax


def manhattan(
    gwas,
    chrom="#CHROM",
    pos="POS",
    p="P",
    snp="ID",
    genomewide=5e-8,
    suggestive=1e-5,
    annotate_top=0,
    figsize=(12, 5),
):
    """
    Manhattan plot from a PLINK2 --glm result.

    Required columns:
        #CHROM
        POS
        P

    Optionally labels the N most significant variants.
    """

    df = gwas[[chrom, pos, p, snp]].copy()

    df = df.dropna(subset=[chrom, pos, p])

    df[p] = pd.to_numeric(df[p], errors="coerce")
    df[pos] = pd.to_numeric(df[pos], errors="coerce")
    df[chrom] = pd.to_numeric(df[chrom], errors="coerce")

    df = df.dropna()

    df = df[
        (df[p] > 0)
        & (df[p] <= 1)
    ]

    df = df.sort_values([chrom, pos])

    df["minus_log10_p"] = -np.log10(df[p])

    # Construct cumulative genomic coordinates
    chromosome_offsets = {}
    offset = 0
    ticks = []
    labels = []

    xs = np.empty(len(df))

    fig, ax = plt.subplots(figsize=figsize)

    for i, (chr_name, d) in enumerate(df.groupby(chrom, sort=True)):

        idx = d.index

        chromosome_offsets[chr_name] = offset

        x = d[pos].values + offset
        xs[df.index.get_indexer(idx)] = x

        ax.scatter(
            x,
            d["minus_log10_p"],
            s=12,
            alpha=0.8
        )

        ticks.append((x.min() + x.max()) / 2)
        labels.append(str(int(chr_name)))

        offset = x.max() + 5_000_000

    df["x"] = xs

    # Significance thresholds
    if genomewide is not None:
        ax.axhline(
            -np.log10(genomewide),
            linestyle="--",
            linewidth=1
        )

    if suggestive is not None:
        ax.axhline(
            -np.log10(suggestive),
            linestyle=":",
            linewidth=1
        )

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(P)")
    ax.set_title("Manhattan plot")

    # Label strongest signals
    if annotate_top > 0:
        top = df.nsmallest(annotate_top, p)

        for _, row in top.iterrows():
            ax.annotate(
                str(row[snp]),
                (row["x"], row["minus_log10_p"]),
                xytext=(0, 5),
                textcoords="offset points",
                fontsize=8,
                rotation=45
            )

    fig.tight_layout()
    return fig, ax


def regional_plot(
    gwas,
    chromosome,
    start,
    end,
    chrom="#CHROM",
    pos="POS",
    p="P",
    snp="ID",
    annotate_top=3,
    figsize=(9, 5),
):
    """
    Simple regional association plot.

    chromosome:
        chromosome number
    start/end:
        genomic coordinates in base pairs
    """

    df = gwas.copy()

    df[chrom] = pd.to_numeric(df[chrom], errors="coerce")
    df[pos] = pd.to_numeric(df[pos], errors="coerce")
    df[p] = pd.to_numeric(df[p], errors="coerce")

    df = df[
        (df[chrom] == chromosome)
        & (df[pos] >= start)
        & (df[pos] <= end)
        & (df[p] > 0)
    ].copy()

    df["minus_log10_p"] = -np.log10(df[p])

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        df[pos] / 1e6,
        df["minus_log10_p"],
        s=20,
        alpha=0.8
    )

    if annotate_top > 0 and len(df):
        top = df.nsmallest(annotate_top, p)

        for _, row in top.iterrows():
            ax.annotate(
                str(row[snp]),
                (row[pos] / 1e6, row["minus_log10_p"]),
                xytext=(0, 5),
                textcoords="offset points",
                fontsize=8
            )

    ax.set_xlabel(f"Chromosome {chromosome} position (Mb)")
    ax.set_ylabel("-log10(P)")
    ax.set_title(
        f"Chr {chromosome}: "
        f"{start / 1e6:.1f}-{end / 1e6:.1f} Mb"
    )

    fig.tight_layout()
    return fig, ax
