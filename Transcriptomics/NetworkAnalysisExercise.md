[⬇️ Download the Jupyter notebook](./NetworkAnalysisExercise.ipynb?raw=1) (Right Click -> Save Link As..)


# Gene Co-expression Network Analysis with `igraph`

**Workflow:** TPM expression → co-expression network → graph properties → centrality → community detection → functional enrichment

This exercise continues the myocardial infarction (**MI**) transcriptomics analysis. Here, instead of testing individual genes for differential expression, we use **gene co-expression networks** to ask which genes show coordinated expression patterns and which genes occupy important positions in the resulting network.

The dataset contains **MI** and **Healthy** samples. The co-expression network is constructed from **TPM expression values** across the included samples.

> **Important distinction:** A co-expression network is not the same as an **MI vs Healthy differential-expression analysis**. An edge indicates that two genes have similar expression patterns across samples; it does not by itself mean that either gene is differentially expressed between MI and Healthy samples.

## Learning objectives

By the end of this notebook, you should be able to:

- understand how a gene co-expression network can be generated from **TPM** data;
- load an edge list into **igraph** and inspect basic network properties;
- identify highly connected **hub genes** and genes with high **betweenness**;
- detect gene modules using the **Leiden** community-detection algorithm; and
- use functional enrichment to help interpret biologically interesting modules.

---

## 1. Setup: check required folders and data

The exercise expects a `data/` folder containing the input files and a `Results/` folder for generated outputs.

If you downloaded only this notebook, the helper script `checkdir.py` & `checknet.py` will:

- create `Results/` if it does not already exist; and
- download the course `data/` folder from GitHub if it is missing.
- download the course `Network/` folder from GitHub if it is missing.

`checkdir.py` & `checknet.py` should be located in the **same folder as this notebook**.

`checkdir.py` should be located in the **same folder as this notebook**.

```python
%run checkdir.py
```

```python
%run checknet.py
```

## 2. Load libraries and define the co-expression helper function

The main packages used in this exercise are:

- **pandas / NumPy** — table and numerical operations;
- **SciPy** — Spearman correlation;
- **statsmodels** — multiple-testing correction;
- **igraph** — network analysis and visualization;
- **GSEApy** — functional enrichment; and

The `coexpression_generation()` function below calculates pairwise Spearman correlations between genes, adjusts the correlation p-values using Benjamini–Hochberg FDR, and returns the significant gene–gene relationships as an edge list.

```python
%matplotlib inline

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import gseapy as gp
import igraph as ig

# Plot settings that make downstream figure editing easier.
sns.set(font="Arial")
sns.set_style("white")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

os.environ["KMP_DUPLICATE_LIB_OK"] = "True"


def coexpression_generation(tpm, padj_thr=0.05):
    """Generate a significant gene co-expression edge list from TPM data."""
    print("Calculating correlations...")
    temp = spearmanr(tpm.T)

    corr = pd.DataFrame(
        temp[0],
        columns=list(tpm.index),
        index=list(tpm.index),
    )
    pval = pd.DataFrame(
        temp[1],
        columns=list(tpm.index),
        index=list(tpm.index),
    )

    print("Keeping one copy of each gene pair...")
    upper_triangle = np.triu(np.ones(corr.shape, dtype=bool))
    corr = corr.where(upper_triangle)
    pval = pval.where(upper_triangle)

    print("Converting correlation matrices to an edge list...")
    corr2 = corr.unstack().reset_index(name="weight")
    pval2 = pval.unstack().reset_index(name="pval")
    res = corr2.merge(pval2, on=["level_0", "level_1"])
    res = res[res["level_0"] != res["level_1"]].dropna()

    print("Adjusting p-values...")
    res["padj"] = multipletests(res["pval"], method="fdr_bh")[1]
    res = res[res["padj"] < padj_thr].reset_index(drop=True)

    print("Done!")
    return res[["level_0", "level_1", "weight"]]
```

---

## 3. Optional — Generate the co-expression network from TPM  (DO NOT DO THIS ON POSIT CLOUD)

A **pre-computed co-expression network** is provided with the exercise so that we can spend most of the lab on downstream network analysis. You can therefore skip this section during the main practical.

It is still useful to understand how the network was generated.

### Why use TPM here?

For **DESeq2**, we used raw counts because the statistical model requires count data. Co-expression analysis asks a different question: whether genes show similar expression patterns across samples. For this exercise, we use **TPM** as the expression measure for calculating correlations.

The network is generated across the available **MI and Healthy samples**. Therefore, the resulting edges represent co-expression across the dataset as a whole; they are **not** direct MI-vs-Healthy tests.

### Optional Step 3.1 — Load TPM expression

```python
tpm = pd.read_csv("data/tpm.txt", sep="\t", index_col=0)

# Remove genes with zero TPM in every sample.
tpm = tpm[tpm.sum(axis=1) > 0]
```

### Optional Step 3.2 — Generate significant positive correlations

```python
coexp_network = coexpression_generation(tpm)

# Keep only positively correlated gene pairs.
coexp_network = coexp_network[coexp_network["weight"] > 0] # We use :99 cut off for this exercise to make it smaller
```

The resulting edge list has three columns:

- `level_0` — first gene in the pair;
- `level_1` — second gene in the pair; and
- `weight` — Spearman correlation coefficient.

> **Why only one triangle of the matrix?** Correlation is symmetric: the correlation of gene A with gene B is the same as gene B with gene A. Keeping only one triangle avoids storing every edge twice.

The helper function is adapted from a co-expression-network workflow used in the following study:  
https://pubmed.ncbi.nlm.nih.gov/37038090/

---

## 4. Load the co-expression network into `igraph`

For the remainder of the exercise, we use the **pre-computed network** supplied in the `data/` folder.

### Step 4.1 — Load the edge list

```python
coexp_network = pd.read_csv("Network/coexpression_network.txt", sep="\t")
coexp_network.head()
```

### Step 4.2 — Create the graph

The edge list contains a correlation `weight`, but in this introductory exercise we will analyze the network as **unweighted**. In other words, an edge is treated simply as present or absent.

`Graph.TupleList()` converts each gene pair into an `igraph` edge.

```python
g = ig.Graph.TupleList(
    zip(coexp_network["level_0"], coexp_network["level_1"])
)
```

### Observe Network Properties

Useful definitions:

- **Vertex / node:** here, a gene.
- **Edge / link:** a significant positive co-expression relationship between two genes.
- **Diameter:** the length of the longest shortest path between connected nodes in the network.
- **Density:** the fraction of all possible gene–gene edges that are actually present.

```python
print("Nodes:", g.vcount())
print("Edges:", g.ecount())
```

```python
print("Diameter:", g.diameter())
print("Density:", g.density())
```

### Step 4.3 — Choose a network layout (SKIP in Posit Cloud or if it takes too long on your own computer)

A network layout controls where nodes are placed for visualization. It does **not** change the network itself.

Here we let `igraph` choose an appropriate layout automatically.

**Reference:** [igraph visualization tutorial](https://igraph.org/python/tutorial/latest/visualisation.html)

```python
if g.vcount() <= 500:
    layout = g.layout("circle")
else:
    layout = None
    print("Network is large — skipping full-network layout.")
```

### Step 4.4 — Visualize the network (SKIP in Posit Cloud or if it takes too long on your own computer)

Large biological networks can look dense or tangled. The purpose of this first plot is mainly to inspect the overall structure before calculating quantitative network properties.

```python
fig, ax = plt.subplots()

ig.plot(
    g,
    layout=layout,
    target=ax,
    vertex_size=5,
    vertex_label=None,
    vertex_color="red",
    edge_color="grey",
    edge_width=0.5,
)

fig.set_size_inches(20, 20)
```

### Checkpoint 1 — Basic network properties

> **Q1.1.** How many **vertices (nodes)** and **edges (links)** are present in the network?  
> **Q1.2.** What are the network **diameter** and **density**? What do these two quantities represent biologically or structurally?

---

## 5. Centrality analysis

Centrality measures help identify genes that occupy important positions in a co-expression network. Different measures capture different meanings of **importance**:

- **Degree:** number of direct connections a gene has. Genes with unusually high degree are often called **hubs**.
- **Betweenness centrality:** how often a gene lies on shortest paths between other genes. High-betweenness genes can act as **bridges** between regions or modules.
- **Closeness centrality:** how close a gene is, on average, to all other reachable genes.
- **Eigenvector centrality:** gives higher scores to genes that are connected to other highly connected genes.

These measures describe **network topology**; they do not by themselves establish biological causality.

**Reference:** [python-igraph documentation](https://python.igraph.org/en/stable/)

### Checkpoint 2 — Identify important genes

A convenient way to work with `igraph` centrality outputs is to convert the returned list to a pandas `Series` using the gene names as the index.

> **Q2.1.** Which genes have the highest **degree**? Report the top 10.  
> **Q2.2.** Which genes have the highest **betweenness centrality**? Report the top 10.  
> **Q2.3.** Calculate one additional centrality measure and compare its ranking with degree and betweenness.  
> **Q2.4.** Do the highest-degree genes also rank highly by betweenness or your other centrality measure? Why might the rankings differ?  
> **Q2.5.** (Optional) Compare the biological functions of hub genes and bridge genes. Do they appear to represent different processes?

For Q2.5, you can use [Enrichr](https://maayanlab.cloud/Enrichr/) or the **GSEApy** workflow introduced in the transcriptomics exercise.

```python
degree_scores = pd.Series(
    g.degree(),
    index=g.vs["name"],
).sort_values(ascending=False)

degree_scores.head(10)
```

In this reference network, several genes may have exactly the same maximum degree. Instead of hard-coding a particular number of genes, we can define the **hub set** as all genes tied for the highest degree.

```python
hubs = degree_scores[degree_scores == degree_scores.max()]
hubs
```

```python
betweenness_scores = pd.Series(
    g.betweenness(),
    index=g.vs["name"],
).sort_values(ascending=False)

betweenness_scores.head(10)
```

### Optional: Use degree as node size (SKIP in Posit Cloud or if it takes too long on your own computer)

One simple way to highlight hub genes visually is to scale node size by degree.

```python
fig, ax = plt.subplots()

ig.plot(
    g,
    layout=layout,
    target=ax,
    vertex_size=g.degree(),
)

fig.set_size_inches(20, 20)
```

---

## 6. Clustering and module detection

Gene co-expression networks often contain **modules**: groups of genes that are more densely connected to each other than to the rest of the network. Such modules can represent shared biological pathways, regulatory programs, or coordinated responses to a biological condition such as MI.

### Common approaches

- **Community-detection algorithms**
  - *Louvain / Leiden* — optimize a community-quality objective such as modularity.
  - *Walktrap* — uses random walks to identify densely connected subgraphs.
  - *Fastgreedy* — greedily optimizes modularity.
- **Hierarchical clustering** — can be applied to adjacency or topological-overlap matrices.
- **WGCNA** — identifies modules using hierarchical clustering and dynamic tree cutting.

### Key concepts

- **Module:** a set of tightly connected/co-expressed genes.
- **Module size:** number of genes assigned to the module.
- **Module eigengene:** in WGCNA-style analyses, the first principal component of a module's expression matrix, used as a summary expression profile.

### Leiden algorithm

In this exercise we use **Leiden community detection** with a modularity objective.

> **Note:** Community labels such as `0`, `1`, `2`, ... are arbitrary identifiers. The numeric label of the largest module can change, so downstream code should identify the largest module programmatically rather than hard-coding a cluster number.

**Reference:** [igraph `community_leiden`](https://python.igraph.org/en/stable/api/igraph.Graph.html#community_leiden)

```python
clustering = g.community_leiden(
    objective_function="modularity",
    n_iterations=-1,
    resolution=0.2,
)
```

### Step 6.1 — Convert module membership to a pandas Series

`clustering.membership` contains one module label for each network node. Converting it to a pandas `Series` makes it easier to count, filter, and match modules to gene names.

```python
memberships = pd.Series(
    clustering.membership,
    index=g.vs["name"],
)

memberships.head()
```

### Checkpoint 3 — Explore the detected modules

> **Q3.1.** How many modules are detected?  
> **Q3.2.** Which module is the largest? Which is the smallest?  
> **Q3.3.** Which module contains the largest number of hub genes identified above?  
> **Q3.4.** For the largest module, what biological functions might it represent? Visualize the top 10 enriched biological processes using Enrichr or GSEApy.  
> **Q3.5 — Optional.** Can you identify the hub gene(s) within the two largest modules? *Hint: create module-specific subgraphs with `.subgraph()`.*

*Try to answer these questions first, then use the short walkthrough below to check your approach.*

### Step 6.2 — Count the detected modules

Each Leiden module has a unique label. `len(clustering)` therefore gives the total number of detected modules and answers **Q3.1**.

```python
# Q3.1 — Number of modules
len(clustering)
```

### Step 6.3 — Compare module sizes

Count how many genes belong to each module. Sorting the counts from largest to smallest lets us identify the **largest** and **smallest** modules for **Q3.2**.

```python
# Q3.2 — Module sizes, ranked from largest to smallest
module_sizes = memberships.value_counts()
module_sizes
```

### Step 6.4 — Find which modules contain the hub genes

Match the hub-gene names to their Leiden memberships, then count how many hubs occur in each module. The module with the largest count answers **Q3.3**.

```python
# Q3.3 — Count high-degree hub genes in each module
hub_modules = memberships.reindex(hubs.index)
hub_modules.value_counts()
```

### Step 6.5 — Functional enrichment of the largest module

To address **Q3.4**, first identify the largest module. Leiden module numbers are arbitrary, so use `.idxmax()` instead of assuming a particular module number will always be the largest.

```python
largest_cluster_id = module_sizes.idxmax()
largest_cluster = memberships[memberships == largest_cluster_id]

print(f"Largest module: {largest_cluster_id}")
print(f"Number of genes: {largest_cluster.shape[0]}")
```

```python
enr = gp.enrichr(
    gene_list=largest_cluster.index.tolist(),
    gene_sets="GO_Biological_Process_2026",
    organism="human",
    outdir=None,
)
```

```python
ax = gp.dotplot(
    enr.res2d,
    title="GO Biological Process",
    cmap="viridis_r",
    size=6,
    figsize=(3, 5),
    cutoff=1,
)
```

### Interpretation prompt

When interpreting enriched terms, ask whether the module is consistent with processes expected in **myocardial infarction**, such as tissue injury, inflammation, extracellular-matrix remodeling, metabolism, vascular responses, or repair.

> Enrichment provides a functional summary of the genes in a module; it does not prove that the module causes the MI phenotype.

---

## 7. Optional — Visualize the detected modules (SKIP in Posit Cloud or if it takes too long on your own computer)

Instead of plotting the original graph object `g`, we can plot the clustering object. Setting `mark_groups=True` adds a separate background region for each detected module.

> **Optional task:** Can you identify which modules contain the high-degree hub genes from the earlier centrality analysis?

```python
fig, ax = plt.subplots()

ig.plot(
    clustering,
    target=ax,
    layout=layout,
    mark_groups=True,
    vertex_size=g.degree(),
)

fig.set_size_inches(20, 20)
```

---

## Workflow summary

You have now moved through a basic gene-network workflow:

**TPM expression → gene–gene correlations → co-expression network → network properties → centrality → Leiden modules → functional enrichment**

### Biological interpretation

The transcriptomics exercise asked which genes differ between **MI and Healthy samples**. This notebook asks a complementary question: **which genes vary together and how are they organized into a network?**

A gene can therefore be biologically interesting for different reasons:

- a **differentially expressed gene** changes in abundance between MI and Healthy samples;
- a **hub gene** has many co-expression partners;
- a **bridge gene** connects otherwise separated network regions; and
- a **module gene** belongs to a coordinated expression program.

These categories can overlap, but they do not have to. Combining differential-expression and network information can provide a richer view of the transcriptomic response to MI.
