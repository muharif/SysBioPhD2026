[⬇️ Download the Jupyter notebook](./TranscriptomicsExercise.ipynb?raw=1)

---

# Differential Expression Analysis with PyDESeq2

**RNA-seq workflow:** data inspection → PCA → sample correlation → differential expression → GO enrichment

**Source Data**: This exercise is based on published data "HIF1 mediates a switch in pyruvate kinase isoforms after myocardial infarction" by Williams et al.
> Manuscript Link: https://pubmed.ncbi.nlm.nih.gov/29652636/
> 
> Data Accession: GSE104187 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104187)

## Learning objectives

By the end of this notebook, you should be able to:

- inspect count, TPM, and metadata tables;
- use PCA and correlation heatmaps to assess sample structure and potential outliers;
- run differential expression analysis with **PyDESeq2**;
- identify significant, up-regulated, and down-regulated genes; and
- perform **GO Biological Process** enrichment with **GSEApy/Enrichr**.

---

## 1. Setup and libraries

Import the Python packages used throughout the notebook. The analysis uses **pandas/numpy** for data handling, **PCA and seaborn** for exploratory analysis, **PyDESeq2** for differential expression, and **GSEApy** for enrichment analysis.

```python
%matplotlib inline
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns
import matplotlib
## Adding these 4 lines will make our life easier if we need 
## to do figure post-processing in other software (like Adobe Illustrator)
#sns.set(font="Arial")
#sns.set_style("white")
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from pca import pca

from sklearn.decomposition import PCA as sklearnPCA
import os, math
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import pickle

from scipy.stats import zscore, spearmanr
import gseapy as gp
```

## 2. Load and inspect the input data

We begin by loading the three tables used in the analysis:

| File | Purpose | Expected layout |
|---|---|---|
| `count.txt` | Raw RNA-seq counts for DESeq2 | genes × samples |
| `tpm.txt` | TPM-normalized expression for exploratory analysis | genes × samples |
| `metadata.txt` | Sample-level experimental information | samples × metadata fields |

The sample identifiers in the expression matrices should match the index of the metadata table.

> **Important:** DESeq2 should be run on **raw counts**, not TPM values. TPM is used here for PCA/exploratory visualization.

```python
count = pd.read_csv("data/count.txt", sep = "\t", index_col = 0) ## Assuming your data is tab-separated (if comma-separated, change sep = ",")
tpm = pd.read_csv("data/tpm.txt", sep = "\t", index_col = 0) ## Assuming your data is tab-separated (if comma-separated, change sep = ",")

metadata = pd.read_csv("data/metadata.txt", sep = "\t", index_col = 0) ## Assuming your metadata data is tab-separated (if comma-separated, change sep = ",")
```

```python
## .head(): showing top 5 rows in a dataframe
count.head()
```

```python
metadata.head()
```

```python
## .shape: showing the number of (rows, columns) in a dataframe
count.shape
```

```python
metadata.shape
```

### Checkpoint 1 — Matrix dimensions

> **Q1.** How many genes and samples are present in the count matrix?  
> Use `.head()` to understand the orientation of the table and `.shape` to identify the number of rows and columns.

## 3. Exploratory data analysis: PCA

**Principal Component Analysis (PCA)** reduces high-dimensional expression data to a small number of components that capture the major sources of variation.

We use PCA here to:

1. assess whether samples cluster by experimental condition;
2. identify major sources of variation; and
3. look for potential outliers or batch effects.

For this exploratory analysis, PCA is performed on **log-transformed TPM values**, with **samples as observations** and **genes as features**.

**Reference:** [pca package documentation](https://erdogant.github.io/pca/pages/html/index.html)

### Choosing the number of principal components

The `pca` package can be initialized in two common ways:

1. **Retain enough PCs to explain a target fraction of variance**  
   `model = pca(n_components=0.95)` keeps enough components to explain 95% of the variance.

2. **Keep a fixed number of PCs**  
   `model = pca(n_components=3)` keeps three components and allows optional 3D visualization.

For this exercise, we will keep **3 principal components**.

### Step 3.1 — Initialize the PCA model

```python
model = pca(n_components=3)
```

### Step 3.2 — Prepare the expression matrix

Before PCA:

1. remove genes with zero expression across all samples;
2. align expression-matrix columns to the metadata sample order;
3. use the metadata `condition` labels for the PCA display; and
4. apply `np.log1p()` to reduce the influence of very large expression values.

> **What does `log1p` do?** It calculates `log(1 + x)`, which safely handles zero values.

```python
# Setting up the dataframes for exploratory analysis
## 1. Remove all-zero genes from the TPM and count matrices.
TPM_PCA = tpm[tpm.sum(1) > 0]
count_PCA = count[count.sum(1) > 0]

## 2. Align sample order with the metadata.
TPM_PCA = TPM_PCA[metadata.index]
count_PCA = count_PCA[metadata.index]

## 3. Replace TPM sample IDs with their condition labels for the PCA display.
TPM_PCA.columns = metadata["condition"]

## PS: Another pythonic way to replace the column names is:
# TPM_PCA = TPM_PCA.rename(columns=metadata["condition"])

## 4. Log-transform the TPM matrix.
TPM_PCA = np.log1p(TPM_PCA)
```

```python
## Observe the TPM_PCA dataframe. In some cases, we may encounter issues, including missing values.
TPM_PCA.head()
```

### Step 3.3 — Fit the PCA model

The `pca` package expects **samples as rows** and **genes as columns**, so the expression matrix is transposed with `.T` before fitting.

```python
results = model.fit_transform(TPM_PCA.T)
```

### Step 3.4 — Plot the PCA

Use `model.scatter()` to visualize the first principal components. Adjust `figsize` as needed for your display.

```python
fig, ax = model.scatter(figsize = (6,5))
```

### Answer key — PCA interpretation

In the expected analysis, the first two principal components explain **90.9% of the total variance**:

- **PC1:** 61.9%
- **PC2:** 29.0%

Observed sample structure:

- **SHAM samples** cluster together near the center of PC1 and lower on PC2.
- **MI_1D samples** cluster on the left side of PC1.
- **MI_3D samples** cluster on the right side of PC1.
- The dominant source of variation is therefore strongly associated with the experimental condition/time after MI.

No major outliers are apparent in the expected PCA result. The strong clustering of biological groups suggests a consistent transcriptomic response across replicates.

> **Interpretation note:** Explained variance can be especially high in tightly controlled experimental systems. Human datasets often contain additional biological and environmental sources of variation.

### Checkpoint 2 — PCA

> **Q2.1.** How would you interpret the PCA? Are there any apparent outliers? Justify your answer.  
> **Q2.2.** Is the total variance shown in the plot title necessarily equal to PC1 + PC2? If not, why?

### PCA quick reference

| Command | Purpose |
|---|---|
| `model.scatter3d()` | 3D PCA plot |
| `model.biplot()` | PCA biplot with influential/loading features |
| `model.biplot3d()` | 3D biplot |
| `model.results` | Detailed PCA output, including explained variance, loadings, and related results |

**Reference:** [pca package documentation](https://erdogant.github.io/pca/pages/html/index.html)

## 4. Exploratory data analysis: sample correlation

A sample-to-sample correlation heatmap provides a second view of sample similarity.

We use it to assess whether:

- biological replicates are highly correlated;
- samples from the same condition cluster together; and
- any sample behaves like a potential outlier.

Here, **Spearman correlation** is calculated between samples and visualized with a clustered heatmap.

### Step 4.1 — Calculate sample correlations

Use pandas `.corr()` on the filtered count matrix. Because samples are columns, the resulting matrix contains **sample-to-sample correlations**.

**Reference:** [pandas `DataFrame.corr`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.corr.html)

```python
corr = count_PCA.corr(method = "spearman")
```

### Step 4.2 — Plot the clustered correlation heatmap

Use `sns.clustermap()` to cluster samples based on their correlation profiles. The `RdYlBu` colormap is used here for readability.

```python
sns.clustermap(corr, cmap = "RdYlBu")
```

### Optional — Replace sample IDs with condition labels

For a more compact visualization, the correlation-matrix row and column labels can be renamed using the metadata `condition` values.

```python
corr = count_PCA.corr(method = "spearman")
corr = corr.rename(columns = metadata["condition"], index = metadata["condition"])
sns.clustermap(corr, cmap = "RdYlBu")
```

### Checkpoint 3 — Sample correlation

> **Q3.1.** What do you conclude from the clustering? Does the heatmap support the PCA result? Are there any apparent sample outliers?  
> **Q3.2.** Which correlation method is being used here, and what does that method measure?

**Reference:** [pandas `DataFrame.corr`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.corr.html)

## 5. Differential expression with PyDESeq2

We next test for **differentially expressed genes (DEGs)** using **PyDESeq2**, a Python implementation of the DESeq2 workflow.

DESeq2 models raw count data with a **negative binomial distribution** and accounts for differences in sequencing depth and biological variability.

The main outputs include:

- **log2 fold change** — direction and magnitude of expression change;
- **p-value** — evidence against the null hypothesis; and
- **adjusted p-value (`padj`)** — multiple-testing corrected significance measure.

> **Important:** Use the original **raw count matrix** for DESeq2, not TPM-normalized values.

**Reference:** [PyDESeq2 documentation](https://pydeseq2.readthedocs.io/en/stable/)

### What happens during DESeq2 analysis?

At a high level, `dds.deseq2()` performs the key steps required for differential expression testing:

1. estimate sample-specific size factors to account for sequencing-depth differences;
2. estimate gene-wise dispersion/variance; and
3. fit the statistical model used to test for differential expression.

We first create a `DeseqDataSet` (`dds`) containing the raw counts, metadata, and experimental design.

### Step 5.1 — Create the DESeq2 dataset

Use the original, unfiltered `count` matrix and `metadata`. The design formula uses the metadata column **`condition`**.

```python
dds = DeseqDataSet(
    counts=np.round(count.T),
    metadata=metadata,
    design="~condition",
    refit_cooks=True,
)
```

### Step 5.2 — Fit the DESeq2 model

Run `dds.deseq2()` to estimate the model parameters. This step may take a few minutes depending on dataset size.

```python
dds.deseq2()
```

> **Optional exploration:** Inspect the attributes stored inside `dds` to see the fitted model components and intermediate results.

### Step 5.3 — Define the comparison

Create a `DeseqStats` object comparing **MI_1D** against **SHAM_1D** using the `condition` variable.

The contrast is interpreted as:

**MI_1D / SHAM_1D**

Therefore, positive log2 fold changes indicate higher expression in **MI_1D**, while negative values indicate lower expression in **MI_1D** relative to **SHAM_1D**.

```python
ds = DeseqStats(dds, contrast=["condition", "MI_1D", "SHAM_1D"])
```

### Step 5.4 — Retrieve the statistical results

Run `ds.summary()` and save `ds.results_df` as `deseq_results`.

Because thousands of genes are tested simultaneously, the raw p-values are corrected for **multiple testing**. We will inspect the adjusted p-values (`padj`) and the MA plot before filtering significant genes.

```python
ds.summary()
```

```python
deseq_results = ds.results_df
```

### Step 5.5 — Inspect adjusted p-values

Plot a histogram of the `padj` column. A larger number of small adjusted p-values indicates stronger evidence for differential expression across the dataset.

Try changing the number of bins or adding labels/title if you want to customize the figure.

**Reference:** [pandas histogram plotting](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.hist.html)

```python
deseq_results["padj"].hist(bins = 100)
```

### Step 5.6 — Inspect the MA plot

Use `ds.plot_MA()` to visualize expression change across the abundance range.

- Each point represents a gene.
- Significant genes (`padj < 0.05`, by default) are highlighted.
- Points shown as triangles lie beyond the current y-axis limits.

The MA plot is useful for checking whether strong fold changes are concentrated among low- or high-expression genes.

```python
ds.plot_MA()
```

### Step 5.7 — Remove genes without an adjusted p-value

Some genes cannot be assigned a valid adjusted p-value, often because their counts are too low to provide sufficient information for testing.

Use `.notna()` on the `padj` column to retain only genes with a valid adjusted p-value.

**Expected answer-key value:** **15,968 genes** with non-missing `padj` values.

**Reference:** [pandas `DataFrame.notna`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.notna.html)

```python
notna_deseq = deseq_results[deseq_results["padj"].notna()]
notna_deseq.shape
```

### Step 5.8 — Identify significant DEGs

For this exercise, genes with **`padj < 0.05`** are considered significantly differentially expressed.

```python
sig_deseq = notna_deseq[notna_deseq["padj"] < 0.05]
sig_deseq.head()
```

```python
sig_deseq.shape
```

### Step 5.9 — Rank significant genes by fold change

We can inspect DEGs in three useful ways:

- **most down-regulated:** smallest log2 fold change;
- **most up-regulated:** largest log2 fold change; and
- **largest absolute change:** largest `|log2FoldChange|`, regardless of direction.

```python
## Ascending order (Top down-regulated genes on the top)
top_down = sig_deseq.sort_values("log2FoldChange", ascending = True)
top_down.head(10) # you can adjust the number based on how many rows you want to show. .head(15) means top 15, etc.
```

```python
## Descending order (Top up-regulated genes on the top)
top_up = sig_deseq.sort_values("log2FoldChange", ascending = False)
top_up.head(10) # you can adjust the number based on how many rows you want to show. .head(15) means top 15, etc.
```

```python
## Sorting it based on the absolute log2FC (mixed of top up and down-regulated genes, based on their log2FC)
top_mixed = sig_deseq.sort_values("log2FoldChange", ascending = False, key = abs)
top_mixed.head(10) # you can adjust the number based on how many rows you want to show. .head(15) means top 15, etc.
```

### Checkpoint 4 — Differential expression

> **Q4.1.** How many significantly differentially expressed genes are detected at `padj < 0.05`?  
> **Q4.2.** *(Optional)* What does the adjusted-p-value histogram suggest? From the MA plot, approximately what are the largest positive and negative log2 fold changes?  
> **Q4.3.** How many genes have `padj = NA`? Compare the total number of genes from Q1 with the number retained in Step 5.7.  
> **Q4.4.** Which 10 genes have the largest **absolute** log2 fold changes, regardless of direction?

## Optional: Redo the analysis to compare MI_3D and SHAM_3d

### Step 5.10 — Save the DESeq2 results

Save all genes with valid adjusted p-values, including non-significant genes, as a tab-separated file:

`Results/Deseq2_Results.txt`

This full table can be reused in later network and visualization exercises and can also be opened in Excel.

```python
notna_deseq.to_csv("Results/Deseq2_Results.txt", sep = "\t")
```

## 6. GO Biological Process enrichment

Differential expression tells us **which genes change**. Enrichment analysis helps us ask **which biological processes are overrepresented** among those genes.

We will perform separate enrichment analyses for:

- significantly **up-regulated** genes; and
- significantly **down-regulated** genes.

The analysis uses **GSEApy** to query the **Enrichr** web service with the **GO Biological Process** gene-set library.

The result is a ranked table of enriched GO terms with statistical significance and enrichment scores.

**References:**  
- [Enrichr](https://maayanlab.cloud/Enrichr/)  
- [GSEApy examples](https://gseapy.readthedocs.io/en/latest/gseapy_example.html)

### Enrichment workflow

### Step 6.1 — Split significant DEGs by direction

Create two gene lists from `sig_deseq`:

- **up-regulated:** `log2FoldChange > 0` and `padj < 0.05`
- **down-regulated:** `log2FoldChange < 0` and `padj < 0.05`

**Expected answer-key values:**

| Direction | Number of genes |
|---|---:|
| Up-regulated | X |
| Down-regulated | Y |

Use `.to_list()` to convert the pandas index/series to a standard Python list.

**Reference:** [pandas `Series.to_list`](https://pandas.pydata.org/docs/reference/api/pandas.Series.to_list.html)

```python
### You can also reload the deseq results
# notna_deseq = pd.read_csv("Results/Deseq2_Results.txt", sep = "\t", index_col = "Gene.name")
# sig_deseq = notna_deseq[notna_deseq["padj"] < 0.05]
# sig_deseq.head()
```

```python
up = sig_deseq[sig_deseq["log2FoldChange"] > 0].index.to_list() 
down = sig_deseq[sig_deseq["log2FoldChange"] < 0].index.to_list()

## you can adjust the log2FoldChange threshold too! Example below: We want log2FC to be above 3 (or below -3 for the down-regulated genes)
# up = sig_deseq[sig_deseq["log2FoldChange"] > 3].index.to_list() 
# down = sig_deseq[sig_deseq["log2FoldChange"] < -3].index.to_list()
```

```python
len(up) # numbers of up-regulated genes
```

```python
len(down) # numbers of down-regulated genes
```

### Step 6.2 — Run Enrichr separately for each direction

Run one enrichment analysis for the up-regulated genes (`enr_up`) and another for the down-regulated genes (`enr_down`).

For consistency in this exercise, use the **GO Biological Process 2026** library and set the organism to **mouse**.

**Reference:** [GSEApy Enrichr web-service example](https://gseapy.readthedocs.io/en/latest/gseapy_example.html#Enrichr-Web-Serives-(without-a-backgound-input))

```python
enr_up = gp.enrichr(gene_list=up, 
                 gene_sets='GO_Biological_Process_2026',
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast, human
                 outdir=None, # don't write to disk
                )
```

```python
enr_down = gp.enrichr(gene_list=down, # or "./tests/data/gene_list.txt",
                 gene_sets='GO_Biological_Process_2026',
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
```

### Step 6.3 — Inspect enrichment results

The detailed Enrichr output is available as a pandas DataFrame in the `.res2d` attribute of each result object.

Inspect the first rows and count how many terms satisfy **Adjusted P-value < 0.05**.

```python
enr_up.res2d.head(10)
```

```python
enr_up.res2d[enr_up.res2d["Adjusted P-value"] < 0.05].shape
```

```python
enr_down.res2d.head(10)
```

```python
enr_down.res2d[enr_down.res2d["Adjusted P-value"] < 0.05].shape
```

### Step 6.4 — Visualize the top enriched processes

Use GSEApy's `dotplot()` to display the strongest GO Biological Process enrichments for each direction. Explore the plotting options to improve readability.

**Reference:** [GSEApy plotting examples](https://gseapy.readthedocs.io/en/latest/gseapy_example.html#Plotting)

```python
ax = gp.dotplot(enr_up.res2d, title='GO BP',cmap='viridis_r', size=3, figsize=(3,5))
```

```python
ax = gp.dotplot(enr_down.res2d, title='GO BP',cmap='viridis_r', size=3, figsize=(3,5))
```

### Checkpoint 5 — Enrichment analysis

> **Q5.1.** How many significantly enriched GO processes are found for the up- and down-regulated gene lists at `Adjusted P-value < 0.05`?  
> **Q5.2.** What does the Enrichr **Combined Score** represent? See the [Enrichr help page](https://maayanlab.cloud/Enrichr/help#basics).  
> **Q5.3.** Which genes contribute to the top enriched process in each direction?  
> **Q5.4.** What are the top 10 enriched processes for the up- and down-regulated genes?  
> **Q5.5.** Replot the results using `barplot()` instead of `dotplot()`. How does the default x-axis differ?

> **Optional extension:** Repeat the enrichment analysis using the KEGG gene-set library (`gene_sets="KEGG_2026"`).

### Step 6.5 — Save the enrichment results

Save the complete enrichment tables as tab-separated files:

- `Results/GO_up.txt`
- `Results/GO_down.txt`

These files can also be opened directly in Excel.

```python
enr_up.res2d.to_csv("Results/GO_up.txt", sep = "\t")
enr_down.res2d.to_csv("Results/GO_down.txt", sep = "\t")
```

---

## Workflow summary

You have now moved through a complete introductory RNA-seq analysis workflow:

**input QC → PCA → sample correlation → DESeq2 → DEG ranking → GO enrichment**

> **Good practice:** Before reusing this workflow with another dataset, confirm the sample IDs, metadata design, comparison/contrast, organism, and gene identifiers.

---

## Optional Exercise — Adjusting for `day` in the DESeq2 model

So far, the main analysis used a combined `condition` variable and tested **MI_1D vs SHAM_1D**. In this optional exercise, we will use the additional metadata variable **`day`** and ask a slightly different question:

> **Which genes differ between MI and SHAM after accounting for the overall effect of day?**

The multifactor metadata contains two variables:

- **`MI`** — experimental group (`0 = SHAM`, `1 = MI`)
- **`day`** — collection time (`1` or `3` days)

We will fit the additive design:

```text
~ day + MI
```

Here, `day` is included as an **adjustment variable (fixed effect)**. The model estimates the MI effect after accounting for systematic differences between day 1 and day 3.

> **Important assumption:** `~ day + MI` assumes that the MI effect is reasonably consistent across days. It does **not** model an `MI × day` interaction. If the biological effect of MI is expected to change substantially over time, an interaction model would be more appropriate.

### Optional Step 1 — Load the multifactor metadata

Reload the raw count matrix and use `metadata_multi.txt`, which contains separate columns for **MI status** and **day**.

As before, DESeq2 should be run on **raw counts**, not TPM values.

```python
count = pd.read_csv("data/count.txt", sep = "\t", index_col = 0)
metadata_multi = pd.read_csv("data/metadata_multi.txt", sep = "\t", index_col = 0)
```

```python
metadata_multi
```

### Optional Step 2 — Prepare the metadata variables

We treat  `MI` and `day` as a categorical variable rather than a continuous numerical measurement. Finally, align the metadata rows to the sample order in the count matrix.

```python
metadata_multi["MI"] = metadata_multi["MI"].astype(str)
metadata_multi["day"] = metadata_multi["day"].astype(str)

metadata_multi = metadata_multi.loc[count.columns]
metadata_multi
```

### Optional Step 3 — Create and fit the multifactor DESeq2 model

Create a new `DeseqDataSet` using:

```text
~ day + MI
```

This design separates the overall **day effect** from the overall **MI effect**. The MI coefficient therefore represents the difference between MI and SHAM **after adjusting for day**.

```python
dds_multi = DeseqDataSet(
    counts=np.round(count.T),
    metadata=metadata_multi,
    design="~ day + MI",
    refit_cooks=True,
)
```

```python
dds_multi.deseq2()
```

### Optional Step 4 — Define the comparison: MI vs SHAM

Now test **MI vs SHAM** using the `MI` variable while retaining `day` in the fitted model.

The contrast is interpreted as:

**MI / SHAM, adjusted for day**

Therefore:

- positive `log2FoldChange` → higher expression in **MI**;
- negative `log2FoldChange` → lower expression in **MI**.

Unlike the earlier **MI_1D vs SHAM_1D** comparison, this is not restricted to a day-1-specific effect. It estimates one common MI effect using information from both days while controlling for the overall difference between day 1 and day 3.

```python
ds_multi = DeseqStats(
    dds_multi,
    contrast=["MI", "MI", "SHAM"]
)
```

```python
ds_multi.summary()
deseq_results_multi = ds_multi.results_df
```

### Optional Step 5 — Identify all significant DEGs

Inspect the adjusted p-values and MA plot as in the main analysis, then retain genes with a valid adjusted p-value and apply the same significance threshold:

```text
padj < 0.05
```

The resulting `sig_deseq_multi` table contains **all significantly differentially expressed genes for MI vs SHAM after adjusting for day**.

```python
deseq_results_multi["padj"].hist(bins = 100)
```

```python
ds_multi.plot_MA()
```

```python
notna_deseq_multi = deseq_results_multi[deseq_results_multi["padj"].notna()]
notna_deseq_multi.shape
```

```python
sig_deseq_multi = notna_deseq_multi[notna_deseq_multi["padj"] < 0.05]

print(f"Number of significant DEGs (MI vs SHAM, adjusted for day): {sig_deseq_multi.shape[0]}")
sig_deseq_multi.head()
```

### Why is this different from `MI_1D vs SHAM_1D`?

The two analyses answer **different biological questions**, so they do not have to produce the same DEG list.

| Analysis | Design / contrast | Question being asked |
|---|---|---|
| **Main analysis** | `~ condition`; `MI_1D vs SHAM_1D` | At **day 1**, which genes differ between MI and SHAM? |
| **Optional analysis** | `~ day + MI`; `MI vs SHAM` | Across the experiment, which genes differ between MI and SHAM **after accounting for the overall effect of day**? |

In the main analysis, `condition` combines treatment and time into separate groups such as `MI_1D`, `SHAM_1D`, `MI_3D`, and `SHAM_3D`. The requested contrast is therefore **day-1 specific**. Samples from the other groups can still contribute to model fitting and dispersion estimation, but the tested difference itself is specifically **MI_1D vs SHAM_1D**.

In the optional model, `day` and `MI` are represented separately. The model uses the samples from both days to estimate a **single MI effect**, while correcting for the overall shift between day 1 and day 3.

This can change the DEG list for several reasons:

- genes may have a strong **day effect** that is separated from the MI effect in the multifactor model;
- combining information across both days can improve precision when the MI effect is consistent over time;
- a gene may respond strongly to MI at day 1 but weakly or differently at day 3, so its overall adjusted MI effect can be smaller than its day-1-specific effect.

> **Key interpretation:** `MI vs SHAM` with `~ day + MI` should not be described as the same comparison with an extra covariate. It asks for a **common MI effect across days, adjusted for day**, whereas `MI_1D vs SHAM_1D` asks specifically for the **MI effect at day 1**.

> **If the MI effect changes with time:** consider an interaction model such as `~ day * MI`. That model can explicitly test whether the effect of MI differs between day 1 and day 3.

### Optional checkpoint

> **QO.1.** How many significant DEGs are detected for **MI vs SHAM after adjusting for day** at `padj < 0.05`?  
> **QO.2.** Is this number larger or smaller than the number obtained for **MI_1D vs SHAM_1D**?  
> **QO.3.** Why should the two analyses not necessarily identify the same genes?  
> **QO.4.** What assumption about the MI effect is made by the additive design `~ day + MI`?
