[⬇️ Download the Jupyter notebook](MetabolomicsExercise.ipynb?raw=1) (Right Click -> Save Link As..)

---

# Introductory Metabolomics Analysis

**Metabolomics workflow:** load data → check samples → quality control → preprocess → PCA/correlation → differential abundance → volcano/heatmap → biological interpretation

**Source Data**: This exercise is based on published data "An Integrative Multiomics Framework for Identification of Therapeutic Targets in Pulmonary Fibrosis" by Arif et al.
> Manuscript Link: https://pubmed.ncbi.nlm.nih.gov/37038090/
> 
> Data is available on Supplementary Material 3 (Raw data, unmapped)

This workshop introduces a **basic downstream metabolomics analysis** in Python. It is designed for participants with little or no computational background and follows the same step-by-step style as the genomics and transcriptomics exercises.

## What are we starting from?

We assume that the raw instrument data have already been processed into a table of metabolite abundances or peak intensities.

In a real LC-MS/GC-MS workflow, important upstream steps can include peak detection, alignment, blank filtering, drift/batch correction, internal-standard normalization, and metabolite annotation. Those steps are **assay- and platform-specific** and are outside the scope of this introductory exercise.

## Learning objectives

By the end of this notebook, you should be able to:

- inspect a metabolomics abundance table and sample metadata;
- verify that sample identifiers match between the two tables;
- assess missing values and simple data-quality metrics;
- understand why transformation and scaling matter in metabolomics;
- use PCA and sample correlation to explore sample structure and possible outliers;
- test metabolites for differences between two experimental groups;
- interpret fold change, p-values, and FDR-adjusted p-values;
- make a volcano plot and a heatmap of changing metabolites; and
- understand why metabolite identification and batch/QC information are important for biological interpretation.

---

## 1. Setup and libraries

We use:

- `pandas` for tables;
- `numpy` for numerical operations;
- `matplotlib` and `seaborn` for plotting;
- `scikit-learn` for PCA and scaling; and
- `scipy` for statistical testing and multiple-testing correction.

```python
%matplotlib inline

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pca import pca
#from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind, false_discovery_control

os.makedirs("Results", exist_ok=True)
```

## 2. Load and inspect the input data

For the whole workshop we will use two main variables:

- **`data`** — the metabolomics abundance table;
- **`metadata`** — sample-level information such as experimental group, batch, sex, time point, or sample type.

For this notebook, the expected layout is:

| Object | Rows | Columns |
|---|---|---|
| `data` | samples | metabolites/features |
| `metadata` | samples | sample information |

The **first column of each file should contain sample IDs** and is loaded as the row index.

> If your files are tab-separated, use `sep="\t"`. If they are comma-separated, use `sep=","` or omit the `sep` argument.

```python
# Change the filenames if needed.
data = pd.read_csv("data/metabolomics_data.txt", index_col=0, sep = "\t")
metadata = pd.read_csv("data/metabolomics_metadata.txt", index_col=0, sep = "\t")
```

### Step 2.1 — Inspect the metabolomics table

Use `.head()` to look at the first few samples and metabolites.

```python
data.head()
```

Check its dimensions. Because rows are samples and columns are metabolites/features, `.shape` reports:

```text
(number of samples, number of metabolites)
```

```python
print("data shape:", data.shape)
```

### Step 2.2 — Inspect the metadata

```python
metadata.head()
```

```python
print("metadata shape:", metadata.shape)
print()
print(metadata.dtypes)
```

### Checkpoint 1 — Getting to know the data

> **Q1.1.** How many samples are present in `data`?  
> **Q1.2.** How many metabolites/features are measured?  
> **Q1.3.** Which variables are available in `metadata`? Which one represents the biological group you want to compare?

## 3. Do the metabolomics data and metadata contain the same samples?

Never assume that two files contain the same samples simply because they come from the same experiment.

We will explicitly check:

1. duplicate sample IDs;
2. samples found only in the metabolomics table;
3. samples found only in the metadata; and
4. the number of samples shared by both tables.

```python
print("Duplicate sample IDs in data:", data.index.duplicated().sum())
print("Duplicate sample IDs in metadata:", metadata.index.duplicated().sum())

samples_data = set(data.index)
samples_metadata = set(metadata.index)

print("Samples in data:", len(samples_data))
print("Samples in metadata:", len(samples_metadata))
print("Samples shared by both:", len(samples_data & samples_metadata))
print("Only in data:", len(samples_data - samples_metadata))
print("Only in metadata:", len(samples_metadata - samples_data))
```

### Step 3.1 — Align the two tables

We keep only shared samples and place the metadata in exactly the same order as `data`.

This is one of the most important bookkeeping steps in any omics analysis: **the sample labels must refer to the correct measurements.**

```python
common_samples = data.index.intersection(metadata.index)

data = data.loc[common_samples].copy()
metadata = metadata.loc[common_samples].copy()

# Make sure the order is identical.
metadata = metadata.loc[data.index]

print("Aligned data shape:", data.shape)
print("Aligned metadata shape:", metadata.shape)
print("Sample order identical:", data.index.equals(metadata.index))
```

### Checkpoint 2 — Sample matching

> **Q2.1.** Were any samples present in only one of the two files?  
> **Q2.2.** Why is it dangerous to compare groups before checking that sample IDs and sample order are correct?

## 4. Basic quality control

Metabolomics data require quality control before statistical testing.

Important questions include:

- Are there missing values?
- Are some metabolites missing in many samples?
- Are some samples missing many metabolites?
- Are values numeric and non-negative where expected?
- Are there pooled QC samples or blanks?
- Are there analytical batches that could explain sample differences?

There is **no single universal QC threshold**. The correct decisions depend on the platform, study design, preprocessing pipeline, and whether the data represent targeted concentrations or untargeted peak intensities.

### Step 4.1 — Make sure the abundance columns are numeric

```python
# Convert every metabolite column to numeric.
# Values that cannot be interpreted as numbers become NaN.
data = data.apply(pd.to_numeric, errors="coerce")

data.head()
```

### Step 4.2 — Count missing values

First inspect missingness across the entire table.

```python
missing_total = data.isna().sum().sum()
missing_fraction = data.isna().mean().mean()

print("Total missing values:", missing_total)
print("Overall missing fraction:", round(missing_fraction, 4))
```

Missingness can also be inspected separately for each **metabolite** and each **sample**.

```python
feature_missing = data.isna().mean(axis=0)
sample_missing = data.isna().mean(axis=1)

print("Feature missingness:")
print(feature_missing.describe())
print()
print("Sample missingness:")
print(sample_missing.describe())
```

```python
plt.figure(figsize=(6, 4))
plt.hist(feature_missing, bins=20)
plt.xlabel("Fraction missing per metabolite")
plt.ylabel("Number of metabolites")
plt.title("Metabolite missingness")
plt.show()
```

### Step 4.3 — Remove metabolites with too much missing data

For this teaching exercise we will remove metabolites missing in more than **20% of samples**.

> **Important:** 20% is an illustrative workshop threshold, not a universal metabolomics rule. In a real project, the threshold should be justified based on the assay and study design.

```python
max_missing = 0.20

keep_features = feature_missing <= max_missing
print("Metabolites before filtering:", data.shape[1])
print("Metabolites retained:", keep_features.sum())
print("Metabolites removed:", (~keep_features).sum())

data_filtered = data.loc[:, keep_features].copy()
```

### Step 4.4 — What about zero values?

A zero can mean different things depending on the dataset:

- a true measured zero;
- a signal below the detection limit;
- a missing peak encoded as zero; or
- a value introduced during preprocessing.

Therefore we **do not automatically convert zeros to missing values**.

Always check how your metabolomics provider or preprocessing pipeline encoded non-detected features.

```python
zero_fraction = (data_filtered == 0).mean().mean()
print("Fraction of values equal to zero:", round(zero_fraction, 4))
```

### Checkpoint 3 — Quality control

> **Q3.1.** What fraction of the dataset is missing?  
> **Q3.2.** How many metabolites were removed using the 20% missingness threshold?  
> **Q3.3.** Why should zeros not automatically be treated as missing values in every metabolomics dataset?  
> **Q3.4.** If pooled QC samples are available, what would you hope to see when you plot them with the biological samples?

## 5. Missing-value handling and transformation

Metabolomics measurements often span several orders of magnitude and can have strongly right-skewed distributions.

For this workshop we will do two simple preprocessing steps:

1. replace remaining missing values with **half of the smallest positive value for that metabolite**;
2. apply a **log2 transformation**.

This is a simple teaching strategy for abundance/intensity data. It is **not universally appropriate**. Missing-value methods should reflect why values are missing, and targeted absolute concentration data may require a different approach.

### Step 5.1 — Impute the remaining missing values

```python
data_imputed = data_filtered.copy()

for metabolite in data_imputed.columns:
    positive = data_imputed.loc[data_imputed[metabolite] > 0, metabolite]
    if len(positive) > 0:
        replacement = positive.min() / 2
        data_imputed[metabolite] = data_imputed[metabolite].fillna(replacement)

# Remove any feature that still contains missing values.
# This can happen if a feature had no positive observations.
data_imputed = data_imputed.dropna(axis=1)

print("Remaining missing values:", data_imputed.isna().sum().sum())
print("Shape after imputation:", data_imputed.shape)
```

### Step 5.2 — Inspect distributions before transformation

A boxplot gives a quick overview of the scale of each sample.

```python
plt.figure(figsize=(10, 4))
sns.boxplot(data=data_imputed.T, orient="v", showfliers=False)
plt.xticks([])
plt.ylabel("Abundance / intensity")
plt.xlabel("Samples")
plt.title("Data before log transformation")
plt.show()
```

### Step 5.3 — Log2-transform the data

A log transformation compresses very large values and often makes metabolomics data more symmetric.

Because some datasets contain zeros, we add a very small constant before taking the logarithm.

```python
positive_values = data_imputed.to_numpy()[data_imputed.to_numpy() > 0]
pseudocount = positive_values.min() / 2 if len(positive_values) else 1e-9

data_log = np.log2(data_imputed + pseudocount)

data_log.head()
```

```python
plt.figure(figsize=(10, 4))
sns.boxplot(data=data_log.T, orient="v", showfliers=False)
plt.xticks([])
plt.ylabel("log2 abundance")
plt.xlabel("Samples")
plt.title("Data after log2 transformation")
plt.show()
```

### What about normalization?

**Normalization, transformation, and scaling are different operations.**

- **Sample normalization** attempts to correct systematic differences in the overall signal between samples. Examples include normalization to sample amount, an internal standard, median normalization, or PQN.
- **Transformation** changes the mathematical scale, for example log transformation.
- **Scaling** changes how much each metabolite contributes to a multivariate analysis such as PCA.

For this workshop, we assume that any assay-specific sample normalization required by the experiment was already performed before the abundance table was provided. We therefore do **not** apply an additional sample-normalization method automatically.

This is an important real-world question to ask before analysing metabolomics data: **What preprocessing and normalization have already been performed?**

## 6. PCA: what are the major sources of variation? (Similar to our Transcriptomics Analysis, refer to that exercise for more details in PCA)

**Principal Component Analysis (PCA)** reduces hundreds or thousands of metabolite measurements into a small number of components that capture the strongest patterns in the dataset.

We use PCA to ask:

1. Do samples cluster by biological group?
2. Are there possible outliers?
3. Do pooled QC samples cluster tightly?
4. Could batch or another technical variable explain the strongest pattern?

```python
group_col = "condition"
```

```python
data_log_PCA = data_log.rename(index = metadata[group_col])
```

```python
model = pca(n_components=3)
```

```python
results = model.fit_transform(data_log_PCA)
```

```python
fig, ax = model.scatter(figsize = (6,5))
```

### Checkpoint 4 — PCA

> **Q4.1.** Do samples separate by biological group?  
> **Q4.2.** Are there any apparent outliers?  
> **Q4.3.** How much variance is explained by PC1 and PC2?  
> **Q4.4.** If you have batch or pooled-QC information, does the PCA suggest a technical effect?

## 8. Differential metabolite abundance

We now ask which metabolites differ between two experimental groups.

For this introductory exercise we will use a **Welch two-sample t-test** on the log2-transformed abundances.

For every metabolite we calculate:

- **log2 fold change** — direction and approximate magnitude of the difference;
- **p-value** — evidence against the null hypothesis; and
- **FDR-adjusted p-value (`padj`)** — correction for testing many metabolites at once.

> **Important:** This simple test is appropriate only for a basic two-group comparison. Paired samples, repeated measures, batches, covariates, and more complex designs require a model that matches the experimental design.

### Step 8.1 — Inspect the group labels and choose two groups

The code below uses the first two group labels by default. If your metadata contain more than two groups, change `group_1` and `group_2` to the comparison you want.

```python
metadata
```

```python
groups = metadata[group_col].dropna().unique().tolist()
print("Available groups:", groups)

group_1 = "Control"
group_2 = "Bleo14D"

print("Comparison:", group_2, "vs", group_1)
```

The contrast is interpreted as:

```text
Bleo14D / Control
```

Therefore:

- positive log2 fold change → higher in `Bleo14D`;
- negative log2 fold change → lower in `Bleo14D`.

### Step 8.2 — Run one test per metabolite

```python
samples_1 = metadata.index[metadata[group_col] == group_1]
samples_2 = metadata.index[metadata[group_col] == group_2]

results = []

for metabolite in data_log.columns:
    values_1 = data_log.loc[samples_1, metabolite]
    values_2 = data_log.loc[samples_2, metabolite]

    test = ttest_ind(values_2, values_1, equal_var=False, nan_policy="omit")
    log2fc = values_2.mean() - values_1.mean()

    results.append({
        "metabolite": metabolite,
        "mean_group_1_log2": values_1.mean(),
        "mean_group_2_log2": values_2.mean(),
        "log2FoldChange": log2fc,
        "pvalue": test.pvalue,
    })

results = pd.DataFrame(results).set_index("metabolite")
results["padj"] = false_discovery_control(results["pvalue"].fillna(1).to_numpy(), method="bh")
results = results.sort_values("padj")

results.head(10)
```

### Step 8.3 — Identify significant metabolites

For this workshop we will use:

```text
padj < 0.05
```

as the statistical significance threshold.

We will also inspect effect size separately instead of treating statistical significance and biological importance as the same thing.

```python
sig = results[results["padj"] < 0.05].copy()

print("Total metabolites tested:", results.shape[0])
print("Significant metabolites (padj < 0.05):", sig.shape[0])
```

### Step 8.4 — Inspect the strongest changes

```python
print("Most increased in", group_2)
display(sig.sort_values("log2FoldChange", ascending=False).head(10))

print("Most decreased in", group_2)
display(sig.sort_values("log2FoldChange", ascending=True).head(10))
```

### Checkpoint 6 — Differential abundance

> **Q6.1.** How many metabolites have `padj < 0.05`?  
> **Q6.2.** Which metabolites have the largest positive and negative log2 fold changes?  
> **Q6.3.** Why do we adjust p-values when testing many metabolites?  
> **Q6.4.** Does a very small p-value necessarily mean that the metabolite has a large biological effect?

## 9. Volcano plot (This can also be used in Transcriptomics)

A volcano plot combines **effect size** and **statistical evidence**.

- x-axis: log2 fold change;
- y-axis: `-log10(adjusted p-value)`;
- farther left/right: larger abundance difference;
- higher: stronger statistical evidence.

For visualization we will highlight metabolites with:

- `padj < 0.05`; and
- `|log2 fold change| >= 1`.

The fold-change threshold is illustrative and should be adapted to the biology and measurement scale of the experiment.

```python
plot_results = results.copy()
plot_results["minus_log10_padj"] = -np.log10(plot_results["padj"].clip(lower=1e-300))
plot_results["highlight"] = (
    (plot_results["padj"] < 0.05) &
    (plot_results["log2FoldChange"].abs() >= 1)
)

plt.figure(figsize=(7, 5))
plt.scatter(
    plot_results.loc[~plot_results["highlight"], "log2FoldChange"],
    plot_results.loc[~plot_results["highlight"], "minus_log10_padj"],
    alpha=0.5,
    s=25,
    label="Other metabolites"
)
plt.scatter(
    plot_results.loc[plot_results["highlight"], "log2FoldChange"],
    plot_results.loc[plot_results["highlight"], "minus_log10_padj"],
    alpha=0.8,
    s=35,
    label="padj < 0.05 and |log2FC| >= 1"
)
plt.axhline(-np.log10(0.05), linestyle="--", linewidth=1)
plt.axvline(-1, linestyle="--", linewidth=1)
plt.axvline(1, linestyle="--", linewidth=1)
plt.xlabel(f"log2 fold change: {group_2} / {group_1}")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano plot")
plt.legend(fontsize=8)
plt.tight_layout()
plt.show()
```

### Optional — Label the strongest metabolites

```python
top_labels = plot_results.sort_values("padj").head(10)

plt.figure(figsize=(7, 5))
plt.scatter(plot_results["log2FoldChange"], plot_results["minus_log10_padj"], alpha=0.5, s=25)
plt.axhline(-np.log10(0.05), linestyle="--", linewidth=1)
plt.axvline(-1, linestyle="--", linewidth=1)
plt.axvline(1, linestyle="--", linewidth=1)

for metabolite, row in top_labels.iterrows():
    plt.text(row["log2FoldChange"], row["minus_log10_padj"], str(metabolite), fontsize=8)

plt.xlabel(f"log2 fold change: {group_2} / {group_1}")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano plot: top metabolites labelled")
plt.tight_layout()
plt.show()
```

### Checkpoint 7 — Volcano plot

> **Q7.1.** Are most changing metabolites increased or decreased in `group_2`?  
> **Q7.2.** Are there metabolites with a large fold change but weak statistical evidence?  
> **Q7.3.** Why might a large fold change fail to reach statistical significance?

## 10. Heatmap of the strongest metabolites

A heatmap lets us look at the abundance pattern of multiple metabolites across individual samples.

We will select the metabolites with the smallest adjusted p-values and display **within-metabolite z-scores**.

This means the colour represents whether a metabolite is relatively high or low **compared with its own average**. Colours should not be used to compare absolute concentrations between different metabolites.

```python
number_to_plot = min(10, results.shape[0])
top_metabolites = results.head(number_to_plot).index

heatmap_data = data_log.loc[:, top_metabolites].copy()
heatmap_data.index = [f"{sample} | {metadata.loc[sample, group_col]}" for sample in heatmap_data.index]

sns.clustermap(
    heatmap_data.T,
    z_score=0,
    cmap="vlag",
    figsize=(10, 8)
)
plt.show()
```

### Checkpoint 8 — Heatmap

> **Q8.1.** Do samples cluster according to experimental group?  
> **Q8.2.** Can you identify groups of metabolites with similar abundance patterns?  
> **Q8.3.** Why does the heatmap use z-scores instead of raw abundance values?

## 11. Save the statistical results

Save the complete table, not only the significant metabolites. Keeping the full table makes it possible to revisit thresholds or use ranked results later.

```python
results.to_csv("Results/metabolite_statistics.csv")
sig.to_csv("Results/significant_metabolites.csv")
```

## 12. From significant metabolites to biology

Statistical analysis tells us **which measured features differ**, but metabolomics has an additional challenge: **what exactly is each feature?**

Before making strong biological claims, ask:

1. Is the feature a confidently identified metabolite or only an m/z/retention-time feature?
2. Was identification supported by an authentic standard, MS/MS spectrum, or database match?
3. Does the metabolite have a stable identifier such as HMDB, KEGG, PubChem, or ChEBI?
4. Could multiple features correspond to the same compound, isotope, adduct, or fragment?
5. Is the direction of change consistent with the underlying biology and sample type?

This is one reason metabolomics interpretation is not simply "take the significant list and run pathway enrichment."

### Optional — Prepare a metabolite list for pathway/enrichment analysis

If your column names are **confident metabolite names or database identifiers**, export the significant list and use an appropriate metabolomics pathway tool such as MetaboAnalyst.

If the features are still unidentified LC-MS peaks, pathway analysis should wait until annotation has been addressed.

```python
significant_names = sig.index.to_series(name="metabolite")
significant_names.to_csv("Results/significant_metabolite_names.txt", index=False, header=False)

print("Saved Results/significant_metabolite_names.txt")
```

## 13. Metabolite-set enrichment with MetaboAnalyst

So far we have identified individual metabolites that differ between the two groups. The next question is:

> **Do several of these metabolites point to the same biological pathway or biological process?**

This is the idea behind **metabolite-set enrichment analysis**.

For this introductory workshop, we will use the **MetaboAnalyst web interface** rather than adding another Python package. This makes it easier to focus on the biological idea rather than software syntax.

### Important before enrichment

Enrichment is most meaningful when the features have reliable metabolite names or identifiers such as **HMDB, KEGG, PubChem, or ChEBI**.

(Not in our case, since all is mapped) If your columns are only LC-MS features such as `m/z 301.123 @ 4.6 min`, do **not** pretend that they are confidently identified metabolites. Untargeted feature-level pathway analysis requires a different strategy, discussed briefly below.

### Step 13.1 — Export the metabolite lists

We will create two simple files:

1. **significant metabolites** — the metabolites passing our statistical threshold;
2. **reference/background metabolites** — all metabolites that were actually measured and tested in this experiment.

Why keep a background list?

A metabolomics platform measures only part of the metabolome. If we compare our significant metabolites against *all metabolites known to biology*, we may obtain enrichment simply because our platform preferentially measures certain metabolite classes.

Using the metabolites that were **actually measurable in this experiment** gives a fairer reference set.

```python
# Significant metabolites
sig.index.to_series(name="metabolite").to_csv(
    "Results/metaboanalyst_significant_metabolites.txt",
    index=False,
    header=False
)

# Background: every metabolite that was included in the statistical test
data_log.columns.to_series(index=None, name="metabolite").to_csv(
    "Results/metaboanalyst_background_metabolites.txt",
    index=False,
    header=False
)

print("Saved files for enrichment:")
print("Results/metaboanalyst_significant_metabolites.txt")
print("Results/metaboanalyst_background_metabolites.txt")
```

### Step 13.2 — Run over-representation analysis in MetaboAnalyst

Open **MetaboAnalyst** in your browser:

https://www.metaboanalyst.ca/

Then choose:

**Enrichment Analysis → Over Representation Analysis**

Use the file:

```text
Results/metaboanalyst_significant_metabolites.txt
```

### **NOTE: open your "metaboanalyst_significant_metabolites.txt" file and remove ALL the double quote `"` mark.**

Then:

1. choose the identifier type that matches your metabolite names/IDs; (Compound names for our current case)
2. check how many metabolites are successfully recognized;
3. use the **measured-metabolite background/reference list** when the interface gives you the option; (Not needed for now.)
4. choose a pathway/metabolite-set library, for example **KEGG** or another library appropriate for your organism and question;
5. inspect the enriched metabolite sets and their FDR-adjusted significance.

### Do not skip the mapping step

If many metabolites are not recognized, the enrichment result may be misleading. Always look at which metabolites successfully mapped before interpreting pathways.

### How to read the enrichment result

A pathway appearing near the top does **not** mean that the entire pathway is activated.

It means that **more metabolites from that set/pathway appeared in your selected list than expected by chance**, given the chosen reference/background.

Look at:

- the pathway/metabolite-set name;
- the number of your metabolites that matched the set;
- the enrichment p-value;
- the FDR-adjusted p-value; and
- **which actual metabolites are responsible for the signal**.

Always return from the pathway result to the individual metabolites and ask whether the direction of change makes biological sense.

### Checkpoint 9 — Enrichment

> **Q9.1.** Which pathway or metabolite set is most strongly enriched?  
> **Q9.2.** How many of your metabolites contribute to that result?  
> **Q9.3.** Are those metabolites increased, decreased, or mixed?

(OPTIONAL)
> **Q9.4.** How many metabolites failed to map in MetaboAnalyst?  
> **Q9.5.** Why is the measured-metabolite background preferable to using the entire known metabolome?

## 15. Bonus — Connect metabolomics and transcriptomics

This course also includes transcriptomics, so there is a natural **multi-omics extension**.

Since we are using a different metabolomics dataset from the one used in our Transcriptomics workshop, the DESeq2 output for the matching transcriptomics dataset is included in the **`Results_Transcriptomics`** folder.

Remember that, by default, we are using the **Bleo14D vs Control** comparison.

If you would like to redo the Transcriptomics analysis using the matching dataset from the same paper, please download **Supplementary Table 1** here:

https://advanced.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fadvs.202207454&file=advs5487-sup-0002-TableS1.xlsx

MetaboAnalyst provides **Joint Pathway Analysis**, where a list of significant genes and a list of significant metabolites can be mapped together onto metabolic pathways.

Conceptually:

```text
Transcriptomics
significant genes
        \
         \
          → shared metabolic pathways
         /
        /
Metabolomics
significant metabolites
```

This can be more informative than asking whether a pathway is supported by only one omics layer.

For example, a pathway becomes especially interesting when:

- several genes encoding enzymes in the pathway change; **and**
- metabolites belonging to the same pathway also change.

### Important limitation

Transcriptomics often measures most genes, while metabolomics typically measures only a fraction of the metabolome. Therefore, absence of a metabolite from your dataset does **not** mean that its pathway is inactive.

### Optional class exercise — Joint pathway analysis

If you still have the significant gene list from the transcriptomics workshop:

1. export the significant genes from the transcriptomics analysis; ("Results_Transcriptomics" folder, filter it based on "padj")
2. use `Results/metaboanalyst_significant_metabolites.txt` from this notebook;
3. open **MetaboAnalyst → Pathway Analysis → Joint Pathway Analysis**; (ID Type: "Compound Name" and "Official Gene Symbol")
4. upload/paste the gene and metabolite lists;
5. select the appropriate organism; (**Mus Musculus** for this example)
6. Choose whether you want to use just Metabolic Pathways or All pathways. Make sure you select the "Integrated" pathways. Explore the options, but you can use the defaults for now.
7. inspect pathways supported by both molecular layers.

> **Question:** Does integrating the two omics layers strengthen the biological story, or do the transcriptomics and metabolomics results point to different processes?

## 16. What if the metabolites are not identified?

For **untargeted LC-MS**, you may have hundreds or thousands of features represented only by **m/z** values (and sometimes retention times), rather than confident metabolite identities.

In that situation, standard name-based enrichment is not appropriate.

MetaboAnalyst includes functional-analysis approaches based on algorithms such as **mummichog**, which try to infer pathway activity directly from ranked mass-spectrometry features without requiring every feature to be identified first.

That is very useful, but it is a more advanced topic because interpretation depends on ion mode, mass accuracy, adducts, feature redundancy, organism, and other MS-specific details.

For this introductory workshop:

> **identified metabolites → enrichment/pathway analysis**  
> **unidentified m/z features → mention mummichog as an advanced extension**

### Checkpoint 9 — Biological interpretation

> **Q9.1.** Are the metabolite names in your dataset confirmed identities, putative annotations, or unknown features?  
> **Q9.2.** Why does annotation confidence matter before pathway analysis?  
> **Q9.3.** Which two or three metabolites would you investigate first, and why? Consider both statistical evidence and effect size.

---

## Workflow summary

You have now moved through a basic downstream metabolomics workflow:

```text
Load abundance data + metadata
        |
        v
Check sample IDs
        |
        v
Missingness / basic QC
        |
        v
Imputation + log transformation
        |
        v
Autoscaling + PCA
        |
        v
Sample correlation
        |
        v
Two-group differential abundance
        |
        v
Volcano plot + heatmap
        |
        v
Metabolite identification + biological interpretation
```

### Main conceptual lessons

1. **Sample IDs must be checked explicitly.**
2. **Metabolomics preprocessing choices depend on how the data were generated.**
3. **Missing values and zeros need biological/technical interpretation, not blind replacement.**
4. **Normalization, transformation, and scaling solve different problems.**
5. **PCA is useful for biological structure, outliers, pooled QC behaviour, and batch effects.**
6. **Multiple-testing correction matters when many metabolites are tested.**
7. **Effect size and statistical significance are not the same thing.**
8. **A statistically significant feature is only biologically interpretable if its metabolite identity is sufficiently reliable.**

---

## Optional extension — What would change in a more advanced analysis?

Depending on the experiment, a more complete metabolomics analysis might include:

- pooled-QC coefficient-of-variation filtering;
- blank/contaminant filtering;
- signal-drift correction using injection order;
- batch correction;
- normalization to sample amount or internal standards;
- paired or repeated-measures models;
- adjustment for age, sex, batch, or other covariates;
- targeted pathway analysis using validated metabolite identifiers; or
- integration with transcriptomics, proteomics, or other omics layers.

These are important topics, but the correct method depends on the experimental design and measurement platform.

---

## References and further reading

- MetaboAnalyst tutorials: https://www.metaboanalyst.ca/docs/Tutorials.xhtml
- MetaboAnalyst normalization overview: https://www.metaboanalyst.ca/MetaboAnalyst/Secure/process/NormalizationView.xhtml
- Alseekh et al. *Mass spectrometry-based metabolomics: a guide for annotation, quantification and best reporting practices*. Nature Methods (2021): https://www.nature.com/articles/s41592-021-01197-1

```python

```

```python

```
