# Pre-course setup

Below are useful web-based resources for working with **transcriptomics, proteomics, metabolomics, and genomics** data.

You do **not** need to learn all of these before the course. This page is intended as a collection of tools and resources that you can explore during and after the course.

## Table of contents

- [Public omics data repositories](#public-omics-data-repositories)
- [Omics data analysis](#omics-data-analysis)
- [Pathway and functional interpretation](#pathway-and-functional-interpretation)
- [Networks and systems biology](#networks-and-systems-biology)
- [Genomics and variant interpretation](#genomics-and-variant-interpretation)
- [General bioinformatics platforms](#general-bioinformatics-platforms)
- [Cloud coding environments](#cloud-coding-environments)
- [Where do these tools fit?](#where-do-these-tools-fit)

---

## Public omics data repositories

### [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/)

GEO is one of the major public repositories for functional genomics data.

It contains data from:

- RNA-seq
- microarrays
- single-cell studies
- epigenomics
- other high-throughput functional genomics experiments

You can search GEO for published datasets, inspect study designs, and download processed or raw data for your own analyses.

**Mainly useful for:** Transcriptomics and functional genomics

---

### [PRIDE](https://www.ebi.ac.uk/pride/)

PRIDE is a major public repository for mass-spectrometry-based proteomics data.

It contains:

- protein-identification datasets
- quantitative proteomics data
- peptide and protein measurements
- raw mass-spectrometry files

It is particularly useful for finding published proteomics datasets that can be reanalysed.

**Mainly useful for:** Proteomics

---

### [MetaboLights](https://www.ebi.ac.uk/metabolights/)

MetaboLights is a public repository for metabolomics studies.

It contains experimental metadata, metabolite measurements, and raw data from a wide range of metabolomics experiments.

**Mainly useful for:** Metabolomics

---

## Omics data analysis

### [BioJupies](https://maayanlab.cloud/biojupies/)

A web-based platform for analysing and visualising RNA-seq data without writing code.

BioJupies can be used for:

- differential gene-expression analysis
- exploratory visualisation
- clustering
- enrichment analysis
- generating reproducible Jupyter notebooks

**Mainly useful for:** Transcriptomics

---

### [GEO2R](https://www.ncbi.nlm.nih.gov/geo/geo2r/)

GEO2R is an analysis interface connected directly to GEO.

It allows you to compare groups of samples from a GEO study and identify differentially expressed genes without downloading and analysing the dataset manually.

**Mainly useful for:** Transcriptomics

---

### [MetaboAnalyst](https://www.metaboanalyst.ca/)

A comprehensive web platform for metabolomics and multi-omics analysis.

It includes tools for:

- data processing and normalisation
- PCA
- clustering
- statistical testing
- biomarker analysis
- pathway analysis
- metabolite-set enrichment
- multi-omics integration

**Mainly useful for:** Metabolomics and multi-omics

---

## Pathway and functional interpretation

### [Enrichr](https://maayanlab.cloud/Enrichr/)

A gene-set enrichment analysis tool.

Provide a list of genes and Enrichr can identify associations with:

- biological pathways
- Gene Ontology terms
- transcription factors
- diseases
- cell types
- protein interactions
- many other gene-set libraries

It is particularly useful after obtaining a list of differentially expressed genes or proteins.

**Useful for:** Transcriptomics, Proteomics, and Genomics

---

### [KEGG](https://www.kegg.jp/)

KEGG — Kyoto Encyclopedia of Genes and Genomes — is a major resource linking genes, proteins, metabolites, reactions, diseases, and biological pathways.

Useful KEGG resources include:

- **KEGG PATHWAY** — biological pathway diagrams
- **KEGG GENES** — genes and proteins
- **KEGG COMPOUND** — metabolites and small molecules
- **KEGG REACTION** — biochemical reactions
- **KEGG ORTHOLOGY (KO)** — functional ortholog groups
- **KEGG DISEASE** — molecular information related to diseases

KEGG is especially useful for connecting molecular measurements to known biological processes.

**Useful for:** Transcriptomics, Proteomics, Metabolomics, and Genomics

---

### [KEGG Mapper Color](https://www.genome.jp/kegg/mapper/color.html)

KEGG Mapper Color allows you to place your own experimental results directly onto KEGG pathway diagrams.

For example, you can:

- highlight genes that are up- or down-regulated
- color proteins according to abundance
- highlight metabolites
- display different experimental groups using different colors

The tool accepts KEGG identifiers and optional colors, then highlights the corresponding objects on KEGG pathways.

This is particularly useful for visualising **omics results in biological pathway context**.

**Useful for:** Transcriptomics, Proteomics, Metabolomics, and multi-omics

---

### [Reactome](https://reactome.org/)

Reactome is a curated database of biological pathways.

You can submit lists of genes or proteins and identify pathways that are over-represented in your dataset.

It is useful for:

- pathway enrichment
- pathway browsing
- visualising molecules within pathways
- exploring relationships between biological processes

**Useful for:** Transcriptomics, Proteomics, and Genomics

---

### [g:Profiler](https://biit.cs.ut.ee/gprofiler/)

A web tool for functional enrichment analysis.

It can analyse lists of genes against resources such as:

- Gene Ontology
- biological pathways
- regulatory motifs
- protein databases

It also provides tools for converting between different gene identifiers.

**Useful for:** Transcriptomics, Proteomics, and Genomics

---

## Networks and systems biology

### [STRING](https://string-db.org/)

STRING is a database and analysis platform for protein-protein association networks.

Provide a list of genes or proteins and STRING can:

- construct interaction networks
- identify connected groups
- perform functional enrichment
- suggest related proteins
- visualize biological relationships

It is particularly useful for moving from a list of significant molecules toward a **systems-level interpretation**.

**Useful for:** Proteomics, Transcriptomics, and Systems Biology

---

### [Cytoscape](https://cytoscape.org/)

Cytoscape is a widely used platform for visualising and analysing biological networks.

It can be used to explore:

- protein-protein interaction networks
- gene-regulatory networks
- pathway networks
- metabolite networks
- integrated multi-omics networks

Cytoscape itself is primarily a desktop application, but it is worth knowing because it is commonly used in systems biology.

**Useful for:** Systems Biology and multi-omics

---

### [NDEx](https://www.ndexbio.org/)

NDEx — Network Data Exchange — is a public platform for storing, sharing, and exploring biological networks.

It integrates particularly well with Cytoscape.

**Useful for:** Network biology and Systems Biology

---

## Genomics and variant interpretation

### [UCSC Genome Browser](https://genome.ucsc.edu/)

An interactive genome browser for exploring genomic regions.

You can inspect:

- genes
- variants
- regulatory regions
- conservation
- chromatin annotations
- sequencing tracks
- many other genomic features

A useful way to think about the Genome Browser is:

> **Where in the genome is my gene or variant, and what is around it?**

**Mainly useful for:** Genomics

---

### [Ensembl](https://www.ensembl.org/)

Ensembl provides genome annotation and comparative genomics resources for many organisms.

You can explore:

- genes and transcripts
- genomic regions
- sequence
- genetic variants
- orthologues
- regulatory features
- comparative genomics

**Mainly useful for:** Genomics

---

### [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/Tools/VEP)

VEP predicts the possible functional consequences of genetic variants.

For example, it can determine whether a variant:

- occurs inside a gene
- changes an amino acid
- creates a stop codon
- occurs in an intron
- occurs in a regulatory region
- has existing population or clinical annotations

**Mainly useful for:** Genomics and variant interpretation

---

### [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

BLAST compares DNA or protein sequences against sequence databases.

Typical questions include:

- What gene does this sequence resemble?
- Which organisms contain similar sequences?
- Does this protein have homologues?
- How conserved is this sequence?

**Useful for:** Genomics and Proteomics

---

## General bioinformatics platforms

### [Galaxy](https://usegalaxy.org/)

Galaxy is a web-based environment for running bioinformatics workflows without installing individual command-line programs.

Galaxy contains thousands of tools and supports workflows for:

- RNA-seq
- genomics
- variant calling
- proteomics
- metabolomics
- sequence analysis
- statistics
- visualisation

Analyses are stored as histories and can be assembled into reproducible workflows.

**Useful for:** Transcriptomics, Proteomics, Metabolomics, and Genomics

---

## Cloud coding environments

These platforms allow you to run Python or R analyses without configuring a full programming environment on your computer.

### [Posit Cloud](https://posit.cloud/)

A browser-based computational environment supporting Jupyter, Python, and R.

For this course, **Posit Cloud is the recommended option** if you do not want to configure a local environment.

It provides a consistent environment for working with notebooks and course materials.

---

### [Google Colab](https://colab.research.google.com/)

Google Colab is a cloud-based Jupyter notebook environment.

It is useful for:

- Python analysis
- interactive notebooks
- sharing analyses
- machine learning
- temporary GPU access

No local Python installation is required.

---

## Where do these tools fit?

A simplified omics workflow might look like:

**Public data → Quality control → Statistical analysis → Significant features → Functional interpretation → Networks and pathways → Biological conclusions**

| Task | Example tools |
|---|---|
| Find transcriptomics data | GEO |
| Find proteomics data | PRIDE |
| Find metabolomics data | MetaboLights |
| Analyse public expression data | GEO2R, BioJupies |
| Analyse metabolomics data | MetaboAnalyst |
| Run general bioinformatics workflows | Galaxy |
| Functional enrichment | Enrichr, g:Profiler, Reactome |
| Explore biological pathways | KEGG, Reactome |
| Overlay omics data onto pathways | KEGG Mapper Color |
| Protein interaction networks | STRING |
| Biological network analysis | Cytoscape, NDEx |
| Explore genomic regions | UCSC Genome Browser, Ensembl |
| Predict variant consequences | Ensembl VEP |
| Compare DNA/protein sequences | NCBI BLAST |
| Write your own analysis | Posit Cloud, Google Colab |

The goal is **not to learn every tool**. Instead, you should begin to recognise what kinds of biological questions different tools can help answer.