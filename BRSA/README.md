<div align="center">

```
██████╗ ██████╗ ███████╗ █████╗ 
██╔══██╗██╔══██╗██╔════╝██╔══██╗
██████╔╝██████╔╝███████╗███████║
██╔══██╗██╔══██╗╚════██║██╔══██║
██████╔╝██║  ██║███████║██║  ██║
╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝
```

# 🕷️ BRSA — Basic RNA-Seq Analysis 🕷️
### *"From the Graveyard of Raw Reads, Life Shall Rise."*

[![Python](https://img.shields.io/badge/Python-3.9%2B-8a2be2?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-ff6600?style=for-the-badge)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Alive%20(barely)-darkred?style=for-the-badge)]()
[![Os](https://img.shields.io/badge/Plataforma-Windows%20%7C%20Linux%20%7C%20macOS-8a2be2?style=for-the-badge)]()

*A complete, unholy RNA-Seq analysis pipeline conjured from the depths of bioinformatics arcana — featuring differential expression, machine learning, protein interactomes, and survival analysis.*

---

</div>

## 🪦 Table of Contents

- [The Curse (Overview)](#-the-curse-overview)
- [The Coven (Features)](#-the-coven-features)
- [Summoning Ingredients (Installation)](#-summoning-ingredients-installation)
- [Raising the Dead (Usage)](#-raising-the-dead-usage)
  - [The Dark CLI Ritual](#the-dark-cli-ritual)
  - [The Graphical Séance (GUI)](#the-graphical-séance-gui)
- [The Sacred Architecture](#-the-sacred-architecture)
- [The Grimoire of Methods](#-the-grimoire-of-methods)
  - [Quality Control & Preprocessing](#i-quality-control--preprocessing)
  - [Alignment & Quantification](#ii-alignment--quantification)
  - [Normalization](#iii-normalization)
  - [Differential Expression Analysis](#iv-differential-expression-analysis)
  - [Dimensionality Reduction](#v-dimensionality-reduction)
  - [Machine Learning Golems](#vi-machine-learning-golems)
  - [Co-expression Network Analysis](#vii-co-expression-network-analysis)
  - [Protein-Protein Interactome](#viii-protein-protein-interactome)
  - [Survival Analysis](#ix-survival-analysis)
  - [Gene Set Enrichment Analysis](#x-gene-set-enrichment-analysis)
- [The Eldritch Parameters (CLI Reference)](#-the-eldritch-parameters-cli-reference)
- [The Cauldron Outputs](#-the-cauldron-outputs)
- [Reproducibility Scrolls](#-reproducibility-scrolls)
- [References & Ancient Tomes](#-references--ancient-tomes)
- [License & Dark Pact](#-license--dark-pact)
- [The Necromancer (Author)](#-the-necromancer-author)

---

## 🕸️ The Curse (Overview)

**BRSA** is a comprehensive, end-to-end RNA sequencing analysis pipeline developed at the **Universidade Federal de Minas Gerais (UFMG)**. It was conjured to guide researchers through the entire transcriptomic dark ritual — from raw FASTQ files all the way to publication-ready figures, interactive HTML artifacts, and reproducibility logs.

The pipeline integrates classical statistical methods with modern machine learning and network biology, allowing mortal bioinformaticians to wield an arsenal of analytical weapons through a single, unified interface — be it a terminal incantation or a full graphical séance.

> *"That which does not kill the pipeline only makes it more statistically significant."*

---

## 🦇 The Coven (Features)

| Module | Power Granted |
|--------|--------------|
| 🧟 **Quality Control** | FastQC interrogation + MultiQC all-seeing eye |
| ⚗️ **Adapter Trimming** | fastp guillotine for adapter & quality purging |
| 🗺️ **Alignment** | HISAT2 splice-aware alignment vortex |
| 💀 **Pseudo-alignment** | Kallisto ultra-fast quantification ghost |
| 🔮 **Normalization** | CPM, TPM, and DESeq2 variance stabilization |
| 🌋 **Differential Expression** | DESeq2-inspired negative binomial testing |
| 🕳️ **Dimensionality Reduction** | PCA, t-SNE, and UMAP projections |
| 🤖 **Machine Learning** | Random Forest + SVM for biomarker discovery |
| 🕸️ **PPI Interactome** | STRING DB protein network summoning |
| 🧬 **Isoform Analysis** | Alternative splicing Frankenstein limb detection |
| 💉 **Survival Analysis** | Kaplan-Meier mortality doom curves |
| 📜 **GSEA** | Gene Set Enrichment with Enrichr/GSEApy |
| 🗺️ **Interactive Plots** | Plotly-forged HTML artifacts |
| 📖 **Reproducibility** | Complete parameter + environment logging |

---

## 🧪 Summoning Ingredients (Installation)

### Prerequisites — The External Daemons

The following command-line spirits must be installed and accessible in your system's `PATH` before invoking BRSA:

```bash
# Quality Control Daemons
fastqc --version
multiqc --version

# Trimming & Alignment Spirits
fastp --version
hisat2 --version
samtools --version
featureCounts -v      # Subread package

# Pseudo-alignment Ghost
kallisto version
```

### Python Alchemical Environment

```bash
# Clone the forbidden repository
git clone https://github.com/VictorCaricatte/BasicBioinfo/tree/main/BRSA/src.git
cd BRSA

# Conjure a virtual sanctuary (recommended)
python -m venv brsa_coven
source brsa_coven/bin/activate          # Linux/macOS
# brsa_coven\Scripts\activate.bat       # Windows (if you dare)

# Install the required arcane libraries
pip install -r requirements.txt
```

### `requirements.txt` — The Forbidden Ingredients

```
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.6.0
seaborn>=0.12.0
scipy>=1.9.0
scikit-learn>=1.1.0
pydeseq2>=0.3.0
gseapy>=1.0.0
plotly>=5.10.0
mplcursors>=0.5.0
networkx>=2.8.0
requests>=2.28.0
lifelines>=0.27.0
umap-learn>=0.5.3
PyQt5>=5.15.0        # For the GUI séance only
```

---

## 💀 Raising the Dead (Usage)

### The Dark CLI Ritual

For those who prefer to work in the shadows of the terminal, BRSA can be invoked directly without loading any graphical spirits:

```bash
python BRSA.py \
  --counts data/counts_matrix.csv \
  --meta data/metadata.csv \
  --cond1 Tumor \
  --cond2 Normal \
  --norm DESeq2 \
  --dim_red UMAP \
  --pval 0.05 \
  --min_counts 10 \
  --batch batch_column \
  --ml \
  --ppi
```

Upon completion, two dark artifacts shall materialize:
- `CLI_differential_results.csv` — The tome of significant genes
- `CLI_reproducibility_log.txt` — The eternal reproducibility scroll

### The Graphical Séance (GUI)

For those who prefer their rituals with buttons and sliders:

```bash
python interface.py
```

The graphical portal shall materialize, allowing you to load your matrices, configure your parameters, and summon all analytical modules with a click.

---

## 🏚️ The Sacred Architecture

```
BRSA/
│
├── 🧠 BRSA.py              — Main entry point: CLI vs GUI dispatcher
├── 🎛️ interface.py         — PyQt5 graphical séance (GUI)
├── ⚙️ job.py               — FrankensteinBioinformaticsOverlordThread
│                              (Worker thread: runs all analyses)
├── 🔬 biotools.py          — External tool wrappers & API callers
│                              (FastQC, HISAT2, Kallisto, STRING DB...)
├── 📊 plot.py              — All visualization functions
│                              (Volcano, MA, UMAP, Clustermap, GSEA...)
└── 📜 args.py              — CLI parser & SacredSessionScroll dataclass
```

**Data flow through the catacombs:**

```
FASTQ Files ──► FastQC/fastp ──► HISAT2/Kallisto ──► Counts Matrix
                                                           │
Metadata CSV ─────────────────────────────────────────────┤
                                                           ▼
                                                   Normalization
                                                  (CPM/TPM/DESeq2)
                                                           │
                                        ┌──────────────────┤
                                        ▼                  ▼
                                   Dim. Reduction     DESeq2 DE Test
                                  (PCA/t-SNE/UMAP)         │
                                        │                  ▼
                                        │         Volcano / MA Plot
                                        │         GSEA / ORA
                                        │         ML Golems
                                        │         PPI Network
                                        │         Survival Analysis
                                        └──────► Interactive HTML Artifacts
```

---

## 📚 The Grimoire of Methods

> *Every dark ritual has its justification. Here lie the academic tomes that birthed these methods.*

---

### I. Quality Control & Preprocessing

**FastQC** is employed as the first inquisition upon raw sequencing reads, evaluating per-base quality scores, GC content distribution, adapter contamination, and duplication levels. Results from multiple samples are aggregated by **MultiQC** into a unified surveillance report.

Adapter trimming and quality filtering are executed via **fastp**, which simultaneously performs sliding-window quality trimming (Phred score ≥ 20), adapter auto-detection, and generation of QC reports in HTML and JSON format.

> Andrews, S. (2010). *FastQC: A Quality Control Tool for High Throughput Sequence Data*. Bioinformatics Group, Babraham Institute. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
>
> Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
>
> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

---

### II. Alignment & Quantification

Reads are aligned to the reference genome using **HISAT2**, a graph-based alignment algorithm that efficiently handles spliced alignments across exon-exon junctions using a Hierarchical Graph FM Index (HGFM). Aligned reads are converted and coordinate-sorted using **SAMtools**.

Gene-level read counts are harvested via **featureCounts** (Subread package), which assigns aligned reads to genomic features based on GTF/GFF3 annotation using a union-overlap strategy.

As an alternative, **Kallisto** provides ultra-fast transcript-level abundance quantification via pseudoalignment — mapping reads to a k-mer index of the reference transcriptome without full alignment, yielding TPM estimates with bootstrap uncertainty quantification.

> Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. *Nature Methods*, 12(4), 357–360. https://doi.org/10.1038/nmeth.3317
>
> Li, H., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
>
> Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656
>
> Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

---

### III. Normalization

BRSA supports three normalization exorcisms to remove sequencing depth demons and compositional biases:

- **CPM (Counts Per Million):** Scales raw counts by total library size, correcting for sequencing depth. Suitable for exploratory analysis and visualization.
- **TPM (Transcripts Per Million):** Normalizes first by gene length (RPK), then by total RPK sum. Preserves cross-sample comparability for transcript abundance.
- **DESeq2 Variance Stabilizing Transformation (VST):** Applies a negative binomial model-based variance stabilization, shrinking the mean-variance dependence. Recommended for downstream statistical tests and clustering.

> Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
>
> Wagner, G. P., Kin, K., & Lynch, V. J. (2012). Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory in Biosciences*, 131(4), 281–285. https://doi.org/10.1007/s12064-012-0162-3

---

### IV. Differential Expression Analysis

Differential expression between two conditions is tested using the **negative binomial generalized linear model** framework as implemented in the DESeq2 paradigm. The pipeline:

1. Filters lowly expressed genes (default: < 10 reads across all samples)
2. Estimates size factors using the **median-of-ratios** normalization
3. Estimates gene-wise dispersion using empirical Bayes shrinkage
4. Tests the null hypothesis (Log₂FC = 0) via a **Wald test**
5. Adjusts p-values using the **Benjamini-Hochberg** procedure to control the False Discovery Rate (FDR)

Results are visualized as **Volcano Plots** (Log₂FC vs −log₁₀ FDR) and **MA Constellation Plots** (mean expression vs Log₂FC), with interactive hover annotations powered by `mplcursors`.

> Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
>
> Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

---

### V. Dimensionality Reduction

Three spectral weapons are available for projecting high-dimensional expression space into visible 2D void:

- **PCA (Principal Component Analysis):** Linear decomposition via singular value decomposition (SVD), preserving global variance structure. Explained variance per component is reported.
- **t-SNE (t-Distributed Stochastic Neighbor Embedding):** Non-linear dimensionality reduction that preserves local neighborhood structure by minimizing KL-divergence between high- and low-dimensional probability distributions.
- **UMAP (Uniform Manifold Approximation and Projection):** Topology-preserving manifold learning based on Riemannian geometry and fuzzy simplicial set theory, balancing local and global structure with superior computational speed.

> Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. *The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science*, 2(11), 559–572.
>
> van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE. *Journal of Machine Learning Research*, 9, 2579–2605. https://jmlr.org/papers/v9/vandermaaten08a.html
>
> McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv preprint*. https://arxiv.org/abs/1802.03426

---

### VI. Machine Learning Golems

When invoked with `--ml`, BRSA trains two predictive golems on the normalized expression matrix:

- **Random Forest:** An ensemble of decorrelated decision trees trained via bootstrap aggregation (bagging). Feature importances are computed as mean decrease in impurity (Gini importance) across all trees, ranking genes by their predictive power for condition classification.
- **Support Vector Machine (SVM):** A maximum-margin classifier using a radial basis function (RBF) kernel, projecting the data into a higher-dimensional feature space to find the optimal separating hyperplane.

Feature importance scores from Random Forest are visualized as horizontal bar charts, enabling biomarker discovery.

> Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
>
> Cortes, C., & Vapnik, V. (1995). Support-vector networks. *Machine Learning*, 20(3), 273–297. https://doi.org/10.1007/BF00994018

---

### VII. Co-expression Network Analysis

Hierarchical clustering of top-variable genes is performed and visualized as a **clustermap** (heatmap with dual dendrograms) using Ward linkage on Euclidean distances. Sample-level clustering reveals batch effects and biological groupings, while gene-level clustering exposes co-regulated transcriptional modules.

This approach follows the principles of **Weighted Gene Co-expression Network Analysis (WGCNA)**, wherein genes with correlated expression profiles across samples are grouped into functional modules.

> Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559
>
> Ward, J. H. (1963). Hierarchical Grouping to Optimize an Objective Function. *Journal of the American Statistical Association*, 58(301), 236–244. https://doi.org/10.2307/2282967

---

### VIII. Protein-Protein Interactome

When `--ppi` is enabled, BRSA queries the **STRING database** REST API to retrieve protein-protein interaction data for the list of significant differentially expressed genes. Interactions are filtered by a confidence score threshold (default: 0.4) and assembled into a **NetworkX** undirected graph, rendered as a force-directed spring layout network — the Hairball of Chaos.

> Szklarczyk, D., et al. (2023). The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any of 14000+ organisms. *Nucleic Acids Research*, 51(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000
>
> Hagberg, A. A., Schult, D. A., & Swart, P. J. (2008). Exploring Network Structure, Dynamics, and Function using NetworkX. In *Proceedings of the 7th Python in Science Conference (SciPy2008)*, pp. 11–15.

---

### IX. Survival Analysis

BRSA implements **Kaplan-Meier survival curve estimation** to correlate individual gene expression levels with patient outcome data. Samples are stratified into "High" and "Low" expression groups based on median expression. Survival functions are estimated non-parametrically for each group and plotted as step functions with confidence intervals.

This analysis requires clinical metadata columns for survival time and event occurrence (censoring indicator).

> Kaplan, E. L., & Meier, P. (1958). Nonparametric Estimation from Incomplete Observations. *Journal of the American Statistical Association*, 53(282), 457–481. https://doi.org/10.2307/2281868
>
> Davidson-Pilon, C. (2019). lifelines: survival analysis in Python. *Journal of Open Source Software*, 4(40), 1317. https://doi.org/10.21105/joss.01317

---

### X. Gene Set Enrichment Analysis

BRSA employs **GSEApy** to perform both:

- **GSEA (Gene Set Enrichment Analysis):** Ranks all detected genes by their signed Log₂ fold change and tests for coordinated enrichment of curated gene sets at the extremes of the ranked list, computing a Running Enrichment Score (RES) and normalized enrichment score (NES).
- **Over-Representation Analysis (ORA):** Tests for significant overlap between the significant gene list and annotated gene sets using Fisher's exact test with FDR correction.

Multiple gene set databases are supported, including GO Biological Process, KEGG, Reactome, and MSigDB Hallmarks.

> Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102
>
> Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*, 39(1), btac757. https://doi.org/10.1093/bioinformatics/btac757

---

## ⚰️ The Eldritch Parameters (CLI Reference)

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--counts` | `str` | *required* | Path to counts matrix CSV (genes × samples) |
| `--meta` | `str` | *required* | Path to metadata CSV (samples × covariates) |
| `--cond1` | `str` | *required* | Reference condition (Group A / Alpha) |
| `--cond2` | `str` | *required* | Test condition (Group B / Omega) |
| `--min_counts` | `int` | `10` | Minimum read threshold for gene filtering |
| `--pval` | `float` | `0.05` | FDR cutoff for significance |
| `--norm` | `str` | `DESeq2` | Normalization: `CPM`, `TPM`, or `DESeq2` |
| `--batch` | `str` | `None` | Metadata column for batch effect correction |
| `--dim_red` | `str` | `PCA` | Dimensionality reduction: `PCA`, `t-SNE`, `UMAP` |
| `--ml` | flag | `False` | Invoke ML golems (Random Forest + SVM) |
| `--ppi` | flag | `False` | Summon STRING DB protein network |
| `--isoforms` | flag | `False` | Hunt alternative splicing events |

### Expected Input Formats

**Counts Matrix (`--counts`):**
```
gene_id,Sample_1,Sample_2,Sample_3,...
ENSG00000000003,1500,2300,800,...
ENSG00000000005,45,12,67,...
```

**Metadata (`--meta`):**
```
sample,condition,batch
Sample_1,Tumor,batch_A
Sample_2,Normal,batch_A
Sample_3,Tumor,batch_B
```

---

## 🏺 The Cauldron Outputs

Upon successful completion, the following dark artifacts shall be forged:

```
outputs/
├── 📄 CLI_differential_results.csv       — Full DE gene table (LFC, p-adj, significance)
├── 📄 CLI_reproducibility_log.txt        — Parameter registry + pip freeze snapshot
│
└── 📁 interactive_artifacts/
    ├── 🌐 interactive_pca.html           — Plotly interactive PCA/t-SNE/UMAP
    ├── 🌐 interactive_volcano.html       — Interactive volcano with gene annotations
    └── 🌐 interactive_ma_constellation.html — Interactive MA plot
```

**DE Results Table columns:**

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier |
| `log2fc` | Log₂ Fold Change (cond2 vs cond1) |
| `p-adj` | Benjamini-Hochberg adjusted p-value |
| `mean_1` | Mean normalized expression in cond1 |
| `mean_2` | Mean normalized expression in cond2 |
| `significant` | Boolean: True if FDR < threshold |

---

## 📜 Reproducibility Scrolls

Every BRSA ritual automatically carves a reproducibility log containing:

- All applied parameters and thresholds (complete `SacredSessionScroll` state)
- The full output of `pip freeze` capturing exact library versions
- Timestamp and execution context

This ensures that the exact analytical conditions may be reconstructed by any future necromancer — or peer reviewer — attempting to replicate the findings.

---

## 📖 References & Ancient Tomes

1. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Ewels, P., et al. (2016). MultiQC. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
3. Chen, S., et al. (2018). fastp. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560
4. Kim, D., et al. (2015). HISAT. *Nature Methods*, 12(4), 357–360. https://doi.org/10.1038/nmeth.3317
5. Li, H., et al. (2009). SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
6. Liao, Y., et al. (2014). featureCounts. *Bioinformatics*, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656
7. Bray, N. L., et al. (2016). Kallisto. *Nature Biotechnology*, 34(5), 525–527. https://doi.org/10.1038/nbt.3519
8. Love, M. I., et al. (2014). DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
9. Benjamini, Y., & Hochberg, Y. (1995). FDR. *JRSS-B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
10. van der Maaten, L., & Hinton, G. (2008). t-SNE. *JMLR*, 9, 2579–2605. https://jmlr.org/papers/v9/vandermaaten08a.html
11. McInnes, L., et al. (2018). UMAP. *arXiv*. https://arxiv.org/abs/1802.03426
12. Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
13. Cortes, C., & Vapnik, V. (1995). SVM. *Machine Learning*, 20(3), 273–297. https://doi.org/10.1007/BF00994018
14. Langfelder, P., & Horvath, S. (2008). WGCNA. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559
15. Szklarczyk, D., et al. (2023). STRING 2023. *Nucleic Acids Research*, 51(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000
16. Kaplan, E. L., & Meier, P. (1958). Kaplan-Meier. *JASA*, 53(282), 457–481. https://doi.org/10.2307/2281868
17. Davidson-Pilon, C. (2019). lifelines. *JOSS*, 4(40), 1317. https://doi.org/10.21105/joss.01317
18. Subramanian, A., et al. (2005). GSEA. *PNAS*, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102
19. Fang, Z., et al. (2023). GSEApy. *Bioinformatics*, 39(1), btac757. https://doi.org/10.1093/bioinformatics/btac757

---

## ⚖️ License & Dark Pact

This software is distributed under the **MIT License**. See `LICENSE` for details.

By using BRSA, you implicitly agree to cite the relevant methodological references in any resulting publications, lest your manuscripts be cursed with peer reviewer #3.

---

## 🧙 The Necromancer (Author)

**Victor S. Caricatte De Araújo (Frankestein)**  
*Universidade Federal de Minas Gerais (UFMG)*  
Version: `0.9.0` 

---

<div align="center">

*"In the graveyard of sequencing errors and batch effects, BRSA lights the lantern of differential expression."*

🕯️ *May your p-values be low and your sample sizes be high.* 🕯️

</div>
