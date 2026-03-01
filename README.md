# Basic Tools for Bioinformaticians

Hello friends,

This repository contains some tools I built while working with bioinformatics. They were created out of research demands I had and still have. Along the way I ended up giving them fun themes, and I improved a few things to make them more complete — but they are still basic bioinformatics analysis tools at their core. Feel free to use them, fork them, and build your own from them. You are also welcome to report bugs, make suggestions, and request improvements or new features (new analyses and functionalities), as I may add them over time.


Best regards,
VictorSC

---

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen)]()
[![Institution](https://img.shields.io/badge/Institution-UFMG-blue)]()

**Repository:** [https://github.com/VictorCaricatte/BasicBioinfo](https://github.com/VictorCaricatte/BasicBioinfo)

---

## Table of Contents

- [Tools Overview](#tools-overview)
- [BRSA — Basic RNA-Seq Analysis](#brsa--basic-rna-seq-analysis)
- [DELIMITER — Delimited File Filter](#delimiter--delimited-file-filter)
- [LenghtCount — Genomic & Proteomic Analysis Suite](#lenghtcount--genomic--proteomic-analysis-suite)
- [PALINDROME — DNA Palindrome Analyzer](#palindrome--dna-palindrome-analyzer)
- [Installation](#installation)
- [License](#license)
- [Author](#author)

---

## Tools Overview

| Tool | Description |
|------|-------------|
| [BRSA](#brsa--basic-rna-seq-analysis) | End-to-end RNA-Seq analysis pipeline |
| [DELIMITER](#delimiter--delimited-file-filter) | Desktop GUI for filtering and transforming tabular/genomic data files |
| [LenghtCount](#lenghtcount--genomic--proteomic-analysis-suite) | Multi-module genomic and proteomic characterization suite |
| [PALINDROME](#palindrome--dna-palindrome-analyzer) | DNA palindromic sequence detection and analysis |

---

## BRSA — Basic RNA-Seq Analysis

**Location:** `BRSA/`

BRSA is a comprehensive, end-to-end RNA sequencing analysis pipeline developed at the Universidade Federal de Minas Gerais (UFMG). It guides researchers through the complete transcriptomic workflow — from raw count matrices all the way to publication-ready figures, interactive HTML reports, and reproducibility logs.

The pipeline integrates classical statistical methods with machine learning and network biology, exposing all analytical modules through a unified CLI or a PyQt5-based graphical interface.

### Key Features

- **Differential Expression Analysis** using a DESeq2-inspired negative binomial model with Benjamini-Hochberg FDR correction
- **Normalization** via CPM, TPM, or DESeq2 variance stabilization
- **Dimensionality Reduction** with PCA, t-SNE, and UMAP
- **Machine Learning** classifiers (Random Forest + SVM) for biomarker discovery
- **Co-expression Network Analysis** based on WGCNA-style hierarchical clustering
- **Protein-Protein Interactome** queried from the STRING database
- **Survival Analysis** using Kaplan-Meier curves stratified by gene expression
- **Gene Set Enrichment Analysis (GSEA and ORA)** via GSEApy, supporting GO, KEGG, Reactome, and MSigDB Hallmarks
- **Interactive Plots** (Plotly HTML artifacts) for volcano, MA, and PCA visualizations
- **Reproducibility Logging** with full parameter registry and `pip freeze` snapshot
- Support for **batch effect correction** via metadata covariates

### Quick Start

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/BRSA
pip install -r requirements.txt

# CLI usage
python BRSA.py --counts data/counts_matrix.csv --meta data/metadata.csv \
  --cond1 Tumor --cond2 Normal --norm DESeq2 --dim_red UMAP --ml --ppi

# GUI usage
python interface.py
```

### External Dependencies

Requires the following tools available in `PATH` for the upstream preprocessing steps: `fastqc`, `multiqc`, `fastp`, `hisat2`, `samtools`, `featureCounts`, and `kallisto`.

---

## DELIMITER — Delimited File Filter

**Location:** `Delimiter/`

DELIMITER is a desktop application for the inspection, filtering, transformation, and export of delimited data files. Built with PyQt5 and Pandas, it provides a multi-threaded data processing interface designed primarily for bioinformatics tabular data, but applicable to any structured file.

### Key Features

- **File Loading** for CSV, TSV, and Excel (`.xlsx`) with automatic delimiter detection; file I/O runs on dedicated threads to prevent UI freezing
- **Multi-condition Filtering** with stacked filter pods supporting operators such as equals, contains, starts/ends with, regex, and numeric comparisons
- **Filter Preset Management** — save and reload filter configurations as JSON
- **Column Operations** — create derived columns via concatenation, addition, subtraction, multiplication, or division
- **NaN Handling** — drop, fill with mean/median, or fill with a custom value
- **Cell Editing** with full undo/redo history
- **Conditional Cell Highlighting** — define visual rules to flag values in the table
- **Export** to CSV, TSV, Excel, and FASTA formats, all written in chunked threads with progress reporting
- **Descriptive Statistics** panel for any selected column
- **Internationalization** with built-in support for English, Brazilian Portuguese, and Spanish
- **Dark and Light Themes** switchable at runtime

### Quick Start

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Delimiter
pip install PyQt5 pandas openpyxl psutil
python delimeter.py
```

---

## LenghtCount — Genomic & Proteomic Analysis Suite

**Location:** `LengthCounter/`

LenghtCount is a multi-modal bioinformatics platform for comprehensive genomic and proteomic characterization. It consolidates a wide range of sequence analysis tasks into a single tool accessible via CLI or a PyQt6 graphical interface.

### Key Features

**Assembly & Sequence Statistics**
- N50, N90, L50, L90, NG50 and full descriptive statistics for FASTA assemblies
- Histogram, boxplot, violin plot, ECDF, and ridgeline plots of contig length distributions

**Nucleotide & Composition Analysis**
- GC content, GC skew (oriC prediction), CpG island detection, k-mer frequency profiling
- Shannon entropy, DUST low-complexity masking, dinucleotide Rho odds ratio
- Tandem repeat and microsatellite detection, pathogenicity island signatures, telomere mapping

**FASTQ Quality Control**
- Per-cycle Phred quality profiles, base composition, GC distribution
- Sequencing depth estimation (Lander-Waterman), rarefaction/saturation curves, PCR duplicate estimation
- Quality-based FASTQ trimming

**Protein Analysis**
- Physicochemical properties: molecular weight, isoelectric point, instability index, GRAVY score, aliphatic index, molar extinction coefficient
- Kyte-Doolittle hydropathy plots
- Secondary structure prediction (Chou-Fasman)

**Comparative & Structural Genomics**
- Global pairwise alignment (Needleman-Wunsch), Ti/Tv ratio, synonymous/non-synonymous substitution analysis
- Dot plot / synteny matrix
- Restriction enzyme mapping and virtual agarose gel simulation
- CRISPR/SpCas9 NGG PAM site scanning
- ORF prediction (6 reading frames)
- UPGMA phylogenetic tree construction from multiple sequence alignments
- Local BLAST+ integration

**Variant & Annotation Analysis**
- VCF file parsing and interactive variant table viewer
- GFF3-guided gene sequence extraction

**Format Conversion**
- Interconversion between FASTA, FASTQ, GenBank, SAM, and BAM (via samtools)

**NCBI Integration**
- Sequence download by accession via NCBI Entrez

### Quick Start

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/LengthCounter
pip install biopython numpy pandas matplotlib seaborn scipy PyQt6

# Assembly statistics
python Lenght.py assembly.fasta --stats output.txt --histogram

# Launch GUI
python Lenght.py --gui
```

---

## PALINDROME — DNA Palindrome Analyzer

**Location:** `Palindrome/`

PALINDROME is a Python-based command-line tool for the detection and characterization of palindromic sequences in DNA. Palindromic regions — where double-stranded DNA reads identically in both 5'→3' directions — are biologically significant as restriction enzyme recognition sites, transcription factor binding sites, and structural regulators of gene expression.

### Key Features

- **Palindrome Detection** with configurable length range and mismatch tolerance
- **Thermodynamic Profiling** — melting temperature (Tm) and Gibbs free energy (ΔG) via the nearest-neighbor model
- **GC Content** and **Shannon Entropy** per palindrome
- **Restriction Enzyme Matching** against 15 canonical recognition sites (EcoRI, BamHI, HindIII, NotI, and others)
- **Methylation Motif Detection** — Dam (GATC), Dcm (CCWGG), and CpG sites
- **CRISPR/Cas9 PAM Site Mapping** (NGG on both strands)
- **CpG Island Identification** via sliding-window GC%/O/E ratio analysis
- **ORF Prediction** in 6 reading frames
- **Tandem Repeat Detection** and **Low-Complexity Masking**
- **PCR Primer Design** with GC%, Tm, and hairpin checks flanking each palindrome
- **Codon Translation** in three reading frames
- **Shine-Dalgarno and Transcription Terminator Prediction**
- **Ti/Tv Ratio** and **UPGMA Phylogenetic Tree** for comparative analysis
- **Multiple Input Sources** — raw sequence string, FASTA/FASTQ/GenBank files, NCBI accession ID, mock sequence generation, or batch directory processing
- **Output Formats** — CSV table, SVG linear map, SVG circular map, and interactive HTML report

### Quick Start

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Palindrome
pip install matplotlib numpy biopython

# Basic scan
python Palindrome.py -s "ATGCGCAT"

# Full analysis from a FASTA file
python Palindrome.py -f genome.fasta -min 4 -max 12 -mis 1 --crispr --report

# Generate a mock sequence and test
python Palindrome.py --mock 2000 --crispr --report
```

---

## Installation

Each tool has its own dependencies. Refer to the respective section above for specific installation instructions. All tools share the following base requirement:

```
Python >= 3.8
```

Clone the full repository to access all tools:

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
```

---

## License

All tools in this repository are distributed under the **MIT License**.  
See the full license text at: [LICENSE](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)

---

## Author

**Victor S. Caricatte De Araújo**  
Universidade Federal de Minas Gerais (UFMG)  
GitHub: [@VictorCaricatte](https://github.com/VictorCaricatte)  
Repository: [BasicBioinfo](https://github.com/VictorCaricatte/BasicBioinfo)
