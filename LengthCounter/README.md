<div align="center">

# 🍃 LenghtCount  🍃

> *"The Will of Fire applied to Bioinformatics."*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python)](https://www.python.org/)
[![BioPython](https://img.shields.io/badge/BioPython-1.80%2B-green?style=for-the-badge)](https://biopython.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20Windows%20%7C%20macOS-lightgrey?style=for-the-badge)]()
[![Repository](https://img.shields.io/badge/GitHub-LengthCounter-orange?style=for-the-badge&logo=github)](https://github.com/VictorCaricatte/BasicBioinfo/tree/main/LengthCounter)

**A full-spectrum genomic and proteomic analysis suite with both a GUI (PyQt6) and a powerful CLI — forged in the Hidden Leaf Village of Bioinformatics.**

</div>

---

## 📜 Table of Contents

1. [Overview](#-overview)
2. [Architecture](#-architecture)
3. [Requirements & Installation](#-requirements--installation)
4. [Quick Start](#-quick-start)
5. [Module I — Assembly Statistics (Sharingan)](#-module-i--assembly-statistics-sharingan)
6. [Module II — Base Composition (Mokuton)](#-module-ii--base-composition-mokuton)
7. [Module III — Protein & Primers (Chidori)](#-module-iii--protein--primers-chidori)
8. [Module IV — FastQ Metrics (Rasengan)](#-module-iv--fastq-metrics-rasengan)
9. [Module V — Advanced Sequence Analysis (Rinnegan)](#-module-v--advanced-sequence-analysis-rinnegan)
10. [Module VI — Variant Viewer (Neji Eye)](#-module-vi--variant-viewer-neji-eye)
11. [Module VII — Anbu Floating HQ (Embedded Shell)](#-module-vii--anbu-floating-hq-embedded-shell)
12. [CLI Reference](#-cli-reference)
13. [Scientific Methods & References](#-scientific-methods--references)
14. [File Format Support](#-file-format-support)
15. [Contributing](#-contributing)

---

## 🌀 Overview

**LenghtCount** is an open-source, multi-modal bioinformatics platform designed for comprehensive genomic and proteomic characterization. It integrates statistical assembly analysis, nucleotide composition profiling, protein physicochemical characterization, sequencing quality control, phylogenetics, variant calling, and CRISPR target prediction into a single cohesive tool.

The project is structured around three files:

| File | Role |
|---|---|
| `Lenght.py` | Entry point — CLI argument parsing and execution orchestration |
| `logic.py` | Core computational engine — all biological algorithms |
| `Interface.py` | PyQt6 graphical interface layer |
| `config.py` | Styling constants and documentation HTML |

The tool is entirely **self-contained**: it requires no external web services beyond NCBI's Entrez API (for the accession downloader) and an optional local BLAST+ installation.

---

## 🏗 Architecture

```
LenghtCount/
│
├── Lenght.py             # CLI entry point (argparse)
├── logic.py              # Core bioinformatics engine
├── Interface.py          # PyQt6 GUI
├── config.py             # Theming, colors, doc strings
│
└── Analyzer Classes (logic.py):
    ├── UchihaSharinganAnalyzer     → Assembly Stats (FASTA/GBK)
    ├── SenjuMokutonAnalyzer        → Base Composition & Motifs
    ├── UzumakiChakraFastqAnalyzer  → FASTQ QC Metrics
    ├── RinneganAdvancedBioTools    → Alignment, ORF, Phylo, CRISPR
    ├── KatsuyuSlugFormatParser     → VCF / GFF / SAM file parsing
    ├── NejiVariantEyeParser        → VCF Variant Extraction
    ├── YamatoWoodFeatureExtractor  → GFF-guided gene extraction
    ├── PainBanshoTeninNCBIFetcher  → NCBI Entrez sequence retrieval
    ├── GaaraSandTrimmer            → FASTQ quality trimming
    ├── HengeFormatShifter          → Format conversion (BAM/SAM/FASTQ/GBK)
    └── AnbuBlackOpsExternalTools   → BLAST+ subprocess wrapper
```

---

## ⚙️ Requirements & Installation

### System Dependencies

- Python **3.8** or higher
- `samtools` (optional, required only for BAM → SAM conversion)
- NCBI BLAST+ (optional, required only for local BLAST alignment)

### Python Dependencies

```bash
pip install biopython numpy pandas matplotlib seaborn scipy PyQt6
```

Full requirements list:

| Package | Purpose |
|---|---|
| `biopython >= 1.80` | Core sequence I/O, alignment, restriction analysis, phylogenetics, Entrez |
| `numpy` | Vectorized array operations (GC skew, CpG detection, dot plots) |
| `pandas` | Tabular data handling for VCF, length distributions, coverage |
| `matplotlib` | Figure rendering (all visualizations) |
| `seaborn` | Statistical plotting (boxplot, violin, KDE, ECDF, histplot) |
| `scipy` | Statistical functions |
| `PyQt6` | GUI framework |

### Installation

```bash
# Clone the repository
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/LengthCounter

# Install Python dependencies
pip install -r requirements.txt

# Launch GUI
python Lenght.py --gui

# Or use the CLI directly
python Lenght.py genome.fasta --histogram --stats report.txt
```

---

## 🚀 Quick Start

```bash
# Basic assembly statistics
python Lenght.py assembly.fasta --stats output.txt

# Generate histogram and boxplot
python Lenght.py assembly.fasta --histogram --boxplot --out-prefix my_analysis

# Download from NCBI and analyze
python Lenght.py --ncbi NC_045512

# FASTQ quality profiling
python Lenght.py reads.fastq.gz --phred --rarefaction

# Protein physicochemical analysis
python Lenght.py --protein-file proteins.faa

# Launch the GUI
python Lenght.py --gui
```

---

## 🔴 Module I — Assembly Statistics (Sharingan)

> *"The Sharingan sees everything. So does this module."*

This module provides a comprehensive statistical framework for evaluating genomic assembly quality through multiple distribution visualizations and summary metrics.

### Assembly Quality Metrics

The core statistical engine (`shikamaru_shadow_stats`) computes the following metrics over a set of contig/scaffold lengths:

| Metric | Description |
|---|---|
| **N50** | The length such that contigs ≥ N50 cover at least 50% of the total assembly length |
| **N90** | The length such that contigs ≥ N90 cover at least 90% of the total assembly length |
| **L50 / L90** | Number of contigs required to reach N50/N90, respectively |
| **NG50** | Reference-aware N50 — computed against an expected genome size *G* |
| **LG50** | Reference-aware L50 |
| **Min / Max / Mean / Median / Std** | Descriptive distribution statistics |
| **Q1 / Q3 (IQR)** | First and third quartiles |

**Bin width selection** for the histogram uses the Freedman-Diaconis rule:

```
bin_width = 2 × IQR × n^(-1/3)
```

This rule is preferred over Sturges' rule for large, skewed genomic length distributions (Freedman & Diaconis, 1981).

### Visualizations

| Plot | Renderer | Use Case |
|---|---|---|
| **Histogram** | Freedman-Diaconis binning | Overall length distribution per dataset |
| **Boxplot** | IQR-based, outliers suppressed | Inter-dataset comparison |
| **Violin Plot** | KDE mirrored | Full distribution shape |
| **KDE (Kernel Density Estimate)** | Gaussian kernel | Smooth density function per dataset |
| **ECDF (Empirical CDF)** | Step function, `seaborn.ecdfplot` | Cumulative proportion analysis |

### File Operations

- **Multi-file and directory loading**: All common genomic formats (FASTA, FASTQ, GBK, VCF, GFF, BAM) are auto-detected by the format sniffer (`neji_byakugan_format_sniffer`).
- **Concatenation (Chibaku Tensei)**: Merges multiple selected FASTA files into a single output.
- **Virtual Trimming**: Filters sequences by a configurable Min/Max bp window before plotting.
- **CSV export**: Exports raw length-per-contig data for downstream analysis in R or Excel.
- **HTML report**: Embeds all generated figures as base64-encoded PNGs into a standalone HTML file.

### NCBI Entrez Downloader

The `PainBanshoTeninNCBIFetcher` wraps BioPython's `Entrez.efetch` to retrieve sequences in FASTA format directly from NCBI Nucleotide by accession number (e.g., `NC_045512` for SARS-CoV-2).

```bash
python Lenght.py --ncbi NC_045512 --out-prefix sarscov2
```

---

## 🟢 Module II — Base Composition (Mokuton)

> *"From a single nucleotide, a forest grows."*

This module focuses on the intrinsic compositional and thermodynamic properties of nucleotide sequences.

### Nucleotide & GC Content

Raw base frequencies (A, T, C, G, U, N) are counted and reported with percentage composition. GC content is used throughout this module as a primary discriminatory feature.

### Transcription, Translation & Reverse Complement

Powered by `Bio.Seq`:

- **Transcription**: DNA → RNA (T→U substitution)
- **Translation**: DNA/RNA → Protein, using NCBI standard genetic code (`translate(to_stop=False)`)
- **Reverse Complement**: Watson-Crick complementarity + strand reversal

### Regex-based Motif Search

The `Sharingan Regex Eye` allows full Python `re` pattern matching over the loaded sequence. IUPAC ambiguity codes and complex patterns (e.g., `AT[GC]T`, `CG{3,}`) are supported. Matches are highlighted positionally.

### K-mer Frequency Analysis

K-mer counting uses a sliding-window generator with Python's `collections.Counter`, reporting the top 15 most frequent k-mers. This approach is equivalent to the frequency-based genomic signature method described by Karlin et al. (1997), used for species-level taxonomic profiling and horizontal gene transfer detection.

```bash
python Lenght.py genome.fasta --kmer 6
```

### Shannon Entropy & Complexity Masking

**Shannon Entropy** is computed over a sliding window of configurable size:

```
H(w) = - Σ p_i × log2(p_i)
```

where *p_i* is the frequency of base *i* in window *w*. High entropy indicates diverse, information-rich sequence; low entropy signals repetitive or low-complexity regions (Shannon, 1948).

**DUST masking** (`gaara_dust_complexity_shield`): Windows with H < 1.0 bit are replaced with `N`, analogous to the DUST algorithm implemented in BLAST+ (Morgulis et al., 2006).

```bash
python Lenght.py genome.fasta --entropy 50
```

### CpG Island Detection (Vectorized)

CpG islands are identified using the Gardiner-Garden & Frommer (1987) criteria:

- Observed/Expected CpG ratio ≥ 0.60
- GC content ≥ 50% (implicit via sliding-window enrichment)

The computation is **fully vectorized with NumPy** using `np.convolve` on binary base masks, providing O(N) performance on chromosome-scale sequences.

```bash
python Lenght.py genome.fasta --cpg
```

### GC Skew

GC skew is calculated with a rolling window of 1,000 bp (default):

```
GC Skew = (G - C) / (G + C)
```

Sign changes in GC skew are classical indicators of the **origin of replication (oriC)** and the **replication terminus (ter)** in bacterial chromosomes (Grigoriev, 1998).

```bash
python Lenght.py genome.fasta --gc-skew
```

### Restriction Enzyme Mapping & Virtual Gel

The `zabuza_executioner_restriction_map` uses `Bio.Restriction` to identify cleavage sites for EcoRI, BamHI, HindIII, XhoI, and NotI. Fragment sizes are plotted on a **log-scale strip plot** to simulate an agarose gel image.

### Pathogenicity Island (PAI) Detection

The `orochimaru_cursed_islands_pai` detects genomic windows with GC content significantly below the global genomic average (> 5% drop). This heuristic mimics the GC deviation method for identifying horizontally acquired DNA, such as pathogenicity islands (Hacker & Kaper, 2000).

### Telomere Mapping

TTAGGG/CCCTAA hexamer repeats are searched within the first and last 5,000 bp of each contig using a regex requiring ≥ 3 consecutive tandem copies — consistent with canonical vertebrate telomeric repeat unit structure (Meyne et al., 1989).

### Dinucleotide Rho Odds Ratio

The Rho (ρ) statistic measures dinucleotide over- or under-representation:

```
ρ(XY) = f(XY) / (f(X) × f(Y))
```

Values > 1.2 indicate over-representation; values < 0.8 indicate suppression. CpG suppression (ρ(CpG) << 1) is a hallmark of vertebrate genomes (Karlin & Ladunga, 1994).

### Kraken-lite Taxonomic Inference

A heuristic GC-content-based organism classifier matching against a small reference database. Intended for educational demonstration of the GC content profiling principle; not a replacement for KRAKEN2 (Wood & Salzberg, 2014) in production workflows.

---

## 🟡 Module III — Protein & Primers (Chidori)

> *"One thousand birds singing physicochemical truths."*

### ORF Finder — All 6 Frames

`jiraiya_sage_mode_orf_finder` scans both strands in all three reading frames (6 total) for ATG-initiated, stop-codon-terminated ORFs using `Bio.Seq.translate(to_stop=True)`. Only ORFs ≥ 100 amino acids (default) are reported, sorted by length descending.

```bash
python Lenght.py genome.fasta --orf
```

### Protein Physicochemical Properties (Gai Eight Gates)

Powered by `Bio.SeqUtils.ProtParam.ProteinAnalysis`:

| Property | Algorithm / Reference |
|---|---|
| **Molecular Weight** | Sum of residue masses + water (Gasteiger et al., 2005) |
| **Isoelectric Point (pI)** | Iterative pH search, Henderson-Hasselbalch equations |
| **Instability Index** | DIWV matrix method (Guruprasad et al., 1990). Score < 40 = stable |
| **GRAVY Score** | Grand Average of Hydropathicity (Kyte & Doolittle, 1982) |
| **Aliphatic Index** | `100 × (Ala + 2.9Val + 3.9(Ile + Leu)) / N` (Ikai, 1980) |
| **Molar Extinction Coefficient** | Pace et al. (1995); reduced and oxidized (Cys-Cys disulfide) forms |
| **2D Secondary Structure Fraction** | Chou-Fasman propensity statistics: α-helix, β-sheet, turn (Chou & Fasman, 1978) |

```bash
python Lenght.py --protein-file proteins.faa
```

### Kyte-Doolittle Hydropathy Plot

A sliding-window (default: 9 residues) mean hydrophobicity is computed using the Kyte-Doolittle scale (Kyte & Doolittle, 1982). Segments exceeding **+1.6** are annotated as putative transmembrane domains, consistent with the TMbase criterion (Hofmann & Stoffel, 1993).

```bash
python Lenght.py --hydro-file proteins.faa
```

### Oligonucleotide Thermodynamics — Primer Wizard

`minato_rasengan_primer_wizard` computes for any given primer sequence:

| Parameter | Method |
|---|---|
| **Tm (Melting Temperature)** | Nearest-Neighbor thermodynamic model (`Bio.SeqUtils.MeltingTemp.Tm_NN`) — SantaLucia, 1998 |
| **ΔG (Gibbs Free Energy)** | `ΔG = ΔH - T × ΔS` at 37°C (310.15 K) |
| **GC%** | Direct count; optimal range: 40–60% |
| **Self-Dimer Risk** | Reverse complement prefix overlap detection |
| **Hairpin Risk** | Intra-sequence reverse complement match detection (≥ 4 bp stem) |

```bash
python Lenght.py --primer-file primers.fasta
```

### CRISPR SpCas9 PAM Scanner

`kakashi_crispr_copy_sgrna` scans the sequence for all `NGG` PAM sites using a lookahead regex. For each valid site, it extracts:

- The 20 nt protospacer (sgRNA)
- GC% of the guide (optimal range: 40–70%)
- Off-target risk estimate based on exact 20-mer frequency in the provided sequence

The algorithm follows the canonical CRISPR-Cas9 targeting rules of Doench et al. (2014) and Hsu et al. (2013).

```bash
python Lenght.py --crispr genome.fasta
```

---

## 🔵 Module IV — FastQ Metrics (Rasengan)

> *"Quality is everything. Even Naruto learned that."*

### FASTQ Parsing (Block/Generator-based)

`naruto_rasengan_read_fastq_blocks` uses `Bio.SeqIO.QualityIO.FastqGeneralIterator` — a generator-based reader that processes FASTQ files line-by-line without loading the entire file into memory. Supports both plain `.fastq` and gzip-compressed `.fastq.gz`.

### Per-Cycle Phred Quality Profile

Phred quality scores (Q) are decoded from ASCII encoding:

```
Q = ord(char) - 33  [Phred+33, Illumina 1.8+]
```

Average quality per cycle position is computed and plotted, revealing the characteristic **3' quality drop-off** common in Illumina short-read sequencing (Ewing & Green, 1998).

```bash
python Lenght.py reads.fastq.gz --phred
```

### Genome Coverage Estimation

Lander-Waterman coverage equation (Lander & Waterman, 1988):

```
Coverage (C) = (N_reads × L_read × multiplier) / G
```

Where:
- `N_reads` = total read count
- `L_read` = average read length (configurable, default 150 bp)
- `multiplier` = 2 for paired-end, 1 for single-end
- `G` = expected genome size (configurable, default 1 Mb)

```bash
python Lenght.py reads.fastq --coverage --gsize 3000000000 --rlen 150
```

### Rarefaction / Saturation Curve

Simulates progressive sampling of reads and plots cumulative unique k-mer (representative) coverage against read count. A plateau in the curve indicates sequencing saturation — i.e., additional sequencing will yield diminishing novel information (Lander & Waterman, 1988; Good, 1953).

```bash
python Lenght.py reads.fastq --rarefaction
```

### PCR Duplication Estimator

Analyzes the first 50 bp of each read. A high proportion (> 30%) of 100%-identical 50-mer fragments is a strong indicator of **PCR amplification bias** — a common artifact in library preparation (Kozarewa et al., 2009). This heuristic is analogous to the deduplication flag in FastQC.

### FASTQ Trimmer (Gaara Sand Trimmer)

Quality-based 3'-end trimming: reads are written to a new FASTQ only if their mean Phred score exceeds a configurable minimum threshold (default Q20, equivalent to 99% base call accuracy).

```bash
python Lenght.py --trim-fastq reads.fastq --min-phred 20
```

---

## 🟣 Module V — Advanced Sequence Analysis (Rinnegan)

> *"The eyes that see all paths."*

### Global Alignment — Needleman-Wunsch

`sasuke_sharingan_snp_deducer_global` implements pairwise global alignment using BioPython's `Align.PairwiseAligner` in **global mode** with the following scoring scheme:

| Parameter | Value |
|---|---|
| Match score | +1 |
| Mismatch penalty | -2 |
| Gap open penalty | -5 |
| Gap extend penalty | -1 |

This is equivalent to the Needleman-Wunsch algorithm (Needleman & Wunsch, 1970) with affine gap penalties.

**Post-alignment SNP analysis:**

- **Ti/Tv Ratio**: Transitions (A↔G, C↔T) vs Transversions (A/G↔C/T). A ratio of ~2.0 is expected for biological mutation; lower values suggest sequencing error artefacts (Collins & Jukes, 1994).
- **Synonymous/Non-Synonymous classification**: Each SNP is classified by comparing the amino acid encoded by the reference codon vs the mutant codon.

```bash
python Lenght.py --snp-files reference.fasta query.fasta
```

### Dot Plot Matrix

`madara_rinnegan_dotplot_matrix_numpy` generates a **sequence similarity dot matrix** using NumPy vectorized broadcasting for O(N × M) comparison, followed by a diagonal convolution with `scipy.signal.convolve2d` to filter matches by a sliding window threshold. Limited to sequences ≤ 4,000 bp to prevent memory overflow. This approach mirrors the classic dotplot method of Gibbs & McIntyre (1970).

```bash
python Lenght.py --dotplot-files seq1.fasta seq2.fasta
```

### Tandem Repeat Detection

`choji_expansion_jutsu_tandem_repeats` uses regex patterns to detect microsatellites and tandem repeats of unit lengths 2–9 bp with ≥ 3 copies, equivalent to the approach used by TRF (Tandem Repeats Finder; Benson, 1999) in its basic motif detection mode.

```bash
python Lenght.py genome.fasta --tandem
```

### UPGMA Phylogenetic Tree

`hashirama_phylo_tree_builder` builds a phylogenetic dendrogram from a **multiple sequence alignment (FASTA)** using:

1. **Identity distance matrix** (`Bio.Phylo.TreeConstruction.DistanceCalculator('identity')`)
2. **UPGMA (Unweighted Pair Group Method with Arithmetic Mean)** clustering (`DistanceTreeConstructor.upgma()`)

UPGMA assumes a molecular clock and uniform substitution rates (Sokal & Michener, 1958). It is appropriate for closely related sequences or educational demonstrations.

```bash
python Lenght.py --phylo alignment.fasta
```

### Local BLAST+ Integration

`deidara_explosive_blast_art` wraps the NCBI BLAST+ `blastn` binary through a subprocess call, writing temporary FASTA files, executing the alignment, and parsing the resulting XML via `Bio.Blast.NCBIXML`. Requires a local BLAST+ installation (Altschul et al., 1990).

```bash
python Lenght.py --blast-path /usr/bin/blastn --blast-query query.fasta --blast-subject subject.fasta
```

### Gene Extraction from GFF

`YamatoWoodFeatureExtractor.mokuton_extract_feature` parses a GFF3 annotation file to retrieve chromosomal coordinates of a named feature, then surgically extracts the corresponding sequence from the reference FASTA using BioPython's `SeqIO`. Strand orientation is respected (reverse complement applied for `-` strand features).

```bash
python Lenght.py --gff annotation.gff3 --fasta-ref genome.fasta --gene "BRCA1"
```

### Gene Density Plot

Uses `seaborn.histplot` with 10,000 bp bins and KDE overlay to visualize the distribution of gene start positions along the chromosome — a common genome annotation quality metric.

---

## 🔷 Module VI — Variant Viewer (Neji Eye)

> *"The Byakugan sees every SNP."*

`NejiVariantEyeParser.eight_trigrams_parse_vcf` parses standard **VCF (Variant Call Format)** files (v4.x), including gzip-compressed `.vcf.gz`. It extracts the mandatory columns:

```
CHROM | POS | ID | REF | ALT | QUAL | INFO
```

The GUI builds a fully interactive and searchable `QTableWidget` from the resulting records, with filtering by chromosome, position, variant ID, and REF/ALT alleles. The parser correctly skips `##` metadata header lines and `#CHROM` column header lines.

---

## ⚔️ Module VII — Anbu Floating HQ (Embedded Shell)

The integrated terminal (`AnbuBlackOpsExternalTools.itachi_mangekyou_run_command`) provides a dockable embedded shell widget within the GUI. It executes arbitrary shell commands via Python's `subprocess.run` with `shell=True`, capturing both stdout and stderr streams.

**Recommended uses:**
- Launching Docker containers with bioinformatics tools (e.g., `docker run broadinstitute/gatk ...`)
- Running third-party scripts while the main GUI remains active
- Triggering samtools, BWA, STAR, or other external tools in parallel

> ⚠️ **Warning:** The shell executes with the same permissions as the running Python process. Use responsibly.

---

## 📟 CLI Reference

```
usage: Lenght.py [-h] [--gui] [--out-prefix OUT_PREFIX] [--label LABEL]
                 [--histogram] [--boxplot] [--stats STATS] [--csv CSV]
                 [--kmer K] [--entropy WINDOW] [--cpg] [--gc-skew]
                 [--orf] [--tandem]
                 [--coverage] [--gsize GSIZE] [--rlen RLEN]
                 [--phred] [--rarefaction]
                 [--primer-file PRIMER_FILE] [--protein-file PROTEIN_FILE]
                 [--hydro-file HYDRO_FILE]
                 [--snp-files REF QUERY] [--dotplot-files SEQ1 SEQ2]
                 [--ncbi NCBI] [--phylo PHYLO]
                 [--restriction RESTRICTION] [--crispr CRISPR]
                 [--vcf-table VCF_TABLE]
                 [--blast-path BLAST_PATH] [--blast-query BLAST_QUERY]
                 [--blast-subject BLAST_SUBJECT]
                 [--gff GFF] [--fasta-ref FASTA_REF] [--gene GENE]
                 [--trim-fastq TRIM_FASTQ] [--min-phred MIN_PHRED]
                 [--convert-in CONVERT_IN] [--convert-out CONVERT_OUT]
                 [--fmt-in FMT_IN] [--fmt-out FMT_OUT]
                 [files ...]
```

### Key Arguments

| Flag | Description |
|---|---|
| `--gui` | Launch the PyQt6 graphical interface |
| `--histogram` | Generate Freedman-Diaconis length histogram |
| `--boxplot` | Generate sequence length boxplot |
| `--stats FILE` | Write assembly statistics report to FILE |
| `--csv FILE` | Export raw length matrix to CSV |
| `--kmer K` | Compute and plot K-mer frequency (e.g., `--kmer 6`) |
| `--entropy W` | Plot Shannon entropy with window size W bp |
| `--cpg` | Detect and plot CpG islands |
| `--gc-skew` | Plot rolling GC skew (oriC prediction) |
| `--orf` | Find ORFs in all 6 reading frames |
| `--tandem` | Detect microsatellites and tandem repeats |
| `--coverage` | Estimate genome coverage (use with `--gsize`, `--rlen`) |
| `--phred` | Plot per-cycle Phred quality profile |
| `--rarefaction` | Plot sequencing saturation/rarefaction curve |
| `--primer-file` | Analyze primer thermodynamics from FASTA |
| `--protein-file` | Protein physicochemical properties from FASTA |
| `--hydro-file` | Kyte-Doolittle hydropathy plot |
| `--snp-files REF QUERY` | Global alignment + Ti/Tv + synonymy analysis |
| `--dotplot-files S1 S2` | Synteny dot matrix |
| `--ncbi ACCESSION` | Download sequence from NCBI by accession |
| `--phylo MSA.fasta` | Build UPGMA phylogenetic tree from MSA |
| `--restriction SEQ` | Restriction digest + virtual agarose gel |
| `--crispr SEQ` | SpCas9 PAM scan + sgRNA list |
| `--vcf-table VCF` | Parse and display VCF variant records |
| `--blast-query / --blast-subject` | Run local BLAST+ alignment |
| `--gff + --fasta-ref + --gene` | Extract gene sequence from GFF annotation |
| `--trim-fastq + --min-phred` | Quality-trim a FASTQ file |
| `--convert-in/out --fmt-in/out` | Convert between BAM/SAM/FASTQ/GBK/FASTA |

---

## 📚 Scientific Methods & References

| Method | Implementation | Reference |
|---|---|---|
| N50 / NG50 assembly metrics | `shikamaru_shadow_stats` | Schatz et al. (2010); Miller et al. (2010) |
| Freedman-Diaconis binning | `sasuke_amaterasu_histogram` | Freedman & Diaconis (1981) |
| Shannon Entropy | `orochimaru_curse_mark_entropy_plot` | Shannon (1948) |
| DUST complexity masking | `gaara_dust_complexity_shield` | Morgulis et al. (2006) |
| CpG Island detection (O/E ratio) | `hidan_jashin_cpg_islands_vectorized` | Gardiner-Garden & Frommer (1987) |
| GC Skew / oriC prediction | `kakashi_kamui_gc_skew_rolling` | Grigoriev (1998) |
| Dinucleotide Rho odds ratio | `neji_trigram_rho_odds` | Karlin & Ladunga (1994) |
| K-mer genomic signature | `shino_kikaichu_kmer_swarm` | Karlin et al. (1997) |
| Pathogenicity Island detection | `orochimaru_cursed_islands_pai` | Hacker & Kaper (2000) |
| Telomere mapping | `kimimaro_bone_telomere_mapper` | Meyne et al. (1989) |
| Restriction enzyme analysis | `zabuza_executioner_restriction_map` | BioPython: Cock et al. (2009) |
| ORF prediction (6 frames) | `jiraiya_sage_mode_orf_finder` | Standard genetic code; BioPython |
| Protein pI, MW, instability | `raikage_lightning_armor_protein_properties` | Gasteiger et al. (2005); Guruprasad et al. (1990) |
| GRAVY / Aliphatic Index | `raikage_lightning_armor_protein_properties` | Kyte & Doolittle (1982); Ikai (1980) |
| Molar extinction coefficient | `raikage_lightning_armor_protein_properties` | Pace et al. (1995) |
| Chou-Fasman 2D structure | `raikage_lightning_armor_protein_properties` | Chou & Fasman (1978) |
| Kyte-Doolittle hydropathy | `kisame_water_prison_hydrophobicity_plot` | Kyte & Doolittle (1982) |
| Nearest-Neighbor Tm (ΔG, ΔH, ΔS) | `minato_rasengan_primer_wizard` | SantaLucia (1998) |
| CRISPR PAM scanning | `kakashi_crispr_copy_sgrna` | Hsu et al. (2013); Doench et al. (2014) |
| Phred quality scoring | `naruto_rasengan_read_fastq_blocks` | Ewing & Green (1998) |
| Lander-Waterman coverage | `naruto_oodama_rasengan_calculate_genomes` | Lander & Waterman (1988) |
| PCR duplication estimation | FASTQ 50-mer fingerprint analysis | Kozarewa et al. (2009) |
| Needleman-Wunsch alignment | `sasuke_sharingan_snp_deducer_global` | Needleman & Wunsch (1970) |
| Ti/Tv ratio | post-alignment SNP parser | Collins & Jukes (1994) |
| Dot plot / synteny matrix | `madara_rinnegan_dotplot_matrix_numpy` | Gibbs & McIntyre (1970) |
| Tandem repeat detection | `choji_expansion_jutsu_tandem_repeats` | Benson (1999) |
| UPGMA phylogenetic tree | `hashirama_phylo_tree_builder` | Sokal & Michener (1958) |
| Local BLAST+ | `deidara_explosive_blast_art` | Altschul et al. (1990) |

---

## 📂 File Format Support

| Format | Extension(s) | Operations |
|---|---|---|
| FASTA | `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz` | Read, write, merge, extract, slice |
| FASTQ | `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` | Read, QC, trim, convert |
| GenBank | `.gb`, `.gbk`, `.gbff` | Read, parse features, convert to FASTA |
| VCF | `.vcf`, `.vcf.gz` | Parse, count, table view |
| GFF3 | `.gff`, `.gff3` | Parse, feature extraction, density plot |
| BED | `.bed` | Feature count |
| SAM | `.sam` | Alignment count |
| BAM | `.bam` | Detect; convert to SAM via `samtools` |

---

## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/new-jutsu`)
3. Commit your changes (`git commit -m 'Add: new bioinformatics jutsu'`)
4. Push to the branch (`git push origin feature/new-jutsu`)
5. Open a Pull Request

---

<div align="center">

**Developed by VictorSC**  
*"The Will of Fire applied to Bioinformatics."*

🍃 *Believe it.* 🍃

</div>
