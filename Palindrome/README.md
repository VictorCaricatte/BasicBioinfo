# ⚔️ PALINDROME 

> *"In the double helix, every strand mirrors its twin — just as every hero is defined by their shadow."*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen)]()
[![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-purple)]()

---

## 📜 The Legend — What Is This Artifact?

**Palindrome** is a Python-based command-line bioinformatics tool designed for the deep analysis of DNA palindromic sequences. Palindromic sequences — regions where the double-stranded DNA reads the same in both 5'→3' directions — are of fundamental biological importance: they serve as recognition sites for restriction endonucleases, transcription factor binding sites, and structural regulators of gene expression.

This tool provides a unified pipeline for sequence input, purification, palindrome detection, thermodynamic characterization, CRISPR site mapping, CpG island identification, repeat masking, ORF prediction, and rich report generation — all from a single command.

---

## 🏛️ The Five Sacred Tomes — Project Architecture

| Module | Codename | Purpose |
|---|---|---|
| `Palindrome.py` | *The Herald's Scroll* | CLI entry point and argument dispatcher |
| `Logic.py` | *The Grand Alchemist* | Core bioinformatics algorithms and analyses |
| `DB.py` | *The Guardian Dragon* | Output management: CSV, SVG, HTML reports |
| `Config.py` | *The Scroll of Eternal Truths* | Persistent JSON configuration management |
| `Interface.py` | *The Oracle's Chamber* | Graphical user interface (GUI) |
| `Dependences.py` | *The Arsenal Summoner* | External tool management (BLAST, DIAMOND) |

---

## ⚗️ Magical Arsenal — Features

### Core Analyses
- **Palindrome Detection** with configurable length range and mismatch tolerance
- **GC Content** calculation across full sequences and individual palindromes
- **Nucleotide Frequency** tallying (A, T, C, G, N)
- **Shannon Entropy** for sequence complexity assessment
- **Thermodynamic Profiling**: Melting Temperature (Tm) and Gibbs Free Energy (ΔG) via Nearest-Neighbor model

### Molecular & Structural
- **Restriction Enzyme Recognition**: Matches palindromes against 15 canonical restriction sites (EcoRI, BamHI, HindIII, NotI, and more)
- **Codon Translation**: Three-frame in-silico protein translation
- **Methylation Motif Detection**: Dam (GATC), Dcm (CCWGG), and CpG site identification
- **ORF Detection**: Six-frame Open Reading Frame prediction (ATG → stop codon) with configurable minimum size
- **CpG Island Mapping**: Sliding-window analysis based on GC% ≥ 50% and O/E ratio ≥ 0.6
- **Tandem Repeat Detection**: Identifies sequential motif echoes (e.g., ATATAT...)
- **Shine-Dalgarno Identification**: Ribosomal binding site prediction upstream of ATG portals
- **Transcription Terminator Prediction**: Intrinsic GC-rich hairpin + poly-T tail detection

### Advanced Utilities
- **CRISPR/Cas9 Target Mapping**: NGG PAM site identification on both strands
- **Regex Pattern Search**: Custom motif hunting within sequences
- **Synteny Block Detection**: Identifies conserved k-mer blocks between two sequences
- **Ti/Tv Ratio Calculation**: Transition/transversion ratio for evolutionary analysis
- **UPGMA Phylogenetic Tree**: Simple Newick tree generation from multiple sequences
- **PCR Primer Design**: Heuristic primer design flanking each palindrome (GC%, Tm, hairpin check)
- **Pathogen Signature Scanning**: Matches sequences against a virulence motif database
- **Molecular Weight Calculation**: Exact double-stranded DNA weight in Daltons
- **Random Occurrence Probability**: Statistical expectation of a palindrome occurring by chance

### Pre-Processing
- **PDF Demon Exorcism**: Strips formatting artifacts from copy-pasted PDF/text sequences
- **Low-Complexity Masking** (`--mask`): Transmutes tandem repeat regions to 'N'
- **Terminal Trimming** (`--trim`): Strips dangling terminal 'N' bases and adapters

### Output Formats
- **CSV Report**: Full tabular output of all palindromes and attributes
- **SVG Linear Map**: Vectorial visualization of palindrome positions along the sequence
- **SVG Circular Map (Ouroboros)**: Polar coordinate genome map highlighting palindromic loci
- **HTML Interactive Report**: Dark-themed full genomic analysis report

### Input Sources
- Raw sequence string (`-s`)
- FASTA / multi-FASTA / FASTQ / GenBank files (`-f`)
- NCBI E-utilities remote fetch by accession ID (`--ncbi`)
- Random mock sequence generation (`--mock`)
- Batch directory processing of all FASTA files (`--batch`)

---

## 🗡️ Installation — Forging the Weapon

### Prerequisites

```bash
Python >= 3.8
```

### Clone the Repository

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Palindrome
```

### Install Dependencies

```bash
pip install matplotlib numpy biopython
```

> **Note:** `matplotlib` and `numpy` are required for SVG/circular map generation. `biopython` is required only for GenBank (`.gbk`) file parsing. All other features run with the Python standard library.

---

## 📖 Usage — Casting the Spells

### Basic Palindrome Scan

```bash
python Palindrome.py -s "ATGCGCAT"
```

### From a FASTA File

```bash
python Palindrome.py -f genome.fasta
```

### Multi-FASTA File

```bash
python Palindrome.py -f multi_record.fasta
```

### From NCBI by Accession

```bash
python Palindrome.py --ncbi NC_000913.3
```

### Generate a Random Mock Sequence

```bash
python Palindrome.py --mock 5000
```

### Full-Featured Run

```bash
python Palindrome.py -f genome.fasta -min 4 -max 12 -mis 1 --crispr --mask --trim --report
```

### Batch Directory Processing

```bash
python Palindrome.py --batch ./my_fasta_folder/
```

### Clean All Output Artifacts

```bash
python Palindrome.py --clean
```

---

## 🧙 CLI Argument Reference

| Argument | Type | Default | Description |
|---|---|---|---|
| `-s`, `--sequence` | `str` | `None` | Raw DNA string input |
| `-f`, `--file` | `str` | `None` | Path to genomic file (FASTA, FASTQ, GenBank) |
| `-n`, `--ncbi` | `str` | `None` | NCBI accession ID for remote fetching |
| `--mock` | `int` | `None` | Length of random DNA sequence for testing |
| `--batch` | `str` | `None` | Directory path for batch FASTA processing |
| `-min`, `--minimum` | `int` | `4` | Minimum palindrome length |
| `-max`, `--maximum` | `int` | `12` | Maximum palindrome length |
| `-mis`, `--mismatches` | `int` | `0` | Maximum allowed mismatches in palindrome |
| `--regex` | `str` | `None` | Regex pattern to search within the sequence |
| `--crispr` | flag | `False` | Enable CRISPR/Cas9 NGG PAM site detection |
| `--mask` | flag | `False` | Mask low-complexity tandem repeat regions with 'N' |
| `--trim` | flag | `False` | Trim terminal 'N' bases from sequence edges |
| `--report` | flag | `False` | Generate full interactive HTML report |
| `--clean` | flag | `False` | Remove all SVG, CSV, and HTML output artifacts |

---

## 📊 Output — The Sacred Columns

Each detected palindrome is reported with the following attributes:

| Column | Description |
|---|---|
| `COORD.` | Genomic coordinates (start–end) |
| `SEQ` | Palindromic sequence |
| `MIS` | Number of mismatches |
| `GC%` | GC content percentage |
| `TM(°C)` | Melting temperature |
| `ΔG (kcal)` | Gibbs free energy of duplex |
| `ENTROPY` | Shannon entropy (sequence complexity) |
| `METHYL.` | Detected methylation motifs |
| `PATHOGEN SIG.` | Virulence-associated motif matches |

---

## 🔬 Scientific Methods & References

The algorithms implemented in this tool are grounded in established bioinformatics and molecular biology methods:

### Palindrome Detection
Reverse complement comparison adapted from foundational sequence analysis principles.
> Watson, J.D., Crick, F.H.C. (1953). Molecular structure of nucleic acids. *Nature*, 171, 737–738. https://doi.org/10.1038/171737a0

### GC Content Calculation
> Marmur, J., Doty, P. (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. *Journal of Molecular Biology*, 5(1), 109–118. https://doi.org/10.1016/S0022-2836(62)80066-7

### Melting Temperature — Wallace Rule (short sequences)
> Wallace, R.B., Shaffer, J., Murphy, R.F., et al. (1979). Hybridization of synthetic oligodeoxyribonucleotides to φX174 DNA: the effect of single base pair mismatch. *Nucleic Acids Research*, 6(11), 3543–3557. https://doi.org/10.1093/nar/6.11.3543

### Melting Temperature — Long Sequence Formula
> Baldino, F. Jr., Chesselet, M.F., Lewis, M.E. (1989). High-resolution in situ hybridization histochemistry. *Methods in Enzymology*, 168, 761–777. https://doi.org/10.1016/0076-6879(89)68055-3

### Thermodynamics — Nearest-Neighbor ΔG Model
> SantaLucia, J. Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the National Academy of Sciences*, 95(4), 1460–1465. https://doi.org/10.1073/pnas.95.4.1460

### Shannon Entropy for Sequence Complexity
> Shannon, C.E. (1948). A mathematical theory of communication. *The Bell System Technical Journal*, 27(3), 379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x

### CpG Island Identification
Sliding window algorithm using GC% ≥ 50% and observed/expected CpG ratio ≥ 0.6.
> Gardiner-Garden, M., Frommer, M. (1987). CpG islands in vertebrate genomes. *Journal of Molecular Biology*, 196(2), 261–282. https://doi.org/10.1016/0022-2836(87)90689-9

### Restriction Enzyme Recognition Sites
> Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2015). REBASE — a database for DNA restriction and modification: enzymes, genes and genomes. *Nucleic Acids Research*, 43(D1), D298–D299. https://doi.org/10.1093/nar/gku1046

### Codon Translation Table
> Nomenclature Committee of the International Union of Biochemistry (NC-IUB). (1985). Nomenclature for incompletely specified bases in nucleic acid sequences. *European Journal of Biochemistry*, 150(1), 1–5. https://doi.org/10.1111/j.1432-1033.1985.tb08977.x

### DNA Methylation Motifs (Dam, Dcm, CpG)
> Casadesús, J., Low, D.A. (2013). Programmed heterogeneity: epigenetic mechanisms in bacteria. *Journal of Biological Chemistry*, 288(20), 13929–13935. https://doi.org/10.1074/jbc.R113.472274

### CRISPR/Cas9 PAM Site Detection (NGG)
> Jinek, M., Chylinski, K., Fonfara, I., et al. (2012). A programmable dual-RNA–guided DNA endonuclease in adaptive bacterial immunity. *Science*, 337(6096), 816–821. https://doi.org/10.1126/science.1225829

### Open Reading Frame (ORF) Prediction
> Zhu, H., Hu, G.Q., Yang, Y.F., Wang, J., She, Z.S. (2007). MED: a new non-supervised gene prediction algorithm for bacterial and archaeal genomes. *BMC Bioinformatics*, 8, 97. https://doi.org/10.1186/1471-2105-8-97

### Shine-Dalgarno Sequence Recognition
> Shine, J., Dalgarno, L. (1974). The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. *Proceedings of the National Academy of Sciences*, 71(4), 1342–1346. https://doi.org/10.1073/pnas.71.4.1342

### Intrinsic Transcription Terminator Prediction
> d'Aubenton Carafa, Y., Brody, E., Thermes, C. (1990). Prediction of rho-independent Escherichia coli transcription terminators. *Journal of Molecular Biology*, 216(4), 835–858. https://doi.org/10.1016/S0022-2836(99)80005-9

### Transition/Transversion Ratio (Ti/Tv)
> Freese, E. (1962). On the evolution of the base composition of DNA. *Journal of Theoretical Biology*, 3(1), 82–101. https://doi.org/10.1016/0022-5193(62)90071-8

### UPGMA Phylogenetic Tree Construction
> Sokal, R.R., Michener, C.D. (1958). A statistical method for evaluating systematic relationships. *University of Kansas Science Bulletin*, 38, 1409–1438.

### Tandem Repeat Detection
> Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research*, 27(2), 573–580. https://doi.org/10.1093/nar/27.2.573

### Molecular Weight Calculation
> Thermo Fisher Scientific. (2020). Calculating the molecular weight of double-stranded DNA. *Technical Reference*. Retrieved from https://www.thermofisher.com

### Smith-Waterman Local Alignment
> Smith, T.F., Waterman, M.S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195–197. https://doi.org/10.1016/0022-2836(81)90087-5

### NCBI E-Utilities Sequence Fetching
> Sayers, E.W., Beck, J., Bolton, E.E., et al. (2022). Database resources of the National Center for Biotechnology Information. *Nucleic Acids Research*, 50(D1), D20–D26. https://doi.org/10.1093/nar/gkab1112

---

## 📂 Output File Structure

After execution, the following files are generated in the `local_db_dungeon/` directory:

```
local_db_dungeon/
├── cli_alchemical_treasures.csv    # Full palindrome table (CSV)
├── constellation_map.svg           # Linear DNA map (SVG)
├── ouroboros_mirror.svg            # Circular genome map (SVG)
└── cli_ultimate_parchment.html     # Interactive HTML report (with --report)

batch_results_dungeon/
└── hellfire_batch_summary.csv      # Batch processing summary (with --batch)
```

---

## 🔧 External Tools Integration

Palindrome supports optional integration with:

- **BLAST** (`makeblastdb`): For building local nucleotide databases
- **DIAMOND** (`diamond makedb`): For building protein/DNA comparison libraries

Configure paths via the GUI (`Interface.py`) or by editing `secret_configurations.json`.

---

## 🧪 Quick Test Run

```bash
# Generate a 2000 bp mock sequence and run full analysis with CRISPR mapping and HTML report
python Palindrome.py --mock 2000 -min 4 -max 12 --crispr --report
```

---

## 📜 License

This project is distributed under the **MIT License**.
See the full license text at: [LICENSE](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)

---

## 🧙‍♂️ Author

**Victor Caricatte**
GitHub: [@VictorCaricatte](https://github.com/VictorCaricatte)
Repository: [BasicBioinfo/Palindrome](https://github.com/VictorCaricatte/BasicBioinfo/tree/main/Palindrome)

---

> *"The sequence speaks. The palindrome whispers back. The alchemist listens."*
