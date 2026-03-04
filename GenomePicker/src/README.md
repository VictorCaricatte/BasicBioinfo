#  GenomePicker

**GenomePicker** is a command-line tool written in Python for automated downloading and processing of genomic data from NCBI and protein structures from UniProt/AlphaFold. It supports parallelized downloads, metadata extraction, and geographic visualization of isolates.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
- [Usage](#usage)
- [Arguments](#arguments)
- [Features](#features)
- [Output Files](#output-files)
- [Examples](#examples)
- [Notes and Limitations](#notes-and-limitations)

---

## Overview

GenomePicker automates three main workflows:

1. **Genome/Proteome Download from NCBI** — Downloads genomic files (`.fna`, `.faa`, `.gbff`, `.gff`) directly from the NCBI FTP server, using a metadata table (TSV/CSV) exported from [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/).
2. **Metadata Extraction and Geographic Mapping** — Scrapes BioSample pages on NCBI to extract epidemiological metadata (host, disease, isolation source, geographic location) and generates geographic distribution maps.
3. **Protein Structure Download** — Downloads protein sequences (FASTA) and 3D structures (PDB/CIF) from UniProt/AlphaFold or NCBI Protein databases.

---

## Requirements

| Dependency | Purpose |
|---|---|
| `Python >= 3.8` | Runtime |
| `pandas` | Table parsing |
| `matplotlib` | Static map generation |
| `mpl_toolkits.basemap` | World map rendering |
| `plotly` | Interactive map generation |

Install all dependencies with:

```bash
pip install pandas matplotlib plotly basemap
```

> **Note:** `basemap` may require additional system-level dependencies. See [Basemap installation guide](https://matplotlib.org/basemap/stable/users/installation.html) for platform-specific instructions.

---

## Installation

No installation is required. Simply clone or download the script and run it directly:

```bash
git clone <repository-url>
cd genomepicker
python genome_picker.py --help
```

---

## Input Files

### NCBI Genome Table (TSV/CSV)

The primary input is a metadata table downloaded from the [NCBI Genome portal](https://www.ncbi.nlm.nih.gov/datasets/genome/). This file must contain the following columns:

| Column | Description |
|---|---|
| `Assembly Accession` | Unique genome accession (used as index) |
| `Assembly Name` | Name of the assembly |
| `Organism Name` | Full organism name (genus + species) |
| `Organism Infraspecific Names Strain` | Strain identifier (primary) |
| `Organism Infraspecific Names Isolate` | Isolate identifier (fallback) |
| `Organism Infraspecific Names Breed` | Breed identifier (fallback) |
| `Organism Infraspecific Names Cultivar` | Cultivar identifier (fallback) |
| `Organism Infraspecific Names Ecotype` | Ecotype identifier (fallback) |
| `Assembly BioSample Accession` | BioSample ID (required for `--metadata`) |

> Supported formats: `.tsv` (tab-separated) and `.csv` (comma-separated).

### Protein ID Lists (TXT)

For protein/structure downloads, provide a plain text file with one identifier per line:

```
P12345
Q9UNQ0
A0A000
```

---

## Usage

```bash
python genome_picker.py [OPTIONS]
```

If no arguments are provided, the help menu is displayed automatically.

---

## Arguments

### Genome Download (`--input` required)

| Argument | Description |
|---|---|
| `--input FILE` | Path to the TSV or CSV table from NCBI Genome |
| `--genome` | Download genome sequences in FASTA nucleotide format (`.fna`) |
| `--fasta` | Same as `--genome`, but saves files with `.fasta` extension |
| `--proteome` | Download protein sequences (`.faa`) |
| `--genbank` | Download GenBank flat files (`.gbff`) |
| `--gbk` | Same as `--genbank`, but saves files with `.gbk` extension |
| `--gff` | Download genome annotation files (`.gff`) |
| `--metadata` | Extract BioSample metadata and generate geographic maps |
| `--cpu N` | Number of parallel threads for downloading (default: `1`) |
| `--outdir DIR` | Custom output directory (auto-generated if not specified) |

### Protein/Structure Download

| Argument | Description |
|---|---|
| `--uniprot_fasta FILE` | Text file with UniProt IDs to download FASTA sequences |
| `--uniprot_pdb FILE` | Text file with UniProt IDs to download PDB structures via AlphaFold |
| `--uniprot_cif FILE` | Text file with UniProt IDs to download CIF structures via AlphaFold |
| `--ncbi_protein FILE` | Text file with NCBI Protein IDs to download FASTA sequences |

---

## Features

###  Parallel Downloads
Downloads are parallelized using Python's `multiprocessing.Pool`, significantly reducing runtime for large datasets. The number of workers is capped at `total_cpus - 1` to avoid system overload.

###  Retry Logic
Each download attempt is retried up to **5 times** with a 3-second interval between attempts, ensuring robustness against transient network failures.

###  Duplicate Strain Handling
If two entries share the same strain identifier, a random 5-character alphanumeric suffix (`_dup_XXXXX`) is appended to avoid filename collisions.

###  PROKKA Script Generation
When downloading genome files (`--genome` or `--fasta`), a `PROKKA.sh` shell script is automatically generated in the output directory. This script contains ready-to-use [Prokka](https://github.com/tseemann/prokka) annotation commands for each downloaded genome, pre-filled with genus, species, strain, and locus tag information.

###  Geographic Mapping
When `--metadata` is used, GenomePicker:
1. Scrapes each BioSample page on NCBI to extract host, disease, isolation source, and geographic location.
2. Downloads a reference coordinate table (`latlon.csv`) to resolve country names to latitude/longitude.
3. Generates two maps:
   - **Static PNG** — High-resolution world map (600 DPI) rendered with Basemap (`matplotlib`).
   - **Interactive HTML** — Zoomable, hoverable map rendered with Plotly Express.

###  Protein Structure Sources
- **UniProt FASTA**: fetched from `rest.uniprot.org`
- **AlphaFold PDB/CIF**: fetched from `alphafold.ebi.ac.uk` (model v4)
- **NCBI Protein FASTA**: fetched via NCBI E-utilities (`efetch`)

---

## Output Files

### Genome Downloads

Files are saved in the output directory with standardized names:

```
{Genus_initial}{species}_{strain}.{extension}
```

Example: `Ecoli_O157-H7.fna`

### Metadata

| File | Description |
|---|---|
| `Metadados_DD-MM-YYYY_HH-MM-SS.tsv` | Tab-separated metadata table with host, disease, isolation source, and geographic location per genome |
| `Metadados_contagem_paises_*.tsv` | Country-level count of isolates used for mapping |
| `Mapa_Basemap_*.png` | High-resolution static world map (600 DPI) |
| `Mapa_Interativo_Plotly_*.html` | Interactive geographic distribution map |

### Prokka

| File | Description |
|---|---|
| `PROKKA.sh` | Shell script with pre-configured Prokka commands for each downloaded genome |

### Proteins and Structures

Files are saved as `{ID}.{format}` (e.g., `P12345.fasta`, `Q9UNQ0.pdb`).

---

## Examples

### Download genome sequences using 8 threads
```bash
python genome_picker.py --input assembly_table.tsv --genome --cpu 8 --outdir MyGenomes/
```

### Download genomes and GFF annotation files simultaneously
```bash
python genome_picker.py --input assembly_table.tsv --genome --gff --cpu 4
```

### Extract metadata and generate geographic maps
```bash
python genome_picker.py --input assembly_table.tsv --metadata --cpu 4 --outdir MyMetadata/
```

### Download proteomes and GenBank files
```bash
python genome_picker.py --input assembly_table.tsv --proteome --genbank --cpu 6
```

### Download protein FASTA from UniProt
```bash
python genome_picker.py --uniprot_fasta uniprot_ids.txt --outdir UniProtFASTA/
```

### Download AlphaFold PDB structures
```bash
python genome_picker.py --uniprot_pdb uniprot_ids.txt --outdir Structures/
```

### Download protein FASTA from NCBI Protein
```bash
python genome_picker.py --ncbi_protein ncbi_protein_ids.txt --outdir NCBIProteins/
```

### Full combined run
```bash
python genome_picker.py \
  --input assembly_table.tsv \
  --genome --gff --metadata \
  --cpu 8 \
  --outdir FullAnalysis/
```

---

## Notes and Limitations

- **Entries without a strain/isolate/breed/cultivar/ecotype identifier are automatically skipped.** Only records with at least one infraspecific name are processed.
- The `--metadata` flag scrapes HTML pages from NCBI BioSample and may be slow for large datasets. Using `--cpu` with multiple threads is strongly recommended.
- The geographic map feature depends on the availability of the external coordinate reference file (`latlon.csv`) hosted on GitHub. Countries not present in this reference file will not appear on the map.
- File downloads are skipped if the output file already exists, enabling safe resumption of interrupted runs.
- CPU count is automatically capped at `total_CPUs - 1` regardless of the value passed to `--cpu`, to preserve system stability.

---
