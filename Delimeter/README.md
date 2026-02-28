# 🧬 DELIMITER — Delimited File Filter Quantum
### *A Particle Accelerator for Genomic Data. A Reactor. A Beast.*

> **"We do not process files here. We detonate them."**
> — *Supreme Genomic Overlord, Session Initiated*

---

```
╔══════════════════════════════════════════════════════════════════════╗
║   ██████  ███████ ██      ██ ███    ███ ██ ████████ ███████ ██████  ║
║   ██   ██ ██      ██      ██ ████  ████ ██    ██    ██      ██   ██ ║
║   ██   ██ █████   ██      ██ ██ ████ ██ ██    ██    █████   ██████  ║
║   ██   ██ ██      ██      ██ ██  ██  ██ ██    ██    ██      ██   ██ ║
║   ██████  ███████ ███████ ██ ██      ██ ██    ██    ███████ ██   ██ ║
╚══════════════════════════════════════════════════════════════════════╝
          QUANTUM FILTRATION REACTOR — CONTAINMENT ACTIVE
```

[![Version](https://img.shields.io/badge/version-0.3.1-brightgreen?style=for-the-badge)](#)
[![Python](https://img.shields.io/badge/python-3.8+-blue?style=for-the-badge&logo=python)](#)
[![PyQt5](https://img.shields.io/badge/UI-PyQt5-orange?style=for-the-badge)](#)
[![License](https://img.shields.io/badge/license-MIT-purple?style=for-the-badge)](#)
[![Reactor](https://img.shields.io/badge/reactor-STABLE-success?style=for-the-badge)](#)
[![Biosafety](https://img.shields.io/badge/NSF%2FANSI%2049-COMPLIANT-red?style=for-the-badge)](#)

---

## ☢️ WARNING: REACTOR AREA — AUTHORIZED PERSONNEL ONLY

**Delimiter** is a desktop application for the **inspection, filtering, transformation, and synthesis** of delimited data files (CSV, TSV, Excel, FASTA). Built on top of PyQt5 and Pandas, it provides a high-performance, multi-threaded data processing pipeline wrapped in a theme that treats every button press like a nuclear event.

Designed at the **Universidade Federal de Minas Gerais (UFMG)** by **Victor S. Caricatte De Araújo** (`victorsc@ufmg.br`), Delimiter was conceived with bioinformatics pipelines in mind — but its reactor is powerful enough to process *any* structured tabular data.

---

## 🗺️ Table of Contents

1. [Core Architecture](#-core-architecture)
2. [Prerequisites](#-prerequisites)
3. [Installation](#-installation--reactor-ignition-sequence)
4. [Project Structure](#-project-structure)
5. [Features — Reactor Subsystems](#-features--reactor-subsystems)
6. [The Filtration Lab — Detailed Walkthrough](#-the-filtration-lab--detailed-walkthrough)
7. [Export Synthesis Protocols](#-export-synthesis-protocols)
8. [Internationalization — The Rosetta Stone Ribosome](#-internationalization--the-rosetta-stone-ribosome)
9. [Theming — Quantum Spectrum Shift](#-theming--quantum-spectrum-shift)
10. [Logging — Biocontainment Runes](#-logging--biocontainment-runes)
11. [Known Anomalies & Troubleshooting](#-known-anomalies--troubleshooting)
12. [Biosafety Guidelines](#-biosafety-guidelines)
13. [Author & Attribution](#-author--attribution)

---

## 🔬 Core Architecture

The application is organized around four modules — each fulfilling a specific role in the reactor pipeline:

| Module | Codename | Role |
|---|---|---|
| `delimeter.py` | **Particle Accelerator** | Main entry point. Orchestrates all threads, signals, and UI events. Instantiates the `Supreme_Genomic_Overlord`. |
| `Interface.py` | **Genomic Holodeck UI** | Full PyQt5 graphical interface. All visual widgets, dialogs, and layout definitions. |
| `logic.py` | **Sterilized Data Engine** | All data processing logic. Filtering, merging, column transmutation, NaN irradiation, deduplication, and export. |
| `xenoglossy_codex.py` | **Rosetta Stone Ribosome** | Internationalization engine. Stores and retrieves UI strings in EN, PT-BR, and ES. |

The `Supreme_Genomic_Overlord` class in `delimeter.py` acts as the **central control rod** of the reactor — it connects all PyQt5 signals from the UI to their corresponding methods in the data engine, spawns and manages `QThread` workers, and prevents the application from going critical.

---

## 🛠️ Prerequisites

Before igniting the reactor, ensure the following particles are available in your environment:

- **Python** ≥ 3.8
- **PyQt5** — UI framework
- **pandas** — Tabular data engine
- **openpyxl** — Excel file synthesis
- **psutil** *(optional but recommended)* — RAM overload detection before injection

---

## 🚀 Installation & Reactor Ignition Sequence

**Step 1 — Clone the genetic repository:**
```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd delimiter
```

**Step 2 — Install dependencies:**
```bash
pip install PyQt5 pandas openpyxl psutil
```

**Step 3 — Verify the `themes/` directory exists** with `dark_theme.qss` and `light_theme.qss`. The reactor defaults to **Dark Matter** mode on startup.

**Step 4 — Initiate the Supreme Overlord:**
```bash
python delimeter.py
```

If a catastrophic failure occurs on startup, a `CRASH_REPORT_BIOLAB.txt` file will be synthesized in the working directory. Consult it for the full traceback before filing a containment breach report.

---

## 📁 Project Structure

```
delimiter/
│
├── delimeter.py              # ⚛️  Main reactor core — entry point
├── Interface.py              # 🖥️  Genomic Holodeck UI (PyQt5)
├── logic.py                  # 🧠  Sterilized Data Engine & worker threads
├── xenoglossy_codex.py       # 🌐  Internationalization — Rosetta Stone Ribosome
│
├── themes/
│   ├── dark_theme.qss        # 🌙  Dark Matter stylesheet
│   └── light_theme.qss       # ☀️  Solar Photon stylesheet
│
├── panvita_quarantine.log    # 📋  Runtime log — Biocontainment Runes
└── CRASH_REPORT_BIOLAB.txt   # 💥  Generated on catastrophic failure only
```

---

## ⚡ Features — Reactor Subsystems

### 🧬 Biomass Injection (File Loading)
- Supports **CSV**, **TSV**, and **Excel** (`.xlsx`) formats.
- **Auto-detection** of delimiters via `csv.Sniffer` — the `DarkMatter_Sniffer` probes the void of the first 2048 bytes of a file and determines the correct separator automatically.
- Files can be loaded by dragging into the path field or using the **🧬 Inject** button.
- Header row detection is configurable per injection.
- Before loading, a **Mass Overload Alert** is triggered if `psutil` detects that the file size × 3 exceeds available system RAM — giving the operator a chance to abort before the reactor crashes.
- File loading runs on a **dedicated `QThread`** (`Hyperdrive_Loader_Thread`) — the UI never freezes during injection, no matter how large the biological material.

---

### 🎛️ Matryoshka Filter Reactor
The signature feature of Delimiter. Filters are structured as **nested, stacked condition pods** — like a Russian Matryoshka doll — each one narrowing down the dataset further.

Each filter pod exposes the following configuration:

| Parameter | Description |
|---|---|
| **Column** | The genomic strand (column) to target |
| **Operation** | The comparison operator (see table below) |
| **Value** | The biological target value |
| **Case Sensitive** | Toggle for case-sensitive string matching |
| **Ignore NaN** | Automatically exclude null/NaN rows from the filter pass |
| **Only NaN** | Isolate and return only rows where the column is null |

**Supported Quantum Operators:**

| Operator | Description |
|---|---|
| `Equals` | Strict string equality |
| `Contains` | Substring presence check |
| `Starts With` | Prefix matching |
| `Ends With` | Suffix matching |
| `Regex` | Full regular expression matching (validated before detonation) |
| `Greater Than (>)` | Numeric comparison |
| `Less Than (<)` | Numeric comparison |
| `Numeric Equals (=)` | Exact numeric equality |

Multiple pods can be added dynamically via **➕ Add Rule**. All active filters are applied sequentially when **⚡ Detonate Collider** is pressed. The collision runs on a separate `QThread` (`Quantum_Collider_Thread`) and reports survival statistics in the Biological X-Ray panel upon completion.

---

### 💾 Preset Crystallization & Resurrection
Filter configurations can be saved and reloaded between sessions:

- **💾 Crystallize Memory** — Serializes all active Matryoshka pods to a `.json` file via `Elephant_Brain_Storage`.
- **📂 Resurrect Filters** — Deserializes and re-populates the UI with a previously saved filter set. All existing pods are vaporized before resurrection to prevent contamination.

---

### 💀 Purge Clones (Deduplication)
The **💀 Purge Clones** operation removes all duplicate rows from the active dataset across all columns. Executed synchronously against the `Sterilized_Data_Engine` and triggers a full UI refresh.

---

### 🧬 Alien Splicing (Dataset Merge)
Inject foreign genetic material — load a secondary file and merge it with the active dataset using a key column. This is a standard `pandas` inner join operation, wrapped in a dedicated dialog (`Xenomorph_Splicer_Dialog`). The fused abomination replaces the current working dataset. A `gc.collect()` is performed post-merge to reclaim memory.

---

### ☢️ N-Terminal Radiation (NaN Imputation)
Fill missing values (`NaN`) in a selected column using one of three irradiation methods:

| Method | Description |
|---|---|
| `Mean` | Replace NaN with the column's arithmetic mean |
| `Median` | Replace NaN with the column's median value |
| `Fixed Value` | Replace NaN with a user-specified constant |

---

### ⚗️ Column Transcription (Column Transmutation)
Create a **new derived column** from two existing ones using arithmetic or string operations:

| Spell | Description |
|---|---|
| `Concat` | String concatenation of two columns |
| `Add (+)` | Numeric addition |
| `Subtract (-)` | Numeric subtraction |
| `Multiply (*)` | Numeric multiplication |
| `Divide (/)` | Numeric division |

The resulting column is appended to the dataset and given a user-defined name.

---

### 📊 Spectrometry (Column Statistics)
Opens a dialog displaying descriptive statistics (`pandas .describe()`) for the currently selected column in the deep scan panel. Includes count, unique values, mean, standard deviation, min/max, and quartile distribution.

---

### 🚨 Illuminate Anomalies (Conditional Cell Highlighting)
Define visual rules to highlight cells in the data table based on their values. Each rule targets a column (or `ALL_COLUMNS`), specifies an operator and threshold value, and assigns a highlight color (Red, Green, or Yellow). These rules are applied live to the `Infinite_Hologram_Model` — the custom `QAbstractTableModel` that drives the table view.

---

### ⏪ Undo / ⏩ Redo (Temporal Dilation)
The `Sterilized_Data_Engine` maintains a state history stack. Every mutating operation (filter detonation, CRISPR cell edit, clone purge, alien splice, radiation, transmutation) saves a snapshot of the current filtered dataset. The undo (`rewind_time_dilation`) and redo (`temporal_redo`) operations traverse this history, allowing operators to reverse or replay transformations without reloading from disk.

---

### ✏️ CRISPR Cell Editing
Individual cells in the data table are directly editable. Double-clicking a cell allows in-place modification. Changes are routed through `signal_crispr_edit` and committed to the `Sterilized_Data_Engine`, which updates the filtered dataset and pushes a new undo state.

---

### 🔬 Biological X-Ray Panel
A live statistics dashboard, updated after every filter collision and data mutation:

- **Raw Reads** — Total number of rows in the original injected file.
- **Survivors** — Number of rows remaining after the current filter configuration.
- **Rate** — Survival rate as a percentage, computed as `(survivors / raw_reads) × 100`.

---

## 📤 Export Synthesis Protocols

When ready to synthesize output, four formats are available. All exports run on a dedicated `QThread` (`Chunk_Devourer_Exporter` or `Fasta_Devourer_Exporter`) with live progress reporting and a **Temporal Paradox Warning** if the destination file already exists.

| Format | Thread | Notes |
|---|---|---|
| **CSV** | `Chunk_Devourer_Exporter` | Written in 50,000-row chunks for memory efficiency |
| **TSV** | `Chunk_Devourer_Exporter` | Same chunked pipeline, tab-separated |
| **Excel** | `Chunk_Devourer_Exporter` | Single-pass via `pandas.to_excel` (`openpyxl`) |
| **FASTA** | `Fasta_Devourer_Exporter` | Requires mapping a header column and a sequence column. Exports in standard `>header\nsequence` format, batch-reporting progress every 1,000 rows. |

All exported files default to the name `{original_filename}_mutated.{ext}`.

---

## 🌐 Internationalization — The Rosetta Stone Ribosome

All UI strings are managed by `Rosetta_Stone_Ribosome` in `xenoglossy_codex.py`. The codex maps string keys to translations across three supported dialects:

| Code | Language |
|---|---|
| `EN` | English |
| `PT-BR` | Brazilian Portuguese |
| `ES` | Spanish |

**Runtime Usage:**
```python
# Switch the active dialect
Rosetta_Stone_Ribosome.mutate_dialect("PT-BR")

# Retrieve a localized string
label = Rosetta_Stone_Ribosome.extract_peptide("btn_inject")
# → "🧬 Injetar"
```

Dialect switching is exposed in the UI. Selecting a language from the combo box fires `signal_dialect_mutation`, which calls `mutate_dialect()` and then triggers `retranscribe_ui()` on the Holodeck UI — causing all visible labels, buttons, and group boxes to re-render with the new dialect instantly, without restarting the application.

To **add a new language**, extend every entry in `_genetic_codex` with the new language code and its translation string.

---

## 🎨 Theming — Quantum Spectrum Shift

Delimiter supports two visual modes, switched at runtime via the **☀️ / 🌙** toggle button:

| Theme | File | Description |
|---|---|---|
| **Dark Matter** | `themes/dark_theme.qss` | Default. Dark background, high-contrast accent colors. |
| **Solar Photon** | `themes/light_theme.qss` | Light mode for operators who fear the void. |

Themes are loaded by `Chromatic_Mutator.force_mutation()`, which reads the `.qss` file and applies it to the `QApplication` instance. If a theme file is not found in `themes/`, the system falls back to searching the working directory and logs an **Aesthetic Anomaly** warning.

To create a custom theme, author a new `.qss` file following Qt stylesheet syntax and load it by modifying the `Chromatic_Mutator.force_mutation()` call.

---

## 📋 Logging — Biocontainment Runes

All significant operations are engraved into `panvita_quarantine.log` via `Runes_Of_Biocontainment.engrave_biocontainment_runes()`. This includes:

- Session initiation
- Filter collisions (how many Matryoshka pods were applied)
- FASTA / CSV / TSV / Excel synthesis completions
- Alien splice events (which key columns were used)
- Column transmutation events (which new strand was created)
- Error traces on containment breach

Log format:
```
YYYY-MM-DD HH:MM:SS,ms - BIOLOGICAL ALERT - [message]
```

On a **catastrophic startup failure**, a separate `CRASH_REPORT_BIOLAB.txt` is written to the working directory with the complete Python traceback before the process exits with code `1`.

---

## 🐛 Known Anomalies & Troubleshooting

| Symptom | Probable Cause | Resolution |
|---|---|---|
| Application fails to start | Missing dependency | Run `pip install PyQt5 pandas openpyxl` |
| Theme not applied | `themes/` directory missing | Create the `themes/` folder and add `.qss` files |
| `Mass Overload Alert` on small files | `psutil` not installed | `pip install psutil` — or proceed and accept the risk |
| Filter returns 0 rows unexpectedly | Case sensitivity or data type mismatch | Check operator type. Use `Numeric Equals` for numbers, not `Equals`. |
| FASTA export fails | Selected columns contain NaN | Apply N-Terminal Radiation before exporting |
| Crash report generated | Critical import failure or unhandled exception | Read `CRASH_REPORT_BIOLAB.txt` — check the traceback |
| UI freezes during load | Threading issue — should not occur | File a containment breach report with the log |

---

## ⚠️ Biosafety Guidelines

> *Brazil does not have an ABNT standard for Biological Safety Cabinets. In the absence of national regulation, the applicable standard is **NSF/ANSI 49**.*

This application processes real data with real operations. The following protocols are mandatory:

- **Always keep a backup copy** of source files before injection and processing.
- **Verify survival rates** in the X-Ray panel after every filter collision before exporting.
- **Do not force-quit** during an active synthesis thread — partial output files will be produced.
- **Audit the `panvita_quarantine.log`** after processing large datasets to verify that all operations completed as expected.
- When a **Temporal Paradox Warning** is raised for file overwrites, confirm deliberately — this action is irreversible from within the application.

---

## 👤 Author & Attribution

| Field | Detail |
|---|---|
| **Author** | Victor S. Caricatte De Araújo |
| **Email** | victorsc@ufmg.br |
| **Institution** | Universidade Federal de Minas Gerais (UFMG) |
| **Version** | 0.3.1 |
| **Codename** | Particle Accelerator |

---

<p align="center">
  <b>🧬 DELIMITER — Filtration complete. Survivors reported. Reactor cooling down. 🧬</b><br>
  <i>Built at UFMG. Powered by Python, Pandas, and the will to detonate data.</i>
</p>
