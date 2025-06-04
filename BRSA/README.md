# BRSA - Basic RNA-Seq Analysis

**Author**: Victor Silveira Caricatte  
**Institution**: Federal University of Minas Gerais (UFMG)

## Overview

**BRSA (Basic RNA-Seq Analysis)** is a graphical application built with PyQt6 that facilitates accessible, reproducible, and automated RNA-Seq data analysis. It allows users to perform everything from data import and normalization to differential expression analysis, visualizations, and functional enrichment.

---

## Features

### 📂 File Import
- **Counts File**: CSV format with genes as rows and samples as columns.
- **Metadata File**: CSV format with at least `sample` and `condition` columns.

### ⚙️ Analysis Parameters
- Selection of two experimental conditions for comparison.
- Filter for low-count genes.
- p-value cutoff with FDR correction.

### 🔬 Analysis Steps

1. **Preprocessing**:
   - Filters low-expression genes.
   - Validates input data consistency.

2. **Normalization**:
   - Counts Per Million (CPM) method.

3. **Differential Expression**:
   - Student’s t-test.
   - False Discovery Rate (FDR) correction (Benjamini-Hochberg).
   - Calculates log2 Fold Change.

4. **Visualizations**:
   - PCA (Principal Component Analysis).
   - Heatmap of top 50 most variable genes.
   - Interactive Volcano Plot (with `mplcursors`).

5. **Functional Enrichment**:
   - Supported databases:
     - GO (BP, MF, CC)
     - KEGG
     - Reactome
     - Disease Ontology
   - Uses `gprofiler` for querying.

6. **Export Options**:
   - 📄 **PDF report** with graphs and tables (via `reportlab`).
   - 📊 **Excel Workbook** with separate sheets for DE results, normalized counts, and enrichment analysis (via `openpyxl`).

7. **Session Management**:
   - Save and load analysis sessions as `.json` files.

---

## Graphical Interface

BRSA provides a user-friendly interface with multiple tabs:

- **Input**: File selection and parameter configuration.
- **Results**:
  - `Differential Expression`: Summary of significant genes.
  - `Visualization`: Plot selector and display.
  - `Enrichment Analysis`: Run and view functional analysis results.
- **Help**: Instructions and application overview.

---

## Requirements

### 🐍 Python
- Version: Python 3.7 or higher

### 📦 Dependencies

Install all required packages with:

```bash
pip install -r requirements.txt
