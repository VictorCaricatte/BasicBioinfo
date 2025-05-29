# LengthAnaly.py - FASTA/FASTQ Sequence Analyzer

**Author:** Victor S. Caricatte De Araújo
**Email:** victorleniwys@gmail.com or victorsc@ufmg.br
**Institution:** Universidade Federal de Minas Gerais (UFMG)
**Version:** 1.1.9
**Version Date:** April 28

## Overview

`LengthAnaly.py` is a graphical user interface (GUI) and command-line interface (CLI) tool developed in Python for performing analyses on biological sequence files. Its main functionalities include:

1.  **Sequence Length Analysis (FASTA Files):**
    * Loads multiple FASTA files.
    * Calculates descriptive statistics of sequence lengths (minimum, maximum, mean, median, standard deviation, quartiles).
    * Generates visualizations such as histograms (with an option for KDE) and boxplots to compare length distributions.
2.  **Base Composition Analysis (DNA/RNA Sequences):**
    * Allows manual input of a DNA/RNA sequence.
    * Calculates the count and frequency of each base (A, T, C, G, U, N).
    * Generates a bar chart of the base composition.
    * Allows customization of base colors in the chart.
3.  **Genome Coverage Estimation (FASTQ Files):**
    * Loads multiple FASTQ files.
    * Calculates estimated genome coverage based on the number of reads, genome size, read length, and sequencing type (single-end or paired-end).

The application utilizes libraries such as PyQt5 for the graphical interface, Matplotlib and Seaborn for graph generation, Biopython for parsing FASTA/FASTQ files, and NumPy for statistical calculations.

## Key Features

* **Intuitive Graphical Interface (GUI):** User-friendly, with tabs for different analysis types.
* **Sequence Length Analysis (FASTA):**
    * Add and remove multiple FASTA files.
    * Assign custom labels for each dataset.
    * Choose between histogram and boxplot for visualization.
    * Adjust the number of bins and KDE display for histograms.
    * View detailed statistics.
    * Save plots (PNG, JPG, PDF, SVG) and statistics (TXT).
* **Base Composition Analysis:**
    * Input DNA/RNA sequence directly into the interface.
    * Graphical visualization of base composition.
    * Customize base colors in the plot.
    * Display composition statistics (counts and percentages).
    * Save the composition plot.
* **FASTQ Coverage Estimation:**
    * Add and remove multiple FASTQ files.
    * Set parameters: genome size, read length, and sequencing type (paired-end/single-end).
    * View coverage estimation results.
    * Save results (TXT).
* **Integrated Help Tab:** Built-in documentation explaining the use of each feature and the coverage calculation formula.
* **Command-Line Interface (CLI):** For quick FASTA sequence length analyses without needing to open the GUI.

## Prerequisites

* Python 3.x
* Python Libraries:
    * PyQt5
    * Matplotlib
    * Seaborn
    * NumPy
    * Biopython

## Installation

1.  Clone the repository or download the `lenghtanaly.py` file.
2.  Install the necessary dependencies. You can use pip:

    ```bash
    pip install PyQt5 matplotlib seaborn numpy biopython
    ```

## Usage

### Graphical User Interface (GUI)

To start the application in GUI mode, run the script without arguments:

```bash
python lenghtanaly.py
