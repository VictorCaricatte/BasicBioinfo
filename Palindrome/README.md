# DNA Palindrome Finder

**Version:** 1.1
**Author:** Victor S Caricatte De Araújo
**Email:** victorleniwys@gmail.com or victorsc@ufmg.br
**Institution:** Universidade Federal de Minas Gerais (UFMG)
**Date:** Apr, 14

## Table of Contents

1.  [Description](#description)
2.  [What are DNA Palindromes?](#what-are-dna-palindromes)
3.  [Features](#features)
4.  [How to Use](#how-to-use)
5.  [User Interface](#user-interface)
    * [Main Tab (DNA Analysis)](#main-tab-dna-analysis)
    * [Help Tab](#help-tab)
    * [Menu](#menu)
6.  [Requirements](#requirements)
7.  [Installing Dependencies](#installing-dependencies)
8.  [How to Run](#how-to-run)
9.  [Technical Details](#technical-details)
10. [About the Author](#about-the-author)

## Description

The **DNA Palindrome Finder** is a bioinformatics tool designed to identify palindromic sequences in DNA molecules. DNA palindromes are sequences that read the same on one strand as they do on their reverse complementary strand. These sequences are biologically significant as they often constitute recognition and cleavage sites for restriction enzymes.

This application provides a user-friendly graphical interface to facilitate the analysis of DNA sequences, allowing users to input sequences manually, load them from FASTA files, and configure parameters for palindrome searching.

## What are DNA Palindromes?

In molecular biology, a DNA palindrome is a nucleotide sequence that is identical to its reverse complement strand.
For example, the sequence 5'-GAATTC-3' is a palindrome because its complementary strand is 3'-CTTAAG-5'. If we read this complementary strand in reverse (from right to left, in the 5' to 3' direction), we get 5'-GAATTC-3', which is identical to the original sequence.
This specific example (GAATTC) is the recognition site for the EcoRI restriction enzyme.

## Features

* **Flexible Sequence Input:**
    * Manual entry of DNA sequences in the text area.
    * Loading sequences from files in FASTA format (`.fasta`, `.fa`, `.faa`, `.fna`).
* **Configurable Search Parameters:**
    * Definition of minimum and maximum length for the palindromes to be identified.
* **Sequence Validation:**
    * Checks if the input sequence contains only valid bases (A, T, C, G), ignoring spaces and line breaks.
* **Efficient Analysis:**
    * Optimized algorithm for palindrome searching.
    * Uses threading to prevent the graphical interface from freezing during long analyses.
* **Results Visualization:**
    * Presentation of found palindromes in a table with columns: Position (start-end), Palindromic Sequence, and Length.
* **Results Management:**
    * Copy results to the clipboard.
    * Save results to text (`.txt`) or CSV (`.csv`) files.
    * Clear input fields and results for a new analysis.
* **User-Friendly Interface:**
    * Intuitive graphical interface built with Tkinter.
    * "Help" tab with documentation and program information.
    * Status bar for user feedback.

## How to Use

1.  **Insert DNA Sequence:**
    * Type or paste the DNA sequence directly into the "DNA Input" text area.
    * Alternatively, go to `Menu > File > Open FASTA file...` to load a sequence from a FASTA file. Only A, T, C, G characters will be considered (case-insensitive).
2.  **Set Parameters:**
    * Adjust the "Min length" and "Max length" fields for the palindromes you wish to find. The default values are 4 and 12, respectively.
3.  **Analyze:**
    * Click the "Analyze DNA" button.
    * The status bar will indicate the analysis progress.
4.  **View Results:**
    * Found palindromes will be listed in the "Results" table, showing their position (1-based), the sequence, and its length.
5.  **Manage Results (Optional):**
    * **Copy Results:** Copies the table data to the clipboard.
    * **Save Results:** Saves the table data to a text or CSV file.
    * **Clear All:** Clears the DNA input area, the results table, and resets lengths to default values.

## User Interface

The main interface is divided into tabs and a menu.

### Main Tab (DNA Analysis)

* **DNA Input:** A text area to enter or view the DNA sequence.
* **Settings Frame:**
    * **Min length:** Field to set the minimum palindrome length.
    * **Max length:** Field to set the maximum palindrome length.
    * **Analyze DNA:** Button to start the palindrome search.
* **Results Frame:**
    * **Results Table:** Displays found palindromes with "Position", "Sequence", and "Length" columns.
    * **Copy Results:** Button to copy table data.
    * **Save Results:** Button to save table data to a file.
    * **Clear All:** Button to clear all fields and results.
* **Status Bar:** Located at the bottom of the window, it displays the current application status.

### Help Tab

Contains information about the program, how to use it, the technical definition of DNA palindromes, and developer details.

### Menu

* **File:**
    * **Open FASTA file...:** Opens a dialog to select and load a FASTA file.
    * **Save results...:** Opens a dialog to save the current results.
    * **Exit:** Closes the application.
* **Help:**
    * **Documentation:** Selects the Help tab.
    * **About:** Displays a window with information about the software version and purpose.

## Requirements

* Python 3.x
* Tkinter (usually included in the standard Python installation)
* Biopython (the `Bio.Seq` library is imported, although the palindrome checking logic is custom in the script)

## Installing Dependencies

If you don't have Biopython installed, you can install it using pip:

```bash
pip install biopython
