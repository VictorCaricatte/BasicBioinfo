# Delimited File Filter for genetic files

**Version:** 1.1.2
**Date:** April, 11
**Author:** Victor S Caricatte De Araújo
**Email:** victorleniwys@gmail.com or victorsc@ufmg.br
**Institution:** Universidade Federal de Minas Gerais (UFMG)
**GitHub:** [VictorCaricatte](https://github.com/VictorCaricatte)

## Description

A graphical application that reads various file formats (CSV, TSV, Excel) and allows filtering rows based on specific values in a given column. Features a tabbed interface with filtering controls and a comprehensive help section. This tool is designed for data manipulation and filtering, particularly useful in contexts of genetic data analysis or general tabular data.

## Key Features

* **Multiple File Format Support:** Reads CSV, TSV, and Excel (.xlsx, .xls) files.
* **Tabbed Interface:** Organizes functionalities into "Filter" and "Help" tabs.
* **Multiple Filtering Conditions:**
    * Equals
    * Contains
    * Starts With
    * Ends With
    * Regex (Regular Expression)
* **Sensitivity Options:** Allows for case-sensitive or case-insensitive filtering.
* **Invert Match:** Excludes rows that match the filter criteria.
* **Header Row Handling:** Option to treat the first row as a header.
* **Results Preview:** Displays a preview of the filtered data (up to 1000 rows).
* **Flexible Export:** Exports filtered results to CSV, TSV, or Excel.
* **Theme Toggle:** Supports light and dark themes for better visualization.

## Prerequisites

To run this software, you will need Python 3 installed, along with the following libraries:

* pandas
* PyQt5
* openpyxl (for Excel file support)

## Installation

1.  Clone or download this repository.
2.  Navigate to the project directory.
3.  Install the dependencies using pip:

    ```bash
    pip install pandas PyQt5 openpyxl
    ```
   

## How to Use

1.  Run the `delimitergenfile.py` script:

    ```bash
    python delimitergenfile.py
    ```

2.  **In the "Filter" tab:**
    * **File Selection:**
        * Click "Browse..." to select your data file.
        * The file type can be auto-detected ("Auto Detect") or specified (CSV, TSV, Excel).
        * Click "Load File" to load the data.
    * **Filter Settings:**
        * **Column:** Select the column on which the filter will be applied (enabled after loading the file).
        * **Operation:** Choose the type of filtering operation (Equals, Contains, etc.).
        * **Value:** Enter the value to search/compare against.
        * **Options:**
            * `Case Sensitive`: Check to differentiate between uppercase and lowercase.
            * `Invert Match`: Check to get rows that *do not* match the filter.
            * `Has Header`: Check if your file has a header row (checked by default).
    * **Actions:**
        * Click "Apply Filter" to view the results in the "Results Preview" table.
        * Click "Reset Filter" to clear the applied filter and show all original data.
    * **Export Results:**
        * After filtering, click "Export to CSV", "Export to TSV", or "Export to Excel" to save the filtered data.

3.  **In the "Help" tab:**
    * Find a quick start guide, description of filter options, supported formats, tips, and contact information.

4.  **"View" Menu:**
    * Allows toggling between the dark theme (default) and the light theme.

## Tips

* For large files, applying complex filters may take some time.
* Regex operations support full Python regular expression syntax.
* You can preview the results before exporting.

## Contributions

Contributions are welcome! Feel free to open an issue or submit a pull request.

## License

This project have a MIT license at this moment.
