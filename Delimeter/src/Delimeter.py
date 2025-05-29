#!/usr/bin/env python3
"""
# ==============================================================================
# Author:       Victor S Caricatte De Araújo
# Email:        victorleniwys@gmail.com or victorsc@ufmg.br
# Intitution:   Universidade federal de Minas Gerais
# Version:      1.1.2
# Date:         Abr, 11
# ...................................
# ==============================================================================
"""
"""
Delimited File Filter

Description:
    A graphical application that reads various file formats (CSV, TSV, Excel) and allows filtering rows
    based on specific values in a given column. Features a tabbed interface with filtering controls
    and a comprehensive help section.

Features:
    - Supports multiple file formats: CSV, TSV, Excel (xlsx, xls)
    - Tabbed interface with filter controls and help section
    - Multiple filtering conditions (equals, contains, starts/ends with, regex)
    - Case-sensitive and case-insensitive options
    - Invert matching (exclude matching rows)
    - Header row handling
    - Preview of filtered results
    - Export to multiple formats
    - Dark/light theme toggle

Dependencies:
    - pandas
    - PyQt5
    - openpyxl (for Excel support)

To install dependencies:
    pip install pandas PyQt5 openpyxl
"""

import sys
import re
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout,
                             QLabel, QComboBox, QLineEdit, QPushButton, QCheckBox, QTextEdit,
                             QFileDialog, QTableWidget, QTableWidgetItem, QMessageBox, QGroupBox,
                             QMenuBar, QMenu, QAction)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

class FileFilterApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Delimited File Filter for genetic files")
        self.setGeometry(100, 100, 1000, 700)
        
        # Initialize variables
        self.df = None
        self.filtered_df = None
        self.file_path = ""
        self.file_type = ""
        self.dark_theme = True  # Default to dark theme
        
        # Create main widget and layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.main_layout = QVBoxLayout(self.main_widget)
        
        # Create menu bar
        self.create_menu_bar()
        
        # Create tab widget
        self.tabs = QTabWidget()
        self.main_layout.addWidget(self.tabs)
        
        # Create tabs
        self.create_filter_tab()
        self.create_help_tab()
        
        # Status bar
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready")
        
        # Apply dark theme by default
        self.apply_theme()
    
    def create_menu_bar(self):
        """Create the menu bar with theme toggle option"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("File")
        
        open_action = QAction("Open", self)
        open_action.triggered.connect(self.browse_file)
        file_menu.addAction(open_action)
        
        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # View menu
        view_menu = menubar.addMenu("View")
        
        self.theme_action = QAction("Switch to Light Theme", self)
        self.theme_action.triggered.connect(self.toggle_theme)
        view_menu.addAction(self.theme_action)
    
    def toggle_theme(self):
        """Toggle between dark and light theme"""
        self.dark_theme = not self.dark_theme
        self.apply_theme()
        
        # Update menu text
        if self.dark_theme:
            self.theme_action.setText("Switch to Light Theme")
        else:
            self.theme_action.setText("Switch to Dark Theme")
    
    def apply_theme(self):
        """Apply the current theme (dark or light)"""
        if self.dark_theme:
            self.setStyleSheet("""
                QMainWindow {
                    background-color: #2D2D2D;
                }
                QTabWidget::pane {
                    border: 1px solid #444;
                    background: #2D2D2D;
                }
                QTabBar::tab {
                    background: #3D3D3D;
                    color: #DDD;
                    padding: 8px;
                    border: 1px solid #444;
                    border-bottom: none;
                }
                QTabBar::tab:selected {
                    background: #505050;
                    color: white;
                }
                QLabel {
                    color: #DDD;
                }
                QLineEdit, QComboBox, QTextEdit, QTableWidget {
                    background-color: #3D3D3D;
                    color: #DDD;
                    border: 1px solid #555;
                }
                QPushButton {
                    background-color: #505050;
                    color: #DDD;
                    border: 1px solid #555;
                    padding: 5px;
                }
                QPushButton:hover {
                    background-color: #606060;
                }
                QGroupBox {
                    border: 1px solid #555;
                    border-radius: 3px;
                    margin-top: 10px;
                    color: #DDD;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 3px;
                }
                QHeaderView::section {
                    background-color: #505050;
                    color: white;
                    padding: 4px;
                    border: 1px solid #444;
                }
                QTableWidget {
                    gridline-color: #444;
                }
                QMenuBar {
                    background-color: #3D3D3D;
                    color: #DDD;
                }
                QMenuBar::item {
                    background-color: transparent;
                    padding: 5px 10px;
                }
                QMenuBar::item:selected {
                    background-color: #505050;
                }
                QMenu {
                    background-color: #3D3D3D;
                    color: #DDD;
                    border: 1px solid #555;
                }
                QMenu::item:selected {
                    background-color: #505050;
                }
            """)
        else:
            self.setStyleSheet("""
                QMainWindow {
                    background-color: #F0F0F0;
                }
                QTabWidget::pane {
                    border: 1px solid #AAA;
                    background: #F0F0F0;
                }
                QTabBar::tab {
                    background: #E0E0E0;
                    color: #333;
                    padding: 8px;
                    border: 1px solid #AAA;
                    border-bottom: none;
                }
                QTabBar::tab:selected {
                    background: #FFFFFF;
                    color: #000;
                }
                QLabel {
                    color: #333;
                }
                QLineEdit, QComboBox, QTextEdit, QTableWidget {
                    background-color: #FFFFFF;
                    color: #333;
                    border: 1px solid #AAA;
                }
                QPushButton {
                    background-color: #E0E0E0;
                    color: #333;
                    border: 1px solid #AAA;
                    padding: 5px;
                }
                QPushButton:hover {
                    background-color: #D0D0D0;
                }
                QGroupBox {
                    border: 1px solid #AAA;
                    border-radius: 3px;
                    margin-top: 10px;
                    color: #333;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 3px;
                }
                QHeaderView::section {
                    background-color: #E0E0E0;
                    color: #333;
                    padding: 4px;
                    border: 1px solid #AAA;
                }
                QTableWidget {
                    gridline-color: #AAA;
                }
                QMenuBar {
                    background-color: #E0E0E0;
                    color: #333;
                }
                QMenuBar::item {
                    background-color: transparent;
                    padding: 5px 10px;
                }
                QMenuBar::item:selected {
                    background-color: #D0D0D0;
                }
                QMenu {
                    background-color: #FFFFFF;
                    color: #333;
                    border: 1px solid #AAA;
                }
                QMenu::item:selected {
                    background-color: #E0E0E0;
                }
            """)
    
    def create_filter_tab(self):
        """Create the main filter tab with all controls"""
        self.filter_tab = QWidget()
        self.tabs.addTab(self.filter_tab, "Filter")
        
        layout = QVBoxLayout(self.filter_tab)
        
        # File selection group
        file_group = QGroupBox("File Selection")
        file_layout = QHBoxLayout()
        
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select a file...")
        file_layout.addWidget(self.file_path_edit)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_file)
        file_layout.addWidget(browse_btn)
        
        self.file_type_combo = QComboBox()
        self.file_type_combo.addItems(["Auto Detect", "CSV", "TSV", "Excel"])
        file_layout.addWidget(self.file_type_combo)
        
        load_btn = QPushButton("Load File")
        load_btn.clicked.connect(self.load_file)
        file_layout.addWidget(load_btn)
        
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        # Filter controls group
        filter_group = QGroupBox("Filter Settings")
        filter_layout = QVBoxLayout()
        
        # Column selection
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Column:"))
        
        self.column_combo = QComboBox()
        self.column_combo.setEnabled(False)
        row1.addWidget(self.column_combo)
        
        row1.addWidget(QLabel("Operation:"))
        
        self.operation_combo = QComboBox()
        self.operation_combo.addItems(["Equals", "Contains", "Starts With", "Ends With", "Regex"])
        row1.addWidget(self.operation_combo)
        
        filter_layout.addLayout(row1)
        
        # Value input
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Value:"))
        
        self.value_edit = QLineEdit()
        row2.addWidget(self.value_edit)
        
        filter_layout.addLayout(row2)
        
        # Options
        row3 = QHBoxLayout()
        
        self.case_check = QCheckBox("Case Sensitive")
        row3.addWidget(self.case_check)
        
        self.invert_check = QCheckBox("Invert Match")
        row3.addWidget(self.invert_check)
        
        self.header_check = QCheckBox("Has Header", checked=True)
        row3.addWidget(self.header_check)
        
        filter_layout.addLayout(row3)
        
        # Action buttons
        row4 = QHBoxLayout()
        
        apply_btn = QPushButton("Apply Filter")
        apply_btn.clicked.connect(self.apply_filter)
        row4.addWidget(apply_btn)
        
        reset_btn = QPushButton("Reset Filter")
        reset_btn.clicked.connect(self.reset_filter)
        row4.addWidget(reset_btn)
        
        filter_layout.addLayout(row4)
        
        filter_group.setLayout(filter_layout)
        layout.addWidget(filter_group)
        
        # Results preview
        preview_group = QGroupBox("Results Preview")
        preview_layout = QVBoxLayout()
        
        self.results_table = QTableWidget()
        self.results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        preview_layout.addWidget(self.results_table)
        
        # Export buttons
        export_layout = QHBoxLayout()
        
        export_csv_btn = QPushButton("Export to CSV")
        export_csv_btn.clicked.connect(lambda: self.export_results("csv"))
        export_layout.addWidget(export_csv_btn)
        
        export_tsv_btn = QPushButton("Export to TSV")
        export_tsv_btn.clicked.connect(lambda: self.export_results("tsv"))
        export_layout.addWidget(export_tsv_btn)
        
        export_excel_btn = QPushButton("Export to Excel")
        export_excel_btn.clicked.connect(lambda: self.export_results("excel"))
        export_layout.addWidget(export_excel_btn)
        
        preview_layout.addLayout(export_layout)
        preview_group.setLayout(preview_layout)
        layout.addWidget(preview_group)
    
    def create_help_tab(self):
        """Create the help tab with documentation"""
        self.help_tab = QWidget()
        self.tabs.addTab(self.help_tab, "Help")
        
        layout = QVBoxLayout(self.help_tab)
        
        help_text = QTextEdit()
        help_text.setReadOnly(True)
        help_text.setFont(QFont("Arial", 10))
        
        help_content = """
        <h1>Delimited File Filter Help</h1>
        
        <h2>Overview</h2>
        <p>This application allows you to filter rows in various file formats (CSV, TSV, Excel) 
        based on values in a specific column.</p>
        
        <h2>Getting Started</h2>
        <ol>
            <li>Click "Browse" to select your input file</li>
            <li>Select the appropriate file type or leave as "Auto Detect"</li>
            <li>Click "Load File" to load the data</li>
            <li>Configure your filter settings</li>
            <li>Click "Apply Filter" to see the results</li>
            <li>Export the filtered results if needed</li>
        </ol>
        
        <h2>Filter Options</h2>
        <ul>
            <li><b>Column:</b> Select which column to filter on</li>
            <li><b>Operation:</b> Choose the type of matching to perform:
                <ul>
                    <li><b>Equals:</b> Exact match</li>
                    <li><b>Contains:</b> Substring match</li>
                    <li><b>Starts With:</b> Match beginning of string</li>
                    <li><b>Ends With:</b> Match end of string</li>
                    <li><b>Regex:</b> Regular expression match</li>
                </ul>
            </li>
            <li><b>Value:</b> The value to match against</li>
            <li><b>Case Sensitive:</b> Toggle case sensitivity</li>
            <li><b>Invert Match:</b> Exclude rows that match the criteria</li>
            <li><b>Has Header:</b> Specify if the file has a header row</li>
        </ul>
        
        <h2>Supported File Formats</h2>
        <ul>
            <li><b>CSV:</b> Comma-separated values</li>
            <li><b>TSV:</b> Tab-separated values</li>
            <li><b>Excel:</b> .xlsx and .xls files</li>
        </ul>
        
        <h2>Tips</h2>
        <ul>
            <li>For large files, applying complex filters may take some time</li>
            <li>Regex operations support full Python regular expression syntax</li>
            <li>You can preview the results before exporting</li>
            <li>Toggle between light and dark theme from the View menu</li>
        </ul>

        <h2>Credits</h2>
        <p> Victor Silveira Caricatte/ Universidade Federal de Minas Gerais.<p>

        <h2>Contact</h2>
        <p>https://github.com/VictorCaricatte or victorsc@ufmg.br<p>

        """
        
        help_text.setHtml(help_content)
        layout.addWidget(help_text)
    
    def browse_file(self):
        """Open file dialog to select input file"""
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self,
            "Open File",
            "",
            "All Supported Files (*.csv *.tsv *.xlsx *.xls);;CSV Files (*.csv);;TSV Files (*.tsv);;Excel Files (*.xlsx *.xls)"
        )
        
        if file_path:
            self.file_path_edit.setText(file_path)
    
    def load_file(self):
        """Load the selected file into a DataFrame"""
        self.file_path = self.file_path_edit.text()
        if not self.file_path:
            QMessageBox.warning(self, "Warning", "Please select a file first")
            return
        
        # Determine file type
        file_type = self.file_type_combo.currentText()
        if file_type == "Auto Detect":
            if self.file_path.lower().endswith('.csv'):
                file_type = "CSV"
            elif self.file_path.lower().endswith('.tsv'):
                file_type = "TSV"
            elif self.file_path.lower().endswith(('.xlsx', '.xls')):
                file_type = "Excel"
        
        try:
            if file_type == "CSV":
                self.df = pd.read_csv(self.file_path, header=0 if self.header_check.isChecked() else None)
            elif file_type == "TSV":
                self.df = pd.read_csv(self.file_path, sep='\t', header=0 if self.header_check.isChecked() else None)
            elif file_type == "Excel":
                self.df = pd.read_excel(self.file_path, header=0 if self.header_check.isChecked() else None)
            
            # Update column combo
            self.column_combo.clear()
            if self.header_check.isChecked() and isinstance(self.df.columns, pd.core.indexes.base.Index):
                self.column_combo.addItems(self.df.columns.astype(str))
            else:
                self.column_combo.addItems([f"Column {i+1}" for i in range(len(self.df.columns))])
            
            self.column_combo.setEnabled(True)
            self.filtered_df = self.df.copy()
            self.update_results_table()
            self.status_bar.showMessage(f"Loaded {len(self.df)} rows from {self.file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{str(e)}")
    
    def apply_filter(self):
        """Apply the filter to the loaded data"""
        if self.df is None:
            QMessageBox.warning(self, "Warning", "Please load a file first")
            return
        
        if not self.column_combo.currentText():
            QMessageBox.warning(self, "Warning", "Please select a column to filter on")
            return
        
        if not self.value_edit.text():
            QMessageBox.warning(self, "Warning", "Please enter a filter value")
            return
        
        column = self.column_combo.currentText()
        operation = self.operation_combo.currentText()
        value = self.value_edit.text()
        case_sensitive = self.case_check.isChecked()
        invert = self.invert_check.isChecked()
        
        try:
            # Get column name or index
            if self.header_check.isChecked():
                col = column
            else:
                col = int(column.split()[1]) - 1
            
            # Apply filter based on operation
            if operation == "Equals":
                if case_sensitive:
                    mask = self.df[col].astype(str) == value
                else:
                    mask = self.df[col].astype(str).str.lower() == value.lower()
            elif operation == "Contains":
                if case_sensitive:
                    mask = self.df[col].astype(str).str.contains(value, regex=False)
                else:
                    mask = self.df[col].astype(str).str.lower().str.contains(value.lower(), regex=False)
            elif operation == "Starts With":
                if case_sensitive:
                    mask = self.df[col].astype(str).str.startswith(value)
                else:
                    mask = self.df[col].astype(str).str.lower().str.startswith(value.lower())
            elif operation == "Ends With":
                if case_sensitive:
                    mask = self.df[col].astype(str).str.endswith(value)
                else:
                    mask = self.df[col].astype(str).str.lower().str.endswith(value.lower())
            elif operation == "Regex":
                flags = 0 if case_sensitive else re.IGNORECASE
                try:
                    mask = self.df[col].astype(str).str.contains(value, regex=True, flags=flags)
                except re.error:
                    QMessageBox.warning(self, "Warning", "Invalid regular expression pattern")
                    return
            
            # Apply invert if needed
            if invert:
                mask = ~mask
            
            self.filtered_df = self.df[mask]
            self.update_results_table()
            self.status_bar.showMessage(f"Filter applied: {len(self.filtered_df)} rows match criteria")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply filter:\n{str(e)}")
    
    def reset_filter(self):
        """Reset the filter to show all data"""
        if self.df is not None:
            self.filtered_df = self.df.copy()
            self.update_results_table()
            self.status_bar.showMessage("Filter reset - showing all rows")
    
    def update_results_table(self):
        """Update the results table with the current filtered DataFrame"""
        if self.filtered_df is None:
            return
        
        # Limit preview to first 1000 rows
        preview_df = self.filtered_df.head(1000)
        
        self.results_table.setRowCount(preview_df.shape[0])
        self.results_table.setColumnCount(preview_df.shape[1])
        
        # Set headers
        if self.header_check.isChecked() and isinstance(preview_df.columns, pd.core.indexes.base.Index):
            self.results_table.setHorizontalHeaderLabels(preview_df.columns.astype(str))
        else:
            self.results_table.setHorizontalHeaderLabels([f"Column {i+1}" for i in range(preview_df.shape[1])])
        
        # Populate table
        for i in range(preview_df.shape[0]):
            for j in range(preview_df.shape[1]):
                item = QTableWidgetItem(str(preview_df.iloc[i, j]))
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                self.results_table.setItem(i, j, item)
        
        self.results_table.resizeColumnsToContents()
    
    def export_results(self, format_type):
        """Export the filtered results to a file"""
        if self.filtered_df is None or len(self.filtered_df) == 0:
            QMessageBox.warning(self, "Warning", "No data to export")
            return
        
        file_dialog = QFileDialog()
        default_name = self.file_path.split('/')[-1].split('.')[0] + "_filtered"
        
        if format_type == "csv":
            file_path, _ = file_dialog.getSaveFileName(
                self,
                "Export to CSV",
                f"{default_name}.csv",
                "CSV Files (*.csv)"
            )
            if file_path:
                try:
                    self.filtered_df.to_csv(file_path, index=False)
                    QMessageBox.information(self, "Success", f"Data exported to {file_path}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to export:\n{str(e)}")
        
        elif format_type == "tsv":
            file_path, _ = file_dialog.getSaveFileName(
                self,
                "Export to TSV",
                f"{default_name}.tsv",
                "TSV Files (*.tsv)"
            )
            if file_path:
                try:
                    self.filtered_df.to_csv(file_path, sep='\t', index=False)
                    QMessageBox.information(self, "Success", f"Data exported to {file_path}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to export:\n{str(e)}")
        
        elif format_type == "excel":
            file_path, _ = file_dialog.getSaveFileName(
                self,
                "Export to Excel",
                f"{default_name}.xlsx",
                "Excel Files (*.xlsx)"
            )
            if file_path:
                try:
                    self.filtered_df.to_excel(file_path, index=False)
                    QMessageBox.information(self, "Success", f"Data exported to {file_path}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to export:\n{str(e)}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FileFilterApp()
    window.show()
    sys.exit(app.exec_())
