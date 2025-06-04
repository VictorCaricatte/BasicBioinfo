"""
# ==============================================================================
# Author:       Victor S Caricatte De Araújo
# Email:        victorleniwys@gmail.com or victorsc@ufmg.br
# Intitution:   Universidade federal de Minas Gerais
# Version:      1.1.9
# Date:         Abr, 28
# ...................................
# ==============================================================================
"""

import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QLabel, QPushButton, QFileDialog, QComboBox, QSpinBox, 
                             QCheckBox, QTextEdit, QTabWidget, QMessageBox, QGroupBox, QInputDialog, QColorDialog)
from PyQt5.QtCore import Qt
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import seaborn as sns
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import os
from typing import Dict, List, Optional


class SequenceLengthAnalyzer:
    """Core analysis logic (same as before)"""
    def __init__(self):
        self.lengths: Dict[str, List[int]] = defaultdict(list)
        self.stats: Dict[str, Dict[str, float]] = defaultdict(dict)
        
    def read_fasta(self, file_path: str, label: Optional[str] = None) -> None:
        if label is None:
            label = os.path.basename(file_path)
            
        try:
            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    self.lengths[label].append(len(record.seq))
                    
            if not self.lengths[label]:
                raise ValueError(f"No sequences found in {file_path}")
            else:
                self._calculate_stats(label)
                return True
                
        except Exception as e:
            raise ValueError(f"Error reading {file_path}: {str(e)}")
    
    def _calculate_stats(self, label: str) -> None:
        lengths = self.lengths[label]
        self.stats[label] = {
            'count': len(lengths),
            'min': min(lengths),
            'max': max(lengths),
            'mean': np.mean(lengths),
            'median': np.median(lengths),
            'std': np.std(lengths),
            'q1': np.percentile(lengths, 25),
            'q3': np.percentile(lengths, 75),
        }
    
    def get_plot_data(self, plot_type: str = 'histogram') -> Figure:
        if not self.lengths:
            raise ValueError("No data to plot")
            
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sns.set_style("whitegrid")
        
        if plot_type == 'histogram':
            all_lengths = [l for lengths in self.lengths.values() for l in lengths]
            iqr = np.percentile(all_lengths, 75) - np.percentile(all_lengths, 25)
            bin_width = 2 * iqr / (len(all_lengths) ** (1/3))
            bins = max(10, int((max(all_lengths) - min(all_lengths)) / bin_width))
            
            for label, lengths in self.lengths.items():
                sns.histplot(lengths, bins=int(bins), ax=ax,
                             alpha=0.6, 
                             label=f"{label} (n={len(lengths)})",
                             kde=True,
                             edgecolor='none')
            
            ax.set_title("Sequence Length Distribution")
            ax.set_xlabel("Sequence Length (bp)")
            ax.set_ylabel("Frequency")
            ax.legend(title="Dataset")
            
        elif plot_type == 'boxplot':
            data = []
            for label, lengths in self.lengths.items():
                for l in lengths:
                    data.append({'length': l, 'dataset': label})
            
            sns.boxplot(x='dataset', y='length', data=data, ax=ax, showfliers=False)
            ax.set_title("Sequence Length Comparison")
            ax.set_xlabel("Dataset")
            ax.set_ylabel("Sequence Length (bp)")
            ax.tick_params(axis='x', rotation=45)
        
        fig.tight_layout()
        return fig
    
    def get_statistics_text(self) -> str:
        if not self.stats:
            return "No statistics available"
            
        text = []
        for label, stats in self.stats.items():
            text.append(f"Statistics for {label}:")
            text.append("-" * (15 + len(label)))
            for stat, value in stats.items():
                text.append(f"{stat:>8}: {value:12.2f}")
            text.append("")
        
        return "\n".join(text)


class BaseCompositionAnalyzer:
    """Analyze base composition of DNA/RNA sequences"""
    def __init__(self):
        self.sequence = ""
        self.base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}
        self.base_colors = {'A': '#3498db', 'T': '#e74c3c', 'C': '#2ecc71', 
                           'G': '#f39c12', 'U': '#9b59b6', 'N': '#95a5a6'}
    
    def set_sequence(self, sequence: str) -> None:
        """Set the sequence to analyze"""
        self.sequence = sequence.upper()
        self._count_bases()
    
    def _count_bases(self) -> None:
        """Count occurrences of each base"""
        self.base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}
        for base in self.sequence:
            if base in self.base_counts:
                self.base_counts[base] += 1
            else:
                self.base_counts['N'] += 1
    
    def set_base_color(self, base: str, color: str) -> None:
        """Set color for a specific base"""
        if base in self.base_colors:
            self.base_colors[base] = color
    
    def get_plot_data(self) -> Figure:
        """Generate bar plot of base composition"""
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sns.set_style("whitegrid")
        
        # Filter out bases with zero counts
        bases = [b for b in self.base_counts if self.base_counts[b] > 0]
        counts = [self.base_counts[b] for b in bases]
        colors = [self.base_colors[b] for b in bases]
        
        # Calculate frequencies if sequence is not empty
        total = sum(counts) if sum(counts) > 0 else 1
        frequencies = [count/total*100 for count in counts]
        
        # Create bar plot
        bars = ax.bar(bases, frequencies, color=colors)
        
        # Add count labels on top of bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}',
                    ha='center', va='bottom')
        
        ax.set_title("Base Composition")
        ax.set_xlabel("Base")
        ax.set_ylabel("Frequency (%)")
        ax.set_ylim(0, 100)
        
        fig.tight_layout()
        return fig
    
    def get_statistics_text(self) -> str:
        """Return text summary of base composition"""
        if not self.sequence:
            return "No sequence loaded"
            
        total = sum(self.base_counts.values())
        if total == 0:
            return "Empty sequence"
            
        text = []
        text.append("Base Composition Analysis:")
        text.append("-" * 25)
        text.append(f"Sequence length: {len(self.sequence)}")
        text.append(f"Total bases counted: {total}")
        text.append("")
        text.append("Base counts:")
        for base, count in sorted(self.base_counts.items()):
            if count > 0:
                text.append(f"{base}: {count} ({count/total*100:.2f}%)")
        
        return "\n".join(text)


class FastqAnalyzer:
    """Analyze FASTQ files to estimate number of genomes"""
    def __init__(self):
        self.read_counts = {}
        self.genome_size = 0
        self.results = {}
    
    def read_fastq(self, file_path: str, label: Optional[str] = None) -> None:
        """Count reads in a FASTQ file"""
        if label is None:
            label = os.path.basename(file_path)
            
        try:
            count = 0
            with open(file_path, "r") as handle:
                for _ in SeqIO.parse(handle, "fastq"):
                    count += 1
            
            if count == 0:
                raise ValueError(f"No reads found in {file_path}")
            else:
                self.read_counts[label] = count
                return True
                
        except Exception as e:
            raise ValueError(f"Error reading {file_path}: {str(e)}")
    
    def calculate_genomes(self, genome_size: int, read_length: int = 150, paired_end: bool = True) -> None:
        """Calculate number of genomes based on formula"""
        if genome_size <= 0:
            raise ValueError("Genome size must be positive")
            
        self.genome_size = genome_size
        self.results = {}
        
        for label, count in self.read_counts.items():
            multiplier = 2 if paired_end else 1
            coverage = (count * read_length * multiplier) / genome_size
            self.results[label] = {
                'reads': count,
                'genomes': coverage,
                'read_length': read_length,
                'paired_end': paired_end,
                'genome_size': genome_size
            }
    
    def get_results_text(self) -> str:
        """Return formatted results"""
        if not self.results:
            return "No results available"
            
        text = []
        text.append("Genome Coverage Estimation:")
        text.append("-" * 30)
        text.append(f"Genome size used: {self.genome_size:,} bp")
        text.append("")
        
        for label, result in self.results.items():
            text.append(f"Results for {label}:")
            text.append(f"  Number of reads: {result['reads']:,}")
            text.append(f"  Read length: {result['read_length']} bp")
            text.append(f"  Sequencing type: {'Paired-end' if result['paired_end'] else 'Single-end'}")
            text.append(f"  Estimated genome coverage: {result['genomes']:.2f}x")
            text.append("")
        
        return "\n".join(text)


class PlotCanvas(FigureCanvas):
    """Matplotlib canvas for embedding plots in Qt"""
    def __init__(self, parent=None, width=8, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_subplot(111)
        sns.set_style("whitegrid")


class FastaAnalyzerGUI(QMainWindow):
    """Main GUI window"""
    def __init__(self):
        super().__init__()
        self.analyzer = SequenceLengthAnalyzer()
        self.base_analyzer = BaseCompositionAnalyzer()
        self.fastq_analyzer = FastqAnalyzer()
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("FASTA/FASTQ Sequence Analyzer")
        self.setGeometry(100, 100, 1000, 700)
        
        # Central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Left panel - controls
        control_panel = QGroupBox("Controls")
        control_layout = QVBoxLayout()
        
        # Tab widget for different analysis types
        analysis_tabs = QTabWidget()
        
        # Tab 1: Sequence Length Analysis
        length_tab = QWidget()
        length_tab_layout = QVBoxLayout()
        
        # File selection
        file_group = QGroupBox("FASTA Files")
        file_layout = QVBoxLayout()
        
        self.file_list = QComboBox()
        self.file_list.setToolTip("Loaded FASTA files")
        
        self.add_file_btn = QPushButton("Add FASTA File")
        self.add_file_btn.clicked.connect(self.add_file)
        
        self.remove_file_btn = QPushButton("Remove Selected")
        self.remove_file_btn.clicked.connect(self.remove_file)
        
        file_layout.addWidget(QLabel("Loaded Files:"))
        file_layout.addWidget(self.file_list)
        file_layout.addWidget(self.add_file_btn)
        file_layout.addWidget(self.remove_file_btn)
        file_group.setLayout(file_layout)
        
        # Plot options
        plot_group = QGroupBox("Plot Options")
        plot_layout = QVBoxLayout()
        
        self.plot_type = QComboBox()
        self.plot_type.addItems(["Histogram", "Boxplot"])
        
        self.kde_check = QCheckBox("Show KDE (Histogram only)")
        self.kde_check.setChecked(True)
        
        self.bins_spin = QSpinBox()
        self.bins_spin.setRange(5, 500)
        self.bins_spin.setValue(20)
        self.bins_spin.setToolTip("Set to 0 for automatic bin calculation")
        
        plot_layout.addWidget(QLabel("Plot Type:"))
        plot_layout.addWidget(self.plot_type)
        plot_layout.addWidget(self.kde_check)
        plot_layout.addWidget(QLabel("Number of Bins (0=auto):"))
        plot_layout.addWidget(self.bins_spin)
        plot_group.setLayout(plot_layout)
        
        # Actions
        action_group = QGroupBox("Actions")
        action_layout = QVBoxLayout()
        
        self.plot_btn = QPushButton("Generate Plot")
        self.plot_btn.clicked.connect(self.update_plot)
        
        self.save_plot_btn = QPushButton("Save Plot...")
        self.save_plot_btn.clicked.connect(self.save_plot)
        
        self.save_stats_btn = QPushButton("Save Statistics...")
        self.save_stats_btn.clicked.connect(self.save_stats)
        
        action_layout.addWidget(self.plot_btn)
        action_layout.addWidget(self.save_plot_btn)
        action_layout.addWidget(self.save_stats_btn)
        action_group.setLayout(action_layout)
        
        # Assemble length analysis tab
        length_tab_layout.addWidget(file_group)
        length_tab_layout.addWidget(plot_group)
        length_tab_layout.addWidget(action_group)
        length_tab_layout.addStretch()
        length_tab.setLayout(length_tab_layout)
        
        # Tab 2: Base Composition Analysis
        base_tab = QWidget()
        base_tab_layout = QVBoxLayout()
        
        # Sequence input
        seq_group = QGroupBox("Sequence Input")
        seq_layout = QVBoxLayout()
        
        self.seq_input = QTextEdit()
        self.seq_input.setPlaceholderText("Paste DNA/RNA sequence here...")
        
        self.analyze_seq_btn = QPushButton("Analyze Sequence")
        self.analyze_seq_btn.clicked.connect(self.analyze_sequence)
        
        seq_layout.addWidget(QLabel("DNA/RNA Sequence:"))
        seq_layout.addWidget(self.seq_input)
        seq_layout.addWidget(self.analyze_seq_btn)
        seq_group.setLayout(seq_layout)
        
        # Color customization
        color_group = QGroupBox("Base Colors")
        color_layout = QVBoxLayout()
        
        self.color_buttons = {}
        bases = ['A', 'T', 'C', 'G', 'U', 'N']
        for base in bases:
            hbox = QHBoxLayout()
            label = QLabel(f"{base}:")
            btn = QPushButton()
            btn.setFixedSize(30, 30)
            btn.setStyleSheet(f"background-color: {self.base_analyzer.base_colors[base]}")
            btn.clicked.connect(lambda _, b=base: self.change_base_color(b))
            hbox.addWidget(label)
            hbox.addWidget(btn)
            hbox.addStretch()
            color_layout.addLayout(hbox)
            self.color_buttons[base] = btn
        
        color_layout.addStretch()
        color_group.setLayout(color_layout)
        
        # Actions
        base_action_group = QGroupBox("Actions")
        base_action_layout = QVBoxLayout()
        
        self.update_base_plot_btn = QPushButton("Update Plot")
        self.update_base_plot_btn.clicked.connect(self.update_base_plot)
        
        self.save_base_plot_btn = QPushButton("Save Plot...")
        self.save_base_plot_btn.clicked.connect(self.save_base_plot)
        
        base_action_layout.addWidget(self.update_base_plot_btn)
        base_action_layout.addWidget(self.save_base_plot_btn)
        base_action_group.setLayout(base_action_layout)
        
        # Assemble base analysis tab
        base_tab_layout.addWidget(seq_group)
        base_tab_layout.addWidget(color_group)
        base_tab_layout.addWidget(base_action_group)
        base_tab_layout.addStretch()
        base_tab.setLayout(base_tab_layout)
        
        # Tab 3: FASTQ Genome Coverage Analysis
        fastq_tab = QWidget()
        fastq_tab_layout = QVBoxLayout()
        
        # FASTQ file selection
        fastq_file_group = QGroupBox("FASTQ Files")
        fastq_file_layout = QVBoxLayout()
        
        self.fastq_file_list = QComboBox()
        self.fastq_file_list.setToolTip("Loaded FASTQ files")
        
        self.add_fastq_file_btn = QPushButton("Add FASTQ File")
        self.add_fastq_file_btn.clicked.connect(self.add_fastq_file)
        
        self.remove_fastq_file_btn = QPushButton("Remove Selected")
        self.remove_fastq_file_btn.clicked.connect(self.remove_fastq_file)
        
        fastq_file_layout.addWidget(QLabel("Loaded Files:"))
        fastq_file_layout.addWidget(self.fastq_file_list)
        fastq_file_layout.addWidget(self.add_fastq_file_btn)
        fastq_file_layout.addWidget(self.remove_fastq_file_btn)
        fastq_file_group.setLayout(fastq_file_layout)
        
        # Parameters
        param_group = QGroupBox("Parameters")
        param_layout = QVBoxLayout()
        
        # Genome size input
        genome_size_layout = QHBoxLayout()
        genome_size_layout.addWidget(QLabel("Genome Size (bp):"))
        self.genome_size_input = QSpinBox()
        self.genome_size_input.setRange(1, 1000000000)
        self.genome_size_input.setValue(1000000)
        genome_size_layout.addWidget(self.genome_size_input)
        genome_size_layout.addStretch()
        
        # Read length
        read_length_layout = QHBoxLayout()
        read_length_layout.addWidget(QLabel("Read Length (bp):"))
        self.read_length_input = QSpinBox()
        self.read_length_input.setRange(1, 10000)
        self.read_length_input.setValue(150)
        read_length_layout.addWidget(self.read_length_input)
        read_length_layout.addStretch()
        
        # Paired-end checkbox
        self.paired_end_check = QCheckBox("Paired-end sequencing")
        self.paired_end_check.setChecked(True)
        
        param_layout.addLayout(genome_size_layout)
        param_layout.addLayout(read_length_layout)
        param_layout.addWidget(self.paired_end_check)
        param_group.setLayout(param_layout)
        
        # Actions
        fastq_action_group = QGroupBox("Actions")
        fastq_action_layout = QVBoxLayout()
        
        self.calculate_btn = QPushButton("Calculate Genome Coverage")
        self.calculate_btn.clicked.connect(self.calculate_genome_coverage)
        
        self.save_fastq_results_btn = QPushButton("Save Results...")
        self.save_fastq_results_btn.clicked.connect(self.save_fastq_results)
        
        fastq_action_layout.addWidget(self.calculate_btn)
        fastq_action_layout.addWidget(self.save_fastq_results_btn)
        fastq_action_group.setLayout(fastq_action_layout)
        
        # Assemble FASTQ analysis tab
        fastq_tab_layout.addWidget(fastq_file_group)
        fastq_tab_layout.addWidget(param_group)
        fastq_tab_layout.addWidget(fastq_action_group)
        fastq_tab_layout.addStretch()
        fastq_tab.setLayout(fastq_tab_layout)
        
        # Add tabs to analysis tabs widget
        analysis_tabs.addTab(length_tab, "Sequence Length")
        analysis_tabs.addTab(base_tab, "Base Composition")
        analysis_tabs.addTab(fastq_tab, "FASTQ Coverage")
        
        # Add analysis tabs to control panel
        control_layout.addWidget(analysis_tabs)
        control_panel.setLayout(control_layout)
        
        # Right panel - results
        result_panel = QWidget()
        result_layout = QVBoxLayout()
        
        self.tabs = QTabWidget()
        
        # Plot tab
        self.plot_tab = QWidget()
        plot_tab_layout = QVBoxLayout()
        self.canvas = PlotCanvas(self)
        plot_tab_layout.addWidget(self.canvas)
        self.plot_tab.setLayout(plot_tab_layout)
        
        # Statistics tab
        self.stats_tab = QWidget()
        stats_tab_layout = QVBoxLayout()
        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        stats_tab_layout.addWidget(self.stats_text)
        self.stats_tab.setLayout(stats_tab_layout)
        
        # Base composition tab
        self.base_comp_tab = QWidget()
        base_comp_tab_layout = QVBoxLayout()
        self.base_canvas = PlotCanvas(self)
        base_comp_tab_layout.addWidget(self.base_canvas)
        
        self.base_stats_text = QTextEdit()
        self.base_stats_text.setReadOnly(True)
        base_comp_tab_layout.addWidget(self.base_stats_text)
        
        self.base_comp_tab.setLayout(base_comp_tab_layout)
        
        # FASTQ results tab
        self.fastq_results_tab = QWidget()
        fastq_results_tab_layout = QVBoxLayout()
        self.fastq_results_text = QTextEdit()
        self.fastq_results_text.setReadOnly(True)
        fastq_results_tab_layout.addWidget(self.fastq_results_text)
        self.fastq_results_tab.setLayout(fastq_results_tab_layout)
        
        # Help tab
        self.help_tab = QWidget()
        help_tab_layout = QVBoxLayout()
        help_text = QTextEdit()
        help_text.setReadOnly(True)
        
        # Documentation content
        documentation = """
        <h1>FASTA/FASTQ Sequence Analyzer</h1>
        
        <h2>About</h2>
        <p>This application analyzes sequence length distributions in FASTA files, base composition, 
        and estimates genome coverage from FASTQ files.</p>
        
        <h2>Features</h2>
        <ul>
            <li>Load multiple FASTA files for comparison</li>
            <li>Generate histogram or boxplot visualizations</li>
            <li>View detailed statistics for each dataset</li>
            <li>Save plots and statistics to files</li>
            <li>Analyze base composition of DNA/RNA sequences</li>
            <li>Customize colors for base composition plots</li>
            <li>Estimate genome coverage from FASTQ files</li>
            <li>Command-line interface available</li>
        </ul>
        
        <h2>Usage</h2>
        <h3>Sequence Length Analysis</h3>
        <ol>
            <li>Click "Add FASTA File" to load sequence files</li>
            <li>Optionally provide custom labels for each dataset</li>
            <li>Select plot type and options</li>
            <li>Click "Generate Plot" to visualize the data</li>
            <li>Use the tabs to switch between visualization and statistics</li>
            <li>Save results using the appropriate buttons</li>
        </ol>
        
        <h3>Base Composition Analysis</h3>
        <ol>
            <li>Paste a DNA/RNA sequence in the input box</li>
            <li>Click "Analyze Sequence" to calculate base composition</li>
            <li>Optionally customize base colors by clicking the color buttons</li>
            <li>View results in the Base Composition tab</li>
        </ol>
        
        <h3>FASTQ Genome Coverage Analysis</h3>
        <ol>
            <li>Click "Add FASTQ File" to load sequencing files</li>
            <li>Enter the genome size in base pairs</li>
            <li>Specify read length (default 150bp)</li>
            <li>Check "Paired-end sequencing" if applicable</li>
            <li>Click "Calculate Genome Coverage" to run the analysis</li>
            <li>View results in the FASTQ Results tab</li>
        </ol>
        
        <h3>Genome Coverage Formula</h3>
        <p>The genome coverage is calculated using the formula:</p>
        <p>Coverage = (Number of reads × Read length × 2 if paired-end) / Genome size</p>
        <p>This gives the estimated coverage in "X" (times coverage).</p>
        
        <h2>Command Line Usage</h2>
        <p>The application can also be run from command line with arguments:</p>
        <pre>
        python script.py file1.fasta file2.fasta [options]
        
        Options:
        -l, --label      Custom labels for files
        --histogram      Output file for histogram plot
        --boxplot        Output file for boxplot
        --stats          Output file for statistics
        </pre>
        
        <h2>Technical Details</h2>
        <p>This application uses:</p>
        <ul>
            <li>Python 3 and PyQt5 for the GUI</li>
            <li>Matplotlib and Seaborn for visualization</li>
            <li>Biopython for FASTA/FASTQ parsing</li>
            <li>NumPy for statistical calculations</li>
        </ul>
        
        <h2>Author</h2>
        <p>Victor Silveira Caricatte/Universidade Federal de Minas Gerais.</p>
        
        <h2>Contact</h2>
        <p>https://github.com/VictorCaricatte or victorsc@ufmg.br.</p>
        """
        
        help_text.setHtml(documentation)
        help_tab_layout.addWidget(help_text)
        self.help_tab.setLayout(help_tab_layout)
        
        self.tabs.addTab(self.plot_tab, "Length Visualization")
        self.tabs.addTab(self.stats_tab, "Length Statistics")
        self.tabs.addTab(self.base_comp_tab, "Base Composition")
        self.tabs.addTab(self.fastq_results_tab, "FASTQ Results")
        self.tabs.addTab(self.help_tab, "Help")
        
        result_layout.addWidget(self.tabs)
        result_panel.setLayout(result_layout)
        
        # Add panels to main layout
        main_layout.addWidget(control_panel, stretch=1)
        main_layout.addWidget(result_panel, stretch=3)
        
        # Status bar
        self.statusBar().showMessage("Ready")
        
    def add_file(self):
        """Add a FASTA file to analyze"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open FASTA File", "", "FASTA Files (*.fasta *.fa *.fna);;All Files (*)")
        
        if file_path:
            label, ok = QInputDialog.getText(
                self, "Dataset Label", 
                "Enter a label for this dataset (or leave blank to use filename):",
                text=os.path.basename(file_path))
            
            if ok:
                try:
                    success = self.analyzer.read_fasta(file_path, label if label else None)
                    if success:
                        self.file_list.addItem(f"{label if label else os.path.basename(file_path)} ({file_path})", 
                                             userData=file_path)
                        self.update_stats()
                        self.statusBar().showMessage(f"Loaded {file_path}", 3000)
                except ValueError as e:
                    QMessageBox.warning(self, "Error", str(e))
    
    def remove_file(self):
        """Remove selected FASTA file from analysis"""
        if self.file_list.count() > 0 and self.file_list.currentIndex() >= 0:
            index = self.file_list.currentIndex()
            label = self.file_list.itemText(index).split(" (")[0]
            self.file_list.removeItem(index)
            
            if label in self.analyzer.lengths:
                del self.analyzer.lengths[label]
                if label in self.analyzer.stats:
                    del self.analyzer.stats[label]
            
            self.update_stats()
            if self.file_list.count() > 0:
                self.update_plot()
            else:
                self.canvas.ax.clear()
                self.canvas.draw()
    
    def update_plot(self):
        """Update the plot with current data and settings"""
        if not self.analyzer.lengths:
            QMessageBox.information(self, "No Data", "Please load FASTA files first")
            return
            
        try:
            plot_type = self.plot_type.currentText().lower()
            
            # Clear previous plot
            self.canvas.fig.clf()
            self.canvas.ax = self.canvas.fig.add_subplot(111)
            
            if plot_type == 'histogram':
                bins = self.bins_spin.value() if self.bins_spin.value() > 0 else 'auto'
                show_kde = self.kde_check.isChecked()
                
                for label, lengths in self.analyzer.lengths.items():
                    sns.histplot(lengths, bins=bins, ax=self.canvas.ax,
                                alpha=0.6, 
                                label=f"{label} (n={len(lengths)})",
                                kde=show_kde,
                                edgecolor='none')
                
                self.canvas.ax.set_title("Sequence Length Distribution")
                self.canvas.ax.set_xlabel("Sequence Length (bp)")
                self.canvas.ax.set_ylabel("Frequency")
                self.canvas.ax.legend(title="Dataset")
                
            elif plot_type == 'boxplot':
                # Create a DataFrame from the data
                import pandas as pd
                data = []
                for label, lengths in self.analyzer.lengths.items():
                    for l in lengths:
                        data.append({'length': l, 'dataset': label})
                
                df = pd.DataFrame(data)
                sns.boxplot(x='dataset', y='length', data=df, ax=self.canvas.ax, showfliers=False)
                
                self.canvas.ax.set_title("Sequence Length Comparison")
                self.canvas.ax.set_xlabel("Dataset")
                self.canvas.ax.set_ylabel("Sequence Length (bp)")
                self.canvas.ax.tick_params(axis='x', rotation=45)
            
            self.canvas.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            QMessageBox.warning(self, "Plot Error", f"Could not generate plot: {str(e)}")
    
    def update_stats(self):
        """Update statistics display"""
        self.stats_text.setPlainText(self.analyzer.get_statistics_text())
    
    def save_plot(self):
        """Save current plot to file"""
        if not self.analyzer.lengths:
            QMessageBox.information(self, "No Data", "Nothing to save - no data loaded")
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "", 
            "PNG Image (*.png);;JPEG Image (*.jpg);;PDF Document (*.pdf);;SVG Image (*.svg)")
        
        if file_path:
            try:
                self.canvas.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                self.statusBar().showMessage(f"Plot saved to {file_path}", 5000)
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Could not save plot: {str(e)}")
    
    def save_stats(self):
        """Save statistics to text file"""
        if not self.analyzer.stats:
            QMessageBox.information(self, "No Data", "No statistics to save")
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Statistics", "", "Text Files (*.txt);;All Files (*)")
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write(self.analyzer.get_statistics_text())
                self.statusBar().showMessage(f"Statistics saved to {file_path}", 5000)
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Could not save statistics: {str(e)}")
    
    def analyze_sequence(self):
        """Analyze the base composition of the input sequence"""
        sequence = self.seq_input.toPlainText().strip()
        if not sequence:
            QMessageBox.warning(self, "Input Error", "Please enter a DNA/RNA sequence")
            return
            
        # Remove whitespace and numbers from sequence
        sequence = ''.join([c for c in sequence if c.isalpha()])
        self.base_analyzer.set_sequence(sequence)
        self.update_base_plot()
        self.update_base_stats()
        self.statusBar().showMessage("Sequence analyzed", 3000)
    
    def update_base_plot(self):
        """Update the base composition plot"""
        if not self.base_analyzer.sequence:
            return
            
        try:
            # Clear previous plot
            self.base_canvas.fig.clf()
            ax = self.base_canvas.fig.add_subplot(111)
            
            # Get base data
            bases = [b for b in self.base_analyzer.base_counts if self.base_analyzer.base_counts[b] > 0]
            counts = [self.base_analyzer.base_counts[b] for b in bases]
            colors = [self.base_analyzer.base_colors[b] for b in bases]
            
            # Calculate frequencies
            total = sum(counts) if sum(counts) > 0 else 1
            frequencies = [count/total*100 for count in counts]
            
            # Create bar plot directly on our canvas
            bars = ax.bar(bases, frequencies, color=colors)
            
            # Add count labels
            for bar, count in zip(bars, counts):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                        f'{count}',
                        ha='center', va='bottom')
            
            ax.set_title("Base Composition")
            ax.set_xlabel("Base")
            ax.set_ylabel("Frequency (%)")
            ax.set_ylim(0, 100)
            
            self.base_canvas.fig.tight_layout()
            self.base_canvas.draw()
            
        except Exception as e:
            QMessageBox.warning(self, "Plot Error", f"Could not generate base composition plot: {str(e)}")
    
    def update_base_stats(self):
        """Update base composition statistics display"""
        self.base_stats_text.setPlainText(self.base_analyzer.get_statistics_text())
    
    def save_base_plot(self):
        """Save base composition plot to file"""
        if not self.base_analyzer.sequence:
            QMessageBox.information(self, "No Data", "Nothing to save - no sequence analyzed")
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "", 
            "PNG Image (*.png);;JPEG Image (*.jpg);;PDF Document (*.pdf);;SVG Image (*.svg)")
        
        if file_path:
            try:
                fig = self.base_analyzer.get_plot_data()
                fig.savefig(file_path, dpi=300, bbox_inches='tight')
                self.statusBar().showMessage(f"Plot saved to {file_path}", 5000)
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Could not save plot: {str(e)}")
    
    def change_base_color(self, base: str):
        """Change the color for a specific base"""
        color = QColorDialog.getColor()
        if color.isValid():
            hex_color = color.name()
            self.base_analyzer.set_base_color(base, hex_color)
            self.color_buttons[base].setStyleSheet(f"background-color: {hex_color}")
            self.update_base_plot()
    
    def add_fastq_file(self):
        """Add a FASTQ file to analyze"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open FASTQ File", "", "FASTQ Files (*.fastq *.fq);;All Files (*)")
        
        if file_path:
            label, ok = QInputDialog.getText(
                self, "Dataset Label", 
                "Enter a label for this dataset (or leave blank to use filename):",
                text=os.path.basename(file_path))
            
            if ok:
                try:
                    success = self.fastq_analyzer.read_fastq(file_path, label if label else None)
                    if success:
                        self.fastq_file_list.addItem(f"{label if label else os.path.basename(file_path)} ({file_path})", 
                                                   userData=file_path)
                        self.statusBar().showMessage(f"Loaded {file_path}", 3000)
                except ValueError as e:
                    QMessageBox.warning(self, "Error", str(e))
    
    def remove_fastq_file(self):
        """Remove selected FASTQ file from analysis"""
        if self.fastq_file_list.count() > 0 and self.fastq_file_list.currentIndex() >= 0:
            index = self.fastq_file_list.currentIndex()
            label = self.fastq_file_list.itemText(index).split(" (")[0]
            self.fastq_file_list.removeItem(index)
            
            if label in self.fastq_analyzer.read_counts:
                del self.fastq_analyzer.read_counts[label]
    
    def calculate_genome_coverage(self):
        """Calculate genome coverage from FASTQ data"""
        if not self.fastq_analyzer.read_counts:
            QMessageBox.information(self, "No Data", "Please load FASTQ files first")
            return
            
        genome_size = self.genome_size_input.value()
        read_length = self.read_length_input.value()
        paired_end = self.paired_end_check.isChecked()
        
        try:
            self.fastq_analyzer.calculate_genomes(genome_size, read_length, paired_end)
            self.fastq_results_text.setPlainText(self.fastq_analyzer.get_results_text())
            self.statusBar().showMessage("Genome coverage calculated", 3000)
        except ValueError as e:
            QMessageBox.warning(self, "Error", str(e))
    
    def save_fastq_results(self):
        """Save FASTQ analysis results to file"""
        if not self.fastq_analyzer.results:
            QMessageBox.information(self, "No Data", "No results to save")
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "", "Text Files (*.txt);;All Files (*)")
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write(self.fastq_analyzer.get_results_text())
                self.statusBar().showMessage(f"Results saved to {file_path}", 5000)
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Could not save results: {str(e)}")


def main():
    # If command line arguments are present, use CLI mode
    if len(sys.argv) > 1:
        import argparse
        parser = argparse.ArgumentParser(description="FASTA Sequence Length Analyzer")
        parser.add_argument("files", nargs="+", help="FASTA file(s) to analyze")
        parser.add_argument("-l", "--label", action="append", help="Custom labels for files")
        parser.add_argument("--histogram", help="Output file for histogram plot")
        parser.add_argument("--boxplot", help="Output file for boxplot")
        parser.add_argument("--stats", help="Output file for statistics")
        args = parser.parse_args()
        
        analyzer = SequenceLengthAnalyzer()
        
        labels = args.label if args.label else [None] * len(args.files)
        if len(labels) != len(args.files):
            print("Warning: Number of labels doesn't match number of files. Using filenames as labels.")
            labels = [None] * len(args.files)
        
        for file, label in zip(args.files, labels):
            try:
                analyzer.read_fasta(file, label)
            except ValueError as e:
                print(f"Error: {str(e)}")
                continue
        
        if analyzer.lengths:
            if args.stats:
                with open(args.stats, 'w') as f:
                    f.write(analyzer.get_statistics_text())
                print(f"Statistics saved to {args.stats}")
            
            if args.histogram:
                try:
                    fig = analyzer.get_plot_data('histogram')
                    fig.savefig(args.histogram, dpi=300, bbox_inches='tight')
                    print(f"Histogram saved to {args.histogram}")
                except Exception as e:
                    print(f"Error saving histogram: {str(e)}")
            
            if args.boxplot:
                try:
                    fig = analyzer.get_plot_data('boxplot')
                    fig.savefig(args.boxplot, dpi=300, bbox_inches='tight')
                    print(f"Boxplot saved to {args.boxplot}")
                except Exception as e:
                    print(f"Error saving boxplot: {str(e)}")
            
            if not any([args.histogram, args.boxplot, args.stats]):
                print(analyzer.get_statistics_text())
        else:
            print("No valid sequence data found in input files.")
    else:
        # Otherwise launch GUI
        app = QApplication(sys.argv)
        app.setStyle('Fusion')
        
        window = FastaAnalyzerGUI()
        window.show()
        sys.exit(app.exec_())


if __name__ == "__main__":
    main()
  
