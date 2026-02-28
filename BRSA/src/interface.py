import sys
import os
import tempfile
import shutil
import json
import pandas as pd
from collections import defaultdict
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QTabWidget, QLabel, QPushButton, QFileDialog, QTextEdit, 
                             QComboBox, QProgressBar, QMessageBox, QCheckBox, QSpinBox,
                             QGroupBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QListWidget, QStackedWidget,
                             QDockWidget, QTableWidget, QTableWidgetItem, QSplitter, QSizePolicy)
from PyQt6.QtCore import Qt, QObject, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib import colors

from job import FastqDevourerThread as FastqProcessingThread
from job import FrankensteinBioinformaticsOverlordThread as RNASeqAnalysisThread
from plot import (paint_dimensional_void_masterpiece, erupt_volcano_plot, summon_clustermap_demon, 
                  paint_gsea_mountain, forge_interactive_plotly_artifacts, paint_ma_constellation, 
                  plot_dendrogram_of_madness, predict_mortal_doom_curves, 
                  interrogate_feature_importance_spirits, paint_hairball_of_chaos, uncover_hidden_isoform_mutations)
from args import SacredSessionScroll as AnalysisSession
from biotools import scribe_reproducibility_chronicles

# Advanced Halloween / Necromantic Console Theme
PROFESSIONAL_STYLE = """
QMainWindow { background-color: #1a1b26; }
QWidget { color: #a9b1d6; font-family: 'Segoe UI', Helvetica, sans-serif; font-size: 13px; }

QGroupBox { 
    border: 1px solid #292e42; 
    border-radius: 6px; 
    margin-top: 24px;
    padding-top: 15px;
    padding-bottom: 10px;
    background-color: #1f2335; 
}
QGroupBox::title { 
    subcontrol-origin: margin;
    subcontrol-position: top left; 
    left: 10px; 
    top: 0px;
    padding: 0px 5px; 
    background-color: transparent;
    color: #ff8c00; 
    font-weight: bold;
}

QPushButton { 
    background-color: #292e42; 
    border: 1px solid #414868; 
    border-radius: 4px; 
    padding: 8px 16px; 
    color: #c0caf5; 
    font-weight: bold; 
}
QPushButton:hover { background-color: #414868; border-color: #ff8c00; color: #ffffff; }
QPushButton:pressed { background-color: #ff8c00; color: #1a1b26; }
QPushButton#btnPrimary { background-color: #ff4500; border-color: #ff8c00; color: #ffffff; }
QPushButton#btnPrimary:hover { background-color: #ff8c00; color: #1a1b26; }
QPushButton#btnWarning { background-color: #8b0000; border-color: #ff0000; color: #ffffff; }
QPushButton#btnWarning:hover { background-color: #ff0000; color: #1a1b26; }

QComboBox, QSpinBox, QDoubleSpinBox { 
    background-color: #16161e; 
    border: 1px solid #414868; 
    border-radius: 4px; 
    padding: 6px; 
    color: #c0caf5; 
    min-width: 120px;
}
QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover { border-color: #ff8c00; }

QTextEdit, QTableWidget, QListWidget { 
    background-color: #16161e; 
    border: 1px solid #292e42; 
    border-radius: 4px; 
    padding: 5px; 
    gridline-color: #292e42; 
}
QTextEdit#halloweenTerminal {
    background-color: #0d001a;
    color: #39ff14;
    font-family: 'Consolas', 'Courier New', monospace;
    border: 1px solid #ff8c00;
}
QHeaderView::section { 
    background-color: #1f2335; 
    color: #ff8c00; 
    padding: 6px; 
    border: 1px solid #292e42; 
    font-weight: bold;
}

QListWidget::item { padding: 14px; border-bottom: 1px solid #1a1b26; }
QListWidget::item:selected { 
    background-color: #1f2335; 
    border-left: 4px solid #ff8c00; 
    color: #ffffff; 
    font-weight: bold;
}

QProgressBar { 
    border: 1px solid #414868; 
    border-radius: 4px; 
    text-align: center; 
    background-color: #16161e; 
    color: #ffffff;
    font-weight: bold;
}
QProgressBar::chunk { background-color: #8a2be2; border-radius: 3px; }

QDockWidget { font-weight: bold; color: #ff8c00; }
QDockWidget::title { 
    background: #1f2335; 
    padding: 10px; 
    border: 1px solid #ff8c00; 
}

QTabWidget::pane { 
    border: 1px solid #292e42; 
    border-radius: 4px; 
    background-color: #1f2335; 
}
QTabBar::tab { 
    background-color: #16161e; 
    color: #565f89;
    padding: 10px 20px; 
    margin-right: 2px; 
    border-top-left-radius: 4px; 
    border-top-right-radius: 4px; 
    border: 1px solid #292e42;
    border-bottom: none;
    font-weight: bold;
}
QTabBar::tab:selected { 
    background-color: #1f2335; 
    color: #ff8c00; 
}
QSplitter::handle { background-color: #ff8c00; width: 2px; }
"""

class SpectralWhisperCatcher(QObject):
    whisper_heard = pyqtSignal(str)
    def write(self, text):
        if text.strip(): self.whisper_heard.emit(text.strip())
    def flush(self): pass

class FastqProcessingDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("FASTQ Dark Pre-processing Pipeline")
        self.setMinimumWidth(750)
        self.setStyleSheet(PROFESSIONAL_STYLE)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(20)
        
        grp_input = QGroupBox("Raw Sequencing Sacrifices (FASTQ)")
        lay_input = QVBoxLayout(grp_input)
        lay_input.setContentsMargins(15, 25, 15, 15)
        self.txt_box = QTextEdit()
        self.txt_box.setReadOnly(True)
        lay_input.addWidget(self.txt_box)
        
        btn_add = QPushButton("Select FASTQ Files")
        btn_add.setFixedWidth(200)
        lay_input.addWidget(btn_add, alignment=Qt.AlignmentFlag.AlignRight)
        btn_add.clicked.connect(self.capture_files)
        layout.addWidget(grp_input)
        
        grp_tools = QGroupBox("Ritual Parameters")
        lay_tools = QVBoxLayout(grp_tools)
        lay_tools.setContentsMargins(15, 25, 15, 15)
        lay_tools.setSpacing(15)
        
        lay_aligner = QHBoxLayout()
        self.cmb_aligner = QComboBox()
        self.cmb_aligner.addItems(["HISAT2 (Alignment-based)", "Kallisto (Pseudo-alignment)"])
        lay_aligner.addWidget(QLabel("Quantification Demon:"))
        lay_aligner.addWidget(self.cmb_aligner)
        lay_tools.addLayout(lay_aligner)
        
        lay_ref = QHBoxLayout()
        self.lbl_ref = QLabel("Reference Index: Not Selected")
        btn_ref = QPushButton("Browse Index")
        btn_ref.clicked.connect(self.find_reference)
        lay_ref.addWidget(self.lbl_ref)
        lay_ref.addWidget(btn_ref)
        lay_tools.addLayout(lay_ref)
        
        lay_ann = QHBoxLayout()
        self.lbl_ann = QLabel("Annotation File: Not Selected")
        btn_ann = QPushButton("Browse GTF/GFF3")
        btn_ann.clicked.connect(self.find_annotation)
        lay_ann.addWidget(self.lbl_ann)
        lay_ann.addWidget(btn_ann)
        lay_tools.addLayout(lay_ann)
        
        self.chk_gff3 = QCheckBox("Input is GFF3 format (Uncheck for GTF)")
        self.chk_trim = QCheckBox("Summon fastp guillotine (Quality trim)")
        self.chk_trim.setChecked(True)
        lay_tools.addWidget(self.chk_gff3)
        lay_tools.addWidget(self.chk_trim)
        
        lay_out = QHBoxLayout()
        self.lbl_out = QLabel("Output Graveyard: Not Selected")
        btn_out = QPushButton("Select Destination")
        btn_out.clicked.connect(self.find_output)
        lay_out.addWidget(self.lbl_out)
        lay_out.addWidget(btn_out)
        lay_tools.addLayout(lay_out)
        
        layout.addWidget(grp_tools)
        
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
        self.arsenal = defaultdict(list)
        self.ref_index = ""
        self.ann_file = ""
        self.out_dir = ""

    def capture_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select FASTQs", "", "FASTQ (*.fastq *.fq *.fastq.gz *.fq.gz)")
        for f in files:
            bn = os.path.basename(f)
            sample = bn.split('_')[0]
            if '_R1' in bn or '_1' in bn: sample = bn.split('_R1')[0] if '_R1' in bn else bn.split('_1')[0]
            elif '_R2' in bn or '_2' in bn: sample = bn.split('_R2')[0] if '_R2' in bn else bn.split('_2')[0]
            self.arsenal[sample].append(f)
            
        text = "Sample ID\tAssociated Files\n"
        for s, fs in self.arsenal.items(): text += f"{s}\t{', '.join([os.path.basename(x) for x in fs])}\n"
        self.txt_box.setPlainText(text)

    def find_reference(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Index", "", "Index Files (*.ht2 *.ht2l *.idx)")
        if f:
            if f.endswith(tuple(f".{i}.ht2" for i in range(1,9))): f = f[:f.rfind(".", 0, f.rfind("."))]
            self.ref_index = f
            self.lbl_ref.setText(f"Index: {os.path.basename(f)}")

    def find_annotation(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Annotation", "", "Annotations (*.gtf *.gff3)")
        if f:
            self.ann_file = f
            self.lbl_ann.setText(f"Annotation: {os.path.basename(f)}")
            if f.endswith('.gff3'): self.chk_gff3.setChecked(True)

    def find_output(self):
        d = QFileDialog.getExistingDirectory(self, "Select Destination")
        if d:
            self.out_dir = d
            self.lbl_out.setText(f"Output: {d}")

class RNASeqDashboard(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BRSA - Bioinformatics RNA-Seq Necromancy Dashboard")
        self.setGeometry(50, 50, 1600, 950)
        self.setStyleSheet(PROFESSIONAL_STYLE)
        
        self.session_scroll = AnalysisSession()
        self.temp_dir = tempfile.mkdtemp()
        
        self.build_menus()
        self.build_docks()
        self.hijack_terminal()
        
        main_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.setCentralWidget(main_splitter)
        
        self.sidebar = QListWidget()
        self.sidebar.setMaximumWidth(280)
        self.sidebar.addItem("📂 Altar of Data Input")
        self.sidebar.addItem("🔍 Dataset Inspector")
        self.sidebar.addItem("📈 Dark Statistics & Modeling")
        self.sidebar.addItem("📊 Advanced Visualization")
        self.sidebar.addItem("🧬 Functional Enrichment (GSEA)")
        self.sidebar.addItem("📖 Grimoire (Help)")
        self.sidebar.currentRowChanged.connect(self.switch_workspaces)
        main_splitter.addWidget(self.sidebar)
        
        self.workspace_stack = QStackedWidget()
        main_splitter.addWidget(self.workspace_stack)
        main_splitter.setSizes([280, 1320])
        
        self.build_setup_workspace()
        self.build_data_inspector_workspace()
        self.build_stats_workspace()
        self.build_viz_workspace()
        self.build_gsea_workspace()
        self.build_help_workspace()
        
        self.devourer_thread = None
        self.cruncher_thread = None
        self.sidebar.setCurrentRow(0)

    def build_menus(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        
        act_pdf = file_menu.addAction("Export Report (PDF - ABNT standard)")
        act_pdf.triggered.connect(self.export_pdf_report)
        act_html = file_menu.addAction("Export Interactive Plots (HTML)")
        act_html.triggered.connect(self.export_interactive_html)
        act_repro = file_menu.addAction("Generate Reproducibility Log")
        act_repro.triggered.connect(self.export_reproducibility_log)
        file_menu.addSeparator()
        act_exit = file_menu.addAction("Flee the Dashboard")
        act_exit.triggered.connect(self.close)

    def build_docks(self):
        self.dock_logs = QDockWidget("Terminal of Spells", self)
        self.dock_logs.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea | Qt.DockWidgetArea.BottomDockWidgetArea)
        self.txt_logs = QTextEdit()
        self.txt_logs.setObjectName("halloweenTerminal")
        self.txt_logs.setReadOnly(True)
        self.txt_logs.append("[SYSTEM] Necromantic application initialized.")
        self.dock_logs.setWidget(self.txt_logs)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.dock_logs)
        
        self.dock_status = QDockWidget("Analysis State", self)
        self.dock_status.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea)
        status_widget = QWidget()
        status_lay = QVBoxLayout(status_widget)
        status_lay.setContentsMargins(15, 15, 15, 15)
        self.lbl_quick_stat = QLabel("Genes Quantified: 0\nSignificant (FDR): 0\nStatus: Awaiting Sacrifice")
        self.prog_bar = QProgressBar()
        status_lay.addWidget(self.lbl_quick_stat)
        status_lay.addWidget(self.prog_bar)
        status_lay.addStretch()
        self.dock_status.setWidget(status_widget)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.dock_status)

    def hijack_terminal(self):
        self.spectral_catcher = SpectralWhisperCatcher()
        self.spectral_catcher.whisper_heard.connect(self.log_event)
        sys.stdout = self.spectral_catcher
        sys.stderr = self.spectral_catcher

    def log_event(self, msg):
        self.txt_logs.append(f"👻 {msg}")
        self.txt_logs.verticalScrollBar().setValue(self.txt_logs.verticalScrollBar().maximum())

    def switch_workspaces(self, index):
        self.workspace_stack.setCurrentIndex(index)

    def build_setup_workspace(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(20, 20, 20, 20)
        
        tabs = QTabWidget()
        
        tab_files = QWidget()
        lay_files = QVBoxLayout(tab_files)
        lay_files.setContentsMargins(20, 20, 20, 20)
        
        grp_files = QGroupBox("Counts & Metadata Configuration")
        lay_grp = QVBoxLayout(grp_files)
        lay_grp.setContentsMargins(20, 30, 20, 20)
        lay_grp.setSpacing(15)
        
        lay_counts = QHBoxLayout()
        self.lbl_counts = QLabel("Expression Matrix: Not Loaded")
        btn_counts = QPushButton("Load Expression Matrix (CSV)")
        btn_counts.setFixedWidth(250)
        btn_counts.clicked.connect(self.inject_matrix)
        lay_counts.addWidget(self.lbl_counts)
        lay_counts.addWidget(btn_counts)
        
        lay_meta = QHBoxLayout()
        self.lbl_meta = QLabel("Metadata Profile: Not Loaded")
        btn_meta = QPushButton("Load Metadata Profile (CSV)")
        btn_meta.setFixedWidth(250)
        btn_meta.clicked.connect(self.inject_metadata)
        lay_meta.addWidget(self.lbl_meta)
        lay_meta.addWidget(btn_meta)
        
        lay_grp.addLayout(lay_counts)
        lay_grp.addLayout(lay_meta)
        lay_files.addWidget(grp_files, 0)
        lay_files.addStretch(1)
        tabs.addTab(tab_files, "Matrix Input")
        
        tab_raw = QWidget()
        lay_raw = QVBoxLayout(tab_raw)
        lay_raw.setContentsMargins(20, 20, 20, 20)
        
        grp_raw = QGroupBox("Automated Sequence Processing (Includes FastQC & MultiQC)")
        lay_grp_raw = QVBoxLayout(grp_raw)
        lay_grp_raw.setContentsMargins(20, 30, 20, 20)
        
        btn_fastq = QPushButton("Configure FASTQ Pre-processing Ritual")
        btn_fastq.setObjectName("btnWarning")
        btn_fastq.setFixedWidth(300)
        btn_fastq.clicked.connect(self.turn_on_processing)
        
        lay_grp_raw.addWidget(QLabel("Initialize trimming and quantification from raw FASTQ files:"))
        lay_grp_raw.addWidget(btn_fastq, alignment=Qt.AlignmentFlag.AlignLeft)
        lay_raw.addWidget(grp_raw, 0)
        lay_raw.addStretch(1)
        tabs.addTab(tab_raw, "FASTQ Pipeline")
        
        layout.addWidget(tabs)
        self.workspace_stack.addWidget(page)

    def build_data_inspector_workspace(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(20, 20, 20, 20)
        
        tabs = QTabWidget()
        
        self.tab_meta_view = QWidget()
        lay_meta = QVBoxLayout(self.tab_meta_view)
        lay_meta.setContentsMargins(15, 15, 15, 15)
        self.table_meta = QTableWidget()
        lay_meta.addWidget(self.table_meta)
        tabs.addTab(self.tab_meta_view, "Metadata Overview")
        
        self.tab_counts_view = QWidget()
        lay_counts = QVBoxLayout(self.tab_counts_view)
        lay_counts.setContentsMargins(15, 15, 15, 15)
        self.table_counts = QTableWidget()
        lay_counts.addWidget(QLabel("Note: This displays the matrix post-normalization."))
        lay_counts.addWidget(self.table_counts)
        tabs.addTab(self.tab_counts_view, "Normalized Expression Matrix")
        
        layout.addWidget(tabs)
        self.workspace_stack.addWidget(page)

    def populate_table(self, table_widget, dataframe):
        if dataframe is None or dataframe.empty: return
        table_widget.setRowCount(dataframe.shape[0])
        table_widget.setColumnCount(dataframe.shape[1] + 1)
        
        headers = ["Feature/Sample"] + list(dataframe.columns)
        table_widget.setHorizontalHeaderLabels(headers)
        
        for i in range(dataframe.shape[0]):
            table_widget.setItem(i, 0, QTableWidgetItem(str(dataframe.index[i])))
            for j in range(dataframe.shape[1]):
                val = dataframe.iloc[i, j]
                if isinstance(val, float): val = f"{val:.4f}"
                table_widget.setItem(i, j+1, QTableWidgetItem(str(val)))
        table_widget.resizeColumnsToContents()

    def build_stats_workspace(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(20, 20, 20, 20)
        
        grp_params = QGroupBox("PyDESeq2, Batch Exorcism & Advanced Witchcraft")
        lay_params = QVBoxLayout(grp_params)
        lay_params.setContentsMargins(20, 30, 20, 20)
        lay_params.setSpacing(15)
        
        lay_cond = QHBoxLayout()
        self.cmb_cond1 = QComboBox()
        self.cmb_cond2 = QComboBox()
        lay_cond.addWidget(QLabel("Test Condition (Group A):"))
        lay_cond.addWidget(self.cmb_cond1)
        lay_cond.addWidget(QLabel("Reference Condition (Group B):"))
        lay_cond.addWidget(self.cmb_cond2)
        lay_params.addLayout(lay_cond)
        
        lay_norm = QHBoxLayout()
        self.cmb_norm = QComboBox()
        self.cmb_norm.addItems(["PyDESeq2 (Recommended)", "CPM", "TPM"])
        self.cmb_batch = QComboBox()
        self.cmb_batch.addItem("None")
        lay_norm.addWidget(QLabel("Normalization/Stat Engine:"))
        lay_norm.addWidget(self.cmb_norm)
        lay_norm.addWidget(QLabel("PyComBat Covariate:"))
        lay_norm.addWidget(self.cmb_batch)
        lay_params.addLayout(lay_norm)
        
        lay_filt = QHBoxLayout()
        self.spn_min = QSpinBox()
        self.spn_min.setValue(10)
        self.spn_pval = QDoubleSpinBox()
        self.spn_pval.setValue(0.05)
        self.spn_pval.setSingleStep(0.01)
        self.cmb_dimred = QComboBox()
        self.cmb_dimred.addItems(["PCA", "t-SNE", "UMAP"])
        lay_filt.addWidget(QLabel("Min Read Filter:"))
        lay_filt.addWidget(self.spn_min)
        lay_filt.addWidget(QLabel("FDR Threshold (q-value):"))
        lay_filt.addWidget(self.spn_pval)
        lay_filt.addWidget(QLabel("Dimensionality Reducer:"))
        lay_filt.addWidget(self.cmb_dimred)
        lay_params.addLayout(lay_filt)

        # Advanced Rituals Checkboxes
        lay_adv = QHBoxLayout()
        self.chk_ml = QCheckBox("Awaken Machine Learning Golems (RF/SVM)")
        self.chk_ppi = QCheckBox("Summon STRING Interactome Spirits")
        self.chk_iso = QCheckBox("Hunt for Frankenstein Isoform Mutations")
        lay_adv.addWidget(self.chk_ml)
        lay_adv.addWidget(self.chk_ppi)
        lay_adv.addWidget(self.chk_iso)
        lay_params.addLayout(lay_adv)
        
        layout.addWidget(grp_params, 0)
        
        self.btn_unleash = QPushButton("Awaken the Frankenstein Overlord")
        self.btn_unleash.setObjectName("btnPrimary")
        self.btn_unleash.setMinimumHeight(45)
        layout.addWidget(self.btn_unleash, 0)
        self.btn_unleash.clicked.connect(self.execute_analysis)
        
        self.txt_de_res = QTextEdit()
        self.txt_de_res.setReadOnly(True)
        layout.addWidget(QLabel("Top Significant Differentially Expressed Genes (Top 20):"))
        layout.addWidget(self.txt_de_res, 1)
        
        self.workspace_stack.addWidget(page)

    def build_viz_workspace(self):
        page = QWidget()
        lay_viz = QVBoxLayout(page)
        lay_viz.setContentsMargins(20, 20, 20, 20)
        
        tabs = QTabWidget()
        
        # Static Plots Tab
        tab_static = QWidget()
        lay_static = QVBoxLayout(tab_static)
        lay_static.setContentsMargins(15, 15, 15, 15)
        
        lay_ctrl = QHBoxLayout()
        self.cmb_viz = QComboBox()
        self.cmb_viz.addItems([
            "Dimensionality Reduction Projection", 
            "Hierarchical Clustermap", 
            "Volcano Plot", 
            "MA Constellation Plot",
            "WGCNA Co-expression Dark Tapestry",
            "Kaplan-Meier Mortal Doom Curve",
            "Machine Learning: Golem Feature Importance",
            "STRING Interactome: Hairball of Chaos",
            "Splicing: Frankenstein Isoform Variance"
        ])
        self.cmb_viz.currentTextChanged.connect(self.render_graphics)
        lay_ctrl.addWidget(self.cmb_viz)
        
        self.cmb_palette = QComboBox()
        self.cmb_palette.addItems(["mako", "viridis", "magma", "coolwarm", "Spectral"])
        self.cmb_palette.currentTextChanged.connect(self.render_graphics)
        lay_ctrl.addWidget(QLabel("Color Palette:"))
        lay_ctrl.addWidget(self.cmb_palette)
        
        self.cmb_genes = QComboBox()
        self.cmb_genes.currentTextChanged.connect(self.render_graphics)
        lay_ctrl.addWidget(QLabel("Select Gene (For KM Curve):"))
        lay_ctrl.addWidget(self.cmb_genes)
        
        lay_static.addLayout(lay_ctrl)
        
        self.canvas = FigureCanvas(Figure(figsize=(8, 6)))
        lay_static.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self)
        lay_static.addWidget(self.toolbar)
        
        btn_save_plot = QPushButton("Export Graphic (PNG/PDF)")
        btn_save_plot.clicked.connect(self.export_static_graphic)
        lay_static.addWidget(btn_save_plot)
        
        tabs.addTab(tab_static, "Publication Graphics")
        
        # Interactive Plots Tab
        tab_inter = QWidget()
        lay_inter = QVBoxLayout(tab_inter)
        lay_inter.setContentsMargins(15, 15, 15, 15)
        lay_inter.addWidget(QLabel("Generate browser-compatible interactive visualisations utilizing the Plotly engine."))
        btn_inter = QPushButton("Export Interactive Visualizations (HTML)")
        btn_inter.setObjectName("btnPrimary")
        btn_inter.setFixedWidth(300)
        btn_inter.clicked.connect(self.export_interactive_html)
        lay_inter.addWidget(btn_inter, alignment=Qt.AlignmentFlag.AlignLeft)
        lay_inter.addStretch()
        tabs.addTab(tab_inter, "Interactive Web Graphics")
        
        lay_viz.addWidget(tabs)
        self.workspace_stack.addWidget(page)

    def build_gsea_workspace(self):
        page = QWidget()
        lay_gsea = QVBoxLayout(page)
        lay_gsea.setContentsMargins(20, 20, 20, 20)
        
        grp_db = QGroupBox("Gene Set Enrichment Configuration")
        lay_db = QHBoxLayout(grp_db)
        lay_db.setContentsMargins(20, 30, 20, 20)
        
        self.cmb_gsea_db = QComboBox()
        self.cmb_gsea_db.addItems(["KEGG_2021_Human", "GO_Biological_Process_2021", "Reactome_2022", "WikiPathway_2021_Human"])
        lay_db.addWidget(QLabel("Biological Database:"))
        lay_db.addWidget(self.cmb_gsea_db)
        
        self.cmb_gsea_terms = QComboBox()
        self.cmb_gsea_terms.currentTextChanged.connect(self.render_gsea_plot)
        lay_db.addWidget(QLabel("Significant Enriched Pathway:"))
        lay_db.addWidget(self.cmb_gsea_terms)
        
        lay_gsea.addWidget(grp_db, 0)
        
        self.canvas_gsea = FigureCanvas(Figure(figsize=(8, 6)))
        lay_gsea.addWidget(self.canvas_gsea, 1)
        
        self.workspace_stack.addWidget(page)

    def build_help_workspace(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(20, 20, 20, 20)
        text_area = QTextEdit()
        text_area.setReadOnly(True)
        help_content = """
        <h2 style="color: #ff8c00;">BRSA Pipeline - The Grimoire</h2>
        <hr style="border: 1px solid #414868;">
        <h3 style="color: #c0caf5;">1. PyDESeq2 & PyComBat</h3>
        <p>This dashboard has been upgraded with empirical Bayes shrinkage via PyDESeq2 and robust batch exorcism via PyComBat. This ensures your statistical models are publication-ready.</p>
        
        <h3 style="color: #c0caf5;">2. Advanced Arcane Mechanics</h3>
        <ul>
            <li><b>WGCNA Dark Tapestry:</b> Groups co-expressed genes into biological guilds (modules).</li>
            <li><b>Kaplan-Meier Mortal Doom Curve:</b> Requires 'time' and 'event' columns in your metadata CSV to predict survival fate based on gene expression.</li>
            <li><b>Machine Learning Golems:</b> Trains Random Forest and SVM models on expression data to identify predictive biomarkers.</li>
            <li><b>Frankenstein Isoform Mutations:</b> Hunts for aberrant variance in isoform-level matrices.</li>
            <li><b>STRING Interactome:</b> Invokes the spirits of the STRING database via API to map physical protein-protein interactions of significant genes.</li>
        </ul>
        """
        text_area.setHtml(help_content)
        layout.addWidget(text_area)
        self.workspace_stack.addWidget(page)

    def inject_matrix(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Expression Matrix", "", "CSV (*.csv)")
        if f:
            self.session_scroll.counts_matrix_path = f
            self.lbl_counts.setText(f"Matrix: {os.path.basename(f)}")
            print(f"[IO] Matrix loaded: {f}")

    def inject_metadata(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Metadata Profile", "", "CSV (*.csv)")
        if f:
            self.session_scroll.metadata_path = f
            self.lbl_meta.setText(f"Metadata: {os.path.basename(f)}")
            print(f"[IO] Metadata loaded: {f}")
            df = pd.read_csv(f)
            self.populate_table(self.table_meta, df)
            
            if 'condition' in df.columns:
                conds = df['condition'].unique()
                self.cmb_cond1.clear(); self.cmb_cond2.clear()
                self.cmb_cond1.addItems(conds); self.cmb_cond2.addItems(conds)
            
            batch_cols = [c for c in df.columns if c not in ['sample', 'condition']]
            self.cmb_batch.clear()
            self.cmb_batch.addItem("None")
            self.cmb_batch.addItems(batch_cols)

    def turn_on_processing(self):
        dlg = FastqProcessingDialog(self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            if not dlg.arsenal or not dlg.ref_index or not dlg.out_dir:
                QMessageBox.warning(self, "Configuration Error", "Required processing parameters are missing.")
                return
                
            aligner = "HISAT2" if "HISAT2" in dlg.cmb_aligner.currentText() else "Kallisto"
            self.devourer_thread = FastqProcessingThread(
                dlg.arsenal, dlg.out_dir, dlg.ref_index, 
                dlg.ann_file, dlg.chk_gff3.isChecked(), dlg.chk_trim.isChecked(), aligner
            )
            self.devourer_thread.progress_echoed.connect(lambda v,m: self.update_status(v, m))
            self.devourer_thread.digestion_completed.connect(self.fastq_completion)
            self.devourer_thread.fatal_choke.connect(self.critical_error)
            self.btn_unleash.setEnabled(False)
            print("[PROC] Spawning pre-processing pipeline thread...")
            self.devourer_thread.start()

    def update_status(self, val, msg):
        self.prog_bar.setValue(val)
        print(f"[UPDATE] {msg}")

    def fastq_completion(self, output_matrix):
        self.session_scroll.counts_matrix_path = output_matrix
        self.lbl_counts.setText(f"Generated Matrix: {os.path.basename(output_matrix)}")
        self.btn_unleash.setEnabled(True)
        print("[SUCCESS] Pre-processing finalized. FastQC and MultiQC reports generated in output directory.")
        QMessageBox.information(self, "Execution Successful", "FASTQ quantification routine completed successfully.")

    def critical_error(self, msg):
        self.btn_unleash.setEnabled(True)
        print(f"[CRITICAL FAILURE] {msg}")
        QMessageBox.critical(self, "Execution Error", msg)

    def execute_analysis(self):
        if not self.session_scroll.counts_matrix_path or not self.session_scroll.metadata_path:
            QMessageBox.warning(self, "Missing Data", "Please input both expression matrix and metadata files.")
            return
            
        self.session_scroll.alpha_condition = self.cmb_cond1.currentText()
        self.session_scroll.omega_condition = self.cmb_cond2.currentText()
        self.session_scroll.normalization_flavor = self.cmb_norm.currentText().split()[0]
        self.session_scroll.batch_demon = self.cmb_batch.currentText()
        self.session_scroll.min_read_threshold = self.spn_min.value()
        self.session_scroll.fdr_cliff = self.spn_pval.value()
        self.session_scroll.dimensionality_weapon = self.cmb_dimred.currentText()
        self.session_scroll.gsea_database_relic = self.cmb_gsea_db.currentText()
        self.session_scroll.summon_ml_golems = self.chk_ml.isChecked()
        self.session_scroll.summon_ppi_interactome = self.chk_ppi.isChecked()
        self.session_scroll.hunt_isoform_mutations = self.chk_iso.isChecked()
        
        self.btn_unleash.setEnabled(False)
        print("[PROC] Awakening the core Frankenstein analysis thread...")
        self.cruncher_thread = RNASeqAnalysisThread(self.session_scroll)
        self.cruncher_thread.brainwaves_updated.connect(lambda v: self.update_status(v, "Calculating PyDESeq2 statistical models..."))
        self.cruncher_thread.monster_awakened.connect(self.analysis_complete)
        self.cruncher_thread.lab_exploded.connect(self.critical_error)
        self.cruncher_thread.start()

    def analysis_complete(self, treasures):
        self.session_scroll.harvested_treasures = treasures
        self.btn_unleash.setEnabled(True)
        self.update_status(100, "Statistical evaluation complete.")
        
        sig_count = len(treasures['significant_genes'])
        self.lbl_quick_stat.setText(f"Genes Quantified: {len(treasures['de_results'])}\nSignificant (FDR): {sig_count}\nStatus: Analysis Complete")
        
        if 'normalized_counts' in treasures:
            self.populate_table(self.table_counts, treasures['normalized_counts'])
            
        self.display_differential_results()
        
        # Populate gene selector for KM
        self.cmb_genes.clear()
        self.cmb_genes.addItems(treasures['de_results'].index.tolist()[:100])
        
        self.render_graphics()
        
        if 'gsea_results' in treasures and not treasures['gsea_results'].empty:
            sig_terms = treasures['gsea_results'][treasures['gsea_results']['FDR q-val'] < 0.25]['Term'].tolist()
            self.cmb_gsea_terms.clear()
            self.cmb_gsea_terms.addItems(sig_terms)
            
        self.sidebar.setCurrentRow(3)
        QMessageBox.information(self, "Execution Successful", "PyDESeq2 and advanced rituals applied successfully.")

    def display_differential_results(self):
        de_res = self.session_scroll.harvested_treasures['de_results']
        sig_genes = de_res[de_res['significant']]
        
        html_table = sig_genes.sort_values('p-value').head(20)[['mean_1', 'mean_2', 'log2fc', 'p-value', 'p-adj']].to_html(
            float_format=lambda x: f"{x:.4e}",
            formatters={'mean_1': '{:.2f}'.format, 'mean_2': '{:.2f}'.format, 'log2fc': '{:.2f}'.format}
        )
        style = "<style>th, td { padding: 6px; border: 1px solid #292e42; color: #a9b1d6; } th { background-color: #1f2335; color: #ff8c00; text-align: left; }</style>"
        self.txt_de_res.setHtml(style + html_table)

    def render_graphics(self, *args):
        if not self.session_scroll.harvested_treasures: return
        
        opt = self.cmb_viz.currentText()
        palette = self.cmb_palette.currentText()
        self.canvas.figure.clf()
        
        if "Clustermap" in opt:
            fig = summon_clustermap_demon(self.session_scroll.harvested_treasures['heatmap_matrix'], self.session_scroll.harvested_treasures['metadata'], palette)
            lay = self.canvas.parentWidget().layout()
            lay.removeWidget(self.canvas)
            self.canvas.deleteLater()
            self.canvas = FigureCanvas(fig)
            lay.insertWidget(1, self.canvas)
            return

        ax = self.canvas.figure.add_subplot(111)
        self.canvas.figure.patch.set_facecolor('#1a1b26')
        ax.set_facecolor('#1a1b26')
        ax.tick_params(colors='#a9b1d6')
        ax.xaxis.label.set_color('#a9b1d6')
        ax.yaxis.label.set_color('#a9b1d6')
        ax.title.set_color('#ff8c00')
        for spine in ax.spines.values(): spine.set_edgecolor('#292e42')

        if "Dimensionality" in opt:
            paint_dimensional_void_masterpiece(ax, self.session_scroll.harvested_treasures['dim_red'], self.session_scroll.harvested_treasures['metadata'])
        elif "Volcano" in opt:
            erupt_volcano_plot(ax, self.session_scroll.harvested_treasures['de_results'], self.session_scroll.fdr_cliff)
        elif "MA Constellation" in opt:
            paint_ma_constellation(ax, self.session_scroll.harvested_treasures['de_results'])
        elif "WGCNA" in opt:
            if 'wgcna_linkage' in self.session_scroll.harvested_treasures:
                plot_dendrogram_of_madness(ax, self.session_scroll.harvested_treasures['wgcna_linkage'])
            else:
                ax.text(0.5, 0.5, "WGCNA Tapestry not available.", ha='center', va='center', color='white')
        elif "Doom Curve" in opt:
            if self.session_scroll.harvested_treasures.get('has_survival'):
                gene = self.cmb_genes.currentText()
                predict_mortal_doom_curves(ax, self.session_scroll.harvested_treasures['metadata'], self.session_scroll.harvested_treasures['normalized_counts'], gene)
            else:
                ax.text(0.5, 0.5, "Mortal Doom Curve requires 'time' and 'event' columns in metadata.", ha='center', va='center', color='white')
        elif "Machine Learning" in opt:
            if self.session_scroll.harvested_treasures.get('ml_models_ready'):
                interrogate_feature_importance_spirits(ax, self.session_scroll.harvested_treasures['ml_rf_importance'])
            else:
                ax.text(0.5, 0.5, "ML Golems were not awakened during analysis.", ha='center', va='center', color='white')
        elif "Interactome" in opt:
            if 'ppi_network' in self.session_scroll.harvested_treasures:
                paint_hairball_of_chaos(ax, self.session_scroll.harvested_treasures['ppi_network'])
            else:
                ax.text(0.5, 0.5, "STRING Interactome was not summoned or no network found.", ha='center', va='center', color='white')
        elif "Splicing" in opt:
            if 'frankenstein_isoforms' in self.session_scroll.harvested_treasures:
                uncover_hidden_isoform_mutations(ax, self.session_scroll.harvested_treasures['frankenstein_isoforms'], self.session_scroll.harvested_treasures['normalized_counts'])
            else:
                ax.text(0.5, 0.5, "Frankenstein Isoform Mutations not hunted.", ha='center', va='center', color='white')
            
        self.canvas.draw()

    def render_gsea_plot(self, term):
        if not term or 'gsea_obj' not in self.session_scroll.harvested_treasures: return
        self.canvas_gsea.figure.clf()
        ax = self.canvas_gsea.figure.add_subplot(111)
        self.canvas_gsea.figure.patch.set_facecolor('#1a1b26')
        ax.set_facecolor('#1a1b26')
        ax.tick_params(colors='#a9b1d6')
        
        paint_gsea_mountain(ax, self.session_scroll.harvested_treasures['gsea_obj'], term)
        
        ax.xaxis.label.set_color('#a9b1d6')
        ax.yaxis.label.set_color('#a9b1d6')
        ax.title.set_color('#ff8c00')
        for spine in ax.spines.values(): spine.set_edgecolor('#292e42')
        self.canvas_gsea.draw()

    def export_static_graphic(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export Visual Artifact", "", "PNG Format (*.png);;PDF Format (*.pdf)")
        if path:
            try:
                self.canvas.figure.savefig(path, bbox_inches='tight', dpi=300)
                print(f"[IO] Graphic exported to {path}")
            except Exception as e:
                self.critical_error(f"Failed to export graphic: {str(e)}")

    def export_interactive_html(self):
        if not self.session_scroll.harvested_treasures: return
        d = QFileDialog.getExistingDirectory(self, "Select Output Directory for HTML Artifacts")
        if d:
            forge_interactive_plotly_artifacts(self.session_scroll.harvested_treasures, d)
            print(f"[IO] Interactive objects forged in directory: {d}")
            QMessageBox.information(self, "Export Successful", "Interactive HTML artifacts generated successfully.")

    def export_reproducibility_log(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export Session Metadata", "bio_reproducibility.txt", "Text Files (*.txt)")
        if path:
            params = self.session_scroll.__dict__.copy()
            scribe_reproducibility_chronicles(path, params)
            print(f"[IO] Protocol logged to {path}")

    def export_pdf_report(self):
        if not self.session_scroll.harvested_treasures: return
        path, _ = QFileDialog.getSaveFileName(self, "Export Analytical Document (ABNT)", "Relatorio_Biologico.pdf", "PDF Files (*.pdf)")
        if not path: return
            
        try:
            cm = 28.3465
            doc = SimpleDocTemplate(path, pagesize=A4, rightMargin=2*cm, leftMargin=3*cm, topMargin=3*cm, bottomMargin=2*cm)
            styles = getSampleStyleSheet()
            
            abnt_style = ParagraphStyle('ABNT', parent=styles['Normal'], fontName='Times-Roman', fontSize=12, leading=18, alignment=TA_JUSTIFY, spaceAfter=12)
            abnt_title = ParagraphStyle('ABNT_Title', parent=abnt_style, fontSize=14, fontName='Times-Bold', alignment=1)
            
            story = []
            story.append(Paragraph("RELATÓRIO DE ANÁLISE BIOINFORMÁTICA E EXPRESSÃO DIFERENCIAL", abnt_title))
            story.append(Spacer(1, 12))
            
            info = f"""
            A metodologia computacional aplicada neste protocolo foi rigidamente estruturada. 
            O modelo de inferência estatística estabeleceu um contraste entre o grupo <b>{self.session_scroll.alpha_condition}</b> e o grupo <b>{self.session_scroll.omega_condition}</b>. 
            Aplicou-se um limiar quantitativo inferior de reads avaliado em {self.session_scroll.min_read_threshold} por unidade transcrita. O algoritmo de escalonamento/modelo estatístico utilizado foi <b>{self.session_scroll.normalization_flavor}</b>. 
            Para múltiplas hipóteses adotou-se a correção FDR (Benjamini-Hochberg) padronizada em {self.session_scroll.fdr_cliff}. 
            A abstração multivariada empregada baseia-se no modelo topológico {self.session_scroll.dimensionality_weapon}.
            """
            story.append(Paragraph(info, abnt_style))
            
            de_res = self.session_scroll.harvested_treasures['de_results']
            sig_genes = de_res[de_res['significant']]
            story.append(Paragraph(f"O universo paramétrico contou com {len(de_res)} feições biológicas, culminando em {len(sig_genes)} marcadores que superaram o crivo de significância estabelecido.", abnt_style))
            story.append(Spacer(1, 12))
            
            top10 = sig_genes.sort_values('p-value').head(10)
            table_data = [['ID Gênico', 'Média Exp (Grp 1)', 'Média Exp (Grp 2)', 'Log2 FC', 'FDR q-val']]
            for gene, row in top10.iterrows():
                table_data.append([gene, f"{row['mean_1']:.2f}", f"{row['mean_2']:.2f}", f"{row['log2fc']:.2f}", f"{row['p-adj']:.2e}"])
            
            de_table = Table(table_data)
            de_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey), ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'), ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('FONTNAME', (0, 0), (-1, -1), 'Times-Roman')
            ]))
            story.append(de_table)
                
            doc.build(story)
            print("[IO] ABNT PDF Documentation generated successfully.")
            QMessageBox.information(self, "Export Successful", f"ABNT structured PDF allocated to {path}")
        except Exception as e:
            self.critical_error(f"Failure formatting PDF stream: {str(e)}")

    def closeEvent(self, event):
        if hasattr(self, 'temp_dir') and os.path.exists(self.temp_dir): shutil.rmtree(self.temp_dir)
        if self.devourer_thread and self.devourer_thread.isRunning(): self.devourer_thread.terminate()
        if self.cruncher_thread and self.cruncher_thread.isRunning(): self.cruncher_thread.terminate()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    dashboard = RNASeqDashboard()
    dashboard.show()
    sys.exit(app.exec())
