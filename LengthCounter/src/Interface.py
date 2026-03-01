import sys
import os
import json
import gzip
import traceback
import tempfile
import re
import pandas as pd
import seaborn as sns
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QLabel, QPushButton, QFileDialog, QComboBox, QSpinBox, 
                             QCheckBox, QTextEdit, QTabWidget, QMessageBox, QGroupBox, 
                             QColorDialog, QAction, QLineEdit, QListWidget, QStackedWidget, 
                             QSplitter, QTableWidget, QTableWidgetItem, QHeaderView, QDockWidget)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QPalette, QTextCursor, QTextCharFormat

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from logic import (UchihaSharinganAnalyzer, SenjuMokutonAnalyzer, UzumakiChakraFastqAnalyzer, AnbuBlackOpsExternalTools,
                   RinneganAdvancedBioTools, KatsuyuSlugFormatParser, ino_mind_transfer_stats_text, 
                   neji_byakugan_format_sniffer, tenten_scroll_export_csv, pain_chibaku_tensei_merge_filter, 
                   madara_infinite_tsukuyomi_html_report, sai_ninja_art_dark_theme, gamabunta_summon_sequence_file,
                   PainBanshoTeninNCBIFetcher, NejiVariantEyeParser, YamatoWoodFeatureExtractor, GaaraSandTrimmer, HengeFormatShifter)
from config import DOCUMENTATION_HTML, NINJA_SESSION_FILE, NARUTO_HOKAGE_STYLE


class KageBunshinAsyncJutsu(QThread):
    finished_signal = pyqtSignal()
    def __init__(self, analyzer, sequence):
        super().__init__()
        self.analyzer = analyzer
        self.sequence = sequence
        
    def run(self):
        self.analyzer.yamato_mokuton_set_sequence(self.sequence)
        self.finished_signal.emit()


class KonohaNinjaPlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=8, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.patch.set_facecolor('#1a1a1a')
        super().__init__(self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor('#121212')
        self.ax.tick_params(colors='#e0e0e0')
        for spine in self.ax.spines.values():
            spine.set_edgecolor('#333333')


class HokageTowerDashboardGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.analyzer = UchihaSharinganAnalyzer()
        self.base_analyzer = SenjuMokutonAnalyzer()
        self.fastq_analyzer = UzumakiChakraFastqAnalyzer()
        self.multi_parser = KatsuyuSlugFormatParser()
        
        self.current_bars = None 
        self.file_paths_cache = set() 
        self.generated_figures = []
        
        self.kuchiyose_init_ui()
        self.raijin_load_session()
        
    def kuchiyose_init_ui(self):
        self.setWindowTitle("LenghtCount")
        self.setGeometry(100, 100, 1400, 800)
        
        self.gokakyu_build_menu()
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        self.lista_ninja_sidebar = QListWidget()
        self.lista_ninja_sidebar.setObjectName("sidebar_ninja")
        self.lista_ninja_sidebar.setFixedWidth(250)
        self.lista_ninja_sidebar.addItem(" 👁️ Assembly Statistics (Sharingan)")
        self.lista_ninja_sidebar.addItem(" 🧬 Base Composition (Mokuton)")
        self.lista_ninja_sidebar.addItem(" ⚡ Protein & Primers (Chidori)")
        self.lista_ninja_sidebar.addItem(" 🌀 FastQ Metrics (Rasengan)")
        self.lista_ninja_sidebar.addItem(" 🌌 Advanced Sequence (Rinnegan)")
        self.lista_ninja_sidebar.addItem(" 🪬 Variant Viewer (Neji Eye)")
        self.lista_ninja_sidebar.addItem(" 📜 Direct Code (Scrolls)")
        self.lista_ninja_sidebar.addItem(" 📖 Documentation")
        self.lista_ninja_sidebar.setCurrentRow(0)
        self.lista_ninja_sidebar.currentRowChanged.connect(self.tsuyoi_mudar_pagina)
        
        self.arena_batalha_stack = QStackedWidget()
        
        self.page_length = QWidget()
        self.kage_bunshin_build_fasta_page(self.page_length)
        self.arena_batalha_stack.addWidget(self.page_length)
        
        self.page_bases = QWidget()
        self.kage_bunshin_build_base_page(self.page_bases)
        self.arena_batalha_stack.addWidget(self.page_bases)
        
        self.page_protein = QWidget()
        self.kage_bunshin_build_protein_page(self.page_protein)
        self.arena_batalha_stack.addWidget(self.page_protein)
        
        self.page_fastq = QWidget()
        self.kage_bunshin_build_fastq_page(self.page_fastq)
        self.arena_batalha_stack.addWidget(self.page_fastq)
        
        self.page_comparative = QWidget()
        self.kage_bunshin_build_comparative_page(self.page_comparative)
        self.arena_batalha_stack.addWidget(self.page_comparative)

        self.page_vcf = QWidget()
        self.kage_bunshin_build_vcf_page(self.page_vcf)
        self.arena_batalha_stack.addWidget(self.page_vcf)
        
        self.page_editor = QWidget()
        self.kage_bunshin_build_editor_page(self.page_editor)
        self.arena_batalha_stack.addWidget(self.page_editor)
        
        self.page_help = QWidget()
        self.kage_bunshin_build_help_page(self.page_help)
        self.arena_batalha_stack.addWidget(self.page_help)
        
        main_layout.addWidget(self.lista_ninja_sidebar)
        main_layout.addWidget(self.arena_batalha_stack)
        
        self.build_anbu_floating_terminal()
        
        self.statusBar().showMessage("LenghtCount System Initialization Complete.", 5000)

    def tsuyoi_mudar_pagina(self, index):
        self.arena_batalha_stack.setCurrentIndex(index)

    def gokakyu_build_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu('System File')
        save_act = QAction('Save Current Session', self)
        save_act.triggered.connect(self.orochimaru_save_session)
        file_menu.addAction(save_act)
        clear_act = QAction('Purge Active Workspace', self)
        clear_act.triggered.connect(self.shinra_tensei_clear_workspace)
        file_menu.addAction(clear_act)

    def build_anbu_floating_terminal(self):
        self.dock_terminal = QDockWidget("Anbu Floating HQ (Embedded Shell)", self)
        self.dock_terminal.setAllowedAreas(Qt.BottomDockWidgetArea | Qt.RightDockWidgetArea)
        
        term_widget = QWidget()
        tl = QVBoxLayout(term_widget)
        
        h = QHBoxLayout()
        self.term_input = QLineEdit()
        self.term_input.setPlaceholderText("Command (e.g., docker exec -it bio_tools bash)")
        self.term_btn = QPushButton("Execute Shell Jutsu")
        self.term_btn.clicked.connect(self.anbu_terminal_execute)
        h.addWidget(self.term_input)
        h.addWidget(self.term_btn)
        
        self.term_output = QTextEdit()
        self.term_output.setReadOnly(True)
        self.term_output.setFont(QFont("Consolas", 10))
        self.term_output.setStyleSheet("background-color: #000; color: #0f0;")
        
        tl.addLayout(h)
        tl.addWidget(self.term_output)
        
        self.dock_terminal.setWidget(term_widget)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.dock_terminal)

    def anbu_terminal_execute(self):
        cmd = self.term_input.text().strip()
        if cmd:
            self.term_output.append(f"$ {cmd}")
            out = AnbuBlackOpsExternalTools.itachi_mangekyou_run_command(cmd)
            self.term_output.append(out + "\n")

    def _ui_generic_file_loader(self, target_input_widget, title="Open Genomic Data File"):
        path, _ = QFileDialog.getOpenFileName(self, title, "", "Text/FASTA Files (*.txt *.fasta *.fa);;All Files (*)")
        if path:
            seq = gamabunta_summon_sequence_file(path)
            if seq:
                target_input_widget.setText(seq) if hasattr(target_input_widget, 'setText') else target_input_widget.setPlainText(seq)

    def kage_bunshin_build_fasta_page(self, page_widget):
        layout = QHBoxLayout(page_widget)
        controls_panel = QWidget()
        controls_panel.setFixedWidth(380)
        controls_layout = QVBoxLayout(controls_panel)
        
        file_group = QGroupBox("Genomic File Management")
        fl = QVBoxLayout(file_group)
        self.file_list = QComboBox()
        btn_h = QHBoxLayout()
        self.add_file_btn = QPushButton("Import File(s)")
        self.add_file_btn.clicked.connect(self.shuriken_add_mass_files)
        self.add_folder_btn = QPushButton("Import Directory")
        self.add_folder_btn.clicked.connect(self.shuriken_add_folder)
        btn_h.addWidget(self.add_file_btn)
        btn_h.addWidget(self.add_folder_btn)
        self.remove_file_btn = QPushButton("Erase Selected File")
        self.remove_file_btn.clicked.connect(self.kunai_remove_file)
        self.merge_files_btn = QPushButton("Concatenate Selected (Chibaku Tensei)")
        self.merge_files_btn.clicked.connect(self.rinnegan_merge_files)
        fl.addWidget(self.file_list)
        fl.addLayout(btn_h)
        fl.addWidget(self.remove_file_btn)
        fl.addWidget(self.merge_files_btn)

        ncbi_group = QGroupBox("Pain Bansho Tenin (NCBI Downloader)")
        nl = QVBoxLayout(ncbi_group)
        h_ncbi = QHBoxLayout()
        self.ncbi_input = QLineEdit()
        self.ncbi_input.setPlaceholderText("Accession No. (e.g., NC_045512)")
        self.ncbi_btn = QPushButton("Pull Sequence")
        self.ncbi_btn.clicked.connect(self.pain_fetch_ncbi_clicked)
        h_ncbi.addWidget(self.ncbi_input)
        h_ncbi.addWidget(self.ncbi_btn)
        nl.addLayout(h_ncbi)

        henge_group = QGroupBox("Henge Format Shifter (Converter)")
        hl = QVBoxLayout(henge_group)
        self.henge_in_path = QLineEdit()
        self.henge_in_path.setPlaceholderText("Select File to Convert...")
        self.henge_in_btn = QPushButton("Browse")
        self.henge_in_btn.clicked.connect(lambda: self.henge_in_path.setText(QFileDialog.getOpenFileName(self, "Open File")[0]))
        h_h1 = QHBoxLayout()
        h_h1.addWidget(self.henge_in_path)
        h_h1.addWidget(self.henge_in_btn)
        h_h2 = QHBoxLayout()
        self.henge_combo_in = QComboBox()
        self.henge_combo_in.addItems(["bam", "fastq", "gbk"])
        self.henge_combo_out = QComboBox()
        self.henge_combo_out.addItems(["sam", "fasta"])
        h_h2.addWidget(QLabel("From:"))
        h_h2.addWidget(self.henge_combo_in)
        h_h2.addWidget(QLabel("To:"))
        h_h2.addWidget(self.henge_combo_out)
        self.henge_run_btn = QPushButton("Transform Data Format")
        self.henge_run_btn.clicked.connect(self.henge_format_convert_clicked)
        hl.addLayout(h_h1)
        hl.addLayout(h_h2)
        hl.addWidget(self.henge_run_btn)
        
        trim_group = QGroupBox("Virtual Trimming Engine")
        tl = QVBoxLayout(trim_group)
        t_h1 = QHBoxLayout()
        t_h1.addWidget(QLabel("Minimum Threshold (bp):"))
        self.min_len_spin = QSpinBox()
        self.min_len_spin.setRange(0, 999999999)
        self.min_len_spin.setValue(0)
        t_h1.addWidget(self.min_len_spin)
        t_h2 = QHBoxLayout()
        t_h2.addWidget(QLabel("Maximum Threshold (bp):"))
        self.max_len_spin = QSpinBox()
        self.max_len_spin.setRange(0, 999999999)
        self.max_len_spin.setValue(999999999)
        t_h2.addWidget(self.max_len_spin)
        self.apply_filter_btn = QPushButton("Apply Hardware Filters")
        self.apply_filter_btn.clicked.connect(self.byakugan_update_stats)
        tl.addLayout(t_h1)
        tl.addLayout(t_h2)
        tl.addWidget(self.apply_filter_btn)
        
        controls_layout.addWidget(file_group)
        controls_layout.addWidget(ncbi_group)
        controls_layout.addWidget(henge_group)
        controls_layout.addWidget(trim_group)
        controls_layout.addStretch()
        
        results_panel = QTabWidget()
        self.plot_tabs_widget = QWidget()
        plot_tabs_layout = QVBoxLayout(self.plot_tabs_widget)
        self.plot_tabs = QTabWidget()
        self.plot_tabs.setTabsClosable(True)
        self.plot_tabs.tabCloseRequested.connect(self.susanoo_close_plot_tab)
        
        ph_layout = QHBoxLayout()
        self.plot_type = QComboBox()
        self.plot_type.addItems(["Histogram Algorithm", "Boxplot Analysis", "Violin Density", "KDE Density", "Cumulative Process (ECDF)"])
        self.normal_check = QCheckBox("Normal Overlay")
        self.plot_btn = QPushButton("Execute Plot Rendering")
        self.plot_btn.clicked.connect(self.raikiri_update_plot)
        ph_layout.addWidget(self.plot_type)
        ph_layout.addWidget(self.normal_check)
        ph_layout.addWidget(self.plot_btn)
        plot_tabs_layout.addLayout(ph_layout)
        plot_tabs_layout.addWidget(self.plot_tabs)
        
        self.stats_tab_inner = QWidget()
        stats_tab_layout = QVBoxLayout(self.stats_tab_inner)
        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setFont(QFont("Consolas", 11))
        stats_tab_layout.addWidget(self.stats_text)
        
        results_panel.addTab(self.stats_tab_inner, "Terminal Output Statistics")
        results_panel.addTab(self.plot_tabs_widget, "Data Renderers")
        
        layout.addWidget(controls_panel)
        layout.addWidget(results_panel, stretch=1)

    def kage_bunshin_build_base_page(self, page_widget):
        layout = QVBoxLayout(page_widget)
        top_split = QHBoxLayout()
        
        input_panel = QWidget()
        il = QVBoxLayout(input_panel)
        
        sg = QGroupBox("Nucleotide Sequence Source Input")
        sl = QVBoxLayout(sg)
        h_file = QHBoxLayout()
        self.load_seq_file_btn = QPushButton("Import from File")
        self.load_seq_file_btn.clicked.connect(lambda: self._ui_generic_file_loader(self.seq_input))
        h_file.addWidget(self.load_seq_file_btn)
        h_file.addStretch()
        
        self.seq_input = QTextEdit()
        self.seq_input.setPlaceholderText("Direct input or load external DNA/RNA FASTA file...")
        
        btn_h = QHBoxLayout()
        self.analyze_seq_btn = QPushButton("Initiate Sequence Profiling")
        self.analyze_seq_btn.clicked.connect(self.chidori_analyze_sequence_async)
        self.translate_btn = QPushButton("In-Silico Translation")
        self.translate_btn.clicked.connect(self.byakugou_translate)
        btn_h.addWidget(self.analyze_seq_btn)
        btn_h.addWidget(self.translate_btn)

        btn_h2 = QHBoxLayout()
        self.revcomp_btn = QPushButton("Reverse Complement (Shadow Reversal)")
        self.revcomp_btn.clicked.connect(self.shikamaru_revcomp_clicked)
        self.rna_btn = QPushButton("Transcribe to RNA")
        self.rna_btn.clicked.connect(self.jiraiya_transcribe_clicked)
        btn_h2.addWidget(self.revcomp_btn)
        btn_h2.addWidget(self.rna_btn)
        
        sl.addLayout(h_file)
        sl.addWidget(self.seq_input)
        sl.addLayout(btn_h)
        sl.addLayout(btn_h2)

        rx_g = QGroupBox("Sharingan Regex Eye (Conditional Expression)")
        rx_l = QHBoxLayout(rx_g)
        self.regex_input = QLineEdit()
        self.regex_input.setPlaceholderText("Enter RegExp pattern (e.g., AT[GC]T)")
        self.regex_btn = QPushButton("Highlight Matches")
        self.regex_btn.clicked.connect(self.sharingan_regex_highlight)
        rx_l.addWidget(self.regex_input)
        rx_l.addWidget(self.regex_btn)
        
        ag = QGroupBox("Advanced Genomic Characteristics")
        al = QVBoxLayout(ag)
        
        h1 = QHBoxLayout()
        self.kmer_spin = QSpinBox()
        self.kmer_spin.setRange(2, 8)
        self.kmer_spin.setValue(3)
        self.kmer_btn = QPushButton("K-mer Swarm Matrix Analysis")
        self.kmer_btn.clicked.connect(self.suiton_plot_kmers)
        self.kraken_btn = QPushButton("Inuzuka Dog Taxonomy (Kraken-lite)")
        self.kraken_btn.clicked.connect(self.kraken_lite_clicked)
        h1.addWidget(QLabel("K-mer:"))
        h1.addWidget(self.kmer_spin)
        h1.addWidget(self.kmer_btn)
        h1.addWidget(self.kraken_btn)
        
        h2 = QHBoxLayout()
        self.ent_btn = QPushButton("Shannon Entropy")
        self.ent_btn.clicked.connect(self.fuuton_plot_entropy)
        self.dust_btn = QPushButton("Gaara Dust Complexity Shield")
        self.dust_btn.clicked.connect(self.dust_complexity_clicked)
        self.cpg_btn = QPushButton("CpG Loci")
        self.cpg_btn.clicked.connect(self.raiton_plot_cpg)
        h2.addWidget(self.ent_btn)
        h2.addWidget(self.dust_btn)
        h2.addWidget(self.cpg_btn)

        h3 = QHBoxLayout()
        self.restrict_btn = QPushButton("Zabuza Restriction Map & Gel")
        self.restrict_btn.clicked.connect(self.zabuza_restriction_map_clicked)
        self.pai_btn = QPushButton("Orochimaru Cursed PAIs")
        self.pai_btn.clicked.connect(self.pai_islands_clicked)
        h3.addWidget(self.restrict_btn)
        h3.addWidget(self.pai_btn)

        h4 = QHBoxLayout()
        self.telo_btn = QPushButton("Kimimaro Bone Telomere Mapper")
        self.telo_btn.clicked.connect(self.telomere_clicked)
        self.rho_btn = QPushButton("Neji Trigram Rho Odds Ratio")
        self.rho_btn.clicked.connect(self.rho_odds_clicked)
        h4.addWidget(self.telo_btn)
        h4.addWidget(self.rho_btn)
        
        al.addLayout(h1)
        al.addLayout(h2)
        al.addLayout(h3)
        al.addLayout(h4)
        
        il.addWidget(sg)
        il.addWidget(rx_g)
        il.addWidget(ag)
        il.addStretch()
        
        results_panel = QWidget()
        rl = QVBoxLayout(results_panel)
        self.base_canvas = KonohaNinjaPlotCanvas(self)
        self.base_stats_text = QTextEdit()
        self.base_stats_text.setReadOnly(True)
        self.base_stats_text.setFont(QFont("Consolas", 11))
        
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.base_canvas)
        splitter.addWidget(self.base_stats_text)
        rl.addWidget(splitter)
        
        top_split.addWidget(input_panel, stretch=1)
        top_split.addWidget(results_panel, stretch=2)
        layout.addLayout(top_split)

    def kage_bunshin_build_protein_page(self, page_widget):
        layout = QHBoxLayout(page_widget)
        left = QWidget()
        ll = QVBoxLayout(left)
        
        og = QGroupBox("Sage Mode ORF Scanner")
        ol = QVBoxLayout(og)
        h_o = QHBoxLayout()
        self.orf_btn_file = QPushButton("Load File")
        self.orf_btn_file.clicked.connect(lambda: self._ui_generic_file_loader(self.orf_input))
        h_o.addWidget(self.orf_btn_file)
        h_o.addStretch()
        self.orf_input = QTextEdit()
        self.orf_input.setPlaceholderText("Input raw genome segment or FASTA file...")
        self.orf_btn = QPushButton("Engage Heuristic ORF Search")
        self.orf_btn.clicked.connect(self.katon_find_orfs)
        ol.addLayout(h_o)
        ol.addWidget(self.orf_input)
        ol.addWidget(self.orf_btn)
        
        pg = QGroupBox("Peptide/Protein Physicochemical Properties")
        pl = QVBoxLayout(pg)
        h_p = QHBoxLayout()
        self.prot_btn_file = QPushButton("Load File")
        self.prot_btn_file.clicked.connect(lambda: self._ui_generic_file_loader(self.prot_input))
        h_p.addWidget(self.prot_btn_file)
        h_p.addStretch()
        self.prot_input = QTextEdit()
        self.prot_input.setPlaceholderText("Input specific amino acid chain data...")
        
        h = QHBoxLayout()
        self.prot_prop_btn = QPushButton("Gai Eight Gates Properties")
        self.prot_prop_btn.clicked.connect(self.suiton_protein_props)
        self.aa_plot_btn = QPushButton("Chouji AA Frequency Plot")
        self.aa_plot_btn.clicked.connect(self.chouji_aa_plot_clicked)
        self.hydro_btn = QPushButton("Hydropathy Plot")
        self.hydro_btn.clicked.connect(self.doton_hydrophobicity)
        h.addWidget(self.prot_prop_btn)
        h.addWidget(self.aa_plot_btn)
        h.addWidget(self.hydro_btn)
        
        pl.addLayout(h_p)
        pl.addWidget(self.prot_input)
        pl.addLayout(h)
        
        prg = QGroupBox("Oligonucleotide Thermodynamic & CRISPR Inspector")
        prl = QVBoxLayout(prg)
        h_pr = QHBoxLayout()
        self.primer_btn_file = QPushButton("Load File")
        self.primer_btn_file.clicked.connect(lambda: self._ui_generic_file_loader(self.primer_input))
        h_pr.addWidget(self.primer_btn_file)
        h_pr.addStretch()
        self.primer_input = QLineEdit()
        self.primer_input.setPlaceholderText("Direct input IUPAC sequence...")
        
        hp2 = QHBoxLayout()
        self.primer_btn = QPushButton("Structural Dynamics")
        self.primer_btn.clicked.connect(self.raiton_primer_analysis)
        self.crispr_btn = QPushButton("Kakashi Copy CRISPR PAM Scan")
        self.crispr_btn.clicked.connect(self.kakashi_crispr_clicked)
        hp2.addWidget(self.primer_btn)
        hp2.addWidget(self.crispr_btn)

        prl.addLayout(h_pr)
        prl.addWidget(self.primer_input)
        prl.addLayout(hp2)
        
        ll.addWidget(og)
        ll.addWidget(pg)
        ll.addWidget(prg)
        
        right = QVBoxLayout()
        self.prot_canvas = KonohaNinjaPlotCanvas(self)
        self.prot_stats = QTextEdit()
        self.prot_stats.setReadOnly(True)
        self.prot_stats.setFont(QFont("Consolas", 10))
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.prot_canvas)
        splitter.addWidget(self.prot_stats)
        right.addWidget(splitter)
        
        layout.addWidget(left, stretch=1)
        layout.addLayout(right, stretch=1)

    def kage_bunshin_build_fastq_page(self, page_widget):
        layout = QHBoxLayout(page_widget)
        controls = QWidget()
        controls.setFixedWidth(350)
        cl = QVBoxLayout(controls)
        
        fg = QGroupBox("Sequencing Telemetry Files (FASTQ)")
        fl = QVBoxLayout(fg)
        self.fastq_file_list = QComboBox()
        self.add_fastq_btn = QPushButton("Import FastQ Read Blocks")
        self.add_fastq_btn.clicked.connect(self.shuriken_add_mass_files)
        fl.addWidget(self.fastq_file_list)
        fl.addWidget(self.add_fastq_btn)

        trim_group = QGroupBox("Gaara Sand Trimmer (FASTQ Cleaner)")
        t_l = QVBoxLayout(trim_group)
        self.gaara_fastq_in = QLineEdit()
        self.gaara_fastq_in.setPlaceholderText("FASTQ File to Trim...")
        self.gaara_fastq_btn = QPushButton("Browse")
        self.gaara_fastq_btn.clicked.connect(lambda: self.gaara_fastq_in.setText(QFileDialog.getOpenFileName(self, "Open FASTQ")[0]))
        h_t1 = QHBoxLayout()
        h_t1.addWidget(self.gaara_fastq_in)
        h_t1.addWidget(self.gaara_fastq_btn)
        h_t2 = QHBoxLayout()
        h_t2.addWidget(QLabel("Min Phred:"))
        self.gaara_phred_spin = QSpinBox()
        self.gaara_phred_spin.setRange(0, 40)
        self.gaara_phred_spin.setValue(20)
        h_t2.addWidget(self.gaara_phred_spin)
        self.gaara_run_btn = QPushButton("Execute Sand Burial (Trim Ends)")
        self.gaara_run_btn.clicked.connect(self.gaara_trimmer_clicked)
        t_l.addLayout(h_t1)
        t_l.addLayout(h_t2)
        t_l.addWidget(self.gaara_run_btn)
        
        pg = QGroupBox("Coverage Mathematics Parameters")
        pl = QVBoxLayout(pg)
        h1 = QHBoxLayout()
        h1.addWidget(QLabel("Genome Size Specimen (bp):"))
        self.genome_size_input = QSpinBox()
        self.genome_size_input.setRange(1000, 1000000000)
        self.genome_size_input.setValue(1000000)
        self.genome_size_input.valueChanged.connect(self._sync_genome_size)
        h1.addWidget(self.genome_size_input)
        h2 = QHBoxLayout()
        h2.addWidget(QLabel("Platform Read Length (bp):"))
        self.read_length_input = QSpinBox()
        self.read_length_input.setRange(35, 10000)
        self.read_length_input.setValue(150)
        h2.addWidget(self.read_length_input)
        self.paired_end_check = QCheckBox("Paired-End Sequencing Protocol Enabled")
        self.paired_end_check.setChecked(True)
        pl.addLayout(h1)
        pl.addLayout(h2)
        pl.addWidget(self.paired_end_check)
        self.calc_btn = QPushButton("Compute Extrapolated Coverage")
        self.calc_btn.clicked.connect(self.rasenshuriken_calculate_coverage)
        pl.addWidget(self.calc_btn)
        
        ag = QGroupBox("Base Quality Pro Tools")
        al = QVBoxLayout(ag)
        self.phred_btn = QPushButton("Phred Assessment Algorithm")
        self.phred_btn.clicked.connect(self.mangekyou_plot_phred)
        self.rare_btn = QPushButton("Rarefaction Trajectory")
        self.rare_btn.clicked.connect(self.rinnegan_plot_rarefaction)
        self.sat_btn = QPushButton("Kisame Saturation Chart")
        self.sat_btn.clicked.connect(self.kisame_saturation_clicked)
        self.dup_btn = QPushButton("Zetsu PCR Duplication Estimator")
        self.dup_btn.clicked.connect(self.zetsu_duplication_clicked)
        al.addWidget(self.phred_btn)
        al.addWidget(self.rare_btn)
        al.addWidget(self.sat_btn)
        al.addWidget(self.dup_btn)
        
        cl.addWidget(fg)
        cl.addWidget(trim_group)
        cl.addWidget(pg)
        cl.addWidget(ag)
        cl.addStretch()
        
        results = QWidget()
        rl = QVBoxLayout(results)
        self.fastq_canvas = KonohaNinjaPlotCanvas(self)
        self.fastq_results_text = QTextEdit()
        self.fastq_results_text.setReadOnly(True)
        self.fastq_results_text.setFont(QFont("Consolas", 11))
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.fastq_canvas)
        splitter.addWidget(self.fastq_results_text)
        rl.addWidget(splitter)
        
        layout.addWidget(controls)
        layout.addWidget(results, stretch=1)

    def _sync_genome_size(self):
        self.analyzer.expected_genome_size = self.genome_size_input.value()

    def kage_bunshin_build_comparative_page(self, page_widget):
        layout = QHBoxLayout(page_widget)
        left = QWidget()
        ll = QVBoxLayout(left)
        
        sg = QGroupBox("Global Alignment & Synteny (Ti/Tv included)")
        sl = QVBoxLayout(sg)
        h_s_1 = QHBoxLayout()
        self.snp_btn_file_ref = QPushButton("Load Reference Matrix")
        self.snp_btn_file_ref.clicked.connect(lambda: self._ui_generic_file_loader(self.comp_ref))
        h_s_1.addWidget(self.snp_btn_file_ref)
        h_s_1.addStretch()
        self.comp_ref = QTextEdit()
        self.comp_ref.setPlaceholderText("Origin Context Sequence (Reference)...")
        
        h_s_2 = QHBoxLayout()
        self.snp_btn_file_query = QPushButton("Load Query Segment")
        self.snp_btn_file_query.clicked.connect(lambda: self._ui_generic_file_loader(self.comp_query))
        h_s_2.addWidget(self.snp_btn_file_query)
        h_s_2.addStretch()
        self.comp_query = QTextEdit()
        self.comp_query.setPlaceholderText("Target Context Sequence (Query)...")
        
        h = QHBoxLayout()
        self.snp_btn = QPushButton("Global Align & Kabuto Mutation Predictor")
        self.snp_btn.clicked.connect(self.chidori_find_snps)
        self.dot_btn = QPushButton("Synthesize Matrix Dot Plot")
        self.dot_btn.clicked.connect(self.raikiri_dot_plot)
        h.addWidget(self.snp_btn)
        h.addWidget(self.dot_btn)
        
        sl.addLayout(h_s_1)
        sl.addWidget(self.comp_ref)
        sl.addLayout(h_s_2)
        sl.addWidget(self.comp_query)
        sl.addLayout(h)

        yg = QGroupBox("Yamato Wood Extractor & Density (GFF + FASTA)")
        yl = QVBoxLayout(yg)
        self.yamato_gff_in = QLineEdit()
        self.yamato_gff_in.setPlaceholderText("GFF/BED File Path")
        self.yamato_gff_btn = QPushButton("File")
        self.yamato_gff_btn.clicked.connect(lambda: self.yamato_gff_in.setText(QFileDialog.getOpenFileName(self, "Open GFF")[0]))
        h_y1 = QHBoxLayout()
        h_y1.addWidget(self.yamato_gff_in)
        h_y1.addWidget(self.yamato_gff_btn)

        self.yamato_fa_in = QLineEdit()
        self.yamato_fa_in.setPlaceholderText("Reference FASTA File Path")
        self.yamato_fa_btn = QPushButton("File")
        self.yamato_fa_btn.clicked.connect(lambda: self.yamato_fa_in.setText(QFileDialog.getOpenFileName(self, "Open FASTA")[0]))
        h_y2 = QHBoxLayout()
        h_y2.addWidget(self.yamato_fa_in)
        h_y2.addWidget(self.yamato_fa_btn)

        h_y3 = QHBoxLayout()
        self.yamato_gene = QLineEdit()
        self.yamato_gene.setPlaceholderText("Target Gene Name")
        self.yamato_run_btn = QPushButton("Extract Feature")
        self.yamato_run_btn.clicked.connect(self.yamato_extractor_clicked)
        self.yamato_dens_btn = QPushButton("Hashirama Density Plot")
        self.yamato_dens_btn.clicked.connect(self.hashirama_density_plot_clicked)
        h_y3.addWidget(self.yamato_gene)
        h_y3.addWidget(self.yamato_run_btn)
        h_y3.addWidget(self.yamato_dens_btn)

        yl.addLayout(h_y1)
        yl.addLayout(h_y2)
        yl.addLayout(h_y3)

        bg = QGroupBox("External BLAST+ & Evolution Trees")
        bl = QVBoxLayout(bg)
        h_b1 = QHBoxLayout()
        self.blast_exe_input = QLineEdit()
        self.blast_exe_input.setPlaceholderText("Path to blastn.exe or 'blastn'")
        self.blast_exe_btn = QPushButton("Locate Armory (Browse)")
        self.blast_exe_btn.clicked.connect(self.jiraiya_summon_blast_executable)
        h_b1.addWidget(QLabel("BLAST+ Executable:"))
        h_b1.addWidget(self.blast_exe_input)
        h_b1.addWidget(self.blast_exe_btn)
        h_b2 = QHBoxLayout()
        self.blast_btn = QPushButton("Deidara Explosive Local BLAST+")
        self.blast_btn.clicked.connect(self.deidara_blast_clicked)
        self.phylo_btn = QPushButton("Hashirama Mokuton UPGMA Tree")
        self.phylo_btn.clicked.connect(self.hashirama_tree_clicked)
        h_b2.addWidget(self.blast_btn)
        h_b2.addWidget(self.phylo_btn)
        
        bl.addLayout(h_b1)
        bl.addLayout(h_b2)
        
        ll.addWidget(sg)
        ll.addWidget(yg)
        ll.addWidget(bg)
        
        right = QVBoxLayout()
        self.comp_canvas = KonohaNinjaPlotCanvas(self)
        self.comp_stats = QTextEdit()
        self.comp_stats.setReadOnly(True)
        self.comp_stats.setFont(QFont("Consolas", 10))
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.comp_canvas)
        splitter.addWidget(self.comp_stats)
        right.addWidget(splitter)
        
        layout.addWidget(left, stretch=1)
        layout.addLayout(right, stretch=1)

    def kage_bunshin_build_vcf_page(self, page_widget):
        layout = QVBoxLayout(page_widget)
        h = QHBoxLayout()
        self.vcf_load_btn = QPushButton("Activate Byakugan (Load VCF File)")
        self.vcf_load_btn.clicked.connect(self.neji_load_vcf_clicked)
        self.vcf_search = QLineEdit()
        self.vcf_search.setPlaceholderText("Search by CHROM, POS, ID, REF, ALT...")
        self.vcf_search.textChanged.connect(self.neji_search_vcf)
        h.addWidget(self.vcf_load_btn)
        h.addWidget(self.vcf_search)
        
        self.vcf_table = QTableWidget()
        self.vcf_table.setColumnCount(7)
        self.vcf_table.setHorizontalHeaderLabels(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "INFO"])
        self.vcf_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.vcf_table.horizontalHeader().setStretchLastSection(True)
        
        layout.addLayout(h)
        layout.addWidget(self.vcf_table)

    def kage_bunshin_build_editor_page(self, page_widget):
        layout = QVBoxLayout(page_widget)
        btn_h = QHBoxLayout()
        self.ed_load_btn = QPushButton("Load Local Resource File")
        self.ed_load_btn.clicked.connect(self.sasuke_editor_load)
        self.ed_save_btn = QPushButton("Commit Disk Changes")
        self.ed_save_btn.clicked.connect(self.sasuke_editor_save)
        btn_h.addWidget(self.ed_load_btn)
        btn_h.addWidget(self.ed_save_btn)
        btn_h.addStretch()
        self.ed_path = QLineEdit()
        self.ed_path.setReadOnly(True)
        self.ed_text = QTextEdit()
        self.ed_text.setFont(QFont("Consolas", 11))
        layout.addLayout(btn_h)
        layout.addWidget(self.ed_path)
        layout.addWidget(self.ed_text)

    def kage_bunshin_build_help_page(self, page_widget):
        layout = QVBoxLayout(page_widget)
        help_text = QTextEdit()
        help_text.setReadOnly(True)
        help_text.setHtml(DOCUMENTATION_HTML)
        layout.addWidget(help_text)

    def _display_figure(self, canvas_widget, figure):
        canvas_widget.figure = figure
        canvas_widget.draw()
        self.generated_figures.append(figure)

    def shinra_tensei_clear_workspace(self):
        self.analyzer = UchihaSharinganAnalyzer()
        self.base_analyzer = SenjuMokutonAnalyzer()
        self.fastq_analyzer = UzumakiChakraFastqAnalyzer()
        self.file_paths_cache.clear()
        self.generated_figures.clear()
        self.file_list.clear()
        self.fastq_file_list.clear()
        self.seq_input.clear()
        self.stats_text.clear()
        self.base_stats_text.clear()
        self.fastq_results_text.clear()
        self.plot_tabs.clear()
        self.vcf_table.setRowCount(0)
        self.base_canvas.fig.clf()
        self.fastq_canvas.fig.clf()
        self.prot_canvas.fig.clf()
        self.comp_canvas.fig.clf()
        self.current_bars = None
        self.statusBar().showMessage("System memory block successfully purged.", 3000)

    def orochimaru_save_session(self):
        state = {
            "files": list(self.file_paths_cache),
            "plot_type": self.plot_type.currentIndex(),
            "genome_size": self.genome_size_input.value(),
            "read_length": self.read_length_input.value(),
            "paired_end": self.paired_end_check.isChecked()
        }
        try:
            with open(NINJA_SESSION_FILE, 'w') as f:
                json.dump(state, f)
            self.statusBar().showMessage(f"Operating session successfully archived.", 3000)
        except Exception as e: QMessageBox.warning(self, "I/O Access Error", f"Operation failed: {e}")

    def raijin_load_session(self):
        if not os.path.exists(NINJA_SESSION_FILE): return
        try:
            with open(NINJA_SESSION_FILE, 'r') as f:
                state = json.load(f)
            self.plot_type.setCurrentIndex(state.get("plot_type", 0))
            self.genome_size_input.setValue(state.get("genome_size", 1000000))
            self.read_length_input.setValue(state.get("read_length", 150))
            self.paired_end_check.setChecked(state.get("paired_end", True))
            for filepath in state.get("files", []):
                if os.path.exists(filepath): self._internal_load_file(filepath)
        except Exception: pass

    def shuriken_add_mass_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "System File Dialogue", "", "Sequence Files (*.fasta *.fa *.fna *.fastq *.fq *.gz *.gb *.vcf *.gff *.sam *.bed);;All (*)")
        for file_path in file_paths: self._internal_load_file(file_path)

    def shuriken_add_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Directory Targeting Utility")
        if folder:
            for file in os.listdir(folder):
                if file.endswith(('.fasta', '.fa', '.fna', '.fastq', '.fq', '.gz', '.gb', '.vcf', '.gff', '.sam', '.bed')):
                    self._internal_load_file(os.path.join(folder, file))

    def _internal_load_file(self, file_path: str):
        label = os.path.basename(file_path)
        try:
            fmt = neji_byakugan_format_sniffer(file_path)
            if fmt in ['fasta', 'genbank']:
                if self.analyzer.sasuke_sharingan_read_fasta(file_path, label):
                    item_text = f"{label} ({fmt.upper()})"
                    if self.file_list.findText(item_text) == -1:
                        self.file_list.addItem(item_text, userData=file_path)
                        self.file_paths_cache.add(file_path)
                        self.byakugan_update_stats()
            elif fmt == 'fastq':
                if self.fastq_analyzer.naruto_rasengan_read_fastq_blocks(file_path, label):
                    item_text = f"{label} (FASTQ)"
                    if self.fastq_file_list.findText(item_text) == -1:
                        self.fastq_file_list.addItem(item_text, userData=file_path)
                        self.file_paths_cache.add(file_path)
            elif fmt in ['vcf', 'gff', 'bed', 'sam']:
                res = self.multi_parser.tsunade_katsuyu_summarize_file(file_path, fmt)
                self.stats_text.append(res['text'] + "\n")
            self.statusBar().showMessage(f"Processing Complete: {label}", 2000)
        except ValueError as e: QMessageBox.warning(self, "Format Parser Warning", str(e))
            
    def kunai_remove_file(self):
        if self.file_list.count() > 0 and self.file_list.currentIndex() >= 0:
            idx = self.file_list.currentIndex()
            file_path = self.file_list.itemData(idx)
            label = self.file_list.itemText(idx).split(" (")[0]
            self.file_list.removeItem(idx)
            self.file_paths_cache.discard(file_path)
            if label in self.analyzer.raw_lengths:
                del self.analyzer.raw_lengths[label]
                if label in self.analyzer.file_paths: del self.analyzer.file_paths[label]
            self.byakugan_update_stats()

    def rinnegan_merge_files(self):
        if self.file_list.count() < 2: return QMessageBox.warning(self, "Operation Cancelled", "Requires multiple files.")
        files_to_merge = [self.file_list.itemData(i) for i in range(self.file_list.count())]
        out_path, _ = QFileDialog.getSaveFileName(self, "Output File Definition", "", "FASTA Files (*.fasta)")
        if out_path:
            try:
                pain_chibaku_tensei_merge_filter(files_to_merge, out_path, min_phred=20)
                self.statusBar().showMessage("File Assembly Sequence Executed", 4000)
            except Exception as e: QMessageBox.warning(self, "Write Error Exception", str(e))

    def byakugan_update_stats(self):
        min_l = self.min_len_spin.value()
        max_l = self.max_len_spin.value()
        self.analyzer.itachi_susanoo_apply_filters(min_l, max_l)
        self.stats_text.setPlainText(ino_mind_transfer_stats_text(self.analyzer.filtered_stats))

    def raikiri_update_plot(self):
        if not self.analyzer.filtered_lengths: return QMessageBox.information(self, "Empty Dataset", "Current filter yields void dataset block.")
        try:
            ptype = self.plot_type.currentText()
            if "Histogram" in ptype: fig = self.analyzer.sasuke_amaterasu_histogram(self.normal_check.isChecked())
            elif "Boxplot" in ptype: fig = self.analyzer.itachi_tsukuyomi_boxplot()
            elif "Violin" in ptype: fig = self.analyzer.danzo_izanagi_violinplot()
            elif "KDE" in ptype: fig = self.analyzer.madara_susanoo_kde()
            elif "Cumulative" in ptype: fig = self.analyzer.itachi_izanami_cumulative()
            else: return
            self.generated_figures.append(fig)
            canvas = FigureCanvas(fig)
            idx = self.plot_tabs.addTab(canvas, f"Analytical Graph {self.plot_tabs.count() + 1}")
            self.plot_tabs.setCurrentIndex(idx)
        except Exception as e: QMessageBox.warning(self, "Graphing Framework Exception", f"Render pipeline halted: {str(e)}")

    def susanoo_close_plot_tab(self, index):
        widget = self.plot_tabs.widget(index)
        if widget is not None: widget.deleteLater()
        self.plot_tabs.removeTab(index)

    def mangekyou_html_report(self):
        if not self.analyzer.filtered_stats: return QMessageBox.warning(self, "Constraint Error", "Lack of statistical parameters to build web report.")
        out_path, _ = QFileDialog.getSaveFileName(self, "Export Web Architecture Report", "", "HTML Files (*.html)")
        if out_path:
            try:
                stats_text = ino_mind_transfer_stats_text(self.analyzer.filtered_stats)
                figs = self.generated_figures if self.generated_figures else [self.analyzer.sasuke_amaterasu_histogram()]
                madara_infinite_tsukuyomi_html_report(stats_text, figs, out_path)
                self.statusBar().showMessage(f"Document generation complete", 4000)
            except Exception as e: QMessageBox.warning(self, "Internal Output Failure", str(e))

    def doton_export_csv(self):
        if not self.analyzer.filtered_lengths: return
        file_path, _ = QFileDialog.getSaveFileName(self, "Dump Matrix Data", "", "CSV Files (*.csv)")
        if file_path:
            try:
                tenten_scroll_export_csv(self.analyzer.filtered_lengths, file_path)
                self.statusBar().showMessage("Database export executed.", 4000)
            except Exception as e: QMessageBox.warning(self, "Export Engine Fault", str(e))

    def pain_fetch_ncbi_clicked(self):
        acc = self.ncbi_input.text().strip()
        if not acc: return
        self.ncbi_btn.setEnabled(False)
        self.ncbi_btn.setText("Pulling...")
        try:
            seq_fasta = PainBanshoTeninNCBIFetcher.universal_pull_sequence(acc)
            if seq_fasta:
                tmp_path = os.path.join(tempfile.gettempdir(), f"{acc}.fasta")
                with open(tmp_path, 'w') as f: f.write(seq_fasta)
                self._internal_load_file(tmp_path)
                QMessageBox.information(self, "Pain Bansho Tenin", f"Sequence successfully retrieved and loaded.")
        except Exception as e:
            QMessageBox.warning(self, "Network Jutsu Failed", str(e))
        finally:
            self.ncbi_btn.setEnabled(True)
            self.ncbi_btn.setText("Pull Sequence")

    def henge_format_convert_clicked(self):
        f_in = self.henge_in_path.text().strip()
        if not f_in: return
        out_path, _ = QFileDialog.getSaveFileName(self, "Save Transformed Output")
        if out_path:
            res = HengeFormatShifter.transformation_jutsu_convert(
                f_in, out_path, self.henge_combo_in.currentText(), self.henge_combo_out.currentText()
            )
            self.stats_text.append(res + "\n")
            QMessageBox.information(self, "Henge Jutsu", res)

    def chidori_analyze_sequence_async(self):
        sequence = self.seq_input.toPlainText().strip()
        if not sequence: return
        sequence = ''.join([c for c in sequence if c.isalpha()])
        self.analyze_seq_btn.setEnabled(False)
        self.worker_thread = KageBunshinAsyncJutsu(self.base_analyzer, sequence)
        self.worker_thread.finished_signal.connect(self._on_analysis_finished)
        self.worker_thread.start()

    def _on_analysis_finished(self):
        self.analyze_seq_btn.setEnabled(True)
        self.katon_update_base_plot(redraw=True)
        self.base_stats_text.setPlainText(self.base_analyzer.kakashi_doton_get_statistics_text())

    def byakugou_translate(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            prot = self.base_analyzer.tsunade_creation_rebirth_translate()
            self.base_stats_text.setPlainText("Translated Synthesized Protein Segment:\n" + prot)

    def shikamaru_revcomp_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.shikamaru_shadow_reversal_complement()
            self.seq_input.setPlainText(res)

    def jiraiya_transcribe_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.jiraiya_transcription_jutsu()
            self.seq_input.setPlainText(res)

    def sharingan_regex_highlight(self):
        pattern = self.regex_input.text()
        if not pattern: return
        cursor = self.seq_input.textCursor()
        cursor.select(QTextCursor.Document)
        cursor.setCharFormat(QTextCharFormat())
        format_hl = QTextCharFormat()
        format_hl.setBackground(QColor('#e74c3c'))
        format_hl.setForeground(QColor('#ffffff'))
        seq = self.seq_input.toPlainText()
        try:
            for m in re.finditer(pattern, seq, re.IGNORECASE):
                cursor.setPosition(m.start())
                cursor.setPosition(m.end(), QTextCursor.KeepAnchor)
                cursor.setCharFormat(format_hl)
        except Exception as e:
            QMessageBox.warning(self, "Invalid Regex Jutsu", str(e))

    def suiton_plot_kmers(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            try:
                fig, stats = self.base_analyzer.shino_kikaichu_kmer_swarm(self.kmer_spin.value())
                self._display_figure(self.base_canvas, fig)
                self.base_stats_text.setPlainText(stats)
            except Exception as e: QMessageBox.warning(self, "Heuristic Function Timeout", str(e))

    def kraken_lite_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.inuzuka_mind_dog_taxonomy()
            self.base_stats_text.setPlainText(res)

    def fuuton_plot_entropy(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            try:
                fig = self.base_analyzer.orochimaru_curse_mark_entropy_plot()
                self._display_figure(self.base_canvas, fig)
            except Exception as e: QMessageBox.warning(self, "Thermodynamic Model Error", str(e))

    def dust_complexity_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.gaara_dust_complexity_shield()
            self.seq_input.setPlainText(res)

    def raiton_plot_cpg(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            try:
                fig, stats = self.base_analyzer.hidan_jashin_cpg_islands_vectorized()
                self._display_figure(self.base_canvas, fig)
                self.base_stats_text.setPlainText(stats)
            except Exception as e: QMessageBox.warning(self, "Scanning Execution Issue", str(e))

    def zabuza_restriction_map_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            try:
                fig, stats = self.base_analyzer.zabuza_executioner_restriction_map()
                self._display_figure(self.base_canvas, fig)
                self.base_stats_text.setPlainText(stats)
            except Exception as e: QMessageBox.warning(self, "Restriction Engine Failure", str(e))

    def pai_islands_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.orochimaru_cursed_islands_pai()
            self.base_stats_text.setPlainText(res)

    def telomere_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.kimimaro_bone_telomere_mapper()
            self.base_stats_text.setPlainText(res)

    def rho_odds_clicked(self):
        seq = self.seq_input.toPlainText().strip()
        if seq:
            self.base_analyzer.yamato_mokuton_set_sequence(seq)
            res = self.base_analyzer.neji_trigram_rho_odds()
            self.base_stats_text.setPlainText(res)

    def katon_update_base_plot(self, redraw=False):
        if not self.base_analyzer.sequence: return
        try:
            if redraw or self.current_bars is None:
                self.base_canvas.fig.clf()
                ax = self.base_canvas.fig.add_subplot(111)
                sai_ninja_art_dark_theme(self.base_canvas.fig, ax)
                
                bases = [b for b in self.base_analyzer.base_counts if self.base_analyzer.base_counts[b] > 0]
                counts = [self.base_analyzer.base_counts[b] for b in bases]
                colors = [self.base_analyzer.base_colors[b] for b in bases]
                total = sum(counts) if sum(counts) > 0 else 1
                freqs = [count/total*100 for count in counts]
                
                self.current_bars = ax.bar(bases, freqs, color=colors)
                for bar, count in zip(self.current_bars, counts):
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height, f'{count}', ha='center', va='bottom', color='#ffd900')
                
                ax.set_title("Nucleotide Proportions Chart")
                ax.set_ylabel("Occurrence Percentage (%)")
                ax.set_ylim(0, 100)
                self.base_canvas.fig.tight_layout()
            self.base_canvas.draw()
            self.generated_figures.append(self.base_canvas.fig)
        except Exception: pass

    def katon_find_orfs(self):
        seq = self.orf_input.toPlainText().strip()
        if seq:
            res = RinneganAdvancedBioTools.jiraiya_sage_mode_orf_finder(seq)
            self.prot_stats.setPlainText(res)

    def suiton_protein_props(self):
        seq = self.prot_input.toPlainText().strip()
        if seq:
            res = RinneganAdvancedBioTools.raikage_lightning_armor_protein_properties(seq)
            self.prot_stats.setPlainText(res)

    def chouji_aa_plot_clicked(self):
        seq = self.prot_input.toPlainText().strip()
        if seq:
            try:
                fig = RinneganAdvancedBioTools.chouji_calorie_amino_plot(seq)
                self._display_figure(self.prot_canvas, fig)
            except Exception as e: QMessageBox.warning(self, "Plotting Pipeline Defect", str(e))

    def doton_hydrophobicity(self):
        seq = self.prot_input.toPlainText().strip()
        if seq:
            try:
                fig = RinneganAdvancedBioTools.kisame_water_prison_hydrophobicity_plot(seq)
                self._display_figure(self.prot_canvas, fig)
            except Exception as e: QMessageBox.warning(self, "Plotting Pipeline Defect", str(e))

    def raiton_primer_analysis(self):
        seq = self.primer_input.text().strip()
        if seq:
            res = RinneganAdvancedBioTools.minato_rasengan_primer_wizard(seq)
            self.prot_stats.setPlainText(res)

    def kakashi_crispr_clicked(self):
        seq = self.primer_input.text().strip()
        if seq:
            res = RinneganAdvancedBioTools.kakashi_crispr_copy_sgrna(seq)
            self.prot_stats.setPlainText(res)

    def gaara_trimmer_clicked(self):
        f_in = self.gaara_fastq_in.text().strip()
        if not f_in: return
        out_path, _ = QFileDialog.getSaveFileName(self, "Save Trimmed FASTQ", "", "FASTQ Files (*.fastq)")
        if out_path:
            res = GaaraSandTrimmer.sabaku_kyuu_trim_fastq(f_in, out_path, self.gaara_phred_spin.value())
            self.fastq_results_text.setPlainText(res)
            QMessageBox.information(self, "Gaara Sand Burial", res)

    def rasenshuriken_calculate_coverage(self):
        if not self.fastq_analyzer.read_counts: return
        try:
            self.fastq_analyzer.naruto_oodama_rasengan_calculate_genomes(
                self.genome_size_input.value(), self.read_length_input.value(), self.paired_end_check.isChecked())
            self.fastq_results_text.setPlainText(self.fastq_analyzer.uzumaki_fuuinjutsu_get_results_text())
        except ValueError as e: QMessageBox.warning(self, "Data Calculation Interruption", str(e))

    def mangekyou_plot_phred(self):
        try:
            fig = self.fastq_analyzer.itachi_tsukuyomi_phred_plot()
            self._display_figure(self.fastq_canvas, fig)
        except Exception as e: QMessageBox.warning(self, "Subsystem Failure", str(e))

    def rinnegan_plot_rarefaction(self):
        try:
            fig = self.fastq_analyzer.pain_shinra_tensei_rarefaction_curve()
            self._display_figure(self.fastq_canvas, fig)
        except Exception as e: QMessageBox.warning(self, "Simulation Subroutine Disruption", str(e))

    def kisame_saturation_clicked(self):
        try:
            fig = self.fastq_analyzer.kisame_water_shark_saturation()
            self._display_figure(self.fastq_canvas, fig)
        except Exception as e: QMessageBox.warning(self, "Saturation Error", str(e))

    def zetsu_duplication_clicked(self):
        res = self.fastq_analyzer.zetsu_spore_duplication_rate()
        self.fastq_results_text.setPlainText(res)

    def chidori_find_snps(self):
        ref = self.comp_ref.toPlainText().strip()
        query = self.comp_query.toPlainText().strip()
        if ref and query:
            res = RinneganAdvancedBioTools.sasuke_sharingan_snp_deducer_global(ref, query)
            self.comp_stats.setPlainText(res)

    def raikiri_dot_plot(self):
        ref = self.comp_ref.toPlainText().strip()
        query = self.comp_query.toPlainText().strip()
        if ref and query:
            try:
                fig = RinneganAdvancedBioTools.madara_rinnegan_dotplot_matrix_numpy(ref, query)
                self._display_figure(self.comp_canvas, fig)
            except Exception as e: QMessageBox.warning(self, "Matrix Array Fault", str(e))

    def yamato_extractor_clicked(self):
        gff = self.yamato_gff_in.text().strip()
        fasta = self.yamato_fa_in.text().strip()
        gene = self.yamato_gene.text().strip()
        if gff and fasta and gene:
            res = YamatoWoodFeatureExtractor.mokuton_extract_feature(gff, fasta, gene)
            self.comp_stats.setPlainText(res)

    def hashirama_density_plot_clicked(self):
        gff = self.yamato_gff_in.text().strip()
        if gff:
            try:
                fig = YamatoWoodFeatureExtractor.hashirama_wood_density_plot(gff)
                self._display_figure(self.comp_canvas, fig)
            except Exception as e: QMessageBox.warning(self, "Density Plot Error", str(e))

    def jiraiya_summon_blast_executable(self):
        path, _ = QFileDialog.getOpenFileName(self, "Summon BLAST+ Executable", "", "Executables (*.exe);;All Files (*)")
        if path:
            self.blast_exe_input.setText(path)

    def deidara_blast_clicked(self):
        blast_path = self.blast_exe_input.text().strip()
        if not blast_path: blast_path = "blastn"
        q_seq = self.comp_ref.toPlainText().strip()
        s_seq = self.comp_query.toPlainText().strip()
        if q_seq and s_seq:
            self.blast_btn.setEnabled(False)
            res = AnbuBlackOpsExternalTools.deidara_explosive_blast_art(blast_path, q_seq, s_seq)
            self.comp_stats.setPlainText(res)
            self.blast_btn.setEnabled(True)

    def hashirama_tree_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load FASTA Alignment (MSA)", "", "FASTA Files (*.fasta *.fa *.aln)")
        if path:
            try:
                fig = RinneganAdvancedBioTools.hashirama_phylo_tree_builder(path)
                self._display_figure(self.comp_canvas, fig)
            except Exception as e:
                QMessageBox.warning(self, "Mokuton Core Error", str(e))

    def neji_load_vcf_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, "Activate Byakugan (Load VCF)", "", "VCF Files (*.vcf *.vcf.gz)")
        if path:
            try:
                self.vcf_data = NejiVariantEyeParser.eight_trigrams_parse_vcf(path)
                self._populate_vcf_table(self.vcf_data)
            except Exception as e:
                QMessageBox.warning(self, "Byakugan Failed", str(e))

    def _populate_vcf_table(self, data):
        self.vcf_table.setRowCount(0)
        limit = min(1000, len(data))
        for i in range(limit):
            self.vcf_table.insertRow(i)
            v = data[i]
            for j, key in enumerate(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "INFO"]):
                item = QTableWidgetItem(v.get(key, ''))
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                self.vcf_table.setItem(i, j, item)
        if len(data) > 1000:
            self.statusBar().showMessage(f"Loaded {len(data)} variants. Displaying top 1000 to save memory.", 5000)
        else:
            self.statusBar().showMessage(f"Loaded {len(data)} variants successfully.", 5000)

    def neji_search_vcf(self, query):
        if not hasattr(self, 'vcf_data') or not self.vcf_data: return
        if not query:
            self._populate_vcf_table(self.vcf_data)
            return
        query_lower = query.lower()
        filtered = [v for v in self.vcf_data if any(query_lower in str(val).lower() for val in v.values())]
        self._populate_vcf_table(filtered)

    def sasuke_editor_load(self):
        path, _ = QFileDialog.getOpenFileName(self, "Locate Editable Resource")
        if path:
            try:
                open_func = gzip.open if path.endswith('.gz') else open
                mode = 'rt' if path.endswith('.gz') else 'r'
                with open_func(path, mode) as f: self.ed_text.setPlainText(f.read())
                self.ed_path.setText(path)
            except Exception as e: QMessageBox.warning(self, "Disc Sector Reading Issue", str(e))

    def sasuke_editor_save(self):
        path = self.ed_path.text()
        if not path: return
        try:
            open_func = gzip.open if path.endswith('.gz') else open
            mode = 'wt' if path.endswith('.gz') else 'w'
            with open_func(path, mode) as f: f.write(self.ed_text.toPlainText())
            self.statusBar().showMessage("File state properly committed to disk.", 3000)
        except Exception as e: QMessageBox.warning(self, "Write Permissions Denied", str(e))

def konoha_hokage_office_main():
    app = QApplication(sys.argv)
    
    app.setStyle("Fusion")
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.Window, QColor(18, 18, 18))
    dark_palette.setColor(QPalette.WindowText, Qt.white)
    dark_palette.setColor(QPalette.Base, QColor(26, 26, 26))
    dark_palette.setColor(QPalette.AlternateBase, QColor(18, 18, 18))
    dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
    dark_palette.setColor(QPalette.ToolTipText, Qt.white)
    dark_palette.setColor(QPalette.Text, Qt.white)
    dark_palette.setColor(QPalette.Button, QColor(18, 18, 18))
    dark_palette.setColor(QPalette.ButtonText, Qt.white)
    dark_palette.setColor(QPalette.BrightText, Qt.red)
    dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(dark_palette)
    
    window = HokageTowerDashboardGUI()
    window.setStyleSheet(NARUTO_HOKAGE_STYLE)
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    konoha_hokage_office_main()
