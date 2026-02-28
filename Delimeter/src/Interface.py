from PyQt5.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLabel, QComboBox, QLineEdit, QPushButton, QCheckBox, QTextEdit,
                             QFileDialog, QTableView, QMessageBox, QGroupBox, QProgressBar,
                             QScrollArea, QFrame, QStackedWidget, QDialog)
from PyQt5.QtCore import Qt, pyqtSignal, QTimer
from PyQt5.QtGui import QCursor, QColor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import pandas as pd
from typing import List, Dict, Any, Optional

from logic import Regex_Sentinel, QuantumOperators, TransmutationSpells, RadiationMethods
from xenoglossy_codex import Rosetta_Stone_Ribosome

class Spectrometry_Canvas(FigureCanvasQTAgg):
    def __init__(self, parent: Optional[QWidget] = None, width: int = 6, height: int = 4, dpi: int = 100) -> None:
        spectral_figure = Figure(figsize=(width, height), dpi=dpi)
        self.axes = spectral_figure.add_subplot(111)
        super(Spectrometry_Canvas, self).__init__(spectral_figure)

class Spectrometry_Dialog(QDialog):
    def __init__(self, spectral_dataframe: pd.DataFrame, target_column: str) -> None:
        super().__init__()
        self.setWindowTitle(f"Visual Spectrometry: {target_column}")
        self.setMinimumSize(700, 500)
        main_layout = QVBoxLayout(self)
        self.canvas = Spectrometry_Canvas(self)
        main_layout.addWidget(self.canvas)
        
        nucleotide_counts = spectral_dataframe[target_column].value_counts().head(15)
        nucleotide_counts.plot(kind='bar', ax=self.canvas.axes, color='#bd93f9')
        self.canvas.axes.set_title(f"Frequency Distribution: {target_column}")
        self.canvas.axes.set_ylabel("Absolute Count")
        self.canvas.axes.tick_params(axis='x', rotation=30)
        self.canvas.figure.tight_layout()

class Xenomorph_Splicer_Dialog(QDialog):
    def __init__(self, local_columns: List[str]) -> None:
        super().__init__()
        self.setWindowTitle("Alien DNA Splicing")
        self.setMinimumSize(400, 200)
        main_layout = QVBoxLayout(self)
        
        self.path_alien = QLineEdit()
        self.path_alien.setPlaceholderText("Alien archive path...")
        btn_browse = QPushButton("Search DNA")
        btn_browse.clicked.connect(self.browse_alien_genome)
        
        row_first = QHBoxLayout()
        row_first.addWidget(self.path_alien); row_first.addWidget(btn_browse)
        
        self.combo_local = QComboBox(); self.combo_local.addItems(local_columns)
        self.line_alien_col = QLineEdit(); self.line_alien_col.setPlaceholderText("Key-Column in alien archive")
        self.btn_splice = QPushButton("Fuse Genomes"); self.btn_splice.setObjectName("WarningActionBtn")
        
        main_layout.addLayout(row_first)
        main_layout.addWidget(QLabel("Local Column (Key):")); main_layout.addWidget(self.combo_local)
        main_layout.addWidget(QLabel("Alien Column (Key):")); main_layout.addWidget(self.line_alien_col)
        main_layout.addSpacing(10); main_layout.addWidget(self.btn_splice)

    def browse_alien_genome(self) -> None:
        foreign_path, _ = QFileDialog.getOpenFileName(self, "Select Foreign Genome", "", "CSV/TSV (*.csv *.tsv)")
        if foreign_path: self.path_alien.setText(foreign_path)

class Fasta_Export_Dialog(QDialog):
    def __init__(self, available_columns: List[str]) -> None:
        super().__init__()
        self.setWindowTitle("FASTA Synthesizer")
        main_layout = QVBoxLayout(self)
        self.combo_id = QComboBox(); self.combo_id.addItems(available_columns)
        self.combo_seq = QComboBox(); self.combo_seq.addItems(available_columns)
        self.btn_export = QPushButton("Compile FASTA")
        main_layout.addWidget(QLabel("Identifier Column (ID):")); main_layout.addWidget(self.combo_id)
        main_layout.addWidget(QLabel("Sequence Column (SEQ):")); main_layout.addWidget(self.combo_seq)
        main_layout.addWidget(self.btn_export)

class Transmutation_Altar_Dialog(QDialog):
    def __init__(self, available_columns: List[str]) -> None:
        super().__init__()
        self.setWindowTitle("Transmutation Altar")
        main_layout = QVBoxLayout(self)
        self.line_new = QLineEdit(); self.line_new.setPlaceholderText("New Mutation Name")
        self.combo_c1 = QComboBox(); self.combo_c1.addItems(available_columns)
        self.combo_op = QComboBox(); self.combo_op.addItems([spell.value for spell in TransmutationSpells])
        self.combo_c2 = QComboBox(); self.combo_c2.addItems(available_columns)
        self.btn_transmute = QPushButton("Transmute")
        
        main_layout.addWidget(QLabel("New Strand:")); main_layout.addWidget(self.line_new)
        main_layout.addWidget(self.combo_c1); main_layout.addWidget(self.combo_op); main_layout.addWidget(self.combo_c2)
        main_layout.addWidget(self.btn_transmute)

class Radiation_Dialog(QDialog):
    def __init__(self, available_columns: List[str]) -> None:
        super().__init__()
        self.setWindowTitle("N-Terminal Radiation")
        main_layout = QVBoxLayout(self)
        self.combo_col = QComboBox(); self.combo_col.addItems(available_columns)
        self.combo_method = QComboBox(); self.combo_method.addItems([method.value for method in RadiationMethods])
        self.line_fixed = QLineEdit(); self.line_fixed.setPlaceholderText("Fixed Value")
        self.btn_radiate = QPushButton("Radiate")
        main_layout.addWidget(QLabel("Biological Target:")); main_layout.addWidget(self.combo_col)
        main_layout.addWidget(QLabel("Method:")); main_layout.addWidget(self.combo_method)
        main_layout.addWidget(self.line_fixed); main_layout.addWidget(self.btn_radiate)

class Anomaly_Marker_Dialog(QDialog):
    def __init__(self, available_columns: List[str]) -> None:
        super().__init__()
        self.setWindowTitle("Anomaly Marker")
        main_layout = QVBoxLayout(self)
        self.combo_col = QComboBox(); self.combo_col.addItem("ALL_COLUMNS"); self.combo_col.addItems(available_columns)
        self.combo_op = QComboBox(); self.combo_op.addItems([QuantumOperators.EQUALS.value, QuantumOperators.CONTAINS.value, QuantumOperators.GREATER_THAN.value, QuantumOperators.LESS_THAN.value])
        self.line_val = QLineEdit(); self.line_val.setPlaceholderText("Activator Value")
        self.combo_color = QComboBox(); self.combo_color.addItems(["Critical Red", "Phosphorescent Green", "Radiation Yellow"])
        self.btn_mark = QPushButton("Add Visible Marker")
        main_layout.addWidget(self.combo_col); main_layout.addWidget(self.combo_op); main_layout.addWidget(self.line_val)
        main_layout.addWidget(self.combo_color); main_layout.addWidget(self.btn_mark)

class Warp_Gate_Card(QFrame):
    signal_warp_triggered = pyqtSignal()

    def __init__(self, icon: str, title_key: str, subtitle_key: str) -> None:
        super().__init__()
        self.setObjectName("ActionCard")
        self.setCursor(QCursor(Qt.PointingHandCursor))
        self.title_key = title_key
        self.subtitle_key = subtitle_key
        card_layout = QHBoxLayout(self)
        
        self.lbl_icon = QLabel(icon); self.lbl_icon.setObjectName("CardIcon"); self.lbl_icon.setFixedSize(40, 40); self.lbl_icon.setAlignment(Qt.AlignCenter)
        self.lbl_title = QLabel(""); self.lbl_title.setObjectName("CardTitle")
        self.lbl_subtitle = QLabel(""); self.lbl_subtitle.setObjectName("CardSubtitle")
        self.lbl_arrow = QLabel("›"); self.lbl_arrow.setObjectName("CardArrow")
        
        self.lbl_icon.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.lbl_title.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.lbl_subtitle.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.lbl_arrow.setAttribute(Qt.WA_TransparentForMouseEvents)

        text_layout = QVBoxLayout()
        text_layout.addWidget(self.lbl_title); text_layout.addWidget(self.lbl_subtitle); text_layout.setAlignment(Qt.AlignVCenter)
        card_layout.addWidget(self.lbl_icon); card_layout.addLayout(text_layout); card_layout.addStretch(); card_layout.addWidget(self.lbl_arrow)

    def mousePressEvent(self, event: Any) -> None:
        if event.button() == Qt.LeftButton:
            self.signal_warp_triggered.emit()

    def translate_card(self) -> None:
        self.lbl_title.setText(Rosetta_Stone_Ribosome.extract_peptide(self.title_key))
        if self.subtitle_key:
            self.lbl_subtitle.setText(Rosetta_Stone_Ribosome.extract_peptide(self.subtitle_key))

class Matryoshka_Filter_Pod(QFrame):
    signal_self_destruct = pyqtSignal(object)

    def __init__(self, available_columns: List[str]) -> None:
        super().__init__()
        self.setObjectName("FilterPod")
        self.main_layout = QHBoxLayout(self); self.main_layout.setContentsMargins(10, 10, 10, 10)

        self.combo_col = QComboBox(); self.combo_col.addItems(available_columns)
        self.combo_op = QComboBox(); self.combo_op.addItems([op.value for op in QuantumOperators])
        
        self.line_value = QLineEdit()
        
        self.neural_timer = QTimer()
        self.neural_timer.setSingleShot(True)
        self.neural_timer.timeout.connect(self.scan_for_anomalies)
        self.line_value.textChanged.connect(self.intercept_neural_feedback)
        
        self.combo_grimoire = QComboBox()
        self.combo_grimoire.addItem("Regex Grimoire")
        self.combo_grimoire.addItems(list(self.summon_regex_grimoire().keys()))
        self.combo_grimoire.currentIndexChanged.connect(self.cast_regex_spell)
        
        self.check_case = QCheckBox("Case"); self.check_ignore_nan = QCheckBox("Ignore NaN"); self.check_only_nan = QCheckBox("Only NaN")
        self.btn_trash = QPushButton("🗑️"); self.btn_trash.setObjectName("DangerButton"); self.btn_trash.setFixedWidth(45)
        self.btn_trash.clicked.connect(lambda: self.signal_self_destruct.emit(self))

        self.main_layout.addWidget(self.combo_col); self.main_layout.addWidget(self.combo_op); self.main_layout.addWidget(self.combo_grimoire)
        self.main_layout.addWidget(self.line_value); self.main_layout.addWidget(self.check_case); self.main_layout.addWidget(self.check_ignore_nan)
        self.main_layout.addWidget(self.check_only_nan); self.main_layout.addWidget(self.btn_trash)
        self.translate_pod()

    def translate_pod(self) -> None:
        self.line_value.setPlaceholderText(Rosetta_Stone_Ribosome.extract_peptide("target_placeholder"))

    def summon_regex_grimoire(self) -> Dict[str, str]:
        return {
            "Start Codon (ATG)": "^ATG",
            "TATA Box (TATAAA)": "TATAAA",
            "EcoRI (GAATTC)": "GAATTC",
            "Email Address": r"[^@]+@[^@]+\.[^@]+",
            "Patient CPF": r"\d{3}\.\d{3}\.\d{3}-\d{2}"
        }

    def intercept_neural_feedback(self) -> None:
        self.neural_timer.start(500)

    def cast_regex_spell(self) -> None:
        ancient_spell = self.combo_grimoire.currentText()
        mystic_grimoire = self.summon_regex_grimoire()
        if ancient_spell in mystic_grimoire:
            self.combo_op.setCurrentText(QuantumOperators.REGEX.value)
            self.line_value.setText(mystic_grimoire[ancient_spell])
        else:
            self.line_value.clear()

    def scan_for_anomalies(self) -> None:
        if self.combo_op.currentText() == QuantumOperators.REGEX.value:
            if Regex_Sentinel.inspect_anomaly(self.line_value.text()):
                self.line_value.setStyleSheet("border: 2px solid #50fa7b;")
            else:
                self.line_value.setStyleSheet("border: 2px solid #ff5555;")
        else:
            self.line_value.setStyleSheet("")

    def extract_dna_sequence(self) -> Dict[str, Any]:
        return {
            "column": self.combo_col.currentText(),
            "operation": self.combo_op.currentText(),
            "value": self.line_value.text(),
            "case_sensitive": self.check_case.isChecked(),
            "ignore_nan": self.check_ignore_nan.isChecked(),
            "only_nan": self.check_only_nan.isChecked()
        }
    
    def load_dna_sequence(self, stored_condition: Dict[str, Any]) -> None:
        self.combo_col.setCurrentText(stored_condition['column']); self.combo_op.setCurrentText(stored_condition['operation'])
        self.line_value.setText(stored_condition['value']); self.check_case.setChecked(stored_condition['case_sensitive'])
        self.check_ignore_nan.setChecked(stored_condition['ignore_nan']); self.check_only_nan.setChecked(stored_condition['only_nan'])

class Genomic_Holodeck_UI(QMainWindow):
    signal_load_biomass = pyqtSignal(str, str, bool)
    signal_ignite_filters = pyqtSignal(list)
    signal_export_genome = pyqtSignal(str)
    signal_toggle_spectrum = pyqtSignal()
    signal_dialect_mutation = pyqtSignal(str)

    def __init__(self, processing_engine: Any = None) -> None:
        super().__init__()
        self.processing_engine = processing_engine
        self.setWindowTitle("Delimiter File Filter")
        self.setGeometry(100, 100, 1380, 850)
        self.setAcceptDrops(True); self.setMenuBar(None)
        self.active_columns: List[str] = []; self.matryoshka_pods: List[Matryoshka_Filter_Pod] = []

        self.core_widget = QWidget(); self.core_widget.setObjectName("CoreMatrix")
        self.setCentralWidget(self.core_widget)
        self.master_layout = QHBoxLayout(self.core_widget); self.master_layout.setContentsMargins(0, 0, 0, 0); self.master_layout.setSpacing(0)
        self.cards: List[Warp_Gate_Card] = []

        self.forge_quantum_navigation_array(); self.forge_dimensional_stack(); self.construct_anxiety_bar()
        self.retranscribe_ui() 

    def forge_quantum_navigation_array(self) -> None:
        self.sidebar = QFrame(); self.sidebar.setObjectName("SidebarMatrix"); self.sidebar.setFixedWidth(240)
        sidebar_layout = QVBoxLayout(self.sidebar); sidebar_layout.setContentsMargins(10, 20, 10, 20); sidebar_layout.setSpacing(10)

        self.logo_label = QLabel(""); self.logo_label.setObjectName("SidebarLogo"); self.logo_label.setAlignment(Qt.AlignCenter)
        sidebar_layout.addWidget(self.logo_label); sidebar_layout.addSpacing(15)

        self.combo_dialect = QComboBox(); self.combo_dialect.addItems(["EN", "PT-BR", "ES"])
        self.combo_dialect.currentTextChanged.connect(self.signal_dialect_mutation.emit)
        sidebar_layout.addWidget(self.combo_dialect); sidebar_layout.addSpacing(15)

        self.btn_nav_home = QPushButton(""); self.btn_nav_home.setObjectName("SidebarButton"); self.btn_nav_home.clicked.connect(lambda: self.dimension_stack.setCurrentIndex(0))
        self.btn_nav_lab = QPushButton(""); self.btn_nav_lab.setObjectName("SidebarButton"); self.btn_nav_lab.clicked.connect(lambda: self.dimension_stack.setCurrentIndex(1))
        self.btn_nav_manual = QPushButton(""); self.btn_nav_manual.setObjectName("SidebarButton"); self.btn_nav_manual.clicked.connect(lambda: self.dimension_stack.setCurrentIndex(2))
        
        sidebar_layout.addWidget(self.btn_nav_home); sidebar_layout.addWidget(self.btn_nav_lab); sidebar_layout.addWidget(self.btn_nav_manual)
        
        sidebar_layout.addSpacing(30); self.lbl_microscopy = QLabel()
        sidebar_layout.addWidget(self.lbl_microscopy)
        
        self.combo_deep_scan = QComboBox()
        self.text_deep_scan = QTextEdit(); self.text_deep_scan.setReadOnly(True); self.text_deep_scan.setMaximumHeight(200)
        self.text_deep_scan.setObjectName("CleanScroll")
        
        sidebar_layout.addWidget(self.combo_deep_scan); sidebar_layout.addWidget(self.text_deep_scan)
        sidebar_layout.addStretch()
        
        self.btn_theme_toggle = QPushButton(""); self.btn_theme_toggle.setObjectName("SidebarButton"); self.btn_theme_toggle.clicked.connect(self.signal_toggle_spectrum.emit)
        sidebar_layout.addWidget(self.btn_theme_toggle)
        self.master_layout.addWidget(self.sidebar)

    def retranscribe_ui(self) -> None:
        self.logo_label.setText(Rosetta_Stone_Ribosome.extract_peptide("sidebar_logo"))
        self.btn_nav_home.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_nav_home"))
        self.btn_nav_lab.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_nav_lab"))
        self.btn_nav_manual.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_nav_manual"))
        self.lbl_microscopy.setText(Rosetta_Stone_Ribosome.extract_peptide("lbl_microscopy"))
        

        is_dark = "Dark" not in self.btn_theme_toggle.text() and "Escura" not in self.btn_theme_toggle.text() and "Oscura" not in self.btn_theme_toggle.text()
        self.mutate_spectrum_button_icon(is_dark)

        self.lbl_hero_title.setText(Rosetta_Stone_Ribosome.extract_peptide("hero_title"))
        self.lbl_hero_desc.setText(Rosetta_Stone_Ribosome.extract_peptide("hero_desc"))
        self.btn_jump.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_jump"))

        for card in self.cards: card.translate_card()

        self.path_display.setPlaceholderText(Rosetta_Stone_Ribosome.extract_peptide("path_placeholder"))
        self.btn_explore.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_explore"))
        self.combo_format.setItemText(0, Rosetta_Stone_Ribosome.extract_peptide("combo_auto"))
        self.check_header.setText(Rosetta_Stone_Ribosome.extract_peptide("check_header"))
        self.btn_inject.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_inject"))
        self.btn_undo.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_undo"))
        self.btn_redo.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_redo"))

        self.group_reactor.setTitle(Rosetta_Stone_Ribosome.extract_peptide("group_reactor"))
        self.btn_add_mutation.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_add_rule"))
        self.btn_save_presets.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_save_presets"))
        self.btn_load_presets.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_load_presets"))
        self.btn_collide.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_collide"))

        self.btn_purge.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_purge"))
        self.btn_splice.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_splice"))
        self.btn_radiate.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_radiate"))
        self.btn_transmute.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_transmute"))
        self.btn_spectro.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_spectro"))
        self.btn_anomaly.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_anomaly"))

        self.group_xray.setTitle(Rosetta_Stone_Ribosome.extract_peptide("group_xray"))
        self.btn_export_csv.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_export_csv"))
        self.btn_export_tsv.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_export_tsv"))
        self.btn_export_xlsx.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_export_xlsx"))
        self.btn_export_fasta.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_export_fasta"))
        
        self.manual_text.setHtml(Rosetta_Stone_Ribosome.extract_peptide("manual_html"))
        
        for pod in self.matryoshka_pods: pod.translate_pod()
        
        self.status_bar.showMessage(Rosetta_Stone_Ribosome.extract_peptide("status_awaiting"))

    def mutate_spectrum_button_icon(self, is_dark_active: bool) -> None:
        if is_dark_active: self.btn_theme_toggle.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_theme_light"))
        else: self.btn_theme_toggle.setText(Rosetta_Stone_Ribosome.extract_peptide("btn_theme_dark"))

    def forge_dimensional_stack(self) -> None:
        self.dimension_stack = QStackedWidget(); self.dimension_stack.setObjectName("DimensionStack")
        self.page_home = QWidget(); self.page_lab = QWidget(); self.page_manual = QWidget()
        self.build_main_bridge(); self.build_laboratory(); self.build_manual()
        self.dimension_stack.addWidget(self.page_home); self.dimension_stack.addWidget(self.page_lab); self.dimension_stack.addWidget(self.page_manual)
        self.master_layout.addWidget(self.dimension_stack)

    def build_main_bridge(self) -> None:
        main_layout = QVBoxLayout(self.page_home); main_layout.setContentsMargins(40, 40, 40, 40); main_layout.setSpacing(20)

        hero_frame = QFrame(); hero_frame.setObjectName("SupernovaBanner")
        hero_layout = QVBoxLayout(hero_frame); hero_layout.setContentsMargins(40, 40, 40, 40)
        self.lbl_hero_title = QLabel(""); self.lbl_hero_title.setObjectName("HeroTitle")
        self.lbl_hero_desc = QLabel(""); self.lbl_hero_desc.setObjectName("HeroDesc")
        self.btn_jump = QPushButton(""); self.btn_jump.setObjectName("PrimaryActionBtn"); self.btn_jump.setFixedWidth(200); self.btn_jump.clicked.connect(lambda: self.dimension_stack.setCurrentIndex(1))

        hero_layout.addWidget(self.lbl_hero_title); hero_layout.addWidget(self.lbl_hero_desc); hero_layout.addSpacing(20); hero_layout.addWidget(self.btn_jump); hero_layout.addStretch()
        main_layout.addWidget(hero_frame, stretch=3)

        cards_layout = QHBoxLayout(); cards_layout.setSpacing(15)
        card1 = Warp_Gate_Card("📁", "card1_title", None); card1.signal_warp_triggered.connect(lambda: self.dimension_stack.setCurrentIndex(1))
        card2 = Warp_Gate_Card("⚡", "card2_title", "card2_sub"); card2.signal_warp_triggered.connect(lambda: self.dimension_stack.setCurrentIndex(1))
        card3 = Warp_Gate_Card("🎨", "card3_title", "card3_sub"); card3.signal_warp_triggered.connect(lambda: self.dimension_stack.setCurrentIndex(1))
        card4 = Warp_Gate_Card("📖", "card4_title", None); card4.signal_warp_triggered.connect(lambda: self.dimension_stack.setCurrentIndex(2))
        
        self.cards.extend([card1, card2, card3, card4])
        for card in self.cards: cards_layout.addWidget(card)
        main_layout.addLayout(cards_layout, stretch=1)

    def build_laboratory(self) -> None:
        lab_layout = QVBoxLayout(self.page_lab); lab_layout.setContentsMargins(20, 20, 20, 20); lab_layout.setSpacing(10)

        top_bar = QHBoxLayout()
        self.path_display = QLineEdit(); self.path_display.setReadOnly(True)
        self.btn_explore = QPushButton(""); self.btn_explore.clicked.connect(self.summon_explorer)
        self.combo_format = QComboBox(); self.combo_format.addItems(["Auto Detect", "CSV", "TSV", "Excel"])
        self.check_header = QCheckBox(""); self.check_header.setChecked(True)
        self.btn_inject = QPushButton(""); self.btn_inject.setObjectName("PrimaryActionBtn"); self.btn_inject.clicked.connect(self.trigger_upload_sequence)
        
        self.btn_undo = QPushButton(""); self.btn_redo = QPushButton("")
        
        top_bar.addWidget(self.path_display); top_bar.addWidget(self.btn_explore); top_bar.addWidget(self.combo_format)
        top_bar.addWidget(self.check_header); top_bar.addWidget(self.btn_inject)
        top_bar.addSpacing(20); top_bar.addWidget(self.btn_undo); top_bar.addWidget(self.btn_redo)
        lab_layout.addLayout(top_bar)

        self.group_reactor = QGroupBox()
        reactor_layout = QVBoxLayout()
        self.scroll_chamber = QScrollArea(); self.scroll_chamber.setWidgetResizable(True); self.scroll_chamber.setObjectName("CleanScroll"); self.scroll_chamber.setMaximumHeight(180)
        self.pod_container = QWidget(); self.pod_layout = QVBoxLayout(self.pod_container); self.pod_layout.setAlignment(Qt.AlignTop)
        self.scroll_chamber.setWidget(self.pod_container)
        
        reactor_buttons = QHBoxLayout()
        self.btn_add_mutation = QPushButton(""); self.btn_add_mutation.clicked.connect(self.spawn_matryoshka_pod)
        self.btn_save_presets = QPushButton(""); self.btn_load_presets = QPushButton("")
        self.btn_collide = QPushButton(""); self.btn_collide.setObjectName("WarningActionBtn"); self.btn_collide.clicked.connect(self.trigger_collision)

        reactor_buttons.addWidget(self.btn_add_mutation); reactor_buttons.addWidget(self.btn_save_presets); reactor_buttons.addWidget(self.btn_load_presets)
        reactor_buttons.addStretch(); reactor_buttons.addWidget(self.btn_collide)
        reactor_layout.addWidget(self.scroll_chamber); reactor_layout.addLayout(reactor_buttons)
        self.group_reactor.setLayout(reactor_layout)
        lab_layout.addWidget(self.group_reactor)

        advanced_tools_layout = QHBoxLayout()
        self.btn_purge = QPushButton(""); self.btn_splice = QPushButton(""); self.btn_radiate = QPushButton("")
        self.btn_transmute = QPushButton(""); self.btn_spectro = QPushButton(""); self.btn_anomaly = QPushButton("")
        advanced_tools_layout.addWidget(self.btn_purge); advanced_tools_layout.addWidget(self.btn_splice); advanced_tools_layout.addWidget(self.btn_radiate)
        advanced_tools_layout.addWidget(self.btn_transmute); advanced_tools_layout.addWidget(self.btn_spectro); advanced_tools_layout.addWidget(self.btn_anomaly)
        lab_layout.addLayout(advanced_tools_layout)

        hologram_layout = QHBoxLayout()
        self.infinite_grid = QTableView(); self.infinite_grid.setObjectName("DataGrid")
        hologram_layout.addWidget(self.infinite_grid, stretch=5)
        
        self.group_xray = QGroupBox()
        xray_layout = QVBoxLayout()
        self.lbl_raw_reads = QLabel("Raw Reads: 0"); self.lbl_surviving = QLabel("Survivors: 0"); self.lbl_survival_rate = QLabel("Rate: 0%")
        xray_layout.addWidget(self.lbl_raw_reads); xray_layout.addWidget(self.lbl_surviving); xray_layout.addWidget(self.lbl_survival_rate); xray_layout.addStretch()
        
        self.btn_export_csv = QPushButton(""); self.btn_export_csv.clicked.connect(lambda: self.signal_export_genome.emit("csv"))
        self.btn_export_tsv = QPushButton(""); self.btn_export_tsv.clicked.connect(lambda: self.signal_export_genome.emit("tsv"))
        self.btn_export_xlsx = QPushButton(""); self.btn_export_xlsx.clicked.connect(lambda: self.signal_export_genome.emit("excel"))
        self.btn_export_fasta = QPushButton(""); self.btn_export_fasta.setObjectName("WarningActionBtn")
        
        xray_layout.addWidget(self.btn_export_csv); xray_layout.addWidget(self.btn_export_tsv); xray_layout.addWidget(self.btn_export_xlsx); xray_layout.addWidget(self.btn_export_fasta)
        self.group_xray.setLayout(xray_layout)
        hologram_layout.addWidget(self.group_xray, stretch=1)
        
        lab_layout.addLayout(hologram_layout)

    def build_manual(self) -> None:
        main_layout = QVBoxLayout(self.page_manual); main_layout.setContentsMargins(40, 40, 40, 40)
        self.manual_text = QTextEdit(); self.manual_text.setReadOnly(True); self.manual_text.setObjectName("CleanScroll")
        main_layout.addWidget(self.manual_text)

    def construct_anxiety_bar(self) -> None:
        self.status_bar = self.statusBar(); self.status_bar.setObjectName("StatusBarMatrix")
        self.progress_bar = QProgressBar(); self.progress_bar.setMaximumWidth(300); self.progress_bar.setVisible(False)
        self.status_bar.addPermanentWidget(self.progress_bar)

    def summon_explorer(self) -> None:
        genomic_path, _ = QFileDialog.getOpenFileName(self, "Open Sequence", "", "Sequences (*.csv *.tsv *.xlsx *.xls);;CSV (*.csv);;TSV (*.tsv)")
        if genomic_path: self.path_display.setText(genomic_path)

    def dragEnterEvent(self, event: Any) -> None:
        if event.mimeData().hasUrls(): event.accept()
        else: event.ignore()

    def dropEvent(self, event: Any) -> None:
        injected_files = [u.toLocalFile() for u in event.mimeData().urls()]
        if injected_files:
            self.path_display.setText(injected_files[0])
            self.dimension_stack.setCurrentIndex(1)

    def trigger_upload_sequence(self) -> None:
        target_path = self.path_display.text()
        if not target_path:
            QMessageBox.warning(self, "Alert", Rosetta_Stone_Ribosome.extract_peptide("alert_no_material"))
            return
        self.progress_bar.setVisible(True); self.progress_bar.setRange(0, 0)
        self.status_bar.showMessage(Rosetta_Stone_Ribosome.extract_peptide("status_injecting"))
        self.signal_load_biomass.emit(target_path, self.combo_format.currentText(), self.check_header.isChecked())

    def synchronize_available_columns(self, detected_columns: List[str]) -> None:
        self.active_columns = detected_columns
        for mutation_pod in self.matryoshka_pods: 
            self.vaporize_matryoshka_pod(mutation_pod)
        self.matryoshka_pods.clear()
        self.spawn_matryoshka_pod()
        self.combo_deep_scan.clear(); self.combo_deep_scan.addItems(detected_columns)

    def spawn_matryoshka_pod(self) -> None:
        if not self.active_columns: return
        new_mutation_pod = Matryoshka_Filter_Pod(self.active_columns)
        new_mutation_pod.signal_self_destruct.connect(self.vaporize_matryoshka_pod)
        self.matryoshka_pods.append(new_mutation_pod); self.pod_layout.addWidget(new_mutation_pod)

    def vaporize_matryoshka_pod(self, target_pod: Matryoshka_Filter_Pod) -> None:
        if len(self.matryoshka_pods) > 1 or target_pod not in self.matryoshka_pods:
            if target_pod in self.matryoshka_pods:
                self.matryoshka_pods.remove(target_pod)
            try:
                target_pod.signal_self_destruct.disconnect()
            except TypeError:
                pass
            target_pod.setParent(None)
            target_pod.deleteLater()

    def trigger_collision(self) -> None:
        active_mutations = [w.extract_dna_sequence() for w in self.matryoshka_pods if w.line_value.text() or w.check_only_nan.isChecked()]
        if not active_mutations: return
        self.signal_ignite_filters.emit(active_mutations)

    def update_xray_panel(self, stat_block: Dict[str, Any]) -> None:
        lbl1 = Rosetta_Stone_Ribosome.extract_peptide("raw_reads")
        lbl2 = Rosetta_Stone_Ribosome.extract_peptide("survivors")
        lbl3 = Rosetta_Stone_Ribosome.extract_peptide("rate")
        self.lbl_raw_reads.setText(f"{lbl1} {stat_block['total_reads']}")
        self.lbl_surviving.setText(f"{lbl2} {stat_block['surviving_reads']}")
        self.lbl_survival_rate.setText(f"{lbl3} {stat_block['survival_rate']}%")
