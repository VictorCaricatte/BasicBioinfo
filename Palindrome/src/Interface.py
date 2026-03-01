import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QStackedWidget, 
                             QTextEdit, QLineEdit, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QFileDialog, QMessageBox, QFrame, QTabWidget,
                             QScrollArea, QInputDialog)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont, QIcon, QColor

from Logic import GrandAlchemistOfDoubleStrands
from Dependences import ArsenalSummoner
from DB import GuardianDragonOfMatrices

class DwarvenMinerThread(QThread):
    smoke_signal = pyqtSignal(dict)
    scream_of_agony = pyqtSignal(str)

    def __init__(self, mage, db_dragon, sequence, min_len, max_len, mismatches):
        super().__init__()
        self.mage = mage
        self.db_dragon = db_dragon
        self.sequence = sequence
        self.min_len = min_len
        self.max_len = max_len
        self.mismatches = mismatches

    def run(self):
        try:
            # Invokes all magics simultaneously
            global_gc = self.mage.summon_dragon_breath_gc(self.sequence)
            runes_tally = self.mage.tally_the_magical_runes(self.sequence)
            treasures = self.mage.conjure_search_for_palindromic_treasures(self.sequence, self.min_len, self.max_len, self.mismatches)
            portals = self.mage.open_orf_portals(self.sequence)
            traps = self.mage.predict_hairpin_traps(self.sequence)
            archipelagos = self.mage.map_archipelagos_of_power_cpg(self.sequence)
            echoes = self.mage.hear_echos_of_tandem_cave(self.sequence)
            
            # The New Grand Magics for the Thread
            crispr_targets = self.mage.summon_homunculus_cas9_cleavage(self.sequence)
            ribosomal_marks = self.mage.illuminate_ribosomal_marks_sd(self.sequence)
            transcription_barriers = self.mage.erect_transcription_barriers(self.sequence)
            
            # Export the results to CSV and paint the Vectorial Maps
            csv_path = self.db_dragon.devour_csv_tribute_advanced("alchemical_treasures.csv", treasures)
            svg_path = self.mage.paint_runic_constellations_on_parchment(len(self.sequence), treasures, self.db_dragon.bone_lair)
            
            # Espelho do Ouroboros (Circular Map)
            mystical_features_for_circos = [{"start": t["start_idx"], "end": t["start_idx"] + t["size"]} for t in treasures]
            ouroboros_path = self.mage.gaze_into_ouroboros_mirror(len(self.sequence), mystical_features_for_circos, self.db_dragon.bone_lair)
            
            spoils = {
                'global_gc': global_gc,
                'tally': runes_tally,
                'palindromes': treasures,
                'orfs': portals,
                'hairpins': traps,
                'cpg': archipelagos,
                'tandem': echoes,
                'homunculus_cas9': crispr_targets,
                'ribosomes_sd': ribosomal_marks,
                'terminators': transcription_barriers,
                'csv': csv_path,
                'svg': svg_path,
                'ouroboros_svg': ouroboros_path
            }
            
            # Pergaminho Definitivo (Ultimate HTML Report)
            html_parchment_path = self.db_dragon.inscribe_ultimate_parchment_html("ultimate_parchment.html", spoils)
            spoils['html_parchment'] = html_parchment_path
            
            self.smoke_signal.emit(spoils)
        except Exception as e:
            self.scream_of_agony.emit(str(e))

class PalindromeTavernGrandGUI(QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.mage = GrandAlchemistOfDoubleStrands()
        self.armory = ArsenalSummoner()
        self.dragon_db = GuardianDragonOfMatrices()
        self.dwarven_miner = None
        
        self.erect_tavern_walls()
        self.paint_with_dark_magic()

    def erect_tavern_walls(self):
        self.setWindowTitle("⚔️ Palindrome")
        self.resize(1600, 950)
        
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(0)

        # Side Bar setup
        self.side_bar = QFrame()
        self.side_bar.setObjectName("SideBar")
        self.side_bar.setFixedWidth(240)
        self.side_layout = QVBoxLayout(self.side_bar)
        self.side_layout.setContentsMargins(10, 20, 10, 20)

        self.tavern_label = QLabel("🛡️ Palindrome")
        self.tavern_label.setObjectName("TavernTitle")
        self.tavern_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.side_layout.addWidget(self.tavern_label)
        self.side_layout.addSpacing(30)

        self.btn_arena = self.forge_side_button("🧬 Arena of Analysis")
        self.btn_spells = self.forge_side_button("📜 Quick Spells")
        self.btn_armory = self.forge_side_button("⚔️ Armory (Config)")
        self.btn_dungeon = self.forge_side_button("🐉 Dungeon (DB)")
        self.btn_grimoire = self.forge_side_button("📖 Grimoire (Help)")
        
        self.btn_arena.clicked.connect(lambda: self.screen_stack.setCurrentIndex(0))
        self.btn_spells.clicked.connect(lambda: self.screen_stack.setCurrentIndex(1))
        self.btn_armory.clicked.connect(lambda: self.screen_stack.setCurrentIndex(2))
        self.btn_dungeon.clicked.connect(lambda: self.screen_stack.setCurrentIndex(3))
        self.btn_grimoire.clicked.connect(lambda: self.screen_stack.setCurrentIndex(4))

        self.side_layout.addStretch()

        # The cleanup spell button
        self.btn_clean = QPushButton("🧹 Banish Shadows (Clean)")
        self.btn_clean.setObjectName("SideButton")
        self.btn_clean.clicked.connect(self.clean_dungeon_remnants)
        self.side_layout.addWidget(self.btn_clean)

        self.screen_stack = QStackedWidget()
        self.screen_stack.setObjectName("ScreenStack")
        
        self.screen_arena = self.build_analysis_arena()
        self.screen_spells = self.build_quick_spells_sanctuary()
        self.screen_armory = self.build_armory_forge()
        self.screen_dungeon = self.build_dragon_lair()
        self.screen_grimoire = self.build_grimoire_help()
        
        self.screen_stack.addWidget(self.screen_arena)
        self.screen_stack.addWidget(self.screen_spells)
        self.screen_stack.addWidget(self.screen_armory)
        self.screen_stack.addWidget(self.screen_dungeon)
        self.screen_stack.addWidget(self.screen_grimoire)

        self.main_layout.addWidget(self.side_bar)
        self.main_layout.addWidget(self.screen_stack)

    def forge_side_button(self, text):
        btn = QPushButton(text)
        btn.setObjectName("SideButton")
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.side_layout.addWidget(btn)
        return btn

    def build_analysis_arena(self):
        arena = QWidget()
        layout = QVBoxLayout(arena)
        layout.setContentsMargins(30, 20, 30, 20)

        title = QLabel("🔬 Multiphase Analysis Arena")
        title.setObjectName("ArenaTitle")
        layout.addWidget(title)

        # File Load and New Retrieval Methods
        row_load = QHBoxLayout()
        self.btn_load_scroll = QPushButton("📂 Unroll Scroll (File)")
        self.btn_load_scroll.setObjectName("LoadButton")
        self.btn_load_scroll.clicked.connect(lambda: self.summon_sequence_from_file(self.dna_box))
        
        self.btn_ncbi = QPushButton("🛸 Abduct NCBI (Accession)")
        self.btn_ncbi.setObjectName("LoadButton")
        self.btn_ncbi.clicked.connect(self.fetch_from_ncbi)
        
        self.btn_mock = QPushButton("🎲 Conjure Mock DNA")
        self.btn_mock.setObjectName("LoadButton")
        self.btn_mock.clicked.connect(self.generate_mock_dna)

        # NEW: Forja de Lote (Batch Forge) Button
        self.btn_batch_forge = QPushButton("🏭 Ignite Batch Forge (Folder)")
        self.btn_batch_forge.setObjectName("LoadButton")
        self.btn_batch_forge.clicked.connect(self.ignite_the_batch_forge_ritual)

        row_load.addWidget(self.btn_load_scroll)
        row_load.addWidget(self.btn_ncbi)
        row_load.addWidget(self.btn_mock)
        row_load.addWidget(self.btn_batch_forge)
        row_load.addStretch()
        layout.addLayout(row_load)

        self.dna_box = QTextEdit()
        self.dna_box.setMaximumHeight(80)
        self.dna_box.setPlaceholderText("Pour the messy PDF/DNA scroll here... We will purify it.")
        layout.addWidget(self.dna_box)

        control_panel = QHBoxLayout()
        
        control_panel.addWidget(QLabel("Min Power:"))
        self.input_min = QLineEdit("4")
        self.input_min.setFixedWidth(40)
        control_panel.addWidget(self.input_min)

        control_panel.addWidget(QLabel("Max Power:"))
        self.input_max = QLineEdit("12")
        self.input_max.setFixedWidth(40)
        control_panel.addWidget(self.input_max)
        
        control_panel.addWidget(QLabel("Cursed Mismatches:"))
        self.input_mismatch = QLineEdit("0")
        self.input_mismatch.setFixedWidth(40)
        control_panel.addWidget(self.input_mismatch)

        self.btn_conjure = QPushButton("⚡ Conjure Supreme Magic")
        self.btn_conjure.setObjectName("ActionButton")
        self.btn_conjure.clicked.connect(self.start_search_ritual)
        control_panel.addWidget(self.btn_conjure)
        
        layout.addLayout(control_panel)
        
        self.label_global_gc = QLabel("Global Dragon's Breath (GC%): Awaiting sacrifice...")
        self.label_global_gc.setObjectName("HighlightLabel")
        layout.addWidget(self.label_global_gc)

        self.label_tally = QLabel("Nucleotide Tally: None")
        self.label_tally.setStyleSheet("color: #a0a0b5; font-size: 12px; margin-bottom: 10px;")
        layout.addWidget(self.label_tally)

        self.result_tabs = QTabWidget()
        
        # Tab 1: Palindromes Core with Delta G
        self.table_treasures = QTableWidget(0, 10)
        self.table_treasures.setHorizontalHeaderLabels(["Coord.", "Rune", "Size", "Mism.", "GC%", "Blade", "Transmute(Prot)", "Scale(Da)", "Fire(Tm)", "ΔG (kcal/mol)"])
        self.table_treasures.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_treasures, "💎 Core Palindromes")

        # Tab 2: Advanced Arcanum (Now showing Mithril Primers)
        self.table_advanced = QTableWidget(0, 8)
        self.table_advanced.setHorizontalHeaderLabels(["Coord.", "Codon Bias (GC3%)", "Pathogen Sig.", "Gaps (N)", "Entropy", "Methylation", "Chance Prob.", "Fwd/Rev Mithril Primers"])
        self.table_advanced.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_advanced, "🔮 Advanced Arcanum")
        
        # Tab 3: ORFs
        self.table_orfs = QTableWidget(0, 5)
        self.table_orfs.setHorizontalHeaderLabels(["Direction", "Start", "End", "DNA Size", "Amino Acids"])
        self.table_orfs.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_orfs, "🌀 Portals (ORFs)")
        
        # Tab 4: Hairpins
        self.table_hairpins = QTableWidget(0, 5)
        self.table_hairpins.setHorizontalHeaderLabels(["Coord.", "Stem", "Loop", "Loop Size", "Threat"])
        self.table_hairpins.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_hairpins, "🪝 Traps (Hairpins)")

        # Tab 5: CpG Islands
        self.table_cpg = QTableWidget(0, 4)
        self.table_cpg.setHorizontalHeaderLabels(["Ship Start", "End", "Avg GC%", "Factor (Obs/Exp)"])
        self.table_cpg.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_cpg, "🏝️ Archipelagos (CpG)")

        # Tab 6: Tandem Repeats
        self.table_echoes = QTableWidget(0, 4)
        self.table_echoes.setHorizontalHeaderLabels(["Start", "End", "Scream Motif", "Repeated Times"])
        self.table_echoes.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_echoes, "🦇 Echoes (Tandem)")

        # NEW Tab 7: Evolutionary Arts (Cas9, Ribosomes, Terminators)
        self.table_evolutionary_arts = QTableWidget(0, 4)
        self.table_evolutionary_arts.setHorizontalHeaderLabels(["Mystical Entity", "Start Coordinate", "Target/Motif", "Additional Lore"])
        self.table_evolutionary_arts.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.table_evolutionary_arts, "🧪 Evolutionary Arts")

        layout.addWidget(self.result_tabs)
        return arena

    def build_quick_spells_sanctuary(self):
        sanctuary = QWidget()
        layout = QVBoxLayout(sanctuary)
        layout.setContentsMargins(30, 30, 30, 30)
        
        title = QLabel("📜 Sanctuary of Quick Spells & Structural Magic")
        title.setObjectName("ArenaTitle")
        layout.addWidget(title)

        row_load = QHBoxLayout()
        self.btn_load_spell_scroll = QPushButton("📂 Unroll Ancient Scroll (Load File)")
        self.btn_load_spell_scroll.setObjectName("LoadButton")
        self.btn_load_spell_scroll.clicked.connect(lambda: self.summon_sequence_from_file(self.input_spell_dna))
        row_load.addWidget(self.btn_load_spell_scroll)
        row_load.addStretch()
        layout.addLayout(row_load)

        self.input_spell_dna = QTextEdit()
        self.input_spell_dna.setMaximumHeight(100)
        self.input_spell_dna.setPlaceholderText("Enter a fragment of DNA for quick manipulation...\n*For Ti/Tv, Synteny, or Yggdrasil Tree, enter multiple sequences separated by a new line!*")
        layout.addWidget(self.input_spell_dna)

        # First row of spells
        btn_layout_first = QHBoxLayout()
        btn_mirror = QPushButton("🪞 Fast Mirror (RevComp)")
        btn_mirror.setObjectName("ActionButton")
        btn_mirror.clicked.connect(self.cast_mirror_of_illusions)
        
        btn_wind = QPushButton("🌬️ Whisper of Wind (mRNA)")
        btn_wind.setObjectName("ActionButton")
        btn_wind.clicked.connect(self.cast_whisper_of_wind)

        btn_regex = QPushButton("🔍 Regex Beast Hunt")
        btn_regex.setObjectName("ActionButton")
        btn_regex.clicked.connect(self.cast_regex_hunt)
        
        btn_masking = QPushButton("🛡️ Escudo de Poeira (Mask)")
        btn_masking.setObjectName("ActionButton")
        btn_masking.clicked.connect(self.cast_dust_shield_spell)

        btn_layout_first.addWidget(btn_mirror)
        btn_layout_first.addWidget(btn_wind)
        btn_layout_first.addWidget(btn_regex)
        btn_layout_first.addWidget(btn_masking)
        layout.addLayout(btn_layout_first)

        # Second row of spells
        btn_layout_second = QHBoxLayout()
        btn_trimming = QPushButton("✂️ Tesoura do Ferreiro (Trim)")
        btn_trimming.setObjectName("ActionButton")
        btn_trimming.clicked.connect(self.wield_blacksmith_scissors_in_ui)

        btn_titv = QPushButton("⚖️ Balança Evolutiva (Ti/Tv)")
        btn_titv.setObjectName("ActionButton")
        btn_titv.clicked.connect(self.weigh_souls_on_titv_scale)

        btn_synteny = QPushButton("🔗 Ligação de Almas (Synteny)")
        btn_synteny.setObjectName("ActionButton")
        btn_synteny.clicked.connect(self.forge_soul_link_in_ui)

        btn_yggdrasil = QPushButton("🌳 Semente de Yggdrasil (Tree)")
        btn_yggdrasil.setObjectName("ActionButton")
        btn_yggdrasil.clicked.connect(self.plant_yggdrasil_seed_in_ui)

        btn_layout_second.addWidget(btn_trimming)
        btn_layout_second.addWidget(btn_titv)
        btn_layout_second.addWidget(btn_synteny)
        btn_layout_second.addWidget(btn_yggdrasil)
        layout.addLayout(btn_layout_second)

        self.spell_output = QTextEdit()
        self.spell_output.setReadOnly(True)
        self.spell_output.setPlaceholderText("The result of your incantation will appear here...")
        layout.addWidget(self.spell_output)

        return sanctuary

    def build_armory_forge(self):
        forge = QWidget()
        layout = QVBoxLayout(forge)
        layout.setContentsMargins(30, 30, 30, 30)

        title = QLabel("⚔️ Armory of Dependencies")
        title.setObjectName("ArenaTitle")
        layout.addWidget(title)

        layout.addWidget(QLabel("Sacred Path of the BLAST executable (makeblastdb):"))
        row_blast = QHBoxLayout()
        self.input_blast = QLineEdit(self.armory.reveal_current_armory()["blast"])
        btn_search_blast = QPushButton("🔍 Map")
        btn_search_blast.clicked.connect(lambda: self.seek_artifact(self.input_blast))
        row_blast.addWidget(self.input_blast)
        row_blast.addWidget(btn_search_blast)
        layout.addLayout(row_blast)

        layout.addWidget(QLabel("Sacred Path of the DIAMOND executable:"))
        row_diamond = QHBoxLayout()
        self.input_diamond = QLineEdit(self.armory.reveal_current_armory()["diamond"])
        btn_search_diamond = QPushButton("🔍 Map")
        btn_search_diamond.clicked.connect(lambda: self.seek_artifact(self.input_diamond))
        row_diamond.addWidget(self.input_diamond)
        row_diamond.addWidget(btn_search_diamond)
        layout.addLayout(row_diamond)

        btn_forge = QPushButton("🔨 Carve into Runestones (Save)")
        btn_forge.setObjectName("ActionButton")
        btn_forge.clicked.connect(self.save_armory)
        layout.addWidget(btn_forge)
        
        layout.addStretch()
        return forge

    def build_dragon_lair(self):
        lair = QWidget()
        layout = QVBoxLayout(lair)
        layout.setContentsMargins(30, 30, 30, 30)

        title = QLabel("🐉 Dungeon of Data (DB)")
        title.setObjectName("ArenaTitle")
        layout.addWidget(title)
        
        layout.addWidget(QLabel("Here the Guardian Dragon forges Matrices. Select a FASTA scroll:"))
        
        row_fasta = QHBoxLayout()
        self.input_fasta = QLineEdit()
        btn_search_fasta = QPushButton("📜 Select FASTA")
        btn_search_fasta.clicked.connect(lambda: self.seek_artifact(self.input_fasta, is_file=True, db_mode=True))
        row_fasta.addWidget(self.input_fasta)
        row_fasta.addWidget(btn_search_fasta)
        layout.addLayout(row_fasta)

        btn_panel = QHBoxLayout()
        btn_blast = QPushButton("🔥 Spew BLAST Matrix")
        btn_blast.setObjectName("ActionButton")
        btn_blast.clicked.connect(self.summon_blast_db)
        
        btn_diamond = QPushButton("💎 Crystallize DIAMOND Library")
        btn_diamond.setObjectName("ActionButton")
        btn_diamond.clicked.connect(self.summon_diamond_db)
        
        btn_panel.addWidget(btn_blast)
        btn_panel.addWidget(btn_diamond)
        layout.addLayout(btn_panel)

        self.dragon_log = QTextEdit()
        self.dragon_log.setReadOnly(True)
        self.dragon_log.setPlaceholderText("The dragon is sleeping...")
        layout.addWidget(self.dragon_log)

        return lair

    def build_grimoire_help(self):
        grimoire = QWidget()
        layout = QVBoxLayout(grimoire)
        layout.setContentsMargins(30, 30, 30, 30)
        
        title = QLabel("📖 The Grand Grimoire (Help & Documentation)")
        title.setObjectName("ArenaTitle")
        layout.addWidget(title)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)

        help_text = """
        <h2 style='color:#ff7b00;'>Welcome to the Genomic Enchanter</h2>
        <p>This sacred tool is designed for severe bioinformatics processing wrapped in the mystical arts.</p>
        
        <h3 style='color:#00ffcc;'>🧬 Ultimate Grand Magics</h3>
        <ul>
            <li><b>Pergaminho Definitivo:</b> HTML report automatically forged when a ritual succeeds.</li>
            <li><b>Forja de Lote (Batch):</b> Process an entire kingdom of FASTA files from a folder.</li>
            <li><b>Chaves de Mithril (Primer3):</b> Simulated heuristic primer generation in the Advanced Arcanum tab.</li>
            <li><b>Semente de Yggdrasil:</b> Plant sequences in the Quick Spells tab to grow a Newick Phylogenetic Tree.</li>
            <li><b>Homunculus Cas9:</b> CRISPR target identification (NGG PAMs) located in the Evolutionary Arts tab.</li>
            <li><b>Espelho do Ouroboros:</b> A circular plasmid-like SVG map is now conjured alongside the linear map.</li>
            <li><b>Escudo de Poeira & Tesoura do Ferreiro:</b> Quick spells for Masking Repeats (N's) and Trimming sequences.</li>
            <li><b>Ligação de Almas (Synteny) & Balança Evolutiva (Ti/Tv):</b> Compare two sequences directly in the Quick Spells tab.</li>
            <li><b>Marcação de Ribossomos & Barreiras:</b> Shine-Dalgarno motifs and Intrinsic Terminators predicted.</li>
        </ul>
        """
        
        text_label = QLabel(help_text)
        text_label.setWordWrap(True)
        text_label.setTextFormat(Qt.TextFormat.RichText)
        text_label.setStyleSheet("font-size: 14px; line-height: 1.5;")
        
        content_layout.addWidget(text_label)
        content_layout.addStretch()
        scroll_area.setWidget(content_widget)
        
        layout.addWidget(scroll_area)
        return grimoire

    def paint_with_dark_magic(self):
        dark_style = """
            QMainWindow { background-color: #121212; }
            #SideBar { background-color: #1A1A24; border-right: 2px solid #ff7b00; }
            #TavernTitle { color: #ff7b00; font-size: 20px; font-weight: bold; font-family: 'Segoe UI', sans-serif; }
            #SideButton { background-color: transparent; color: #a0a0b5; text-align: left; padding: 12px; font-size: 14px; border: none; border-radius: 5px; }
            #SideButton:hover { background-color: #2b2b3c; color: #ffffff; }
            #ArenaTitle { color: #ffffff; font-size: 24px; font-weight: bold; margin-bottom: 10px; }
            #HighlightLabel { color: #00ffcc; font-size: 16px; font-weight: bold; margin: 5px 0px; }
            QLabel { color: #d1d1e0; font-size: 14px; }
            QTextEdit, QLineEdit { background-color: #232332; border: 1px solid #3c3c54; color: #ffffff; padding: 8px; border-radius: 5px; font-family: 'Courier New', monospace; }
            QPushButton { background-color: #3c3c54; color: white; border: none; padding: 10px; border-radius: 5px; font-weight: bold; }
            QPushButton:hover { background-color: #505070; }
            #ActionButton { background-color: #ff7b00; color: #000000; }
            #ActionButton:hover { background-color: #ff9933; }
            #LoadButton { background-color: #1b5e20; color: white; padding: 6px 12px; border-radius: 4px; }
            #LoadButton:hover { background-color: #2e7d32; }
            QTableWidget { background-color: #1e1e2c; color: white; border: 1px solid #3c3c54; gridline-color: #3c3c54; }
            QHeaderView::section { background-color: #2b2b3c; color: #ff7b00; padding: 5px; border: 1px solid #3c3c54; font-weight: bold; font-size: 11px; }
            QTabWidget::pane { border: 1px solid #3c3c54; background: #1e1e2c; }
            QTabBar::tab { background: #2b2b3c; color: #a0a0b5; padding: 8px 15px; border-top-left-radius: 4px; border-top-right-radius: 4px; font-weight: bold; }
            QTabBar::tab:selected { background: #ff7b00; color: #000000; }
            QScrollArea { border: none; background-color: transparent; }
            QScrollArea > QWidget > QWidget { background-color: transparent; }
        """
        self.setStyleSheet(dark_style)

    
    def generate_mock_dna(self):
        mock = self.mage.conjure_homunculus_mock_dna(1500)
        self.dna_box.setText(mock)
        QMessageBox.information(self, "Homunculus Spawned", "A 100% random DNA sequence of 1500 bp was conjured into the cauldron.")

    def fetch_from_ncbi(self):
        acc, ok = QInputDialog.getText(self, "Abduct NCBI", "Enter Accession Code (e.g., NM_001308):")
        if ok and acc:
            try:
                seq = self.mage.abduct_alien_dna_from_ncbi_dimension(acc)
                self.dna_box.setText(seq)
                QMessageBox.information(self, "Alien Abducted", f"Sequence {acc} ripped from the NCBI dimension and loaded.")
            except Exception as e:
                QMessageBox.critical(self, "Portal Failed", str(e))

    def ignite_the_batch_forge_ritual(self):
        folder_path = QFileDialog.getExistingDirectory(self, "Select Territory for Batch Forge")
        if folder_path:
            try:
                scrolls_forged, summary_tome = self.mage.ignite_the_hellfire_batch_forge(folder_path, self.dragon_db)
                QMessageBox.information(self, "Batch Forge Complete", f"🔥 Successfully forged {scrolls_forged} ancient scrolls!\n\nThe Grand Summary Tome is mapped at:\n{summary_tome}")
            except Exception as e:
                QMessageBox.critical(self, "Batch Forge Failed", str(e))

    def clean_dungeon_remnants(self):
        count = self.dragon_db.banish_shadow_remnants_to_void()
        QMessageBox.information(self, "Purged!", f"{count} shadowy temp files (CSV/SVG/HTML) have been banished to the void.")

    def summon_sequence_from_file(self, target_text_box):
        path, _ = QFileDialog.getOpenFileName(
            self, 
            "Locate the Ancient Genomic Tome", 
            "", 
            "Genomic Scrolls (*.fasta *.fastq *.fastap *.faa *.ffn *.gbk *.gbff *.gbf);;All Artifacts (*)"
        )
        if path:
            try:
                sequence = self.mage.decipher_ancient_genomic_tome(path)
                if len(sequence) > 1000000:
                    QMessageBox.warning(self, "Massive Scroll", "This tome is huge! Visualizing it might stutter the tavern.")
                target_text_box.setText(sequence[:100000]) 
                QMessageBox.information(self, "Tome Deciphered", f"Successfully extracted magical runes. Supports multi-FASTA reading!")
            except Exception as e:
                QMessageBox.critical(self, "Failed to Unroll", str(e))

    def start_search_ritual(self):
        sequence = self.dna_box.toPlainText()
        try:
            min_len = int(self.input_min.text())
            max_len = int(self.input_max.text())
            mismatches = int(self.input_mismatch.text())
            purified_seq = self.mage.exorcise_pdf_demons_and_purify(sequence)
        except ValueError as err:
            QMessageBox.critical(self, "Magic Failed!", str(err))
            return

        self.btn_conjure.setEnabled(False)
        self.btn_conjure.setText("⏳ Channeling universal energies...")
        
        self.table_treasures.setRowCount(0)
        self.table_advanced.setRowCount(0)
        self.table_orfs.setRowCount(0)
        self.table_hairpins.setRowCount(0)
        self.table_cpg.setRowCount(0)
        self.table_echoes.setRowCount(0)
        self.table_evolutionary_arts.setRowCount(0)

        self.dwarven_miner = DwarvenMinerThread(self.mage, self.dragon_db, purified_seq, min_len, max_len, mismatches)
        self.dwarven_miner.smoke_signal.connect(self.display_battle_spoils)
        self.dwarven_miner.scream_of_agony.connect(self.display_fatal_error)
        self.dwarven_miner.start()

    def display_battle_spoils(self, spoils):
        self.label_global_gc.setText(f"🐉 Global Dragon's Breath (GC%): {spoils['global_gc']:.2f}%")
        self.label_tally.setText(f"Nucleotide Tally: {spoils['tally']}")
        
        palindromes = spoils['palindromes']
        self.table_treasures.setRowCount(len(palindromes))
        self.table_advanced.setRowCount(len(palindromes))
        for row, pal in enumerate(palindromes):
            self.table_treasures.setItem(row, 0, QTableWidgetItem(pal['position']))
            self.table_treasures.setItem(row, 1, QTableWidgetItem(pal['sequence']))
            self.table_treasures.setItem(row, 2, QTableWidgetItem(str(pal['size'])))
            self.table_treasures.setItem(row, 3, QTableWidgetItem(str(pal['mismatches'])))
            self.table_treasures.setItem(row, 4, QTableWidgetItem(f"{pal['gc']}%"))
            self.table_treasures.setItem(row, 5, QTableWidgetItem(pal['monster']))
            self.table_treasures.setItem(row, 6, QTableWidgetItem(pal['protein']))
            self.table_treasures.setItem(row, 7, QTableWidgetItem(str(pal['weight_da'])))
            self.table_treasures.setItem(row, 8, QTableWidgetItem(f"{pal['tm_celsius']}°C"))
            self.table_treasures.setItem(row, 9, QTableWidgetItem(str(pal['delta_g'])))
            
            self.table_advanced.setItem(row, 0, QTableWidgetItem(pal['position']))
            self.table_advanced.setItem(row, 1, QTableWidgetItem(f"{pal['codon_bias']}%"))
            self.table_advanced.setItem(row, 2, QTableWidgetItem(pal['pathogen']))
            self.table_advanced.setItem(row, 3, QTableWidgetItem(str(pal['gaps'])))
            self.table_advanced.setItem(row, 4, QTableWidgetItem(str(pal['entropy'])))
            self.table_advanced.setItem(row, 5, QTableWidgetItem(pal['methylation']))
            self.table_advanced.setItem(row, 6, QTableWidgetItem(pal['probability']))
            self.table_advanced.setItem(row, 7, QTableWidgetItem(f"F:{pal['fwd_primer'][:8]}... R:{pal['rev_primer'][:8]}..."))

        orfs = spoils['orfs']
        self.table_orfs.setRowCount(len(orfs))
        for row, orf in enumerate(orfs):
            self.table_orfs.setItem(row, 0, QTableWidgetItem(orf['direction']))
            self.table_orfs.setItem(row, 1, QTableWidgetItem(str(orf['start'])))
            self.table_orfs.setItem(row, 2, QTableWidgetItem(str(orf['end'])))
            self.table_orfs.setItem(row, 3, QTableWidgetItem(str(orf['bp_size'])))
            self.table_orfs.setItem(row, 4, QTableWidgetItem(str(orf['amino_acids'])))
            
        hairpins = spoils['hairpins']
        self.table_hairpins.setRowCount(len(hairpins))
        for row, hp in enumerate(hairpins):
            self.table_hairpins.setItem(row, 0, QTableWidgetItem(hp['position']))
            self.table_hairpins.setItem(row, 1, QTableWidgetItem(hp['stem']))
            self.table_hairpins.setItem(row, 2, QTableWidgetItem(hp['loop']))
            self.table_hairpins.setItem(row, 3, QTableWidgetItem(str(hp['loop_size'])))
            self.table_hairpins.setItem(row, 4, QTableWidgetItem(hp['threat_level']))

        cpgs = spoils['cpg']
        self.table_cpg.setRowCount(len(cpgs))
        for row, island in enumerate(cpgs):
            self.table_cpg.setItem(row, 0, QTableWidgetItem(str(island['start'])))
            self.table_cpg.setItem(row, 1, QTableWidgetItem(str(island['end'])))
            self.table_cpg.setItem(row, 2, QTableWidgetItem(f"{island['avg_gc']}%"))
            self.table_cpg.setItem(row, 3, QTableWidgetItem(str(island['obs_exp'])))

        tandems = spoils['tandem']
        self.table_echoes.setRowCount(len(tandems))
        for row, echo in enumerate(tandems):
            self.table_echoes.setItem(row, 0, QTableWidgetItem(str(echo['start'])))
            self.table_echoes.setItem(row, 1, QTableWidgetItem(str(echo['end'])))
            self.table_echoes.setItem(row, 2, QTableWidgetItem(echo['motif']))
            self.table_echoes.setItem(row, 3, QTableWidgetItem(str(echo['repeats'])))

        # Fill the New Evolutionary Arts Table
        evo_arts_data = []
        for crispr in spoils['homunculus_cas9']: 
            evo_arts_data.append(("Homunculus Cas9 (CRISPR)", str(crispr['start']), crispr['target'], f"Strand {crispr['strand']}"))
        for ribosome in spoils['ribosomes_sd']: 
            evo_arts_data.append(("Ribosome (Shine-Dalgarno)", str(ribosome['orf_start_coordinate']), "AGGAGG / GGAGG Found", "Upstream of ATG Portal"))
        for terminator in spoils['terminators']: 
            evo_arts_data.append(("Intrinsic Terminator", str(terminator['barrier_start']), terminator['barrier_motif'], f"Ends at {terminator['barrier_end']}"))
        
        self.table_evolutionary_arts.setRowCount(len(evo_arts_data))
        for row, evo_data in enumerate(evo_arts_data):
            self.table_evolutionary_arts.setItem(row, 0, QTableWidgetItem(evo_data[0]))
            self.table_evolutionary_arts.setItem(row, 1, QTableWidgetItem(evo_data[1]))
            self.table_evolutionary_arts.setItem(row, 2, QTableWidgetItem(evo_data[2]))
            self.table_evolutionary_arts.setItem(row, 3, QTableWidgetItem(evo_data[3]))

        self.btn_conjure.setEnabled(True)
        self.btn_conjure.setText("⚡ Conjure Supreme Magic")
        QMessageBox.information(self, "Alchemical Discoveries", 
            f"The Ritual Succeeded!\n\n"
            f"📄 HTML Pergaminho Definitivo: {spoils['html_parchment']}\n"
            f"🎨 Ouroboros Circular Map: {spoils.get('ouroboros_svg', 'Not Forged')}\n"
            f"📁 Linear SVG Map: {spoils['svg']}\n"
            f"📂 CSV Treasury: {spoils['csv']}")

    def display_fatal_error(self, error_msg):
        self.btn_conjure.setEnabled(True)
        self.btn_conjure.setText("⚡ Conjure Supreme Magic")
        QMessageBox.critical(self, "Disaster in Battle", error_msg)

    def seek_artifact(self, input_widget, is_file=True, db_mode=False):
        if is_file:
            if db_mode:
                filter_str = "FASTA Scrolls (*.fasta *.fas *.fna *.faa);;All Artifacts (*)"
            else:
                filter_str = "Genomic Scrolls (*.fasta *.fastq *.gbk);;All Artifacts (*)"
            path, _ = QFileDialog.getOpenFileName(self, "Find the Relic", "", filter_str)
        else:
            path = QFileDialog.getExistingDirectory(self, "Find the Territory")
        if path:
            input_widget.setText(path)

    def save_armory(self):
        self.armory.tame_blast_beast(self.input_blast.text())
        self.armory.find_lost_diamond(self.input_diamond.text())
        QMessageBox.information(self, "Armory Updated", "The weapons have been properly cataloged.")

    def summon_blast_db(self):
        fasta = self.input_fasta.text()
        if not fasta:
            self.dragon_log.append("⚠️ The dragon demands a FASTA scroll for the sacrifice!")
            return
        try:
            self.dragon_log.append("🔥 The dragon is summoning BLAST...")
            result = self.dragon_db.spew_blast_matrix(fasta)
            self.dragon_log.append(f"✔️ Matrix forged: {result}")
        except Exception as e:
            self.dragon_log.append(f"❌ Forge failure: {str(e)}")

    def summon_diamond_db(self):
        fasta = self.input_fasta.text()
        if not fasta:
            self.dragon_log.append("⚠️ Missing the FASTA scroll!")
            return
        try:
            self.dragon_log.append("💎 The dragon is crushing DIAMOND...")
            result = self.dragon_db.crystallize_diamond_library(fasta)
            self.dragon_log.append(f"✔️ Crystallization complete: {result}")
        except Exception as e:
            self.dragon_log.append(f"❌ Shattered: {str(e)}")

    
    def extract_multiple_runes_from_sanctuary(self):
        raw_text = self.input_spell_dna.toPlainText()
        return [line.strip() for line in raw_text.split('\n') if line.strip()]

    def cast_mirror_of_illusions(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if not runes: return
        try:
            purified = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            rev_comp = self.mage.gaze_into_mirror_of_illusions_revcomp(purified)
            self.spell_output.setText(f"🪞 Fast Reverse Complement:\n{rev_comp}")
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def cast_whisper_of_wind(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if not runes: return
        try:
            purified = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            mrna = self.mage.whisper_of_the_wind_rna_transcription(purified)
            self.spell_output.setText(f"🌬️ Transcribed mRNA:\n{mrna}")
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def cast_regex_hunt(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if not runes: return
        try:
            purified = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            pat, ok = QInputDialog.getText(self, "Regex Hunt", "Enter Regex Pattern (e.g., AAT.G):")
            if ok and pat:
                res = self.mage.track_aberrations_via_regex(purified, pat)
                output = f"🔍 Found {len(res)} matching beasts:\n"
                for r in res:
                    output += f"Coord {r['start']}-{r['end']}: {r['found_beast']}\n"
                self.spell_output.setText(output)
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def cast_dust_shield_spell(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if not runes: return
        try:
            purified = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            shielded = self.mage.raise_dust_shield_masking(purified)
            self.spell_output.setText(f"🛡️ Masked Sequence (Escudo de Poeira):\n{shielded}")
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def wield_blacksmith_scissors_in_ui(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if not runes: return
        try:
            purified = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            trimmed = self.mage.wield_blacksmith_scissors_trimming(purified)
            self.spell_output.setText(f"✂️ Trimmed Sequence (Tesoura do Ferreiro):\n{trimmed}")
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def weigh_souls_on_titv_scale(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if len(runes) < 2:
            self.spell_output.setText("⚠️ A Balança Evolutiva necessita de duas almas (sequências separadas por quebra de linha) para pesar.")
            return
        try:
            alpha = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            omega = self.mage.exorcise_pdf_demons_and_purify(runes[1])
            results = self.mage.weigh_on_evolutionary_scale_titv(alpha, omega)
            self.spell_output.setText(f"⚖️ Ti/Tv Scale (Balança Evolutiva):\n"
                                      f"Total Transitions: {results['Total_Transitions']}\n"
                                      f"Total Transversions: {results['Total_Transversions']}\n"
                                      f"Ti/Tv Ratio: {results['Ti_Tv_Sacred_Ratio']}")
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def forge_soul_link_in_ui(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if len(runes) < 2:
            self.spell_output.setText("⚠️ A Ligação de Almas (Synteny) precisa de duas almas (sequências em linhas separadas).")
            return
        try:
            alpha = self.mage.exorcise_pdf_demons_and_purify(runes[0])
            omega = self.mage.exorcise_pdf_demons_and_purify(runes[1])
            blocks = self.mage.forge_soul_link_synteny(alpha, omega)
            if not blocks:
                self.spell_output.setText("🔗 Nenhuma ligação de alma (sintenia) foi encontrada entre as sequências.")
                return
            
            out_text = f"🔗 Encontrados {len(blocks)} Blocos de Ligação de Almas:\n\n"
            for b in blocks:
                out_text += f"Bloco K-mer: {b['soul_block']}\nCoord na Seq Alpha: {b['alpha_coordinate']}\nCoord na Seq Omega: {b['omega_coordinate']}\n---\n"
            self.spell_output.setText(out_text)
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

    def plant_yggdrasil_seed_in_ui(self):
        runes = self.extract_multiple_runes_from_sanctuary()
        if len(runes) < 2:
            self.spell_output.setText("⚠️ A Semente de Yggdrasil requer no mínimo duas sequências em linhas separadas para germinar.")
            return
        try:
            dict_of_souls = {}
            for i, raw_rune in enumerate(runes):
                purified = self.mage.exorcise_pdf_demons_and_purify(raw_rune)
                dict_of_souls[f"Alma_{i+1}"] = purified
                
            newick_tree = self.mage.plant_yggdrasil_seed_upgma(dict_of_souls)
            self.spell_output.setText(f"🌳 Árvore Yggdrasil Germinada (Formato Newick):\n\n{newick_tree}")
        except Exception as e:
            self.spell_output.setText(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    summoning_ritual = QApplication(sys.argv)
    main_tavern = PalindromeTavernGrandGUI()
    main_tavern.show()
    sys.exit(summoning_ritual.exec())
