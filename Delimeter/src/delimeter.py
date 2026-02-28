# ==============================================================================
# Author:       Victor S Caricatte De Araújo
# Email:        victorsc@ufmg.br
# Institution:  Universidade Federal de Minas Gerais
# Version:      0.3.1
# Description:  Particle Accelerator
# ==============================================================================

import sys
import os
import traceback
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional

class Spatial_Coordinate_Aligner:
    @staticmethod
    def calibrate_gps() -> None:
        reactor_core_path = Path(__file__).resolve().parent
        sys.path.insert(0, str(reactor_core_path))
        os.chdir(str(reactor_core_path))

Spatial_Coordinate_Aligner.calibrate_gps()

try:
    from PyQt5.QtWidgets import QApplication, QMessageBox, QFileDialog
    from PyQt5.QtGui import QColor
    from PyQt5.QtCore import QThread
    import pandas as pd
    
    from Interface import (Genomic_Holodeck_UI, Spectrometry_Dialog, Xenomorph_Splicer_Dialog, 
                           Fasta_Export_Dialog, Transmutation_Altar_Dialog, Radiation_Dialog, Anomaly_Marker_Dialog)
    from logic import (Sterilized_Data_Engine, Hyperdrive_Loader_Thread, Infinite_Hologram_Model, 
                       Chunk_Devourer_Exporter, Elephant_Brain_Storage, Fasta_Devourer_Exporter,
                       Runes_Of_Biocontainment, Quantum_Collider_Thread)
    from xenoglossy_codex import Rosetta_Stone_Ribosome

    class Chromatic_Mutator:
        @staticmethod
        def force_mutation(app_instance: QApplication, stylesheet_path: Path) -> None:
            if not stylesheet_path.exists() and Path(stylesheet_path.name).exists():
                stylesheet_path = Path(stylesheet_path.name)
                
            if stylesheet_path.exists():
                with open(stylesheet_path, "r", encoding='utf-8') as css_file:
                    app_instance.setStyleSheet(css_file.read())
            else:
                logging.warning(f"Aesthetic Anomaly: Encapsulated CSS {stylesheet_path.name} not found in /themes or root.")

    class Supreme_Genomic_Overlord:
        def __init__(self) -> None:
            self.app = QApplication(sys.argv)
            self.data_engine = Sterilized_Data_Engine()
            
            self.ui = Genomic_Holodeck_UI(self.data_engine)
            
            self.table_model: Optional[Infinite_Hologram_Model] = None
            self.anomaly_rules_global: List[Dict[str, Any]] = []
            
            self.active_synthesis_threads: List[QThread] = []
            
            self.themes_directory = Path("themes")
            self.dark_envelope = self.themes_directory / "dark_theme.qss"
            self.light_envelope = self.themes_directory / "light_theme.qss"
            
            Chromatic_Mutator.force_mutation(self.app, self.dark_envelope)
            self.is_dark_matter_active = True
            self.ui.mutate_spectrum_button_icon(self.is_dark_matter_active)

            self.wire_synapses()

        def wire_synapses(self) -> None:
            self.ui.signal_load_biomass.connect(self.verify_and_inject_dark_matter)
            self.ui.signal_ignite_filters.connect(self.detonate_filter_collider)
            self.ui.signal_export_genome.connect(self.synthesize_output)
            self.ui.signal_toggle_spectrum.connect(self.shift_quantum_spectrum)
            self.ui.signal_dialect_mutation.connect(self.apply_linguistic_mutation)
            
            self.ui.btn_undo.clicked.connect(self.execute_undo)
            self.ui.btn_redo.clicked.connect(self.execute_redo)
            
            self.ui.btn_purge.clicked.connect(self.execute_purge)
            self.ui.btn_spectro.clicked.connect(self.open_spectrometry)
            self.ui.btn_splice.clicked.connect(self.open_splicer)
            self.ui.btn_export_fasta.clicked.connect(self.open_fasta_export)
            self.ui.btn_transmute.clicked.connect(self.open_transmutation)
            self.ui.btn_radiate.clicked.connect(self.open_radiation)
            self.ui.btn_anomaly.clicked.connect(self.open_anomaly_marker)
            
            self.ui.btn_save_presets.clicked.connect(self.crystallize_current_mutations)
            self.ui.btn_load_presets.clicked.connect(self.resurrect_extinct_filters)
            
            self.ui.combo_deep_scan.currentTextChanged.connect(self.calculate_nucleotide_entropy)

        def apply_linguistic_mutation(self, newly_selected_dialect: str) -> None:
            Rosetta_Stone_Ribosome.mutate_dialect(newly_selected_dialect)
            self.ui.retranscribe_ui()

        def shift_quantum_spectrum(self) -> None:
            if self.is_dark_matter_active: 
                Chromatic_Mutator.force_mutation(self.app, self.light_envelope)
            else: 
                Chromatic_Mutator.force_mutation(self.app, self.dark_envelope)
            self.is_dark_matter_active = not self.is_dark_matter_active
            self.ui.mutate_spectrum_button_icon(self.is_dark_matter_active)

        def update_ui_from_engine(self) -> None:
            if self.data_engine.filtered_biomass is not None and self.table_model is not None:
                self.table_model.mutate_hologram(self.data_engine.filtered_biomass)
                self.calculate_nucleotide_entropy(self.ui.combo_deep_scan.currentText())
                raw_length = len(self.data_engine.raw_biomass) if self.data_engine.raw_biomass is not None else 0
                surviving_length = len(self.data_engine.filtered_biomass)
                
                survival_rate = round((surviving_length / raw_length) * 100, 2) if raw_length > 0 else 0.0
                genomic_stats = {
                    "total_reads": raw_length,
                    "surviving_reads": surviving_length,
                    "survival_rate": survival_rate
                }
                self.ui.update_xray_panel(genomic_stats)

        def execute_crispr_edit(self, row_idx: int, col_idx: int, nucleotide_val: str) -> None:
            self.data_engine.crispr_cas9_cell_edit(row_idx, col_idx, nucleotide_val)
            self.calculate_nucleotide_entropy(self.ui.combo_deep_scan.currentText())

        def execute_undo(self) -> None:
            if self.data_engine.rewind_time_dilation(): self.update_ui_from_engine()

        def execute_redo(self) -> None:
            if self.data_engine.temporal_redo(): self.update_ui_from_engine()

        def execute_purge(self) -> None:
            if self.data_engine.filtered_biomass is not None:
                self.data_engine.purge_genetic_clones(["ALL_COLUMNS"])
                self.update_ui_from_engine()

        def calculate_nucleotide_entropy(self, genomic_column: str) -> None:
            if self.data_engine.filtered_biomass is not None and genomic_column in self.data_engine.filtered_biomass.columns:
                descriptive_stats = self.data_engine.filtered_biomass[genomic_column].describe(include='all')
                self.ui.text_deep_scan.setText(descriptive_stats.to_string())

        def open_spectrometry(self) -> None:
            if not self.ui.active_columns: return
            dialog = Spectrometry_Dialog(self.data_engine.filtered_biomass, self.ui.combo_deep_scan.currentText())
            dialog.exec_()

        def open_splicer(self) -> None:
            if not self.ui.active_columns: return
            self.splice_dlg = Xenomorph_Splicer_Dialog(self.ui.active_columns)
            self.splice_dlg.btn_splice.clicked.connect(self.execute_splicing)
            self.splice_dlg.exec_()

        def execute_splicing(self) -> None:
            alien_path = self.splice_dlg.path_alien.text()
            local_key = self.splice_dlg.combo_local.currentText()
            alien_key = self.splice_dlg.line_alien_col.text()
            try:
                foreign_biomass = pd.read_csv(alien_path, engine='pyarrow') if alien_path.endswith('.csv') else pd.read_excel(alien_path)
                self.data_engine.splice_alien_dna_strands(foreign_biomass, alien_key, local_key)
                self.update_ui_from_engine()
                self.splice_dlg.accept()
            except Exception as collision_error:
                QMessageBox.critical(self.ui, "Splicing Error", str(collision_error))

        def open_fasta_export(self) -> None:
            if not self.ui.active_columns: return
            self.fasta_dlg = Fasta_Export_Dialog(self.ui.active_columns)
            self.fasta_dlg.btn_export.clicked.connect(self.synthesize_fasta_envelope)
            self.fasta_dlg.exec_()

        def synthesize_fasta_envelope(self) -> None:
            if self.data_engine.filtered_biomass is None: return
            id_column = self.fasta_dlg.combo_id.currentText()
            seq_column = self.fasta_dlg.combo_seq.currentText()
            target_path, _ = QFileDialog.getSaveFileName(self.ui, "Save FASTA", "mutant.fasta", "FASTA (*.fasta)")
            
            if target_path:
                if self.verify_disk_overwrite(target_path):
                    self.ui.progress_bar.setVisible(True)
                    fasta_thread = Fasta_Devourer_Exporter(self.data_engine.filtered_biomass, target_path, id_column, seq_column)
                    self.active_synthesis_threads.append(fasta_thread)
                    fasta_thread.signal_progress.connect(self.ui.progress_bar.setValue)
                    fasta_thread.signal_done.connect(self.synthesis_complete)
                    fasta_thread.finished.connect(lambda: self.active_synthesis_threads.remove(fasta_thread))
                    fasta_thread.start()
                    self.fasta_dlg.accept()

        def open_transmutation(self) -> None:
            if not self.ui.active_columns: return
            self.trans_dlg = Transmutation_Altar_Dialog(self.ui.active_columns)
            self.trans_dlg.btn_transmute.clicked.connect(self.execute_transmutation)
            self.trans_dlg.exec_()

        def execute_transmutation(self) -> None:
            new_strand = self.trans_dlg.line_new.text()
            if not new_strand or self.data_engine.filtered_biomass is None: return
            self.data_engine.transcribe_rna_to_dna(new_strand, self.trans_dlg.combo_c1.currentText(), self.trans_dlg.combo_c2.currentText(), self.trans_dlg.combo_op.currentText())
            self.ui.synchronize_available_columns(list(self.data_engine.filtered_biomass.columns))
            self.update_ui_from_engine()
            self.trans_dlg.accept()

        def open_radiation(self) -> None:
            if not self.ui.active_columns: return
            self.rad_dlg = Radiation_Dialog(self.ui.active_columns)
            self.rad_dlg.btn_radiate.clicked.connect(self.execute_radiation)
            self.rad_dlg.exec_()

        def execute_radiation(self) -> None:
            self.data_engine.quarantine_nan_ghosts(self.rad_dlg.combo_col.currentText(), self.rad_dlg.line_fixed.text(), self.rad_dlg.combo_method.currentText())
            self.update_ui_from_engine()
            self.rad_dlg.accept()

        def open_anomaly_marker(self) -> None:
            if not self.ui.active_columns: return
            self.anom_dlg = Anomaly_Marker_Dialog(self.ui.active_columns)
            self.anom_dlg.btn_mark.clicked.connect(self.apply_anomaly_rule)
            self.anom_dlg.exec_()

        def apply_anomaly_rule(self) -> None:
            col_target = self.anom_dlg.combo_col.currentText()
            operator = self.anom_dlg.combo_op.currentText()
            activator_value = self.anom_dlg.line_val.text()
            color_name = self.anom_dlg.combo_color.currentText()
            
            glow_color = QColor('#ff5555') if "Red" in color_name else QColor('#50fa7b') if "Green" in color_name else QColor('#f1fa8c')
            
            self.anomaly_rules_global.append({'col': col_target, 'op': operator, 'val': activator_value, 'color': glow_color})
            if self.table_model: self.table_model.illuminate_anomalies(self.anomaly_rules_global)
            self.anom_dlg.accept()

        def crystallize_current_mutations(self) -> None:
            crystal_path, _ = QFileDialog.getSaveFileName(self.ui, "Crystallize Memory", "preset.json", "JSON (*.json)")
            if crystal_path:
                active_conditions = [w.extract_dna_sequence() for w in self.ui.matryoshka_pods if w.line_value.text() or w.check_only_nan.isChecked()]
                Elephant_Brain_Storage.crystallize_memory(crystal_path, active_conditions)

        def resurrect_extinct_filters(self) -> None:
            ancient_path, _ = QFileDialog.getOpenFileName(self.ui, "Resurrect Memory", "", "JSON (*.json)")
            if ancient_path:
                stored_mutations = Elephant_Brain_Storage.revive_memory(ancient_path)
                for pod in self.ui.matryoshka_pods: 
                    self.ui.vaporize_matryoshka_pod(pod)
                self.ui.matryoshka_pods.clear()
                for condition in stored_mutations:
                    self.ui.spawn_matryoshka_pod()
                    self.ui.matryoshka_pods[-1].load_dna_sequence(condition)

        def check_mass_overload(self, file_path: str) -> bool:
            try:
                import psutil
                file_mass_bytes = Path(file_path).stat().st_size
                free_ram = psutil.virtual_memory().available
                
                if file_mass_bytes * 3 > free_ram:
                    reply = QMessageBox.warning(self.ui, "Mass Overload Alert", 
                                              "The selected genetic material exceeds safe RAM limits. The reactor may crash. Proceed anyway?",
                                              QMessageBox.Yes | QMessageBox.No)
                    if reply == QMessageBox.No: return False
            except ImportError:
                logging.warning("Psutil absent. Bypassing mass overload scan.")
            return True

        def verify_and_inject_dark_matter(self, file_path: str, format_type: str, has_header: bool) -> None:
            if not self.check_mass_overload(file_path): 
                self.ui.progress_bar.setVisible(False)
                self.ui.status_bar.showMessage("Injection aborted by operator.")
                return
            
            self.reader_thread = Hyperdrive_Loader_Thread(file_path, format_type, has_header)
            self.reader_thread.signal_done.connect(self.payload_secured)
            self.reader_thread.signal_error.connect(self.containment_breach)
            self.reader_thread.start()

        def payload_secured(self, extracted_biomass: pd.DataFrame) -> None:
            self.data_engine.inject_raw_material(extracted_biomass)
            
            self.table_model = Infinite_Hologram_Model(self.data_engine.filtered_biomass) 
            self.table_model.signal_crispr_edit.connect(self.execute_crispr_edit)
            self.ui.infinite_grid.setModel(self.table_model)
            self.table_model.illuminate_anomalies(self.anomaly_rules_global)
            
            if self.data_engine.raw_biomass is not None:
                detected_columns = [str(c) for c in self.data_engine.raw_biomass.columns]
                self.ui.synchronize_available_columns(detected_columns)
            
            self.ui.progress_bar.setVisible(False)
            self.ui.status_bar.showMessage(f"Quantum Scanner: {len(extracted_biomass)} sequences stabilized in memory.")
            self.update_ui_from_engine()

        def detonate_filter_collider(self, collision_conditions: List[Dict[str, Any]]) -> None:
            self.ui.progress_bar.setVisible(True)
            self.ui.progress_bar.setRange(0, 0)
            self.ui.status_bar.showMessage("Quantum Collider active...")
            
            self.collider_thread = Quantum_Collider_Thread(self.data_engine, collision_conditions)
            self.collider_thread.signal_done.connect(self.collider_finished)
            self.collider_thread.signal_error.connect(self.containment_breach)
            self.collider_thread.start()

        def collider_finished(self, mutant_strain: pd.DataFrame, genomic_stats: Dict[str, Any]) -> None:
            self.ui.progress_bar.setVisible(False)
            if self.table_model: self.table_model.mutate_hologram(mutant_strain)
            self.ui.update_xray_panel(genomic_stats)
            self.ui.status_bar.showMessage(f"Collision successful. {genomic_stats['surviving_reads']} survive.")

        def verify_disk_overwrite(self, dest_path: str) -> bool:
            if Path(dest_path).exists():
                reply = QMessageBox.question(self.ui, "Temporal Paradox Warning", 
                                           f"A timeline already exists at {Path(dest_path).name}. Obliterate it?", 
                                           QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.No: return False
            return True

        def synthesize_output(self, file_format: str) -> None:
            if self.data_engine.filtered_biomass is None or len(self.data_engine.filtered_biomass) == 0:
                QMessageBox.warning(self.ui, "Biological Warning", "Reactor empty.")
                return

            origin_path = self.ui.path_display.text()
            default_synthesis_name = Path(origin_path).stem + "_mutated"
            target_path, _ = QFileDialog.getSaveFileName(self.ui, f"Save Modified ({file_format.upper()})", f"{default_synthesis_name}.{file_format}")
            
            if target_path:
                if self.verify_disk_overwrite(target_path):
                    self.ui.progress_bar.setVisible(True)
                    self.ui.status_bar.showMessage("Synthesizing... Do not pull the plug.")
                    
                    export_thread = Chunk_Devourer_Exporter(self.data_engine.filtered_biomass, target_path, file_format)
                    self.active_synthesis_threads.append(export_thread)
                    export_thread.signal_progress.connect(self.ui.progress_bar.setValue)
                    export_thread.signal_done.connect(self.synthesis_complete)
                    export_thread.signal_error.connect(self.containment_breach)
                    export_thread.finished.connect(lambda: self.active_synthesis_threads.remove(export_thread))
                    export_thread.start()

        def synthesis_complete(self, completion_msg: str) -> None:
            self.ui.progress_bar.setVisible(False)
            self.ui.status_bar.showMessage(completion_msg)
            QMessageBox.information(self.ui, "Synthesis Successful", completion_msg)

        def containment_breach(self, error_string: str) -> None:
            self.ui.progress_bar.setVisible(False)
            self.ui.status_bar.showMessage("Dimensional Anomaly.")
            Runes_Of_Biocontainment.engrave_biocontainment_runes(f"ERROR: {error_string}")
            QMessageBox.critical(self.ui, "Critical Error", error_string)

        def execute_directive(self) -> None:
            self.ui.show()
            sys.exit(self.app.exec_())

    if __name__ == "__main__":
        Runes_Of_Biocontainment.engrave_biocontainment_runes("Session Initiated in Delimiter Laboratory")
        overlord = Supreme_Genomic_Overlord()
        overlord.execute_directive()

except Exception as catastrophic_failure:
    crash_log_location = Path.cwd() / "CRASH_REPORT_BIOLAB.txt"
    with open(crash_log_location, "w", encoding="utf-8") as panic_log:
        panic_log.write("=== GENOMIC REACTOR CONTAINMENT BREACH ===\n\n")
        panic_log.write(traceback.format_exc())
    logging.critical(f"CRITICAL EXPLOSION! Log forged at: {crash_log_location}")
    sys.exit(1)
