import pandas as pd
import re
import csv
import json
import logging
import gc
from enum import Enum
from PyQt5.QtCore import QAbstractTableModel, Qt, QThread, pyqtSignal, QModelIndex
from typing import List, Dict, Any, Optional, Tuple

logging.basicConfig(filename='panvita_quarantine.log', level=logging.INFO,
                    format='%(asctime)s - BIOLOGICAL ALERT - %(message)s')

class QuantumOperators(Enum):
    EQUALS = "Equals"
    CONTAINS = "Contains"
    STARTS_WITH = "Starts With"
    ENDS_WITH = "Ends With"
    REGEX = "Regex"
    GREATER_THAN = "Greater Than (>)"
    LESS_THAN = "Less Than (<)"
    NUMERIC_EQUALS = "Numeric Equals (=)"

class TransmutationSpells(Enum):
    CONCAT = "Concat"
    ADD = "Add (+)"
    SUBTRACT = "Subtract (-)"
    MULTIPLY = "Multiply (*)"
    DIVIDE = "Divide (/)"

class RadiationMethods(Enum):
    MEAN = "Mean"
    MEDIAN = "Median"
    FIXED_VALUE = "Fixed Value"

class Runes_Of_Biocontainment:
    # Engraves the steps in digital stone
    @staticmethod
    def engrave_biocontainment_runes(log_message: str) -> None:
        logging.info(log_message)

class DarkMatter_Sniffer:
    @staticmethod
    def sniff_the_void(file_path: str) -> str:
        try:
            with open(file_path, 'r', encoding='utf-8') as genomic_file:
                sample_chunk = genomic_file.read(2048)
                dialect = csv.Sniffer().sniff(sample_chunk)
                return dialect.delimiter
        except Exception:
            if file_path.endswith('.tsv'): return '\t'
            return ','

class Regex_Sentinel:
    @staticmethod
    def inspect_anomaly(pattern: str) -> bool:
        try:
            re.compile(pattern)
            return True
        except re.error:
            return False

class Elephant_Brain_Storage:
    @staticmethod
    def crystallize_memory(path: str, memory_data: List[Dict[str, Any]]) -> None:
        with open(path, 'w', encoding='utf-8') as structural_file:
            json.dump(memory_data, structural_file, indent=4)
            
    @staticmethod
    def revive_memory(path: str) -> List[Dict[str, Any]]:
        import os
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as structural_file:
                return json.load(structural_file)
        return []

class Infinite_Hologram_Model(QAbstractTableModel):
    signal_crispr_edit = pyqtSignal(int, int, str)

    def __init__(self, biomass_dataframe: pd.DataFrame) -> None:
        super(Infinite_Hologram_Model, self).__init__()
        self._genomic_data = biomass_dataframe
        self._display_limit = 50
        self.anomaly_rules: List[Dict[str, Any]] = []

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return min(len(self._genomic_data), self._display_limit)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return len(self._genomic_data.columns)

    def canFetchMore(self, parent: QModelIndex = QModelIndex()) -> bool:
        return self._display_limit < len(self._genomic_data)

    def fetchMore(self, parent: QModelIndex = QModelIndex()) -> None:
        current_page = self._display_limit // 50
        self.engage_nanopore_pagination(current_page + 1)

    def engage_nanopore_pagination(self, page_number: int) -> None:
        # Transmits new nucleotide packets to the interface
        remainder = len(self._genomic_data) - self._display_limit
        items_to_fetch = min(50, remainder)
        if items_to_fetch > 0:
            self.beginInsertRows(QModelIndex(), self._display_limit, self._display_limit + items_to_fetch - 1)
            self._display_limit += items_to_fetch
            self.endInsertRows()

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:
        if index.isValid():
            nucleotide_val = self._genomic_data.iloc[index.row(), index.column()]
            col_name = self._genomic_data.columns[index.column()]
            
            if role == Qt.DisplayRole or role == Qt.EditRole:
                return str(nucleotide_val) if not pd.isna(nucleotide_val) else "NaN"
                
            if role == Qt.BackgroundRole:
                for rule in self.anomaly_rules:
                    if rule['col'] == col_name or rule['col'] == "ALL_COLUMNS":
                        try:
                            target = rule['val']
                            op = rule['op']
                            if op == QuantumOperators.EQUALS.value and str(nucleotide_val) == target:
                                return rule['color']
                            elif op == QuantumOperators.CONTAINS.value and target in str(nucleotide_val):
                                return rule['color']
                            elif op == QuantumOperators.GREATER_THAN.value and float(nucleotide_val) > float(target):
                                return rule['color']
                            elif op == QuantumOperators.LESS_THAN.value and float(nucleotide_val) < float(target):
                                return rule['color']
                        except Exception:
                            pass
        return None

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.EditRole) -> bool:
        if role == Qt.EditRole:
            self.signal_crispr_edit.emit(index.row(), index.column(), str(value))
            return True
        return False

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

    def headerData(self, col: int, orientation: Qt.Orientation, role: int = Qt.DisplayRole) -> Any:
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return str(self._genomic_data.columns[col])
        if orientation == Qt.Vertical and role == Qt.DisplayRole:
            return str(self._genomic_data.index[col])
        return None
    
    def mutate_hologram(self, new_biomass_data: pd.DataFrame) -> None:
        self.beginResetModel()
        self._genomic_data = new_biomass_data
        self._display_limit = min(50, len(new_biomass_data))
        self.endResetModel()

    def illuminate_anomalies(self, visual_rules: List[Dict[str, Any]]) -> None:

        self.anomaly_rules = visual_rules
        self.layoutChanged.emit()

class Hyperdrive_Loader_Thread(QThread):
    signal_done = pyqtSignal(pd.DataFrame)
    signal_error = pyqtSignal(str)

    def __init__(self, target_path: str, file_type: str, has_header: bool) -> None:
        super().__init__()
        self.target_path = target_path
        self.file_type = file_type
        self.has_header = has_header

    def run(self) -> None:
        try:
            header_idx = 0 if self.has_header else None
            separator = DarkMatter_Sniffer.sniff_the_void(self.target_path)

            if self.file_type == "Excel" or self.target_path.endswith(('.xlsx', '.xls')):
                extracted_biomass = pd.read_excel(self.target_path, header=header_idx)
            else:
                try:

                    extracted_biomass = pd.read_csv(self.target_path, sep=separator, header=header_idx, engine='pyarrow')
                except Exception:

                    extracted_biomass = pd.read_csv(self.target_path, sep=separator, header=header_idx, low_memory=False)
            

            if not extracted_biomass.columns.empty and header_idx is not None:
                extracted_biomass.columns = extracted_biomass.columns.astype(str).str.strip()

            self.signal_done.emit(extracted_biomass)
            
        except pd.errors.EmptyDataError:
            self.signal_error.emit("Biological Vacuum: The file contains no genetic material (Empty).")
        except UnicodeDecodeError:
            self.signal_error.emit("Alien Encoding: Failed to decode nucleotides. Check file encoding.")
        except Exception as catastrophic_anomaly:
            self.signal_error.emit(f"Critical Extraction Anomaly: {str(catastrophic_anomaly)}")

class Sterilized_Data_Engine:
    def __init__(self) -> None:
        self.raw_biomass: Optional[pd.DataFrame] = None
        self.filtered_biomass: Optional[pd.DataFrame] = None
        self.history_stack: List[pd.DataFrame] = []
        self.history_index: int = -1

    def inject_raw_material(self, incoming_biomass: pd.DataFrame) -> None:
        self.raw_biomass = incoming_biomass
        self.filtered_biomass = incoming_biomass.copy()
        self.history_stack = [self.filtered_biomass.copy()]
        self.history_index = 0

    def force_save_state(self) -> None:
        self.history_stack = self.history_stack[:self.history_index + 1]
        self.history_stack.append(self.filtered_biomass.copy())
        self.history_index += 1
        Runes_Of_Biocontainment.engrave_biocontainment_runes("Temporal Paradox Saved: DataFrame state captured.")

    def rewind_time_dilation(self) -> bool:
        if self.history_index > 0:
            self.history_index -= 1
            self.filtered_biomass = self.history_stack[self.history_index].copy()
            Runes_Of_Biocontainment.engrave_biocontainment_runes("Time Dilation: Undid the last mutation (Undo).")
            gc.collect()
            return True
        return False

    def temporal_redo(self) -> bool:
        if self.history_index < len(self.history_stack) - 1:
            self.history_index += 1
            self.filtered_biomass = self.history_stack[self.history_index].copy()
            Runes_Of_Biocontainment.engrave_biocontainment_runes("Time Dilation: Redid the mutation (Redo).")
            gc.collect()
            return True
        return False

    def purge_genetic_clones(self, target_columns: List[str]) -> None:
        initial_len = len(self.filtered_biomass)
        if not target_columns or "ALL_COLUMNS" in target_columns:
            self.filtered_biomass.drop_duplicates(inplace=True)
        else:
            self.filtered_biomass.drop_duplicates(subset=target_columns, inplace=True)
        
        self.force_save_state()
        gc.collect() # Explicit Memory Cleanup
        Runes_Of_Biocontainment.engrave_biocontainment_runes(f"Genetic Clones Purged: {initial_len - len(self.filtered_biomass)} destroyed.")

    def crispr_cas9_cell_edit(self, row_idx: int, col_idx: int, new_sequence: str) -> None:
        col_name = self.filtered_biomass.columns[col_idx]
        self.filtered_biomass.iat[row_idx, col_idx] = new_sequence
        self.force_save_state()
        Runes_Of_Biocontainment.engrave_biocontainment_runes(f"CRISPR-Cas9 Edit: Row {row_idx}, Column '{col_name}' altered to '{new_sequence}'.")

    def quarantine_nan_ghosts(self, col_target: str, filler_value: str, method: str) -> None:
        if method == RadiationMethods.MEAN.value:
            synthesized_val = pd.to_numeric(self.filtered_biomass[col_target], errors='coerce').mean()
        elif method == RadiationMethods.MEDIAN.value:
            synthesized_val = pd.to_numeric(self.filtered_biomass[col_target], errors='coerce').median()
        else:
            synthesized_val = filler_value
        self.filtered_biomass[col_target] = self.filtered_biomass[col_target].fillna(synthesized_val)
        self.force_save_state()
        Runes_Of_Biocontainment.engrave_biocontainment_runes(f"NaN Quarantine: Filled column {col_target} with '{synthesized_val}'.")

    def transcribe_rna_to_dna(self, new_strand_name: str, strand_a: str, strand_b: str, operation: str) -> None:
        if operation == TransmutationSpells.CONCAT.value:
            self.filtered_biomass[new_strand_name] = self.filtered_biomass[strand_a].astype(str) + "_" + self.filtered_biomass[strand_b].astype(str)
        else:
            sequence_1 = pd.to_numeric(self.filtered_biomass[strand_a], errors='coerce')
            sequence_2 = pd.to_numeric(self.filtered_biomass[strand_b], errors='coerce')
            if operation == TransmutationSpells.ADD.value: self.filtered_biomass[new_strand_name] = sequence_1 + sequence_2
            elif operation == TransmutationSpells.SUBTRACT.value: self.filtered_biomass[new_strand_name] = sequence_1 - sequence_2
            elif operation == TransmutationSpells.MULTIPLY.value: self.filtered_biomass[new_strand_name] = sequence_1 * sequence_2
            elif operation == TransmutationSpells.DIVIDE.value: self.filtered_biomass[new_strand_name] = sequence_1 / sequence_2
        self.force_save_state()
        Runes_Of_Biocontainment.engrave_biocontainment_runes(f"RNA->DNA Transcription: Created strand '{new_strand_name}' from '{strand_a}' and '{strand_b}'.")

    def splice_alien_dna_strands(self, foreign_biomass: pd.DataFrame, alien_key: str, local_key: str) -> None:
        self.filtered_biomass = pd.merge(self.filtered_biomass, foreign_biomass, left_on=local_key, right_on=alien_key, how='inner')
        self.force_save_state()
        gc.collect()
        Runes_Of_Biocontainment.engrave_biocontainment_runes(f"Alien Splicing: Fusion complete using local key '{local_key}'.")

    def unleash_matryoshka_filters(self, mutation_conditions: List[Dict[str, Any]]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        if self.filtered_biomass is None:
            raise ValueError("The reactor is empty.")

        mutant_strain = self.filtered_biomass.copy()
        
        for condition_pod in mutation_conditions:
            col = condition_pod['column']
            operation = condition_pod['operation']
            target_val = condition_pod['value']
            is_case_sensitive = condition_pod['case_sensitive']
            ignore_nan = condition_pod['ignore_nan']
            only_nan = condition_pod['only_nan']
            
            if only_nan:
                mutant_strain = mutant_strain[mutant_strain[col].isna()]
                continue
                
            if ignore_nan:
                valid_nucleotide_mask = mutant_strain[col].notna()
                mutant_strain = mutant_strain[valid_nucleotide_mask]

            if operation in [QuantumOperators.GREATER_THAN.value, QuantumOperators.LESS_THAN.value, QuantumOperators.NUMERIC_EQUALS.value]:
                try:
                    num_val = float(target_val)
                    mutant_strain[col] = pd.to_numeric(mutant_strain[col], errors='coerce')
                    if operation == QuantumOperators.GREATER_THAN.value: mutant_strain = mutant_strain[mutant_strain[col] > num_val]
                    elif operation == QuantumOperators.LESS_THAN.value: mutant_strain = mutant_strain[mutant_strain[col] < num_val]
                    elif operation == QuantumOperators.NUMERIC_EQUALS.value: mutant_strain = mutant_strain[mutant_strain[col] == num_val]
                except ValueError:
                    pass 
            else:
                str_sequence = mutant_strain[col].astype(str)
                if not is_case_sensitive and operation != QuantumOperators.REGEX.value:
                    str_sequence = str_sequence.str.lower()
                    target_val = target_val.lower()

                if operation == QuantumOperators.EQUALS.value:
                    mutant_strain = mutant_strain[str_sequence == target_val]
                elif operation == QuantumOperators.CONTAINS.value:
                    mutant_strain = mutant_strain[str_sequence.str.contains(target_val, regex=False)]
                elif operation == QuantumOperators.STARTS_WITH.value:
                    mutant_strain = mutant_strain[str_sequence.str.startswith(target_val)]
                elif operation == QuantumOperators.ENDS_WITH.value:
                    mutant_strain = mutant_strain[str_sequence.str.endswith(target_val)]
                elif operation == QuantumOperators.REGEX.value:
                    
                    flags = 0 if is_case_sensitive else re.IGNORECASE
                    compiled_pattern = re.compile(target_val, flags=flags)
                    
                    mutant_strain = mutant_strain[mutant_strain[col].astype(str).apply(lambda strand: bool(compiled_pattern.search(strand)))]

        self.filtered_biomass = mutant_strain
        self.force_save_state()
        gc.collect() 
        Runes_Of_Biocontainment.engrave_biocontainment_runes(f"Filter Collision: Applied {len(mutation_conditions)} Matryoshka pods.")

        survival_rate = 0.0
        if len(self.raw_biomass) > 0:
            survival_rate = round((len(self.filtered_biomass) / len(self.raw_biomass)) * 100, 2)

        genomic_stats = {
            "total_reads": len(self.raw_biomass),
            "surviving_reads": len(self.filtered_biomass),
            "survival_rate": survival_rate
        }
        return self.filtered_biomass, genomic_stats

class Quantum_Collider_Thread(QThread):
    signal_done = pyqtSignal(pd.DataFrame, dict)
    signal_error = pyqtSignal(str)

    def __init__(self, processing_engine: Sterilized_Data_Engine, mutation_conditions: List[Dict[str, Any]]) -> None:
        super().__init__()
        self.engine = processing_engine
        self.conditions = mutation_conditions

    def run(self) -> None:
        try:
            mutant_strain, genomic_stats = self.engine.unleash_matryoshka_filters(self.conditions)
            self.signal_done.emit(mutant_strain, genomic_stats)
        except Exception as anomaly:
            self.signal_error.emit(str(anomaly))

class Chunk_Devourer_Exporter(QThread):
    signal_progress = pyqtSignal(int)
    signal_done = pyqtSignal(str)
    signal_error = pyqtSignal(str)

    def __init__(self, export_biomass: pd.DataFrame, dest_path: str, format_type: str) -> None:
        super().__init__()
        self.df = export_biomass
        self.dest_path = dest_path
        self.format_type = format_type

    def run(self) -> None:
        try:
            chunk_size = 50000
            total_rows = len(self.df)
            
            if total_rows == 0:
                self.signal_error.emit("Biological Warning: No data to export.")
                return

            if self.format_type == "excel":
                self.df.to_excel(self.dest_path, index=False)
                self.signal_progress.emit(100)
            else:
                separator = ',' if self.format_type == 'csv' else '\t'
                for i in range(0, total_rows, chunk_size):
                    chunk_strain = self.df.iloc[i:i+chunk_size]
                    write_mode = 'w' if i == 0 else 'a'
                    write_header = True if i == 0 else False
                    chunk_strain.to_csv(self.dest_path, sep=separator, index=False, mode=write_mode, header=write_header)
                    
                    progress_val = int(((i + chunk_size) / total_rows) * 100)
                    self.signal_progress.emit(min(progress_val, 100))
            
            Runes_Of_Biocontainment.engrave_biocontainment_runes(f"Exportation: File synthesized as {self.format_type.upper()}.")
            self.signal_done.emit(f"Exportation of {total_rows} sequences complete.")
        except Exception as synthesis_anomaly:
            self.signal_error.emit(f"Output Synthesis Failure: {str(synthesis_anomaly)}")

class Fasta_Devourer_Exporter(QThread):
    signal_progress = pyqtSignal(int)
    signal_done = pyqtSignal(str)
    signal_error = pyqtSignal(str)

    def __init__(self, export_biomass: pd.DataFrame, dest_path: str, id_col: str, seq_col: str) -> None:
        super().__init__()
        self.df = export_biomass
        self.dest_path = dest_path
        self.id_col = id_col
        self.seq_col = seq_col

    def run(self) -> None:
        try:
            total_rows = len(self.df)
            if total_rows == 0:
                self.signal_error.emit("Warning: Empty DataFrame for FASTA conversion.")
                return

            with open(self.dest_path, 'w', encoding='utf-8') as fasta_file:
                for idx, row in enumerate(self.df.itertuples(index=False)):
                    sequence_header = getattr(row, self.id_col)
                    nucleotide_seq = getattr(row, self.seq_col)
                    fasta_file.write(f">{sequence_header}\n{nucleotide_seq}\n")
                    if idx % 1000 == 0:
                        progress_val = int((idx / total_rows) * 100)
                        self.signal_progress.emit(progress_val)
                        
            self.signal_progress.emit(100)
            Runes_Of_Biocontainment.engrave_biocontainment_runes("FASTA envelope successfully synthesized.")
            self.signal_done.emit(f"FASTA compiled with {total_rows} sequences.")
        except Exception as fasta_anomaly:
            self.signal_error.emit(f"FASTA Synthesis Failure: {str(fasta_anomaly)}")
