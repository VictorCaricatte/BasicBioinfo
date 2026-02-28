import argparse
from dataclasses import dataclass, field
from typing import Dict, Any

# The sacred scroll holding the session's soul. Keeps track of all configurations.
@dataclass
class SacredSessionScroll:
    counts_matrix_path: str = ""
    metadata_path: str = ""
    alpha_condition: str = ""
    omega_condition: str = ""
    min_read_threshold: int = 10
    fdr_cliff: float = 0.05
    normalization_flavor: str = "CPM"
    batch_demon: str = "None"
    aligner_choice: str = "HISAT2"
    dimensionality_weapon: str = "PCA" 
    gsea_database_relic: str = "GO_Biological_Process_2021"
    visual_palette: str = "mako" 
    survival_time_rune: str = "time"
    survival_event_rune: str = "event"
    summon_ml_golems: bool = False
    summon_ppi_interactome: bool = False
    hunt_isoform_mutations: bool = False
    harvested_treasures: Dict[str, Any] = field(default_factory=dict)

# Summons the command line interface deity to parse mortal inputs
def summon_cli_oracle() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="BRSA - Basic RNA-Seq Analysis")
    parser.add_argument("--counts", type=str, help="Path to the counts CSV file.")
    parser.add_argument("--meta", type=str, help="Path to the metadata CSV file.")
    parser.add_argument("--cond1", type=str, help="First condition for the test (Group A).")
    parser.add_argument("--cond2", type=str, help="Second condition for the test (Group B).")
    parser.add_argument("--min_counts", type=int, default=10, help="Minimum read count threshold per gene.")
    parser.add_argument("--pval", type=float, default=0.05, help="False Discovery Rate (FDR) cutoff.")
    parser.add_argument("--norm", type=str, choices=["CPM", "TPM", "DESeq2"], default="DESeq2", help="Normalization method.")
    parser.add_argument("--batch", type=str, default="None", help="Covariate column name for Batch Effect correction.")
    parser.add_argument("--dim_red", type=str, choices=["PCA", "t-SNE", "UMAP"], default="PCA", help="Dimensionality reduction algorithm.")
    parser.add_argument("--ml", action="store_true", help="Train Random Forest Golems and SVM Crystal Balls.")
    parser.add_argument("--ppi", action="store_true", help="Invoke protein-protein interactome spirits.")
    parser.add_argument("--isoforms", action="store_true", help="Hunt for alternative splicing Frankenstein limbs.")
    
    return parser.parse_args()
