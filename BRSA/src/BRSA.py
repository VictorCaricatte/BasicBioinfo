# ==============================================================================
# Author:       Victor S Caricatte De Araújo (Frankestein)
# Institution:  Universidade federal de Minas Gerais (UFMG)
# Version:      0.9.0 
# Notes:        You need to have a little faith for it to work. 
# ==============================================================================

import sys
from args import summon_cli_oracle, SacredSessionScroll
from job import FrankensteinBioinformaticsOverlordThread
from biotools import scribe_reproducibility_chronicles

# Runs the routine strictly via CLI without loading any UI graphical modules
def engage_dark_cli_mode(args):
    print("[CLI EXECUTION] Starting the advanced analysis pipeline via terminal...")
    scroll = SacredSessionScroll(
        counts_matrix_path=args.counts,
        metadata_path=args.meta,
        alpha_condition=args.cond1,
        omega_condition=args.cond2,
        min_read_threshold=args.min_counts,
        fdr_cliff=args.pval,
        normalization_flavor=args.norm,
        batch_demon=args.batch,
        dimensionality_weapon=args.dim_red,
        summon_ml_golems=args.ml,
        summon_ppi_interactome=args.ppi,
        hunt_isoform_mutations=args.isoforms
    )
    
    cruncher = FrankensteinBioinformaticsOverlordThread(scroll)
    
    def catch_brainwaves(x): print(f"[PROCESS] Analysis alchemy progress: {x}%...")
    def colossal_failure(err): print(f"[ERROR] Critical failure: {err}"); sys.exit(1)
    def it_lives(treasures): 
        print("\n[SUCCESS] Data analysis complete.")
        sig = len(treasures['significant_genes'])
        print(f"[RESULTS] Number of significant genes (FDR < {scroll.fdr_cliff}): {sig}")
        treasures['de_results'].to_csv("CLI_differential_results.csv")
        print("[INFO] Results saved to: CLI_differential_results.csv")
        
        scribe_reproducibility_chronicles("CLI_reproducibility_log.txt", scroll.__dict__)
        print("[INFO] Reproducibility script saved to: CLI_reproducibility_log.txt")
        sys.exit(0)
        
    cruncher.brainwaves_updated.connect(catch_brainwaves)
    cruncher.lab_exploded.connect(colossal_failure)
    cruncher.monster_awakened.connect(it_lives)
    
    cruncher.run() 
    it_lives(cruncher.results_vault)

def main():
    magic_args = summon_cli_oracle()
    
    if magic_args.counts and magic_args.meta and magic_args.cond1 and magic_args.cond2:
        engage_dark_cli_mode(magic_args)
    else:
        print("\n[INFO] BRSA is currently operating in CLI Mode.")
        print("Required arguments are missing to perform a CLI analysis.")
        print("\nTo use the Command Line Interface, run:")
        print("  python BRSA.py --counts <file.csv> --meta <file.csv> --cond1 <A> --cond2 <B> --norm DESeq2 --dim_red UMAP --ml --ppi")
        print("\nTo launch the Advanced Graphical User Interface (GUI), run:")
        print("  python interface.py\n")

if __name__ == "__main__":
    main()
