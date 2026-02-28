from PyQt6.QtCore import QThread, pyqtSignal
import pandas as pd
import numpy as np
import os
import re
import traceback
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from statsmodels.stats.multitest import multipletests
import gseapy as gp

# Advanced Necromancy Imports
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    print("[SYSTEM] PyDESeq2 not found. The old statistical gods will mourn.")

# CORREÇÃO AQUI: Invocação direta do pycombat
try:
    from pycombat import pycombat
except ImportError:
    print("[SYSTEM] PyComBat not found. Batch exorcism will be primitive.")

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist
except ImportError:
    pass

try:
    from lifelines import KaplanMeierFitter
except ImportError:
    pass

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
except ImportError:
    print("[SYSTEM] Scikit-learn not found. The Machine Learning golems will remain asleep.")

from biotools import (summon_fastqc_inquisitor, execute_fastp_guillotine, cast_hisat2_alignment_spell, 
                      invoke_samtools_sorting_ritual, reap_featurecounts_harvest,
                      unleash_kallisto_quantification, conjure_multiqc_all_seeing_eye,
                      summon_string_protein_web, invoke_interactome_spirits)

# Cleanses the dirty sample names so Pandas doesn't cry during merges
def exorcise_dirty_strings(text: str) -> str:
    return re.sub(r'[^a-zA-Z0-9_]', '_', str(text).strip())

# A ravenous thread that eats FASTQ files in the background so the GUI doesn't freeze to death
class FastqDevourerThread(QThread):
    progress_echoed = pyqtSignal(int, str)
    digestion_completed = pyqtSignal(str)
    fatal_choke = pyqtSignal(str)

    def __init__(self, fastq_arsenal, output_bunker, reference_index, annotation_file, is_gff3, perform_trimming, aligner="HISAT2", threads=4):
        super().__init__()
        self.fastq_arsenal = fastq_arsenal
        self.output_bunker = output_bunker
        self.reference_index = reference_index
        self.annotation_file = annotation_file
        self.is_gff3 = is_gff3
        self.perform_trimming = perform_trimming
        self.aligner = aligner
        self.threads = threads

    def run(self):
        try:
            os.makedirs(self.output_bunker, exist_ok=True)
            sample_counts = {}
            gene_lengths = None
            total_victims = len(self.fastq_arsenal)
            
            for i, (raw_sample_tag, fq_paths) in enumerate(self.fastq_arsenal.items()):
                sample_tag = exorcise_dirty_strings(raw_sample_tag)
                base_progress = int((i/total_victims)*80)
                
                # Throw to FastQC Inquisitor
                self.progress_echoed.emit(base_progress, f"Interrogating raw sequences for {sample_tag} via FastQC...")
                summon_fastqc_inquisitor(fq_paths, self.output_bunker, sample_tag, self.threads)
                
                sorted_bam = os.path.join(self.output_bunker, f"{sample_tag}.sorted.bam")
                kallisto_out = os.path.join(self.output_bunker, f"{sample_tag}_kallisto", "abundance.tsv")
                counts_file = os.path.join(self.output_bunker, f"{sample_tag}.counts.txt")
                
                if self.aligner == "HISAT2" and os.path.exists(counts_file):
                    self.progress_echoed.emit(base_progress + 15, f"Cache found! Skipping alignment for {sample_tag}...")
                    df_counts = pd.read_csv(counts_file, sep='\t', comment='#', index_col=0, header=0)
                    if gene_lengths is None and 'Length' in df_counts.columns:
                        gene_lengths = df_counts['Length']
                    sample_counts[sample_tag] = df_counts.iloc[:, -1]
                    continue
                elif self.aligner == "Kallisto" and os.path.exists(kallisto_out):
                    self.progress_echoed.emit(base_progress + 15, f"Cache found! Skipping pseudo-alignment for {sample_tag}...")
                    df_kal = pd.read_csv(kallisto_out, sep='\t', index_col=0)
                    if gene_lengths is None and 'eff_length' in df_kal.columns:
                        gene_lengths = df_kal['eff_length']
                    sample_counts[sample_tag] = df_kal['est_counts']
                    continue

                active_reads = fq_paths
                if self.perform_trimming:
                    self.progress_echoed.emit(base_progress + 5, f"Running fastp guillotine for {sample_tag}...")
                    active_reads = execute_fastp_guillotine(fq_paths, self.output_bunker, sample_tag, self.threads)
                
                if self.aligner == "HISAT2":
                    sam_file = os.path.join(self.output_bunker, f"{sample_tag}.sam")
                    bam_file = os.path.join(self.output_bunker, f"{sample_tag}.bam")
                    
                    self.progress_echoed.emit(base_progress + 10, f"Casting HISAT2 alignment spell for {sample_tag}...")
                    cast_hisat2_alignment_spell(active_reads, self.reference_index, sam_file, self.threads)
                    
                    self.progress_echoed.emit(base_progress + 15, f"Invoking samtools sorting ritual for {sample_tag}...")
                    invoke_samtools_sorting_ritual(sam_file, bam_file, sorted_bam)
                    
                    self.progress_echoed.emit(base_progress + 20, f"Reaping featureCounts harvest for {sample_tag}...")
                    reap_featurecounts_harvest(sorted_bam, self.annotation_file, counts_file, self.threads, self.is_gff3)
                    
                    df_counts = pd.read_csv(counts_file, sep='\t', comment='#', index_col=0, header=0)
                    if gene_lengths is None and 'Length' in df_counts.columns:
                        gene_lengths = df_counts['Length']
                    sample_counts[sample_tag] = df_counts.iloc[:, -1]
                    
                    if os.path.exists(sam_file): os.remove(sam_file)
                    if os.path.exists(bam_file): os.remove(bam_file)
                    
                elif self.aligner == "Kallisto":
                    self.progress_echoed.emit(base_progress + 10, f"Unleashing Kallisto ghost for {sample_tag}...")
                    abundance_file = unleash_kallisto_quantification(active_reads, self.reference_index, self.output_bunker, sample_tag, self.threads)
                    df_kal = pd.read_csv(abundance_file, sep='\t', index_col=0)
                    if gene_lengths is None and 'eff_length' in df_kal.columns:
                        gene_lengths = df_kal['eff_length']
                    sample_counts[sample_tag] = df_kal['est_counts']
            
            self.progress_echoed.emit(85, "Conjuring the all-seeing eye of MultiQC...")
            conjure_multiqc_all_seeing_eye(self.output_bunker)

            self.progress_echoed.emit(95, "Merging counts into master matrix...")
            final_df = pd.DataFrame(sample_counts)
            if gene_lengths is not None:
                final_df.insert(0, 'GeneLength', gene_lengths)
                
            master_file = os.path.join(self.output_bunker, "master_counts_compiled.csv")
            final_df.to_csv(master_file)
            
            self.progress_echoed.emit(100, "Processing complete.")
            self.digestion_completed.emit(master_file)
            
        except Exception as e:
            self.fatal_choke.emit(f"Critical error during FASTQ processing: {str(e)}\n{traceback.format_exc()}")


# The mad scientist thread that crunches numbers and brings the statistical monster to life
class FrankensteinBioinformaticsOverlordThread(QThread):
    brainwaves_updated = pyqtSignal(int)
    monster_awakened = pyqtSignal(dict)
    lab_exploded = pyqtSignal(str)

    def __init__(self, session_scroll):
        super().__init__()
        self.session_scroll = session_scroll
        self.results_vault = {}

    def perform_advanced_batch_exorcism(self, data: pd.DataFrame, metadata: pd.DataFrame, batch_col: str) -> pd.DataFrame:
        if batch_col not in metadata.columns:
            return data
        
        print(f"[RITUAL] Performing PyComBat advanced batch exorcism on covariate: {batch_col}")
        batch_labels = metadata.set_index('sample').loc[data.columns, batch_col].tolist()
        try:
            purified_data = pycombat(data, batch_labels)
            return purified_data.clip(lower=0)
        except Exception as e:
            print(f"[WARNING] PyComBat failed, falling back to primitive median centering: {e}")
            corrected_data = data.copy()
            batches = metadata[batch_col].unique()
            global_mean = data.mean(axis=1)
            for batch in batches:
                samples_in_batch = metadata[metadata[batch_col] == batch]['sample'].values
                valid_samples = [s for s in samples_in_batch if s in data.columns]
                if not valid_samples: continue
                batch_mean = data[valid_samples].mean(axis=1)
                corrected_data[valid_samples] = data[valid_samples].subtract(batch_mean, axis=0).add(global_mean, axis=0)
            return corrected_data.clip(lower=0)

    def summon_pydeseq2_necromancer(self, counts_df, meta_df):
        print("[MAGIC] Summoning PyDESeq2 Necromancer for empirical Bayes shrinkage...")
        dds = DeseqDataSet(
            counts=counts_df.T,
            metadata=meta_df,
            design_factors="condition"
        )
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=("condition", self.session_scroll.alpha_condition, self.session_scroll.omega_condition))
        stat_res.summary()
        
        res_df = stat_res.results_df
        res_df = res_df.rename(columns={'log2FoldChange': 'log2fc', 'pvalue': 'p-value', 'padj': 'p-adj', 'baseMean': 'mean_exp'})
        res_df['significant'] = res_df['p-adj'] < self.session_scroll.fdr_cliff
        
        # Approximate mean_1 and mean_2 for legacy compatibility in UI
        alpha_samples = meta_df[meta_df['condition'] == self.session_scroll.alpha_condition].index
        omega_samples = meta_df[meta_df['condition'] == self.session_scroll.omega_condition].index
        norm_counts = dds.layers["log1p"] # log1p normalized counts
        
        res_df['mean_1'] = pd.DataFrame(norm_counts).loc[alpha_samples].mean(axis=0).values
        res_df['mean_2'] = pd.DataFrame(norm_counts).loc[omega_samples].mean(axis=0).values
        
        return res_df, pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T

    def weave_wgcna_dark_tapestry(self, normalized_matrix):
        print("[WEAVING] Weaving the dark tapestry of WGCNA Co-expression...")
        corr_matrix = np.corrcoef(normalized_matrix.values)
        power = 6
        adj_matrix = np.power(np.abs(corr_matrix), power)
        dist_matrix = 1 - adj_matrix
        dist_array = dist_matrix[np.triu_indices(dist_matrix.shape[0], k=1)]
        Z = linkage(dist_array, method='average')
        return Z

    def extract_coexpression_souls(self, linkage_matrix, threshold):
        print("[SOUL EXTRACTION] Grouping genes into co-expression modules...")
        from scipy.cluster.hierarchy import fcluster
        clusters = fcluster(linkage_matrix, t=threshold, criterion='distance')
        return clusters

    def consult_kaplan_meier_oracle(self, metadata):
        time_rune = self.session_scroll.survival_time_rune
        event_rune = self.session_scroll.survival_event_rune
        if time_rune in metadata.columns and event_rune in metadata.columns:
            return True
        return False

    def calculate_hazard_ratio_curse(self):
        print("[RITUAL] Calculating the Hazard Ratio Curse is reserved for the plotting dimension.")
        pass

    def train_random_forest_golems(self, training_matrix, labels):
        print("[MACHINE LEARNING] Training Random Forest Golems...")
        rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
        rf.fit(training_matrix, labels)
        return rf

    def forge_svm_crystal_ball(self, training_matrix, labels):
        print("[MACHINE LEARNING] Forging the SVM Crystal Ball...")
        svm = SVC(kernel='linear', probability=True, random_state=42)
        svm.fit(training_matrix, labels)
        return svm

    def splice_rna_frankenstein_limbs(self, expression_matrix):
        print("[SPLICING] Uncovering hidden isoform mutations in the Frankenstein limbs...")
        # Simulates differential isoform usage by measuring extreme variance across samples
        isoform_variance = expression_matrix.var(axis=1)
        suspicious_isoforms = isoform_variance.sort_values(ascending=False).head(50).index.tolist()
        return suspicious_isoforms

    def run(self):
        try:
            self.brainwaves_updated.emit(5)
            raw_counts = pd.read_csv(self.session_scroll.counts_matrix_path, index_col=0)
            metadata_df = pd.read_csv(self.session_scroll.metadata_path)
            
            metadata_df['sample'] = metadata_df['sample'].apply(exorcise_dirty_strings)
            raw_counts.columns = [exorcise_dirty_strings(c) if c != 'GeneLength' else c for c in raw_counts.columns]
            
            gene_lengths = None
            if 'GeneLength' in raw_counts.columns:
                gene_lengths = raw_counts['GeneLength']
                raw_counts = raw_counts.drop(columns=['GeneLength'])
                
            self.results_vault['metadata'] = metadata_df
            
            valid_samples = [s for s in metadata_df['sample'] if s in raw_counts.columns]
            metadata_df = metadata_df[metadata_df['sample'].isin(valid_samples)].set_index('sample')
            raw_counts = raw_counts[valid_samples]
            
            raw_counts = raw_counts.round().astype(int)
            filtered_counts = raw_counts.loc[raw_counts.sum(axis=1) >= self.session_scroll.min_read_threshold]
            self.results_vault['raw_counts'] = filtered_counts
            
            self.brainwaves_updated.emit(20)
            
            try:
                de_results, normalized_matrix = self.summon_pydeseq2_necromancer(filtered_counts, metadata_df)
            except Exception as e:
                print(f"[FALLBACK] PyDESeq2 failed ({e}), using traditional CPM.")
                normalized_matrix = filtered_counts.div(filtered_counts.sum(axis=0), axis=1) * 1e6
                from bioinfokit.analys import stat
                stat_analysis = stat()
                stat_analysis.df = normalized_matrix.T
                stat_analysis.df['group'] = metadata_df.loc[normalized_matrix.columns, 'condition'].values
                stat_analysis.t_test('group', self.session_scroll.alpha_condition, self.session_scroll.omega_condition)
                de_results = stat_analysis.result
                de_results['gene'] = normalized_matrix.index
                de_results = de_results.set_index('gene')
                de_results['log2fc'] = np.log2((de_results['mean_1'] + 1e-5) / (de_results['mean_2'] + 1e-5))
                rejected, corrected_pvals, _, _ = multipletests(de_results['p-value'], alpha=self.session_scroll.fdr_cliff, method='fdr_bh')
                de_results['p-adj'] = corrected_pvals
                de_results['significant'] = rejected
                
            if self.session_scroll.batch_demon and self.session_scroll.batch_demon != "None":
                normalized_matrix = self.perform_advanced_batch_exorcism(normalized_matrix, metadata_df.reset_index(), self.session_scroll.batch_demon)
                
            self.results_vault['normalized_counts'] = normalized_matrix
            self.results_vault['de_results'] = de_results
            self.results_vault['significant_genes'] = de_results[de_results['significant']].index.tolist()
            
            self.brainwaves_updated.emit(50)
            
            # Dimensionality Reduction Oracles
            X = np.log2(normalized_matrix.T + 1)
            n_samples = X.shape[0]
            
            if self.session_scroll.dimensionality_weapon == "PCA":
                from sklearn.preprocessing import StandardScaler
                X_scaled = StandardScaler().fit_transform(X)
                pca = PCA(n_components=min(2, n_samples))
                coords = pca.fit_transform(X_scaled)
                var_ratio = pca.explained_variance_ratio_ if n_samples > 1 else [0,0]
                self.results_vault['dim_red'] = {'type': 'PCA', 'components': coords, 'explained_variance': var_ratio}
                
            elif self.session_scroll.dimensionality_weapon == "t-SNE":
                perp = min(5, n_samples - 1) if n_samples > 1 else 1
                tsne = TSNE(n_components=2, perplexity=perp, random_state=42)
                coords = tsne.fit_transform(X)
                self.results_vault['dim_red'] = {'type': 't-SNE', 'components': coords}
                
            elif self.session_scroll.dimensionality_weapon == "UMAP":
                n_neighbors = min(15, n_samples - 1) if n_samples > 2 else 2
                reducer = umap.UMAP(n_neighbors=n_neighbors, random_state=42)
                coords = reducer.fit_transform(X)
                self.results_vault['dim_red'] = {'type': 'UMAP', 'components': coords}
            
            self.brainwaves_updated.emit(65)
            # Full Hierarchical Matrix Prep
            top_mutant_genes = normalized_matrix.var(axis=1).sort_values(ascending=False).head(100).index
            prepped_heatmap_matrix = np.log2(normalized_matrix.loc[top_mutant_genes] + 1)
            self.results_vault['heatmap_matrix'] = prepped_heatmap_matrix
            
            # WGCNA Dark Tapestry
            try:
                top_wgcna_genes = normalized_matrix.var(axis=1).sort_values(ascending=False).head(500).index
                wgcna_matrix = np.log2(normalized_matrix.loc[top_wgcna_genes] + 1)
                z_linkage = self.weave_wgcna_dark_tapestry(wgcna_matrix)
                modules = self.extract_coexpression_souls(z_linkage, 0.5)
                self.results_vault['wgcna_linkage'] = z_linkage
                self.results_vault['wgcna_genes'] = top_wgcna_genes
                self.results_vault['wgcna_modules'] = modules
            except Exception as e:
                print(f"[WARNING] Dark tapestry failed to weave: {e}")

            # Machine Learning Golems
            if self.session_scroll.summon_ml_golems:
                try:
                    labels = metadata_df.loc[normalized_matrix.columns, 'condition'].values
                    rf_model = self.train_random_forest_golems(X, labels)
                    svm_model = self.forge_svm_crystal_ball(X, labels)
                    
                    # Interrogate feature importances
                    importances = pd.Series(rf_model.feature_importances_, index=normalized_matrix.index)
                    self.results_vault['ml_rf_importance'] = importances.sort_values(ascending=False).head(20)
                    self.results_vault['ml_models_ready'] = True
                except Exception as e:
                    print(f"[WARNING] ML Golems refused to wake up: {e}")

            # Isoform Mutations
            if self.session_scroll.hunt_isoform_mutations:
                suspicious_isoforms = self.splice_rna_frankenstein_limbs(normalized_matrix)
                self.results_vault['frankenstein_isoforms'] = suspicious_isoforms

            # Protein Interactome
            if self.session_scroll.summon_ppi_interactome and len(self.results_vault['significant_genes']) > 0:
                raw_string_data = summon_string_protein_web(self.results_vault['significant_genes'])
                ppi_graph = invoke_interactome_spirits(raw_string_data)
                self.results_vault['ppi_network'] = ppi_graph

            # Survival Doom Curve Check
            self.results_vault['has_survival'] = self.consult_kaplan_meier_oracle(metadata_df)
            
            self.brainwaves_updated.emit(85)
            # The True GSEA Ritual
            try:
                ranked_list = de_results[['log2fc']].dropna().sort_values(by='log2fc', ascending=False)
                gsea_res = gp.prerank(rnk=ranked_list, gene_sets=self.session_scroll.gsea_database_relic, 
                                      threads=4, min_size=5, max_size=1000, permutation_num=1000, 
                                      outdir=None, seed=6, verbose=True)
                self.results_vault['gsea_results'] = gsea_res.res2d
                self.results_vault['gsea_obj'] = gsea_res
            except Exception as e:
                print(f"GSEA failed (probably network or gene name mismatch), skipping: {e}")
                self.results_vault['gsea_results'] = pd.DataFrame()

            self.brainwaves_updated.emit(100)
            self.monster_awakened.emit(self.results_vault)
            
        except Exception as e:
            self.lab_exploded.emit(f"Fatal analysis error: {str(e)}\n{traceback.format_exc()}")
