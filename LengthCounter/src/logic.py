import os
import csv
import gzip
import re
import math
import base64
import subprocess
import tempfile
from io import BytesIO
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.signal import convolve2d

from Bio import SeqIO, Align, AlignIO, Phylo, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII, XhoI, NotI
from Bio.Blast import NCBIXML

from config import SCROLL_OF_COLORS

Entrez.email = "ninja_shinobi@hiddenleaf.konoha"

def shikamaru_shadow_stats(lengths: List[int], expected_genome_size: int = 0) -> Dict[str, float]:
    if not lengths: return {}
    arr = np.array(lengths)
    sorted_arr = np.sort(arr)[::-1]
    total_len = np.sum(sorted_arr)
    cum_sum = np.cumsum(sorted_arr)
    
    n50_idx = np.searchsorted(cum_sum, total_len * 0.5)
    n50 = sorted_arr[n50_idx] if n50_idx < len(sorted_arr) else 0
    l50 = n50_idx + 1

    n90_idx = np.searchsorted(cum_sum, total_len * 0.9)
    n90 = sorted_arr[n90_idx] if n90_idx < len(sorted_arr) else 0
    l90 = n90_idx + 1

    stats_dict = {
        'count': len(lengths), 'min': int(np.min(arr)), 'max': int(np.max(arr)),
        'mean': float(np.mean(arr)), 'median': float(np.median(arr)), 'std': float(np.std(arr)),
        'q1': float(np.percentile(arr, 25)), 'q3': float(np.percentile(arr, 75)),
        'N50': int(n50), 'L50': int(l50), 'N90': int(n90), 'L90': int(l90)
    }

    if expected_genome_size > 0:
        ng50_idx = np.searchsorted(cum_sum, expected_genome_size * 0.5)
        ng50 = sorted_arr[ng50_idx] if ng50_idx < len(sorted_arr) else 0
        stats_dict['NG50'] = int(ng50)
        stats_dict['LG50'] = int(ng50_idx + 1)

    return stats_dict

def ino_mind_transfer_stats_text(stats_dict: Dict[str, Dict[str, float]]) -> str:
    if not stats_dict: return "No statistics data available for report generation."
    text = []
    for label, stats in stats_dict.items():
        text.append(f"Statistical Summary for Dataset: {label}")
        text.append("-" * 45)
        for stat, value in stats.items():
            if stat in ['count', 'min', 'max', 'N50', 'L50', 'N90', 'L90', 'NG50', 'LG50']:
                text.append(f"{stat:>8}: {int(value):12d}")
            else:
                text.append(f"{stat:>8}: {value:12.2f}")
        text.append("")
    return "\n".join(text)

def neji_byakugan_format_sniffer(file_path: str) -> str:
    lower_path = file_path.lower()
    if lower_path.endswith('.vcf') or lower_path.endswith('.vcf.gz'): return 'vcf'
    if lower_path.endswith('.gff') or lower_path.endswith('.gff3'): return 'gff'
    if lower_path.endswith('.bed'): return 'bed'
    if lower_path.endswith('.sam'): return 'sam'
    if lower_path.endswith('.bam'): return 'bam'
    if lower_path.endswith('.gb') or lower_path.endswith('.gbk') or lower_path.endswith('.gbff'): return 'genbank'
        
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    try:
        with open_func(file_path, mode) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                if line.startswith('>'): return 'fasta'
                if line.startswith('@'): return 'fastq'
                raise ValueError("Unrecognized sequence file signature.")
    except Exception as e:
        raise ValueError(f"I/O Error processing file {file_path}: {str(e)}")

def gamabunta_summon_sequence_file(file_path: str) -> str:
    if not os.path.exists(file_path): return ""
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    try:
        with open_func(file_path, mode) as f:
            content = f.read().strip()
            if content.startswith('>'):
                lines = content.split('\n')
                return "".join([line.strip() for line in lines if not line.startswith('>')])
            else:
                return "".join(content.split())
    except:
        return ""

def tenten_scroll_export_csv(lengths_dict: Dict[str, List[int]], output_path: str) -> None:
    if not lengths_dict: raise ValueError("Insufficient data pool for CSV export.")
    with open(output_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Dataset", "Sequence_Length_bp"])
        for label, lengths in lengths_dict.items():
            for length in lengths:
                writer.writerow([label, length])

def pain_chibaku_tensei_merge_filter(files: List[str], output_path: str, min_phred: int = 20) -> None:
    with open(output_path, 'w') as out_f:
        for f in files:
            open_func = gzip.open if f.endswith('.gz') else open
            mode = 'rt' if f.endswith('.gz') else 'r'
            fmt = neji_byakugan_format_sniffer(f)
            with open_func(f, mode) as in_f:
                if fmt == 'fastq':
                    for title, seq, qual in FastqGeneralIterator(in_f):
                        avg_qual = sum(ord(c) - 33 for c in qual) / len(qual)
                        if avg_qual >= min_phred:
                            out_f.write(f">{title}\n{seq}\n")
                else:
                    parse_fmt = 'genbank' if fmt == 'genbank' else 'fasta'
                    for record in SeqIO.parse(in_f, parse_fmt):
                        SeqIO.write(record, out_f, 'fasta')

def madara_infinite_tsukuyomi_html_report(stats_text: str, figures: List[Figure], output_path: str) -> None:
    img_tags = []
    for fig in figures:
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='#1a1a1a')
        img_str = base64.b64encode(buf.getvalue()).decode('utf-8')
        img_tags.append(f'<img src="data:image/png;base64,{img_str}" alt="Distribution Plot" />')
        
    html_content = f"""
    <html>
    <head>
        <title>Genomic Analysis Report</title>
        <style>
            body {{ font-family: 'Segoe UI', Arial, sans-serif; background-color: #121212; color: #e0e0e0; padding: 20px; }}
            h1 {{ color: #ffd900; border-bottom: 1px solid #333; padding-bottom: 10px; }}
            pre {{ background-color: #1a1a1a; padding: 15px; border-radius: 5px; border: 1px solid #333333; }}
            img {{ max-width: 100%; border: 1px solid #333333; border-radius: 5px; margin-top: 20px; page-break-after: always; }}
        </style>
    </head>
    <body>
        <h1>Comprehensive Sequence Analysis Report</h1>
        <h2>Statistical Overview</h2>
        <pre>{stats_text}</pre>
        <h2>Visualizations</h2>
        {"".join(img_tags)}
    </body>
    </html>
    """
    with open(output_path, 'w') as f:
        f.write(html_content)

def sai_ninja_art_dark_theme(fig, ax):
    fig.patch.set_facecolor('#1a1a1a')
    ax.set_facecolor('#121212')
    ax.tick_params(colors='#e0e0e0')
    ax.xaxis.label.set_color('#ffd900')
    ax.yaxis.label.set_color('#ffd900')
    ax.title.set_color('#ffd900')
    for spine in ax.spines.values():
        spine.set_edgecolor('#333333')

def sai_ninja_art_set_color(legend):
    if legend:
        legend.get_title().set_color('#e0e0e0')
        for text in legend.get_texts():
            text.set_color('#e0e0e0')
        legend.get_frame().set_facecolor('#121212')
        legend.get_frame().set_edgecolor('#333333')

class KatsuyuSlugFormatParser:
    def __init__(self): pass
    def tsunade_katsuyu_summarize_file(self, file_path: str, fmt: str) -> Dict[str, Any]:
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'
        result = {"text": "", "count": 0, "format": fmt}
        try:
            if fmt == 'vcf':
                variants = sum(1 for line in open_func(file_path, mode) if not line.startswith('#'))
                result["text"] = f"[FILE PROFILED] VCF File '{os.path.basename(file_path)}' contains {variants:,} variant records."
                result["count"] = variants
            elif fmt in ['gff', 'bed']:
                features = sum(1 for line in open_func(file_path, mode) if line.strip() and not line.startswith('#'))
                result["text"] = f"[FILE PROFILED] Genomic Annotation File '{os.path.basename(file_path)}' contains {features:,} features/regions."
                result["count"] = features
            elif fmt == 'sam':
                alignments = sum(1 for line in open_func(file_path, mode) if not line.startswith('@'))
                result["text"] = f"[FILE PROFILED] SAM Alignment File '{os.path.basename(file_path)}' contains {alignments:,} mapped records."
                result["count"] = alignments
            elif fmt == 'bam':
                result["text"] = f"[FILE INFO] BAM File '{os.path.basename(file_path)}' detected. Use Henge Format Shifter to convert to SAM."
        except Exception as e:
            result["text"] = f"[ERROR] Parser failure on {file_path}: {str(e)}"
        return result

class HengeFormatShifter:
    @staticmethod
    def transformation_jutsu_convert(input_file: str, output_file: str, in_fmt: str, out_fmt: str) -> str:
        try:
            if in_fmt == 'bam' and out_fmt == 'sam':
                cmd = f"samtools view -h {input_file} > {output_file}"
                subprocess.run(cmd, shell=True, check=True)
                return f"Henge Jutsu Successful: BAM to SAM via samtools."
            
            parsed_in = 'fastq' if in_fmt == 'fastq' else 'genbank' if in_fmt == 'gbk' else 'fasta'
            parsed_out = 'fasta'
            count = SeqIO.convert(input_file, parsed_in, output_file, parsed_out)
            return f"Henge Jutsu Successful: Converted {count} records from {in_fmt.upper()} to {out_fmt.upper()}."
        except Exception as e:
            return f"Transformation failed: {str(e)}"

class PainBanshoTeninNCBIFetcher:
    @staticmethod
    def universal_pull_sequence(accession: str) -> str:
        if not accession: return ""
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession.strip(), rettype="fasta", retmode="text")
            seq_data = handle.read()
            handle.close()
            return seq_data
        except Exception as e:
            raise ValueError(f"Pain's Universal Pull failed to fetch from NCBI: {str(e)}")

class NejiVariantEyeParser:
    @staticmethod
    def eight_trigrams_parse_vcf(filepath: str) -> List[Dict[str, str]]:
        variants = []
        open_f = gzip.open if filepath.endswith('.gz') else open
        mode = 'rt' if filepath.endswith('.gz') else 'r'
        try:
            with open_f(filepath, mode) as f:
                for line in f:
                    if line.startswith('#'): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        variants.append({
                            'CHROM': parts[0], 'POS': parts[1], 'ID': parts[2],
                            'REF': parts[3], 'ALT': parts[4], 
                            'QUAL': parts[5] if len(parts)>5 else '',
                            'INFO': parts[7] if len(parts)>7 else ''
                        })
            return variants
        except Exception as e:
            raise ValueError(f"Byakugan failed to decode VCF matrix: {str(e)}")

class UchihaSharinganAnalyzer:
    def __init__(self):
        self.raw_lengths: Dict[str, List[int]] = defaultdict(list)
        self.file_paths: Dict[str, str] = {}
        self.filtered_lengths: Dict[str, List[int]] = defaultdict(list)
        self.filtered_stats: Dict[str, Dict[str, float]] = defaultdict(dict)
        self.expected_genome_size = 0
        
    def sasuke_sharingan_read_fasta(self, file_path: str, label: Optional[str] = None) -> bool:
        if label is None: label = os.path.basename(file_path)
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'
        fmt = neji_byakugan_format_sniffer(file_path)
        parse_fmt = 'genbank' if fmt == 'genbank' else 'fasta'
            
        try:
            with open_func(file_path, mode) as handle:
                self.file_paths[label] = file_path
                for record in SeqIO.parse(handle, parse_fmt):
                    self.raw_lengths[label].append(len(record.seq))
                    
            if not self.raw_lengths[label]:
                raise ValueError(f"No valid sequence records found in {file_path}")
            return True
        except Exception as e:
            raise ValueError(f"File parsing error for {file_path}: {str(e)}")

    def itachi_susanoo_apply_filters(self, min_len: int, max_len: int) -> None:
        self.filtered_lengths.clear()
        self.filtered_stats.clear()
        for label, lengths in self.raw_lengths.items():
            filtered = [l for l in lengths if min_len <= l <= max_len]
            if filtered:
                self.filtered_lengths[label] = filtered
                self.filtered_stats[label] = shikamaru_shadow_stats(filtered, self.expected_genome_size)

    def _prepare_dataframe(self) -> pd.DataFrame:
        data = []
        for label, lengths in self.filtered_lengths.items():
            for l in lengths: data.append({'length': l, 'dataset': label})
        return pd.DataFrame(data)

    def sasuke_amaterasu_histogram(self, show_normal: bool = False) -> Figure:
        if not self.filtered_lengths: raise ValueError("Empty dataset for plotting operation.")
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        all_lengths = [l for lengths in self.filtered_lengths.values() for l in lengths]
        iqr = np.percentile(all_lengths, 75) - np.percentile(all_lengths, 25)
        bin_width = 2 * iqr / (len(all_lengths) ** (1/3)) if iqr > 0 else 10
        bins = max(10, int((max(all_lengths) - min(all_lengths)) / bin_width)) if bin_width > 0 else 10
        for label, lengths in self.filtered_lengths.items():
            sns.histplot(lengths, bins=int(bins), ax=ax, alpha=0.6, label=f"{label} (n={len(lengths)})", kde=True, edgecolor='none')
            if show_normal and len(lengths) > 1:
                mu, std = np.mean(lengths), np.std(lengths)
                xmin, xmax = ax.get_xlim()
                x = np.linspace(xmin, xmax, 100)
                p = stats.norm.pdf(x, mu, std)
                ax.plot(x, p * len(lengths) * bin_width, linewidth=2, linestyle='--', label=f'Theoretical Normal ({label})')
        ax.set_title("Sequence Length Distribution Profile")
        ax.set_xlabel("Sequence Length (bp)")
        ax.set_ylabel("Frequency")
        sai_ninja_art_set_color(ax.legend(title="Dataset Selection"))
        fig.tight_layout()
        return fig

    def itachi_tsukuyomi_boxplot(self) -> Figure:
        if not self.filtered_lengths: raise ValueError("Empty dataset for plotting operation.")
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        df = self._prepare_dataframe()
        sns.boxplot(x='dataset', y='length', data=df, ax=ax, showfliers=False, color='#ffd900')
        ax.set_title("Sequence Length Variability Comparison")
        ax.set_xlabel("Dataset Cohort")
        ax.set_ylabel("Sequence Length (bp)")
        ax.tick_params(axis='x', rotation=15)
        fig.tight_layout()
        return fig

    def danzo_izanagi_violinplot(self) -> Figure:
        if not self.filtered_lengths: raise ValueError("Empty dataset for plotting operation.")
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        df = self._prepare_dataframe()
        sns.violinplot(x='dataset', y='length', data=df, ax=ax, color='#3498db', inner="quartile")
        ax.set_title("Sequence Distribution Kernel Density (Violin Plot)")
        ax.set_xlabel("Dataset Cohort")
        ax.set_ylabel("Sequence Length (bp)")
        ax.tick_params(axis='x', rotation=15)
        fig.tight_layout()
        return fig

    def madara_susanoo_kde(self) -> Figure:
        if not self.filtered_lengths: raise ValueError("Empty dataset for plotting operation.")
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        df = self._prepare_dataframe()
        sns.kdeplot(data=df, x='length', hue='dataset', fill=True, alpha=0.4, ax=ax, common_norm=False)
        ax.set_title("Probability Density Function (KDE)")
        ax.set_xlabel("Sequence Length (bp)")
        ax.set_ylabel("Density Amplitude")
        sai_ninja_art_set_color(ax.get_legend())
        fig.tight_layout()
        return fig

    def itachi_izanami_cumulative(self) -> Figure:
        if not self.filtered_lengths: raise ValueError("Empty dataset for plotting operation.")
        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        df = self._prepare_dataframe()
        sns.ecdfplot(data=df, x='length', hue='dataset', ax=ax, linewidth=2)
        ax.set_title("Empirical Cumulative Distribution Function (ECDF)")
        ax.set_xlabel("Sequence Length (bp)")
        ax.set_ylabel("Cumulative Proportion")
        sai_ninja_art_set_color(ax.get_legend())
        fig.tight_layout()
        return fig


class SenjuMokutonAnalyzer:
    def __init__(self):
        self.sequence = ""
        self.base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}
        self.base_colors = SCROLL_OF_COLORS.copy()
    
    def yamato_mokuton_set_sequence(self, sequence: str) -> None:
        self.sequence = sequence.upper()
        self.naruto_kage_bunshin_count_bases()
    
    def naruto_kage_bunshin_count_bases(self) -> None:
        self.base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}
        for base in self.sequence:
            if base in self.base_counts: self.base_counts[base] += 1
            else: self.base_counts['N'] += 1
                
    def tsunade_creation_rebirth_translate(self) -> str:
        if not self.sequence: return ""
        try: return str(Seq(self.sequence.replace('-', '')).translate(to_stop=False))
        except Exception as e: return f"Translation algorithm failure: {str(e)}"

    def shikamaru_shadow_reversal_complement(self) -> str:
        if not self.sequence: return ""
        return str(Seq(self.sequence.replace('-', '')).reverse_complement())

    def jiraiya_transcription_jutsu(self) -> str:
        if not self.sequence: return ""
        return str(Seq(self.sequence.replace('-', '')).transcribe())
            
    def kakashi_doton_get_statistics_text(self) -> str:
        if not self.sequence: return "Sequence buffer empty."
        total = sum(self.base_counts.values())
        if total == 0: return "Sequence data invalid."
        text = ["Nucleotide Composition Analysis:", "-" * 35, f"Total Sequence Length: {len(self.sequence)} bp", f"Validated Bases Analyzed: {total}", "", "Distribution Metrics:"]
        for base, count in sorted(self.base_counts.items()):
            if count > 0: text.append(f"{base} Nucleotides: {count} ({count/total*100:.2f}%)")
        return "\n".join(text)

    def shino_kikaichu_kmer_swarm(self, k: int) -> Tuple[Figure, str]:
        if not self.sequence or len(self.sequence) < k: raise ValueError("Sequence length insufficient for K-mer.")
        def kmer_generator(seq, k):
            for i in range(len(seq) - k + 1): yield seq[i:i+k]
        kmer_counts = Counter(kmer_generator(self.sequence, k))
        top_kmers = kmer_counts.most_common(15)
        
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        labels, counts = zip(*top_kmers)
        ax.bar(labels, counts, color='#3498db')
        ax.set_title(f"Genomic Signature: Top 15 K-mers (k={k})")
        ax.set_ylabel("Absolute Frequency")
        ax.tick_params(axis='x', rotation=45)
        fig.tight_layout()
        
        stats_txt = f"Total Distinct {k}-mer Profiles: {len(kmer_counts)}\n"
        for label, count in top_kmers: stats_txt += f"Motif [{label}]: {count} occurrences\n"
        return fig, stats_txt

    def orochimaru_curse_mark_entropy_plot(self, window_size: int = 50) -> Figure:
        if len(self.sequence) < window_size: raise ValueError("Sequence length insufficient.")
        entropies = []
        for i in range(len(self.sequence) - window_size + 1):
            window = self.sequence[i:i+window_size]
            counts = Counter(window)
            ent = 0.0
            for count in counts.values():
                p = count / window_size
                ent -= p * math.log2(p)
            entropies.append(ent)
            
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        ax.plot(entropies, color='#9b59b6')
        ax.set_title(f"Shannon Entropy Profiling (Window: {window_size}bp)")
        ax.set_xlabel("Genomic Position (bp)")
        ax.set_ylabel("Entropy Coefficient (bits)")
        fig.tight_layout()
        return fig

    def gaara_dust_complexity_shield(self, window_size: int = 64) -> str:
        if len(self.sequence) < window_size: return self.sequence
        masked_seq = list(self.sequence)
        for i in range(len(self.sequence) - window_size + 1):
            window = self.sequence[i:i+window_size]
            counts = Counter(window)
            ent = sum(- (c/window_size) * math.log2(c/window_size) for c in counts.values())
            if ent < 1.0: # Arbitrary threshold for low complexity
                for j in range(window_size):
                    masked_seq[i+j] = 'N'
        return "".join(masked_seq)

    def neji_trigram_rho_odds(self) -> str:
        if len(self.sequence) < 2: return "Insufficient data."
        nucs = ['A', 'T', 'G', 'C']
        freqs = {n: self.sequence.count(n)/len(self.sequence) for n in nucs}
        di_counts = Counter(self.sequence[i:i+2] for i in range(len(self.sequence)-1))
        total_di = sum(di_counts.values())
        
        res = ["Dinucleotide Rho Odds Ratio (Neji Trigram Analysis):", "-"*50]
        for n1 in nucs:
            for n2 in nucs:
                di = n1 + n2
                obs = di_counts.get(di, 0) / total_di if total_di > 0 else 0
                exp = freqs[n1] * freqs[n2]
                rho = obs / exp if exp > 0 else 0
                if rho > 1.2: flag = "HIGH (Over-represented)"
                elif rho < 0.8: flag = "LOW (Under-represented)"
                else: flag = "NORMAL"
                res.append(f"Pair [{di}]: Rho = {rho:.3f} -> {flag}")
        return "\n".join(res)

    def hidan_jashin_cpg_islands_vectorized(self, window_size: int = 200) -> Tuple[Figure, str]:
        if len(self.sequence) < window_size: raise ValueError("Sequence length insufficient.")
        seq_arr = np.array(list(self.sequence))
        c_mask = (seq_arr == 'C').astype(int)
        g_mask = (seq_arr == 'G').astype(int)
        cg_mask = np.array([1 if self.sequence[i:i+2] == 'CG' else 0 for i in range(len(self.sequence))], dtype=int)
        
        ones = np.ones(window_size, dtype=int)
        c_counts = np.convolve(c_mask, ones, mode='valid')
        g_counts = np.convolve(g_mask, ones, mode='valid')
        cg_counts = np.convolve(cg_mask, ones, mode='valid')
        
        denom = c_counts * g_counts
        with np.errstate(divide='ignore', invalid='ignore'):
            obs_exp_ratios = np.where(denom > 0, (cg_counts * window_size) / denom, 0)
            
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        ax.plot(obs_exp_ratios, color='#2ecc71')
        ax.axhline(y=0.6, color='r', linestyle='--', label='Statutory Threshold (0.6)')
        ax.set_title(f"CpG Island Identification Profile Vectorized (Window={window_size})")
        ax.set_ylabel("Observed/Expected Ratio")
        sai_ninja_art_set_color(ax.legend())
        fig.tight_layout()
        
        high_cpg = np.sum(obs_exp_ratios >= 0.6)
        return fig, f"Identified High-Density CpG Zones (Ratio >= 0.6): {high_cpg} distinct windows."

    def zabuza_executioner_restriction_map(self) -> Tuple[Figure, str]:
        if not self.sequence: raise ValueError("Sequence empty for restriction digest.")
        seq_obj = Seq(self.sequence.upper())
        rb = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI, NotI])
        ana = Analysis(rb, seq_obj, linear=True)
        results = ana.full()
        
        stats = ["Restriction Digest Cleavage Sites:", "-"*45]
        fragments = []
        
        for enz, sites in results.items():
            stats.append(f"Enzyme {enz}: {len(sites)} cuts -> Positions: {sites}")
            if sites:
                cuts = [0] + sites + [len(self.sequence)]
                for i in range(len(cuts)-1):
                    size = cuts[i+1] - cuts[i]
                    fragments.append({'Enzyme': str(enz), 'Size': size})
            else:
                fragments.append({'Enzyme': str(enz), 'Size': len(self.sequence)})

        fig = Figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        
        if fragments:
            df = pd.DataFrame(fragments)
            sns.stripplot(x='Enzyme', y='Size', data=df, ax=ax, color='#2ecc71', size=8, jitter=True)
            ax.set_yscale('log')
            ax.set_title("Virtual Agarose Gel Simulation (Log Scale)")
            ax.set_ylabel("Fragment Size (bp)")
            ax.set_xlabel("Restriction Enzyme")
        else:
            ax.text(0.5, 0.5, "No cleavage sites matched.", ha='center', va='center', color='white')
            
        fig.tight_layout()
        return fig, "\n".join(stats)

    def inuzuka_mind_dog_taxonomy(self) -> str:
        if not self.sequence: return "No sequence provided."
        gc_content = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100
        # Mini fake-database for demonstration (Kraken-lite heuristic)
        db = {
            'Escherichia coli': 50.8,
            'Staphylococcus aureus': 32.7,
            'Mycobacterium tuberculosis': 65.6,
            'Pseudomonas aeruginosa': 66.6,
            'Plasmodium falciparum': 19.3
        }
        closest = min(db.keys(), key=lambda k: abs(db[k] - gc_content))
        return f"Inuzuka Kraken-Lite Match: Closest Taxonomic Profile is {closest} (Seq GC: {gc_content:.1f}%, DB Ref: {db[closest]:.1f}%)"

    def orochimaru_cursed_islands_pai(self, window: int=5000) -> str:
        if len(self.sequence) < window: return "Contig too small to detect Pathogenicity Islands."
        gc_content = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence)
        
        islands = []
        for i in range(0, len(self.sequence) - window, window):
            sub = self.sequence[i:i+window]
            sub_gc = (sub.count('G') + sub.count('C')) / window
            # If GC drops significantly below global average
            if sub_gc < gc_content - 0.05: 
                islands.append(f"PAI Candidate at {i}-{i+window} (Local GC: {sub_gc*100:.1f}%, Global: {gc_content*100:.1f}%)")
                
        if not islands: return "Orochimaru found no anomalies. Sequence appears native."
        return "Orochimaru Cursed Pathogenicity Islands (PAI):\n" + "-"*50 + "\n" + "\n".join(islands)

    def kimimaro_bone_telomere_mapper(self) -> str:
        if len(self.sequence) < 100: return "Sequence too small."
        # Look for TTAGGG repeats at the edges
        motif = r'(TTAGGG|CCCTAA){3,}'
        start_edge = self.sequence[:5000]
        end_edge = self.sequence[-5000:]
        
        res = ["Kimimaro Telomere Edge Mapper:", "-"*40]
        m1 = list(re.finditer(motif, start_edge))
        m2 = list(re.finditer(motif, end_edge))
        
        if m1: res.append(f"5' Telomeric Repeats found: {len(m1)} loci.")
        if m2: res.append(f"3' Telomeric Repeats found: {len(m2)} loci.")
        if not m1 and not m2: res.append("No Telomeric signatures found at contig edges.")
        return "\n".join(res)


class UzumakiChakraFastqAnalyzer:
    def __init__(self):
        self.read_counts = {}
        self.genome_size = 0
        self.results = {}
        self.quality_profiles = {}
        self.first_fragments = defaultdict(list)
    
    def naruto_rasengan_read_fastq_blocks(self, file_path: str, label: Optional[str] = None) -> bool:
        if label is None: label = os.path.basename(file_path)
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'
        try:
            count = 0
            qual_sums = defaultdict(int)
            qual_counts = defaultdict(int)
            with open_func(file_path, mode) as handle:
                for title, seq, qual in FastqGeneralIterator(handle):
                    count += 1
                    for i, char in enumerate(qual):
                        qual_sums[i] += (ord(char) - 33)
                        qual_counts[i] += 1
                    if len(seq) >= 50:
                        self.first_fragments[label].append(seq[:50])
                        
            if count == 0: raise ValueError(f"No sequence reads extracted from {file_path}")
            self.read_counts[label] = count
            max_len = max(qual_counts.keys()) + 1
            avg_quals = [qual_sums[i] / qual_counts[i] if qual_counts[i] > 0 else 0 for i in range(max_len)]
            self.quality_profiles[label] = avg_quals
            return True
        except Exception as e:
            raise ValueError(f"FASTQ extraction failure on {file_path}: {str(e)}")
    
    def naruto_oodama_rasengan_calculate_genomes(self, genome_size: int, read_length: int = 150, paired_end: bool = True) -> None:
        if genome_size <= 0: raise ValueError("Genome dimension parameter must be strictly positive.")
        self.genome_size = genome_size
        self.results = {}
        for label, count in self.read_counts.items():
            multiplier = 2 if paired_end else 1
            coverage = (count * read_length * multiplier) / genome_size
            self.results[label] = {
                'reads': count, 'genomes': coverage, 'read_length': read_length,
                'paired_end': paired_end, 'genome_size': genome_size
            }
    
    def uzumaki_fuuinjutsu_get_results_text(self) -> str:
        if not self.results: return "Insufficient data for coverage extrapolation."
        text = ["Sequencing Depth and Coverage Estimation:", "-" * 45, f"Reference Assembly Size: {self.genome_size:,} bp", ""]
        for label, result in self.results.items():
            text.append(f"Metrics for Protocol {label}:\n  Absolute Read Count: {result['reads']:,}\n  Estimated Coverage: {result['genomes']:.2f}X\n")
        return "\n".join(text)
        
    def itachi_tsukuyomi_phred_plot(self) -> Figure:
        if not self.quality_profiles: raise ValueError("Phred telemetry data unavailable.")
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        for label, avg_quals in self.quality_profiles.items():
            ax.plot(range(1, len(avg_quals)+1), avg_quals, label=label)
        ax.axhline(y=20, color='r', linestyle='--', alpha=0.5, label='Q20 Baseline')
        ax.axhline(y=30, color='g', linestyle='--', alpha=0.5, label='Q30 Optimal')
        ax.set_title("Base Calling Quality Assessment (Phred Scale)")
        ax.set_xlabel("Read Nucleotide Cycle")
        ax.set_ylabel("Mean Phred Q-Score")
        ax.set_ylim(0, 45)
        sai_ninja_art_set_color(ax.legend())
        fig.tight_layout()
        return fig

    def pain_shinra_tensei_rarefaction_curve(self, step: int = 1000) -> Figure:
        if not self.read_counts: raise ValueError("Dataset read count missing for rarefaction mechanics.")
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        for label, count in self.read_counts.items():
            subsamples = list(range(1, count + 1, max(1, count // 20)))
            if subsamples[-1] != count: subsamples.append(count)
            sim_richness = [100 * math.log(n + 1) for n in subsamples] 
            ax.plot(subsamples, sim_richness, marker='o', label=label)
        ax.set_title("Library Complexity Simulation (Rarefaction Analysis)")
        ax.set_xlabel("Random Sequencing Effort (Reads)")
        ax.set_ylabel("Extrapolated Variant Diversity")
        sai_ninja_art_set_color(ax.legend())
        fig.tight_layout()
        return fig

    def kisame_water_shark_saturation(self) -> Figure:
        if not self.results: raise ValueError("Calculate Genome Coverage first.")
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        for label, res in self.results.items():
            covs = [res['genomes'] * (p/100) for p in range(10, 110, 10)]
            reads = [res['reads'] * (p/100) for p in range(10, 110, 10)]
            ax.plot(reads, covs, marker='s', label=label)
        ax.axhline(30.0, color='yellow', linestyle='--', label='Ideal 30X Coverage')
        ax.set_title("Genomic Coverage Saturation Prediction (Kisame Water Shark)")
        ax.set_xlabel("Sequencing Reads Processed")
        ax.set_ylabel("Fold Coverage (X)")
        sai_ninja_art_set_color(ax.legend())
        fig.tight_layout()
        return fig

    def zetsu_spore_duplication_rate(self) -> str:
        if not self.first_fragments: return "No sequences available to measure duplication."
        res = ["Zetsu Spore Estimator (PCR Duplication Profile):", "-"*50]
        for label, frags in self.first_fragments.items():
            total = len(frags)
            unique = len(set(frags))
            dup_rate = (1 - (unique / total)) * 100 if total > 0 else 0
            res.append(f"Dataset [{label}]: {dup_rate:.2f}% identical reads (PCR Artifact Estimation)")
        return "\n".join(res)


class GaaraSandTrimmer:
    @staticmethod
    def sabaku_kyuu_trim_fastq(input_fastq: str, output_fastq: str, min_phred: int = 20) -> str:
        open_in = gzip.open if input_fastq.endswith('.gz') else open
        mode_in = 'rt' if input_fastq.endswith('.gz') else 'r'
        kept = 0
        discarded = 0
        try:
            with open_in(input_fastq, mode_in) as in_handle, open(output_fastq, 'w') as out_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    start, end = 0, len(qual)
                    while start < end and (ord(qual[start]) - 33) < min_phred: start += 1
                    while end > start and (ord(qual[end-1]) - 33) < min_phred: end -= 1
                    
                    if end - start > 30: 
                        out_handle.write(f"@{title}\n{seq[start:end]}\n+\n{qual[start:end]}\n")
                        kept += 1
                    else:
                        discarded += 1
            return f"Sand Trimmer Results: Kept {kept} reads, discarded {discarded} reads."
        except Exception as e:
            return f"Trimming failure: {str(e)}"


class RinneganAdvancedBioTools:
    @staticmethod
    def kakashi_kamui_gc_skew_rolling(sequence: str, window: int = 1000) -> Figure:
        seq = sequence.upper()
        if len(seq) < window: raise ValueError("Contig length below statistical sliding window requirement.")
        seq_arr = np.array(list(seq))
        g_mask = (seq_arr == 'G').astype(float)
        c_mask = (seq_arr == 'C').astype(float)
        g_roll = pd.Series(g_mask).rolling(window=window).sum().values
        c_roll = pd.Series(c_mask).rolling(window=window).sum().values
        skew = (g_roll - c_roll) / (g_roll + c_roll + 1e-9)
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        ax.plot(range(len(skew)), skew, color='#f39c12')
        ax.axhline(0, color='red', linestyle='--')
        ax.set_title(f"Genomic GC Skew Rolling (Window={window}bp)")
        ax.set_xlabel("Chromosomal Coordinate (bp)")
        ax.set_ylabel("GC Skew Index: (G-C)/(G+C)")
        fig.tight_layout()
        return fig

    @staticmethod
    def jiraiya_sage_mode_orf_finder(sequence: str, min_len: int = 100) -> str:
        seq = Seq(sequence.upper())
        orfs = []
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                frame_seq = nuc[frame:]
                for i in range(0, len(frame_seq) - 2, 3):
                    if frame_seq[i:i+3] == "ATG":
                        prot = str(frame_seq[i:].translate(to_stop=True))
                        if len(prot) >= min_len:
                            orfs.append((len(prot), strand, frame+1, prot))
        orfs = list(set(orfs))
        orfs.sort(key=lambda x: x[0], reverse=True)
        if not orfs: return "System failed to identify valid Open Reading Frames meeting threshold constraints."
        res = [f"Computational ORF Prediction (Top Candidates) - Min Threshold: {min_len} AA", "-"*65]
        for idx, (length, strand, frame, prot) in enumerate(orfs[:5]):
            str_dir = "Sense (+)" if strand == 1 else "Antisense (-)"
            res.append(f"Identifier: ORF_{idx+1} | Extent: {length} AA | Frame Phase: {str_dir} {frame}")
            res.append(f"Translation Output: {prot[:50]}..." if length > 50 else f"Translation Output: {prot}\n")
        return "\n".join(res)

    @staticmethod
    def minato_rasengan_primer_wizard(primer: str) -> str:
        p = primer.upper()
        if not all(c in 'ATGC' for c in p): return "Format violation: Oligonucleotide must strictly conform to IUPAC ATGC characters."
        gc_perc = (p.count('G') + p.count('C')) / len(p) * 100
        try:
            tm = mt.Tm_NN(p)
            dH, dS = mt.calc_Tm_NN(p, return_dH_dS=True)
            dG = dH - (310.15 * dS) / 1000 
        except:
            tm = 4*(p.count('G') + p.count('C')) + 2*(p.count('A') + p.count('T'))
            dG = 0.0
        rc = str(Seq(p).reverse_complement())
        dimer_risk = "Elevated" if rc[:4] in p or p[:4] in rc else "Negligible"
        hairpin_risk = "Negligible"
        for i in range(len(p)-4):
            if str(Seq(p[i:i+4]).reverse_complement()) in p[i+4:]:
                hairpin_risk = "Elevated"
                break
        res = [
            f"Oligonucleotide Thermodynamic Analysis Profile", "-" * 45,
            f"Query Sequence: {p}", f"Molecular Length: {len(p)} bp",
            f"G/C Proportion: {gc_perc:.2f}% (Optimum Parameter: 40-60%)",
            f"Calculated Melting Temperature (Tm): {tm:.2f} °C (Optimum Parameter: 50-65°C)",
            f"Gibbs Free Energy (\u0394G): {dG:.2f} kcal/mol",
            f"Structural Self-Dimer Probability: {dimer_risk}",
            f"Secondary Structure (Hairpin) Probability: {hairpin_risk}"
        ]
        return "\n".join(res)

    @staticmethod
    def kakashi_crispr_copy_sgrna(seq_str: str) -> str:
        seq = seq_str.upper()
        pattern = re.compile(r'(?=([ATGC]{20}[ATGC]GG))')
        matches = []
        for m in pattern.finditer(seq):
            target = m.group(1)
            core_20 = target[:20]
            gc = (core_20.count('G') + core_20.count('C')) / 20.0 * 100
            occurrences = seq.count(core_20)
            off_target_risk = "HIGH" if occurrences > 1 else "LOW"
            matches.append(f"Position: {m.start():<6} | Target: {core_20} | PAM: {target[20:]} | GC%: {gc:.1f}% | Off-Target Risk: {off_target_risk}")
        if not matches: return "No SpCas9 CRISPR (NGG) targets found in the given sequence."
        res = ["CRISPR SpCas9 sgRNA Identification Scan (Sharingan Copied):", "-"*80]
        res.extend(matches[:50])
        if len(matches) > 50: res.append(f"...Data block truncated (Total targets: {len(matches)}).")
        return "\n".join(res)

    @staticmethod
    def chouji_calorie_amino_plot(protein_seq: str) -> Figure:
        seq = ''.join([c for c in protein_seq.upper() if c.isalpha() and c not in 'BXZJ'])
        if not seq: raise ValueError("Empty protein sequence.")
        counts = Counter(seq)
        labels, values = zip(*counts.most_common())
        
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        ax.bar(labels, values, color='#e74c3c')
        ax.set_title("Chouji Calorie Plot: Amino Acid Frequency")
        ax.set_ylabel("Absolute Counts")
        fig.tight_layout()
        return fig

    @staticmethod
    def raikage_lightning_armor_protein_properties(protein_seq: str) -> str:
        seq = ''.join([c for c in protein_seq.upper() if c.isalpha() and c not in 'BXZJ'])
        if not seq: return "Invalid biological sequence input."
        try:
            analysed = ProteinAnalysis(seq)
            mw = analysed.molecular_weight()
            pi = analysed.isoelectric_point()
            instab = analysed.instability_index()
            gravy = analysed.gravy()
            extinction = analysed.molar_extinction_coefficient()
            ss = analysed.secondary_structure_fraction() # Helix, Turn, Sheet
            
            aliphatic_idx = 100 * (seq.count('A') + 2.9*seq.count('V') + 3.9*(seq.count('I') + seq.count('L'))) / len(seq)
            
            res = [
                f"Macromolecular Physicochemical Report (Gai Eight Gates):", "-" * 55,
                f"Amino Acid Residue Count: {len(seq)}",
                f"Molecular Mass: {mw:.2f} Daltons",
                f"Theoretical Isoelectric Point (pI): {pi:.2f}",
                f"Structural Instability Index: {instab:.2f} ({'Stable' if instab < 40 else 'Unstable'})",
                f"GRAVY Score: {gravy:.3f} (Positive = Hydrophobic)",
                f"Aliphatic Index: {aliphatic_idx:.2f}",
                f"Molar Extinction Coefficient (Darui Black Panther):",
                f"   - Reduced: {extinction[0]} | Oxidized: {extinction[1]}",
                f"2D Sec Struct (Rock Lee Helix Kick):",
                f"   - Alpha-Helix: {ss[0]*100:.1f}% | Beta-Sheet: {ss[2]*100:.1f}% | Turn: {ss[1]*100:.1f}%"
            ]
            return "\n".join(res)
        except Exception as e: return f"Peptide computation fault: {str(e)}"

    @staticmethod
    def kisame_water_prison_hydrophobicity_plot(protein_seq: str, window: int = 9) -> Figure:
        seq = ''.join([c for c in protein_seq.upper() if c.isalpha()])
        if len(seq) < window: raise ValueError("Peptide chain insufficient for designated structural window.")
        kd_scale = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 
                    'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 
                    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
        scores = []
        for i in range(len(seq) - window + 1):
            chunk = seq[i:i+window]
            score = sum(kd_scale.get(aa, 0) for aa in chunk) / window
            scores.append(score)
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        ax.plot(scores, color='#3498db')
        ax.axhline(0, color='red', linestyle='--')
        ax.axhline(1.6, color='yellow', linestyle=':', label='Putative Transmembrane Boundary (1.6)')
        ax.set_title(f"Kyte-Doolittle Hydropathy Mapping (Window={window})")
        ax.set_xlabel("Linear Residue Coordinate")
        ax.set_ylabel("Computed Hydropathy Index")
        sai_ninja_art_set_color(ax.legend())
        fig.tight_layout()
        return fig

    @staticmethod
    def sasuke_sharingan_snp_deducer_global(ref: str, query: str) -> str:
        ref, query = ref.upper(), query.upper()
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = -2
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        try:
            best_alignment = aligner.align(ref, query)[0]
            fmt_aln = format(best_alignment).split('\n')
            
             
            transitions = 0
            transversions = 0
            mut_effects = []
            
            for i in range(min(len(ref), len(query))):
                r_nt, q_nt = ref[i], query[i]
                if r_nt != q_nt and r_nt in 'ATGC' and q_nt in 'ATGC':
                    # Ti/Tv
                    if {r_nt, q_nt} == {'A', 'G'} or {r_nt, q_nt} == {'C', 'T'}: transitions += 1
                    else: transversions += 1
                    
                    
                    frame_pos = i % 3
                    codon_start = i - frame_pos
                    if codon_start + 3 <= len(ref) and codon_start + 3 <= len(query):
                        r_codon = Seq(ref[codon_start:codon_start+3])
                        q_codon = Seq(query[codon_start:codon_start+3])
                        if len(r_codon) == 3 and len(q_codon) == 3:
                            r_aa = r_codon.translate()
                            q_aa = q_codon.translate()
                            effect = "Synonymous" if r_aa == q_aa else f"Non-Synonymous ({r_aa}->{q_aa})"
                            mut_effects.append(f"Pos {i}: {r_nt}->{q_nt} | {effect}")

            res = [f"Needleman-Wunsch Global Alignment Mapping:", "-"*45]
            res.extend(fmt_aln[:10])
            res.append(f"...Data block truncated. Alignment Score: {best_alignment.score}\n")
            
            titv = transitions / transversions if transversions > 0 else float('inf')
            res.append(f"Shikamaru Ti/Tv Shadow Ratio:")
            res.append(f"  Transitions (A<->G, C<->T): {transitions}")
            res.append(f"  Transversions: {transversions}")
            res.append(f"  Ratio (Ti/Tv): {titv:.2f}")
            
            if mut_effects:
                res.append("\nKabuto Snake Mutation Predictor:")
                res.extend(mut_effects[:10])
                if len(mut_effects) > 10: res.append("...")
                
            return "\n".join(res)
        except Exception as e: return f"Alignment matrix failed to convergence: {str(e)}"

    @staticmethod
    def madara_rinnegan_dotplot_matrix_numpy(seq1: str, seq2: str, window: int = 10, threshold: int = 8) -> Figure:
        s1, s2 = seq1.upper(), seq2.upper()
        if len(s1) > 4000 or len(s2) > 4000: raise ValueError("Computation restriction: Matrices capped at 4000x4000 bp.")
        seq1_arr = np.array(list(s1))
        seq2_arr = np.array(list(s2))
        matrix = (seq1_arr[:, None] == seq2_arr).astype(int)
        eye_mat = np.eye(window)
        smoothed = convolve2d(matrix, eye_mat, mode='valid')
        x_pts, y_pts = np.where(smoothed >= threshold)
        fig = Figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        ax.scatter(x_pts, y_pts, s=1, color='#ffd900')
        ax.set_title(f"Synteny Alignment Matrix (Window={window}, Match={threshold})")
        ax.set_xlabel("Sequence Reference Vector (X)")
        ax.set_ylabel("Sequence Query Vector (Y)")
        fig.tight_layout()
        return fig

    @staticmethod
    def choji_expansion_jutsu_tandem_repeats(seq: str, min_repeat_len: int = 2, min_copies: int = 3) -> str:
        s = seq.upper()
        results = []
        for k in range(min_repeat_len, 10):
            pattern = re.compile(rf"([ATGC]{{{k}}})\1{{{min_copies-1},}}")
            for match in pattern.finditer(s):
                motif = match.group(1)
                full_match = match.group(0)
                copies = len(full_match) // len(motif)
                results.append(f"Motif Element: {motif} | Duplications: {copies}x | Locus: {match.start()}-{match.end()}")
        if not results: return "Algorithm confirms absence of tandem repeat signatures."
        return "Microsatellite / Tandem Repeat Identification:\n" + "-"*45 + "\n" + "\n".join(results[:50])

    @staticmethod
    def hashirama_phylo_tree_builder(alignment_path: str) -> Figure:
        try:
            aln = AlignIO.read(alignment_path, "fasta")
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(aln)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            fig = Figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            sai_ninja_art_dark_theme(fig, ax)
            Phylo.draw(tree, axes=ax, do_show=False, branch_labels=None, label_colors={'*': '#ffd900'})
            ax.set_title("Hashirama Mokuton: UPGMA Phylogenetic Tree")
            fig.tight_layout()
            return fig
        except Exception as e:
            raise ValueError(f"Mokuton Tree generation failed. Ensure input is a valid FASTA Multiple Sequence Alignment. Error: {str(e)}")


class YamatoWoodFeatureExtractor:
    @staticmethod
    def mokuton_extract_feature(gff_path: str, fasta_path: str, feature_name: str) -> str:
        try:
            target_chrom, target_start, target_end, target_strand = None, None, None, None
            with open(gff_path, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 9 and feature_name.lower() in parts[8].lower():
                        target_chrom = parts[0]
                        target_start = int(parts[3])
                        target_end = int(parts[4])
                        target_strand = parts[6]
                        break
            
            if not target_chrom: return f"Feature '{feature_name}' not found in the provided annotation file."
            
            extracted_seq = ""
            for record in SeqIO.parse(fasta_path, "fasta"):
                if record.id == target_chrom:
                    sliced = record.seq[target_start-1 : target_end]
                    extracted_seq = str(sliced.reverse_complement()) if target_strand == '-' else str(sliced)
                    break
                    
            if not extracted_seq: return f"Chromosome '{target_chrom}' not found in the reference FASTA."
            
            res = [f"Yamato Wood Extractor: '{feature_name}'", "-"*50]
            res.append(f"Coordinates: {target_chrom}:{target_start}-{target_end} ({target_strand})")
            res.append(f"Extracted Length: {len(extracted_seq)} bp\n")
            res.append(extracted_seq)
            return "\n".join(res)
            
        except Exception as e:
            return f"Yamato Extraction failed: {str(e)}"

    @staticmethod
    def hashirama_wood_density_plot(gff_path: str) -> Figure:
        starts = []
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'gene':
                    starts.append(int(parts[3]))
        
        if not starts: raise ValueError("No 'gene' features found in GFF.")
        
        fig = Figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        sai_ninja_art_dark_theme(fig, ax)
        sns.histplot(starts, binwidth=10000, ax=ax, color='#2ecc71', kde=True)
        ax.set_title("Hashirama Gene Density Plot (10kb Bins)")
        ax.set_xlabel("Genomic Coordinate (bp)")
        ax.set_ylabel("Gene Frequency")
        fig.tight_layout()
        return fig


class AnbuBlackOpsExternalTools:
    @staticmethod
    def itachi_mangekyou_run_command(command: str) -> str:
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=60)
            output = result.stdout if result.stdout else ""
            if result.stderr: output += "\n[SYSTEM ERROR STREAM]\n" + result.stderr
            return output if output else "Process executed silently. Return code 0."
        except Exception as e: return f"Kernel execution failure: {str(e)}"

    @staticmethod
    def deidara_explosive_blast_art(blast_exe_path: str, query_seq: str, subject_seq: str) -> str:
        if not os.path.exists(blast_exe_path) and not blast_exe_path.lower() in ["blastn", "blastn.exe"]:
            return "Fatal Error: Ensure the path to the BLAST+ executable is correct."
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as q_file, \
             tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as s_file, \
             tempfile.NamedTemporaryFile(mode='r', delete=False, suffix='.xml') as out_file:
             q_file.write(f">query\n{query_seq}\n")
             s_file.write(f">subject\n{subject_seq}\n")
             q_name, s_name, out_name = q_file.name, s_file.name, out_file.name

        cmd = [blast_exe_path, '-query', q_name, '-subject', s_name, '-outfmt', '5', '-out', out_name]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            res_str = ["Deidara's Explosive BLAST+ Alignment Art:", "-"*50]
            with open(out_name) as f:
                records = NCBIXML.parse(f)
                hits_found = False
                for rec in records:
                    for aln in rec.alignments:
                        for hsp in aln.hsps:
                            hits_found = True
                            res_str.append(f"Alignment Score: {hsp.score} | E-value: {hsp.expect}")
                            res_str.append(f"Query: {hsp.query[:60]}...")
                            res_str.append(f"Match: {hsp.match[:60]}...")
                            res_str.append(f"Subj:  {hsp.sbjct[:60]}...\n")
            return "\n".join(res_str) if hits_found else "Deidara's Blast missed: No hits found."
        except Exception as e:
            return f"BLAST+ execution imploded. Make sure BLAST+ is properly installed. Trace: {str(e)}"
        finally:
            if os.path.exists(q_name): os.remove(q_name)
            if os.path.exists(s_name): os.remove(s_name)
            if os.path.exists(out_name): os.remove(out_name)
