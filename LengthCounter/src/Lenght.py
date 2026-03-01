import sys
import argparse
import os
from Bio.Seq import Seq
from Bio import SeqIO
from logic import (UchihaSharinganAnalyzer, SenjuMokutonAnalyzer, UzumakiChakraFastqAnalyzer, 
                   RinneganAdvancedBioTools, KatsuyuSlugFormatParser, neji_byakugan_format_sniffer, 
                   ino_mind_transfer_stats_text, tenten_scroll_export_csv, gamabunta_summon_sequence_file,
                   PainBanshoTeninNCBIFetcher, NejiVariantEyeParser, AnbuBlackOpsExternalTools, 
                   YamatoWoodFeatureExtractor, GaaraSandTrimmer, HengeFormatShifter)

def tenten_kuchiyose_arg_parser():
    parser = argparse.ArgumentParser(description="LenghtCount Ninja Edition - Six Paths Suite")
    
    
    parser.add_argument("files", nargs="*", help="Sequence or annotation file(s) to analyze")
    parser.add_argument("-l", "--label", action="append", help="Custom dataset labels")
    parser.add_argument("--out-prefix", default="lc_analysis", help="Prefix for outputs")
    parser.add_argument("--gui", action="store_true", help="Launch the Graphical User Interface (GUI)")
    
    
    parser.add_argument("--histogram", action="store_true", help="Generate histogram")
    parser.add_argument("--boxplot", action="store_true", help="Generate a sequence length comparison boxplot")
    parser.add_argument("--stats", help="Output file path for basic sequence statistics")
    parser.add_argument("--csv", help="Specify output file path for raw length data export")
    
    
    parser.add_argument("--kmer", type=int, metavar="K", help="Calculate and plot K-mer frequency distribution")
    parser.add_argument("--entropy", type=int, metavar="WINDOW", help="Plot Shannon Entropy sequence complexity using a sliding window")
    parser.add_argument("--cpg", action="store_true", help="Detect and plot CpG islands across the sequence")
    parser.add_argument("--gc-skew", action="store_true", help="Plot GC Skew to predict oriC and replication terminus")
    parser.add_argument("--orf", action="store_true", help="Identify Open Reading Frames (ORFs)")
    parser.add_argument("--tandem", action="store_true", help="Identify tandem repeats within the sequence")
    
    parser.add_argument("--coverage", action="store_true", help="Calculate estimated genome coverage from sequencing reads")
    parser.add_argument("--gsize", type=int, default=1000000, help="Expected genome size in bp")
    parser.add_argument("--rlen", type=int, default=150, help="Average read length in bp")
    parser.add_argument("--phred", action="store_true", help="Generate Phred Quality Profile plot")
    parser.add_argument("--rarefaction", action="store_true", help="Generate sequence rarefaction curve")
    
    parser.add_argument("--primer-file", type=str, help="Path to primer FASTA file for thermodynamic analysis")
    parser.add_argument("--protein-file", type=str, help="Path to protein FASTA file for physicochemical profiling")
    parser.add_argument("--hydro-file", type=str, help="Path to protein FASTA file for Kyte-Doolittle hydrophobicity plot")
    parser.add_argument("--snp-files", nargs=2, metavar=('REF', 'QUERY'), help="Paths to reference and query files for SNP detection")
    parser.add_argument("--dotplot-files", nargs=2, metavar=('SEQ1', 'SEQ2'), help="Paths to two sequence files for Dot Plot")

    parser.add_argument("--ncbi", type=str, help="Download FASTA directly from NCBI (Pain Bansho Tenin)")
    parser.add_argument("--phylo", type=str, help="Path to FASTA alignment for phylogenetic tree (Hashirama Phylo Tree)")
    parser.add_argument("--restriction", type=str, help="Path to FASTA file for Restriction Gel (Zabuza Executioner)")
    parser.add_argument("--crispr", type=str, help="Path to FASTA file for SpCas9 CRISPR sgRNA targets (Kakashi Crispr Copy)")
    parser.add_argument("--vcf-table", type=str, help="Path to VCF file to parse variant information (Neji Variant Eye)")
    parser.add_argument("--blast-path", type=str, default="blastn", help="Local path to blastn executable")
    parser.add_argument("--blast-query", type=str, help="Path to query file for external BLAST+")
    parser.add_argument("--blast-subject", type=str, help="Path to subject file for external BLAST+")
    
    
    parser.add_argument("--gff", type=str, help="GFF/BED file path for Yamato Wood Extractor")
    parser.add_argument("--fasta-ref", type=str, help="Reference FASTA for Yamato Wood Extractor")
    parser.add_argument("--gene", type=str, help="Gene/Feature name to extract")
    parser.add_argument("--trim-fastq", type=str, help="FASTQ to trim (Gaara Sand Trimmer)")
    parser.add_argument("--min-phred", type=int, default=20, help="Minimum Phred score for Trimmer")
    parser.add_argument("--convert-in", type=str, help="Input file for Format Shifter")
    parser.add_argument("--convert-out", type=str, help="Output file for Format Shifter")
    parser.add_argument("--fmt-in", type=str, help="Input format (bam, fastq, gbk)")
    parser.add_argument("--fmt-out", type=str, help="Output format (sam, fasta)")

    return parser.parse_args()

def hiruzen_shiki_fuujin_targeted_analysis(args):
    executed = False
    
    if args.ncbi:
        print(f"\n[INFO] Engaging Pain's Bansho Tenin for NCBI Accession: {args.ncbi}...")
        try:
            seq = PainBanshoTeninNCBIFetcher.universal_pull_sequence(args.ncbi)
            with open(f"{args.out_prefix}_{args.ncbi}.fasta", 'w') as f: f.write(seq)
            print(f"[SUCCESS] Sequence written to {args.out_prefix}_{args.ncbi}.fasta")
        except Exception as e: print(f"[ERROR] {e}")
        executed = True

    if args.gff and args.fasta_ref and args.gene:
        print(f"\n[INFO] Activating Yamato Wood Extractor for gene '{args.gene}'...")
        res = YamatoWoodFeatureExtractor.mokuton_extract_feature(args.gff, args.fasta_ref, args.gene)
        print(res)
        executed = True

    if args.trim_fastq:
        print(f"\n[INFO] Activating Gaara Sand Trimmer (Phred > {args.min_phred})...")
        res = GaaraSandTrimmer.sabaku_kyuu_trim_fastq(args.trim_fastq, f"{args.out_prefix}_trimmed.fastq", args.min_phred)
        print(res)
        executed = True

    if args.convert_in and args.convert_out and args.fmt_in and args.fmt_out:
        print(f"\n[INFO] Using Henge Format Shifter ({args.fmt_in} -> {args.fmt_out})...")
        res = HengeFormatShifter.transformation_jutsu_convert(args.convert_in, args.convert_out, args.fmt_in, args.fmt_out)
        print(res)
        executed = True

    if args.phylo:
        try:
            fig = RinneganAdvancedBioTools.hashirama_phylo_tree_builder(args.phylo)
            out = f"{args.out_prefix}_phylo_tree.png"
            fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
            print(f"[SUCCESS] Phylogenetic tree rendered to {out}")
        except Exception as e: print(f"[ERROR] {e}")
        executed = True
        
    if args.restriction:
        seq = gamabunta_summon_sequence_file(args.restriction)
        if seq:
            analyzer = SenjuMokutonAnalyzer()
            analyzer.yamato_mokuton_set_sequence(seq)
            try:
                fig, stats = analyzer.zabuza_executioner_restriction_map()
                out = f"{args.out_prefix}_virtual_gel.png"
                fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                print(f"\n{stats}\n[SUCCESS] Virtual Agarose Gel saved to {out}")
            except Exception as e: print(f"[ERROR] Restriction map failed: {e}")
            executed = True

    if args.crispr:
        seq = gamabunta_summon_sequence_file(args.crispr)
        if seq:
            print("\n" + RinneganAdvancedBioTools.kakashi_crispr_copy_sgrna(seq))
            executed = True

    if args.vcf_table:
        try:
            vars_list = NejiVariantEyeParser.eight_trigrams_parse_vcf(args.vcf_table)
            print(f"\n[INFO] Neji Byakugan extracted {len(vars_list)} variant records.")
            for i, v in enumerate(vars_list[:10]):
                print(f"CHROM: {v['CHROM']} | POS: {v['POS']} | REF: {v['REF']} -> ALT: {v['ALT']}")
            if len(vars_list) > 10: print("...Output truncated.")
        except Exception as e: print(f"[ERROR] VCF parsing failed: {e}")
        executed = True

    if args.blast_query and args.blast_subject:
        q_seq = gamabunta_summon_sequence_file(args.blast_query)
        s_seq = gamabunta_summon_sequence_file(args.blast_subject)
        if q_seq and s_seq:
            print(f"\n[INFO] Triggering Deidara's Explosive BLAST+ Art at {args.blast_path}...")
            print(AnbuBlackOpsExternalTools.deidara_explosive_blast_art(args.blast_path, q_seq, s_seq))
            executed = True
            
    if args.primer_file:
        seq = gamabunta_summon_sequence_file(args.primer_file)
        if seq:
            print("\n" + RinneganAdvancedBioTools.minato_rasengan_primer_wizard(seq))
            executed = True
        
    if args.protein_file:
        seq = gamabunta_summon_sequence_file(args.protein_file)
        if seq:
            print("\n" + RinneganAdvancedBioTools.raikage_lightning_armor_protein_properties(seq))
            executed = True
        
    if args.hydro_file:
        seq = gamabunta_summon_sequence_file(args.hydro_file)
        if seq:
            try:
                fig = RinneganAdvancedBioTools.kisame_water_prison_hydrophobicity_plot(seq)
                out_file = f"{args.out_prefix}_hydrophobicity.png"
                fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                print(f"[INFO] Hydrophobicity plot successfully saved to {out_file}")
            except Exception as e: print(f"[ERROR] Hydrophobicity plot generation failed: {str(e)}")
            executed = True
            
    if args.snp_files:
        ref_seq = gamabunta_summon_sequence_file(args.snp_files[0])
        query_seq = gamabunta_summon_sequence_file(args.snp_files[1])
        if ref_seq and query_seq:
            print("\n" + RinneganAdvancedBioTools.sasuke_sharingan_snp_deducer_global(ref_seq, query_seq))
            executed = True
        
    if args.dotplot_files:
        seq1 = gamabunta_summon_sequence_file(args.dotplot_files[0])
        seq2 = gamabunta_summon_sequence_file(args.dotplot_files[1])
        if seq1 and seq2:
            try:
                fig = RinneganAdvancedBioTools.madara_rinnegan_dotplot_matrix_numpy(seq1, seq2)
                out_file = f"{args.out_prefix}_dotplot.png"
                fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                print(f"[INFO] Alignment Dot Plot successfully saved to {out_file}")
            except Exception as e: print(f"[ERROR] Dot Plot generation failed: {str(e)}")
            executed = True

    return executed

def orochimaru_edo_tensei_largest_contig(analyzer: UchihaSharinganAnalyzer) -> str:
    max_len = 0
    best_seq = ""
    import gzip
    for label, f_path in analyzer.file_paths.items():
        fmt = 'genbank' if neji_byakugan_format_sniffer(f_path) == 'genbank' else 'fasta'
        open_func = gzip.open if f_path.endswith('.gz') else open
        mode = 'rt' if f_path.endswith('.gz') else 'r'
        try:
            with open_func(f_path, mode) as handle:
                for record in SeqIO.parse(handle, fmt):
                    if len(record.seq) > max_len:
                        max_len = len(record.seq)
                        best_seq = str(record.seq)
        except: pass
    return best_seq

def chunin_exam_main_execution():
    args = tenten_kuchiyose_arg_parser()
    
    if args.gui:
        from Interface import konoha_hokage_office_main
        konoha_hokage_office_main()
        return

    targeted_executed = hiruzen_shiki_fuujin_targeted_analysis(args)
    
    if not args.files:
        if not targeted_executed:
            print("[ERROR] Insufficient input. Use --help for usage information.")
        return

    fasta_analyzer = UchihaSharinganAnalyzer()
    fasta_analyzer.expected_genome_size = args.gsize
    fastq_analyzer = UzumakiChakraFastqAnalyzer()
    multi_parser = KatsuyuSlugFormatParser()
    
    labels = args.label if args.label else [None] * len(args.files)
    if len(labels) != len(args.files): labels = [None] * len(args.files)
    
    fasta_loaded, fastq_loaded, other_loaded = False, False, False
    
    for file, label in zip(args.files, labels):
        try:
            fmt = neji_byakugan_format_sniffer(file)
            if fmt in ['fasta', 'genbank']:
                fasta_analyzer.sasuke_sharingan_read_fasta(file, label)
                fasta_loaded = True
            elif fmt == 'fastq':
                fastq_analyzer.naruto_rasengan_read_fastq_blocks(file, label)
                fastq_loaded = True
            elif fmt in ['vcf', 'gff', 'bed', 'sam']:
                res = multi_parser.tsunade_katsuyu_summarize_file(file, fmt)
                print(res['text'])
                other_loaded = True
        except ValueError as e:
            print(f"[ERROR] Processing failure for {file}: {str(e)}")
            continue
            
    if fasta_loaded:
        fasta_analyzer.itachi_susanoo_apply_filters(0, 999999999) 
        if args.stats:
            with open(args.stats, 'w') as f: f.write(ino_mind_transfer_stats_text(fasta_analyzer.filtered_stats))
            print(f"[INFO] Sequence statistics report written to {args.stats}")
        elif not other_loaded and not fastq_loaded:
            print("\n" + ino_mind_transfer_stats_text(fasta_analyzer.filtered_stats))
            
        if args.csv:
            tenten_scroll_export_csv(fasta_analyzer.filtered_lengths, args.csv)
            print(f"[INFO] Raw length distributions exported to {args.csv}")
            
        if args.histogram:
            fig = fasta_analyzer.sasuke_amaterasu_histogram()
            fig.savefig(f"{args.out_prefix}_histogram.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
            
        if args.boxplot:
            fig = fasta_analyzer.itachi_tsukuyomi_boxplot()
            fig.savefig(f"{args.out_prefix}_boxplot.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')

        target_sequence = orochimaru_edo_tensei_largest_contig(fasta_analyzer)
        if target_sequence:
            base_analyzer = SenjuMokutonAnalyzer()
            base_analyzer.yamato_mokuton_set_sequence(target_sequence)
            
            if args.kmer:
                fig, stats = base_analyzer.shino_kikaichu_kmer_swarm(args.kmer)
                fig.savefig(f"{args.out_prefix}_kmer_{args.kmer}.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                print(f"\n[INFO] K-mer Frequency Analysis:\n{stats}")
                
            if args.entropy:
                fig = base_analyzer.orochimaru_curse_mark_entropy_plot(args.entropy)
                fig.savefig(f"{args.out_prefix}_entropy.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                
            if args.cpg:
                fig, stats = base_analyzer.hidan_jashin_cpg_islands_vectorized()
                fig.savefig(f"{args.out_prefix}_cpg_islands.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                print(f"\n[INFO] {stats}")
                
            if args.gc_skew:
                fig = RinneganAdvancedBioTools.kakashi_kamui_gc_skew_rolling(target_sequence)
                fig.savefig(f"{args.out_prefix}_gc_skew.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
                
            if args.orf: print("\n" + RinneganAdvancedBioTools.jiraiya_sage_mode_orf_finder(target_sequence))
            if args.tandem: print("\n" + RinneganAdvancedBioTools.choji_expansion_jutsu_tandem_repeats(target_sequence))

    if fastq_loaded:
        if args.coverage:
            fastq_analyzer.naruto_oodama_rasengan_calculate_genomes(args.gsize, args.rlen, paired_end=True)
            print("\n" + fastq_analyzer.uzumaki_fuuinjutsu_get_results_text())
            
        if args.phred:
            fig = fastq_analyzer.itachi_tsukuyomi_phred_plot()
            fig.savefig(f"{args.out_prefix}_phred_quality.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')
            
        if args.rarefaction:
            fig = fastq_analyzer.pain_shinra_tensei_rarefaction_curve()
            fig.savefig(f"{args.out_prefix}_rarefaction.png", dpi=300, bbox_inches='tight', facecolor='#1a1a1a')

if __name__ == "__main__":
    chunin_exam_main_execution()
