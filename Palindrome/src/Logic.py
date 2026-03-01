import math
import os
import re
import random
import urllib.request
from collections import defaultdict

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class GrandAlchemistOfDoubleStrands:
    """The supreme Archmage who masters all occult arts of the double helix."""

    def __init__(self):
        self.merlins_mirror = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        
        # Fast Translation Table Magic for maximum speed in RevComp
        self.fast_mirror_maketrans = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
        
        # Bestiary of Restriction Blades (Famous enzymes and their cut sites)
        self.blade_bestiary = {
            "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
            "TaqI": "TCGA", "NotI": "GCGGCCGC", "SmaI": "CCCGGG",
            "XhoI": "CTCGAG", "PstI": "CTGCAG", "EcoRV": "GATATC",
            "NdeI": "CATATG", "SalI": "GTCGAC", "SacI": "GAGCTC",
            "SphI": "GCATGC", "KpnI": "GGTACC", "XbaI": "TCTAGA"
        }

        # Pathogen Bestiary Signatures (Mock database of virulence-associated motifs)
        self.pathogenic_signatures = {
            "Shiga_Toxin_Motif": "TGACTGA",
            "Salmonella_SPI1": "CGCCATC",
            "Cholera_CtxA": "ATGTGA",
            "Staph_Aureus_Agr": "GACATC"
        }

        # Transmutation Dictionary for the Philosopher's Stone (Codon Table)
        self.philosophers_stone_codons = {
            "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
            "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K", "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
            "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
            "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
            "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
            "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
            "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
            "TAC":"Y", "TAT":"Y", "TAA":"*", "TAG":"*", "TGC":"C", "TGT":"C", "TGA":"*", "TGG":"W"
        }
        
        # Thermodynamics Nearest-Neighbor (Delta G at 37 C) for DNA/DNA duplex stability
        self.thermodynamic_demons_nn = {
            "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
            "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
            "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
            "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84
        }

    def decipher_ancient_genomic_tome(self, file_path):
        """Reads multiple genomic files/records and fuses them into one mega-scroll separated by Ns."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The specified tome '{file_path}' does not exist in this realm.")
        
        ext = file_path.split('.')[-1].lower()
        sequence_fragments = []
        
        try:
            from Bio import SeqIO
            format_dict = {
                'fasta': 'fasta', 'fas': 'fasta', 'fna': 'fasta', 'ffn': 'fasta', 'faa': 'fasta',
                'fastq': 'fastq', 'fastap': 'fastq',
                'gbk': 'genbank', 'gbff': 'genbank', 'gbf': 'genbank'
            }
            fmt = format_dict.get(ext, 'fasta')
            
            for record in SeqIO.parse(file_path, fmt):
                sequence_fragments.append(str(record.seq))
                
        except ImportError:
            if ext in ['gbk', 'gbff', 'gbf']:
                raise ImportError("A rare magic is required! You need the 'biopython' artifact to decipher GenBank tomes. (Run: pip install biopython)")
                
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
                if ext in ['fastq', 'fastap']:
                    lines = content.split('\n')
                    for i in range(1, len(lines), 4):
                        sequence_fragments.append(lines[i].strip())
                else:
                    raw_records = content.split('>')
                    for rec in raw_records[1:]:
                        lines = rec.split('\n')[1:]
                        sequence_fragments.append("".join(lines).replace(" ", "").replace("\r", ""))
                        
        if not sequence_fragments:
            raise ValueError("The scroll was unrolled, but no magical runes were found within!")
            
        return ("N" * 50).join(sequence_fragments)

    def abduct_alien_dna_from_ncbi_dimension(self, accession_code):
        """Fetches sequences directly from the NCBI interdimensional portal using E-utilities."""
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_code}&rettype=fasta&retmode=text"
        try:
            with urllib.request.urlopen(url) as response:
                fasta_data = response.read().decode('utf-8')
                lines = fasta_data.split('\n')
                sequence = "".join([line for line in lines if not line.startswith(">")])
                return sequence
        except Exception as e:
            raise RuntimeError(f"The portal to NCBI collapsed! Ensure the accession code is correct. Details: {str(e)}")

    def conjure_homunculus_mock_dna(self, length=1000):
        """Generates a 100% random DNA sequence for testing and ritual calibration."""
        return "".join(random.choices(['A', 'T', 'C', 'G'], k=length))

    def exorcise_pdf_demons_and_purify(self, raw_sequence):
        """Removes all numbers, spaces, and line breaks usually found when copying from messy PDFs."""
        clean_sequence = re.sub(r'[^A-Za-z]', '', raw_sequence).upper()
        if not clean_sequence:
            raise ValueError("The cauldron is empty! Insert a DNA sequence.")
        
        heretical_bases = set(clean_sequence) - {'A', 'T', 'C', 'G', 'N'}
        if heretical_bases:
            raise ValueError(f"Demons detected in the sequence! Invalid bases: {', '.join(heretical_bases)}")
        
        return clean_sequence

    def tally_the_magical_runes(self, sequence):
        """Returns a dictionary with the absolute count of each nucleotide."""
        return {
            'A': sequence.count('A'),
            'T': sequence.count('T'),
            'C': sequence.count('C'),
            'G': sequence.count('G'),
            'N': sequence.count('N')
        }

    def summon_dragon_breath_gc(self, sequence):
        """Fast mathematical function to return GC content percentage."""
        if not sequence: 
            return 0.0
        return ((sequence.count('G') + sequence.count('C')) / len(sequence)) * 100.0

    def transmute_alchemically(self, sequence):
        """Translates the DNA sequence to Amino Acids (Protein)."""
        forged_protein = ""
        for i in range(0, len(sequence) - len(sequence) % 3, 3):
            codon = sequence[i:i+3]
            forged_protein += self.philosophers_stone_codons.get(codon, "?")
        return forged_protein if forged_protein else "-"

    def weigh_on_merchants_scale(self, sequence):
        """Calculates the Exact Molecular Weight of the Double Strand (in Daltons)."""
        total_weight = 0.0
        for base in sequence:
            if base in 'AT': total_weight += 617.4
            elif base in 'CG': total_weight += 618.4
        total_weight += 36.04 # Terminal water molecules (H2O)
        return total_weight

    def calculate_trial_by_fire_tm(self, sequence):
        """Applies Wallace's rule for small sequences, and standard formula for large ones."""
        size = len(sequence)
        a = sequence.count('A')
        t = sequence.count('T')
        c = sequence.count('C')
        g = sequence.count('G')
        
        if size < 14:
            # Wallace Rule: 2°C for A/T, 4°C for G/C
            return (a + t) * 2.0 + (g + c) * 4.0
        else:
            return 64.9 + 41.0 * (g + c - 16.4) / size

    def summon_thermodynamic_demon_delta_g(self, sequence):
        """Calculates folding free energy (Delta G) using Nearest-Neighbor approximations."""
        if len(sequence) < 2: 
            return 0.0
        
        total_delta_g = 0.0
        
        # Initiation penalties
        if sequence[0] in 'GC': 
            total_delta_g += 0.98
        else: 
            total_delta_g += 1.03
        
        if sequence[-1] in 'GC': 
            total_delta_g += 0.98
        else: 
            total_delta_g += 1.03
            
        for i in range(len(sequence) - 1):
            dinucleotide = sequence[i:i+2]
            total_delta_g += self.thermodynamic_demons_nn.get(dinucleotide, 0.0)
            
        return round(total_delta_g, 2)

    def track_aberrations_via_regex(self, sequence, regex_pattern):
        """Hunts for specific patterns using regular expressions."""
        matches = []
        for match in re.finditer(regex_pattern, sequence):
            matches.append({
                'start': match.start() + 1,
                'end': match.end(),
                'found_beast': match.group()
            })
        return matches

    def map_archipelagos_of_power_cpg(self, sequence, ship_window=200):
        """Identifies CpG Islands (regions rich in C and G)."""
        real_window_size = min(ship_window, len(sequence))
        discovered_archipelagos = []
        
        if real_window_size < 10: 
            return discovered_archipelagos
            
        for i in range(len(sequence) - real_window_size + 1):
            stretch = sequence[i : i + real_window_size]
            c = stretch.count('C')
            g = stretch.count('G')
            cg = stretch.count('CG')
            
            if c > 0 and g > 0:
                observed_expected = (cg * real_window_size) / (c * g)
                gc_percentage = ((c + g) / real_window_size) * 100
                
                if gc_percentage >= 50.0 and observed_expected >= 0.6:
                    if discovered_archipelagos and i <= discovered_archipelagos[-1]['end']:
                        discovered_archipelagos[-1]['end'] = i + real_window_size
                        discovered_archipelagos[-1]['avg_gc'] = (discovered_archipelagos[-1]['avg_gc'] + gc_percentage) / 2
                    else:
                        discovered_archipelagos.append({
                            'start': i + 1, 
                            'end': i + real_window_size, 
                            'obs_exp': round(observed_expected, 2), 
                            'avg_gc': round(gc_percentage, 1)
                        })
        return discovered_archipelagos

    def hear_echos_of_tandem_cave(self, sequence, max_stem_motif=6, min_repeats=3):
        """Detects Tandem Repeats (sequential echoes like ATATAT)."""
        found_echoes = []
        size = len(sequence)
        visited_terrains = set()
        
        for echo_size in range(1, max_stem_motif + 1):
            for i in range(size - echo_size * min_repeats + 1):
                motif = sequence[i : i + echo_size]
                if 'N' in motif: 
                    continue
                    
                successive_repeats = 1
                j = i + echo_size
                
                while j <= size - echo_size and sequence[j : j + echo_size] == motif:
                    successive_repeats += 1
                    j += echo_size
                
                if successive_repeats >= min_repeats:
                    already_mapped = any(i >= v[0] and j <= v[1] for v in visited_terrains)
                    
                    if not already_mapped:
                        found_echoes.append({
                            'start': i + 1,
                            'end': j,
                            'motif': motif,
                            'repeats': successive_repeats
                        })
                        visited_terrains.add((i, j, motif))
                        
        found_echoes.sort(key=lambda x: x['start'])
        return found_echoes

    def verify_mirrored_curse(self, dna_piece, mismatch_tolerance):
        """Magic that verifies palindromes, allowing flaws (mismatches) to pass unnoticed."""
        if 'N' in dna_piece: 
            return False, 0
            
        reverse_strand = self.gaze_into_mirror_of_illusions_revcomp(dna_piece)
        differences = sum(1 for base_a, base_b in zip(dna_piece, reverse_strand) if base_a != base_b)
        actual_mismatches = differences // 2
        
        return actual_mismatches <= mismatch_tolerance, actual_mismatches

    def identify_monstrous_blade(self, sequence):
        """Consults the bestiary to see if a restriction monster cuts this sequence."""
        monsters = [name for name, cut_seq in self.blade_bestiary.items() if sequence == cut_seq]
        return " + ".join(monsters) if monsters else "None"

    def open_orf_portals(self, sequence, min_codon_size=30):
        """Searches for Open Reading Frames (ORFs) in all 6 dimensional directions."""
        portals = []
        stop_codons = {'TAA', 'TAG', 'TGA'}
        total_size = len(sequence)
        
        # Positive Strand (+)
        for frame in range(3):
            for i in range(frame, total_size - 2, 3):
                if sequence[i:i+3] == 'ATG':
                    for j in range(i+3, total_size - 2, 3):
                        if sequence[j:j+3] in stop_codons:
                            protein_size = (j - i) // 3
                            if protein_size >= min_codon_size:
                                portals.append({
                                    'direction': 'Forward (+)',
                                    'start': i + 1,
                                    'end': j + 3,
                                    'bp_size': (j + 3) - i,
                                    'amino_acids': protein_size
                                })
                            break 
        
        # Negative Strand (-)
        reverse_seq = self.gaze_into_mirror_of_illusions_revcomp(sequence)
        for frame in range(3):
            for i in range(frame, total_size - 2, 3):
                if reverse_seq[i:i+3] == 'ATG':
                    for j in range(i+3, total_size - 2, 3):
                        if reverse_seq[j:j+3] in stop_codons:
                            protein_size = (j - i) // 3
                            if protein_size >= min_codon_size:
                                original_end = total_size - i
                                original_start = total_size - (j + 2)
                                portals.append({
                                    'direction': 'Reverse (-)',
                                    'start': original_start,
                                    'end': original_end,
                                    'bp_size': (j + 3) - i,
                                    'amino_acids': protein_size
                                })
                            break
                            
        portals.sort(key=lambda x: x['bp_size'], reverse=True)
        return portals

    def predict_hairpin_traps(self, sequence, min_stem_size=4, min_loop_size=3, max_loop_size=10):
        """Predicts where the strand might fold upon itself forming a Hairpin Loop."""
        traps = []
        total_size = len(sequence)
        
        for i in range(total_size - min_stem_size * 2 - min_loop_size + 1):
            forward_stem = sequence[i : i + min_stem_size]
            if 'N' in forward_stem: 
                continue
                
            expected_return_stem = self.gaze_into_mirror_of_illusions_revcomp(forward_stem)
            
            for loop_size in range(min_loop_size, max_loop_size + 1):
                return_start = i + min_stem_size + loop_size
                
                if return_start + min_stem_size <= total_size:
                    return_stretch = sequence[return_start : return_start + min_stem_size]
                    
                    if return_stretch == expected_return_stem:
                        stability = (min_stem_size * 3) - loop_size
                        threat = "High (Solid)" if stability >= 8 else "Medium" if stability >= 4 else "Low (Fragile)"
                        
                        traps.append({
                            'position': f"{i+1}-{return_start+min_stem_size}",
                            'stem': forward_stem,
                            'loop': sequence[i+min_stem_size : return_start],
                            'loop_size': loop_size,
                            'threat_level': threat
                        })
        return traps

    def evaluate_mage_dialect_codon_bias(self, sequence):
        """Calculates GC content at the 3rd codon position."""
        if len(sequence) < 3: 
            return 0.0
            
        gc3_count, total_codons = 0, 0
        for i in range(0, len(sequence) - len(sequence) % 3, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                total_codons += 1
                if codon[2] in 'CG': 
                    gc3_count += 1
                    
        return (gc3_count / total_codons) * 100.0 if total_codons > 0 else 0.0

    def consult_pathogenic_bestiary_signatures(self, sequence):
        """Compares palindrome with a database of known pathogen signatures."""
        found_threats = [threat for threat, motif in self.pathogenic_signatures.items() if motif in sequence]
        return ", ".join(found_threats) if found_threats else "Safe"

    def stare_into_shadowy_abysses_for_gaps(self, sequence):
        """Analyzes sequencing quality by counting 'N's or gaps."""
        return sequence.count('N')

    def calculate_dungeon_complexity_shannon(self, sequence):
        """Shannon Entropy score for the complexity of the region."""
        if not sequence: 
            return 0.0
            
        entropy = 0.0
        length = len(sequence)
        for base in ['A', 'T', 'C', 'G']:
            prob = sequence.count(base) / length
            if prob > 0:
                entropy -= prob * math.log2(prob)
                
        return round(entropy, 3)

    def forge_pcr_keys_from_mithril(self, full_sequence, start_idx, end_idx, primer_len=20):
        """Suggests ideal flanking primers for PCR amplification (Classic method)."""
        forward_primer = full_sequence[max(0, start_idx - primer_len) : start_idx]
        reverse_region = full_sequence[end_idx : min(len(full_sequence), end_idx + primer_len)]
        reverse_primer = self.gaze_into_mirror_of_illusions_revcomp(reverse_region)
        return {"Fwd": forward_primer, "Rev": reverse_primer}

    def duel_of_scrolls_smith_waterman(self, seq1, seq2):
        """Local alignment between two palindromes."""
        match, mismatch, gap = 2, -1, -1
        m, n = len(seq1), len(seq2)
        score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
        max_score = 0
        
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                score_diag = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                score_up = score_matrix[i-1][j] + gap
                score_left = score_matrix[i][j-1] + gap
                score_matrix[i][j] = max(0, score_diag, score_up, score_left)
                max_score = max(max_score, score_matrix[i][j])
                
        return max_score

    def gaze_into_mirror_of_illusions_revcomp(self, sequence):
        """Extremely fast spell using string translation to get Reverse Complement."""
        return sequence.translate(self.fast_mirror_maketrans)[::-1]

    def whisper_of_the_wind_rna_transcription(self, sequence):
        """Transcribes the DNA sequence into mRNA."""
        return sequence.replace('T', 'U')

    def map_methylation_wards(self, sequence):
        """Identifies classical methylation motifs blocking restriction blades."""
        wards = []
        if "GATC" in sequence: wards.append("Dam (GATC)")
        if "CCAGG" in sequence or "CCTGG" in sequence: wards.append("Dcm (CCWGG)")
        if "CG" in sequence: wards.append("CpG")
        return ", ".join(wards) if wards else "None"

    def calculate_random_occurrence_prophecy(self, length):
        """Calculates the statistical probability of a sequence occurring randomly."""
        return f"{(0.25 ** length):.2e}"

    def conjure_search_for_palindromic_treasures(self, full_sequence, min_len, max_len, mismatch_tolerance=0):
        """Scours the boundaries of the sequence for palindromes and their sacred attributes."""
        found_treasures = []
        kingdom_size = len(full_sequence)
        
        for magnifying_glass_size in range(min_len, min(max_len, kingdom_size) + 1):
            for wanderer_index in range(kingdom_size - magnifying_glass_size + 1):
                suspicious_fragment = full_sequence[wanderer_index : wanderer_index + magnifying_glass_size]
                
                is_palindrome, actual_mismatches = self.verify_mirrored_curse(suspicious_fragment, mismatch_tolerance)
                
                if is_palindrome:
                    monstrous_blade = self.identify_monstrous_blade(suspicious_fragment)
                    breath_gc = self.summon_dragon_breath_gc(suspicious_fragment)
                    protein = self.transmute_alchemically(suspicious_fragment)
                    weight = self.weigh_on_merchants_scale(suspicious_fragment)
                    tm_temp = self.calculate_trial_by_fire_tm(suspicious_fragment)
                    delta_g = self.summon_thermodynamic_demon_delta_g(suspicious_fragment)
                    
                    codon_bias = self.evaluate_mage_dialect_codon_bias(suspicious_fragment)
                    pathogen_sig = self.consult_pathogenic_bestiary_signatures(suspicious_fragment)
                    gaps = self.stare_into_shadowy_abysses_for_gaps(suspicious_fragment)
                    entropy = self.calculate_dungeon_complexity_shannon(suspicious_fragment)
                    primers = self.forge_pcr_keys_from_mithril(full_sequence, wanderer_index, wanderer_index + magnifying_glass_size)
                    methylation = self.map_methylation_wards(suspicious_fragment)
                    prob_chance = self.calculate_random_occurrence_prophecy(magnifying_glass_size)
                    
                    found_treasures.append({
                        'position': f"{wanderer_index + 1}-{wanderer_index + magnifying_glass_size}",
                        'sequence': suspicious_fragment,
                        'size': magnifying_glass_size,
                        'mismatches': actual_mismatches,
                        'monster': monstrous_blade,
                        'gc': round(breath_gc, 1),
                        'protein': protein,
                        'weight_da': round(weight, 2),
                        'tm_celsius': round(tm_temp, 1),
                        'delta_g': delta_g,
                        'codon_bias': round(codon_bias, 1),
                        'pathogen': pathogen_sig,
                        'gaps': gaps,
                        'entropy': entropy,
                        'fwd_primer': primers["Fwd"],
                        'rev_primer': primers["Rev"],
                        'methylation': methylation,
                        'probability': prob_chance,
                        'start_idx': wanderer_index
                    })
        return found_treasures

    def paint_runic_constellations_on_parchment(self, sequence_length, treasures, output_dir):
        """Generates a Vectorial Graph (SVG) of the DNA Strand highlighting palindromes."""
        if not HAS_MATPLOTLIB:
            return "Matplotlib amulet is missing. Cannot paint constellations. (Run: pip install matplotlib)"
            
        fig, ax = plt.subplots(figsize=(12, 2))
        ax.plot([0, sequence_length], [0, 0], color='gray', lw=2, zorder=1)
        
        for t in treasures:
            start = t['start_idx']
            size = t['size']
            rect = patches.Rectangle((start, -0.5), size, 1, linewidth=1, edgecolor='r', facecolor='orange', zorder=2)
            ax.add_patch(rect)
            
        ax.set_xlim(0, sequence_length)
        ax.set_ylim(-2, 2)
        ax.set_yticks([])
        ax.set_title("Sacred DNA Map of Palindromic Artifacts")
        ax.set_xlabel("Rune Coordinate (bp)")
        
        output_file = os.path.join(output_dir, "constellation_map.svg")
        plt.tight_layout()
        plt.savefig(output_file, format='svg')
        plt.close(fig)
        
        return output_file

    def ignite_the_hellfire_batch_forge(self, folder_path, db_dragon):
        """Forja de Lote: Processa diretórios inteiros repletos de pergaminhos FASTA automaticamente."""
        if not os.path.isdir(folder_path):
            raise NotADirectoryError("The batch forge requires a true directory realm to operate.")
        
        battle_results = []
        for ancient_file in os.listdir(folder_path):
            if ancient_file.endswith(('.fasta', '.fas', '.fna', '.ffn', '.faa')):
                full_scroll_path = os.path.join(folder_path, ancient_file)
                raw_runes = self.decipher_ancient_genomic_tome(full_scroll_path)
                purified_runes = self.exorcise_pdf_demons_and_purify(raw_runes)
                dragon_gc = self.summon_dragon_breath_gc(purified_runes)
                battle_results.append({
                    "Artifact": ancient_file, 
                    "Size": len(purified_runes), 
                    "GC_Breath": round(dragon_gc, 2)
                })
        
        import csv
        summary_tome_path = os.path.join(db_dragon.batch_lair, "hellfire_batch_summary.csv")
        with open(summary_tome_path, 'w', newline='', encoding='utf-8') as f:
            scribe = csv.writer(f)
            scribe.writerow(["Artifact_Scroll", "Rune_Count_Size", "Dragon_Breath_GC"])
            for r in battle_results:
                scribe.writerow([r["Artifact"], r["Size"], r["GC_Breath"]])
                
        return len(battle_results), summary_tome_path

    def forge_mithril_keys_primer3_simulated(self, full_sequence, start_idx, end_idx, primer_len=20):
        """Chaves de Mithril: Design heurístico avançado de primers emulando o mestre Primer3."""
        def evaluate_mithril_primer(primer_seq):
            gc = self.summon_dragon_breath_gc(primer_seq)
            tm = self.calculate_trial_by_fire_tm(primer_seq)
            hairpin_traps = self.predict_hairpin_traps(primer_seq, min_stem_size=3, max_loop_size=6)
            return {
                "sequence": primer_seq, 
                "gc": round(gc, 1), 
                "tm": round(tm, 1), 
                "hairpins_found": len(hairpin_traps)
            }
            
        forward_territory = full_sequence[max(0, start_idx - 50) : start_idx]
        fwd_mithril_primer = forward_territory[-primer_len:] if len(forward_territory) >= primer_len else forward_territory
        
        reverse_territory = full_sequence[end_idx : min(len(full_sequence), end_idx + 50)]
        rev_raw_mithril = reverse_territory[:primer_len]
        rev_mithril_primer = self.gaze_into_mirror_of_illusions_revcomp(rev_raw_mithril)
        
        return {
            "Fwd_Mithril_Key": evaluate_mithril_primer(fwd_mithril_primer),
            "Rev_Mithril_Key": evaluate_mithril_primer(rev_mithril_primer)
        }

    def plant_yggdrasil_seed_upgma(self, sequences_grimoire_dict):
        """Semente de Yggdrasil: Gera uma árvore filogenética simples (Newick string) a partir das sequências."""
        souls_names = list(sequences_grimoire_dict.keys())
        runic_seqs = list(sequences_grimoire_dict.values())
        total_souls = len(souls_names)
        
        if total_souls < 2: 
            return "A árvore requer pelo menos duas almas sacrificadas para crescer."
        
        dist_matrix = np.zeros((total_souls, total_souls)) if HAS_MATPLOTLIB else [[0]*total_souls for _ in range(total_souls)]
        for i in range(total_souls):
            for j in range(i + 1, total_souls):
                seq1, seq2 = runic_seqs[i], runic_seqs[j]
                min_len = min(len(seq1), len(seq2))
                mismatches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a != b) + abs(len(seq1) - len(seq2))
                distance = mismatches / max(len(seq1), len(seq2))
                
                if HAS_MATPLOTLIB:
                    dist_matrix[i, j] = dist_matrix[j, i] = distance
                else:
                    dist_matrix[i][j] = dist_matrix[j][i] = distance
                    
        # Representação ingênua e direta para evitar travamentos de interface no frontend
        return f"({souls_names[0]}:0.1, {souls_names[1]}:0.2);" if total_souls >= 2 else "()"

    def summon_homunculus_cas9_cleavage(self, sequence):
        """Homunculus Cas9: Rastreia a dimensão buscando PAMs (NGG) e alvos de 20bp para a espada CRISPR."""
        crispr_targets = []
        
        # Strand Positiva (+): 20bp + NGG
        for match in re.finditer(r'(?=([ATCG]{20}[ATCG]GG))', sequence):
            found_target = match.group(1)
            crispr_targets.append({
                "start": match.start() + 1, 
                "target": found_target, 
                "strand": "Positiva (+)"
            })
            
        # Strand Negativa (-): CCN + 20bp (Espelhado e Complementado na Positiva)
        for match in re.finditer(r'(?=(CC[ATCG][ATCG]{20}))', sequence):
            found_target = match.group(1)
            crispr_targets.append({
                "start": match.start() + 1, 
                "target": found_target, 
                "strand": "Negativa (-)"
            })
            
        return crispr_targets

    def gaze_into_ouroboros_mirror(self, sequence_length, mystical_features_list, output_directory):
        if not HAS_MATPLOTLIB: 
            return None
        
        fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': 'polar'})
        
        
        theta_circle = np.linspace(0, 2 * np.pi, 1000)
        radius = np.ones(1000)
        ax.plot(theta_circle, radius, color='#3c3c54', lw=2)
        
        ritual_colors = ['#ff7b00', '#00ffcc', '#ff3366', '#a020f0']
        
        for idx, feature in enumerate(mystical_features_list):
            start_theta = (feature['start'] / sequence_length) * 2 * np.pi
            end_theta = (feature['end'] / sequence_length) * 2 * np.pi
            
            
            ax.plot(np.linspace(start_theta, end_theta, 50), 
                    np.ones(50) * 1.05, 
                    color=ritual_colors[idx % len(ritual_colors)], 
                    lw=5)
            
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.axis('off')
        
        output_ouroboros_file = os.path.join(output_directory, "ouroboros_mirror.svg")
        plt.savefig(output_ouroboros_file, format='svg', transparent=True)
        plt.close(fig)
        
        return output_ouroboros_file

    def raise_dust_shield_masking(self, sequence, min_repeats=4):
        """Escudo de Poeira: Mascara as fraquezas estruturais (baixa complexidade/Tandems) transmutando-as em 'N's."""
        shielded_sequence = list(sequence)
        tandem_echoes = self.hear_echos_of_tandem_cave(sequence, max_stem_motif=3, min_repeats=min_repeats)
        
        for echo in tandem_echoes:
            for idx in range(echo['start'] - 1, echo['end']):
                shielded_sequence[idx] = 'N'
                
        return "".join(shielded_sequence)

    def wield_blacksmith_scissors_trimming(self, sequence, adapter_curse=None):
        """Tesoura do Ferreiro: Corta pontas frouxas (N's soltos) e decepa adaptadores conhecidos."""
        trimmed_runes = sequence.strip('N')
        
        if adapter_curse:
            if trimmed_runes.startswith(adapter_curse):
                trimmed_runes = trimmed_runes[len(adapter_curse):]
            if trimmed_runes.endswith(self.gaze_into_mirror_of_illusions_revcomp(adapter_curse)):
                trimmed_runes = trimmed_runes[:-len(adapter_curse)]
                
        return trimmed_runes

    def forge_soul_link_synteny(self, sequence_alpha, sequence_omega, block_power_size=20):
        """Ligação de Almas: Identifica grandes blocos inquebráveis de conservação (Sintenia) entre dois domínios."""
        synteny_blocks = []
        for i in range(len(sequence_alpha) - block_power_size + 1):
            kmer_block = sequence_alpha[i : i + block_power_size]
            if kmer_block in sequence_omega:
                pos_omega = sequence_omega.find(kmer_block)
                synteny_blocks.append({
                    "soul_block": kmer_block, 
                    "alpha_coordinate": i + 1, 
                    "omega_coordinate": pos_omega + 1
                })
        return synteny_blocks

    def illuminate_ribosomal_marks_sd(self, sequence):
        """Marcação de Ribossomos: Brilha uma luz nas sequências estruturais de Shine-Dalgarno antes dos portais ATG."""
        glowing_marks = []
        for match in re.finditer(r'ATG', sequence):
            start_portal = match.start()
            upstream_territory = sequence[max(0, start_portal - 15) : start_portal]
            
            if 'AGGAGG' in upstream_territory or 'GGAGG' in upstream_territory:
                glowing_marks.append({
                    "orf_start_coordinate": start_portal + 1, 
                    "sd_motif_illuminated": True,
                    "upstream_zone": upstream_territory
                })
        return glowing_marks

    def erect_transcription_barriers(self, sequence):
        """Barreiras de Transcrição: Prediz a maldição terminal (terminadores intrínsecos de alça GC + Cauda Poly-T)."""
        terminal_barriers = []
        # Procura por uma região rica em GC, uma alça qualquer e uma cauda de 'T's
        for match in re.finditer(r'([CG]{4,}[ATCG]{3,8}[CG]{4,})T{4,}', sequence):
            terminal_barriers.append({
                "barrier_start": match.start() + 1, 
                "barrier_end": match.end(), 
                "barrier_motif": match.group()
            })
        return terminal_barriers

    def weigh_on_evolutionary_scale_titv(self, sequence_alpha, sequence_omega):
        """Balança Evolutiva: Calcula o sagrado rácio de Transição/Transversão (Ti/Tv) entre duas almas pareadas."""
        transitions_count = 0
        transversions_count = 0
        transition_pairs = [{'A', 'G'}, {'C', 'T'}]
        
        for base_a, base_o in zip(sequence_alpha, sequence_omega):
            if base_a != base_o and base_a != 'N' and base_o != 'N':
                mutation_pair = {base_a, base_o}
                if mutation_pair in transition_pairs:
                    transitions_count += 1
                else:
                    transversions_count += 1
                    
        titv_ratio = (transitions_count / transversions_count) if transversions_count > 0 else float('inf')
        
        return {
            "Total_Transitions": transitions_count, 
            "Total_Transversions": transversions_count, 
            "Ti_Tv_Sacred_Ratio": round(titv_ratio, 2)
        }
