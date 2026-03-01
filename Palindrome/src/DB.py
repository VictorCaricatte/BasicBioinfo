import os
import csv
import subprocess
from Config import ScrollOfEternalTruths

class GuardianDragonOfMatrices:
    """Ancestral beast that spews local matrices, devours CSV data, and cleans up the mess."""
    
    def __init__(self):
        self.scroll = ScrollOfEternalTruths()
        self.bone_lair = self.scroll.ancestral_memories["database_cavern"]
        self.batch_lair = self.scroll.ancestral_memories["batch_forge_cavern"]
        self.carve_lair_in_stone()

    def carve_lair_in_stone(self):
        """Hollows out the directory for the databases and batch processing."""
        if not os.path.exists(self.bone_lair):
            os.makedirs(self.bone_lair)
        if not os.path.exists(self.batch_lair):
            os.makedirs(self.batch_lair)

    def devour_csv_tribute_advanced(self, tome_name, discovered_treasures):
        """Saves the palindrome results in a highly detailed CSV scroll."""
        full_path = os.path.join(self.bone_lair, tome_name)
        with open(full_path, mode='w', newline='', encoding='utf-8') as cursed_file:
            scribe = csv.writer(cursed_file)
            scribe.writerow([
                "Map_Coordinates", "DNA_Spell", "Size_Power", "Mismatches", 
                "GC_Percent", "Transmuted_Protein", "Tm_Celsius", "Delta_G_Energy", 
                "Pathogen_Signature", "Methylation_Wards", "Entropy_Complexity"
            ])
            for t in discovered_treasures:
                scribe.writerow([
                    t['position'], t['sequence'], t['size'], t['mismatches'],
                    t['gc'], t['protein'], t['tm_celsius'], t.get('delta_g', 0.0),
                    t['pathogen'], t['methylation'], t['entropy']
                ])
        return full_path

    def inscribe_ultimate_parchment_html(self, tome_name, spoils):
        """Generates an interactive HTML final report (Pergaminho Definitivo)."""
        full_path = os.path.join(self.bone_lair, tome_name)
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Ultimate Parchment - Genomic Analysis</title>
            <style>
                body {{ background-color: #121212; color: #d1d1e0; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 20px; }}
                h1, h2, h3 {{ color: #ff7b00; }}
                table {{ border-collapse: collapse; width: 100%; margin-top: 20px; background-color: #1e1e2c; }}
                th, td {{ border: 1px solid #3c3c54; padding: 10px; text-align: left; }}
                th {{ background-color: #2b2b3c; color: #00ffcc; }}
                .highlight {{ color: #00ffcc; font-weight: bold; }}
            </style>
        </head>
        <body>
            <h1>📜 The Ultimate Parchment of Genomic Discovery</h1>
            <h2>Global Dragon's Breath (GC%): <span class="highlight">{spoils.get('global_gc', 0):.2f}%</span></h2>
            <h2>Nucleotide Tally: <span class="highlight">{spoils.get('tally', {{}})}</span></h2>
            
            <h3>💎 Unearthed Palindromic Relics ({len(spoils.get('palindromes', []))})</h3>
            <table>
                <tr>
                    <th>Coordinates</th><th>Rune Sequence</th><th>Size</th><th>GC%</th><th>ΔG (kcal/mol)</th><th>Pathogen Sig.</th>
                </tr>
        """
        for pal in spoils.get('palindromes', []):
            html_content += f"""
                <tr>
                    <td>{pal['position']}</td><td>{pal['sequence']}</td><td>{pal['size']}</td>
                    <td>{pal['gc']}%</td><td>{pal['delta_g']}</td><td>{pal['pathogen']}</td>
                </tr>
            """
            
        html_content += """
            </table>
        </body>
        </html>
        """
        
        with open(full_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        return full_path

    def spew_blast_matrix(self, input_fasta_file):
        """Uses BLAST magic to forge a nucleotide database."""
        blast_path = self.scroll.ancestral_memories.get("blast_path", "")
        if not blast_path or not os.path.exists(blast_path):
            raise FileNotFoundError("The dragon searched, but the BLAST sword is not in the armory!")
        
        shadow_database = os.path.join(self.bone_lair, "monstrous_blast_db")
        elven_command = f'"{blast_path}" -in "{input_fasta_file}" -dbtype nucl -out "{shadow_database}"'
        
        battle_result = subprocess.run(elven_command, shell=True, capture_output=True, text=True)
        if battle_result.returncode != 0:
            raise RuntimeError(f"The BLAST forge failed: {battle_result.stderr}")
        return shadow_database

    def crystallize_diamond_library(self, input_fasta_file):
        """Uses DIAMOND to crystallize a protein/DNA library."""
        diamond_path = self.scroll.ancestral_memories.get("diamond_path", "")
        if not diamond_path or not os.path.exists(diamond_path):
            raise FileNotFoundError("The DIAMOND staff was not found. The magic failed!")
        
        crystal_database = os.path.join(self.bone_lair, "crystal_diamond_db")
        dwarven_command = f'"{diamond_path}" makedb --in "{input_fasta_file}" -d "{crystal_database}"'
        
        battle_result = subprocess.run(dwarven_command, shell=True, capture_output=True, text=True)
        if battle_result.returncode != 0:
            raise RuntimeError(f"The DIAMOND spell shattered: {battle_result.stderr}")
        return crystal_database

    def banish_shadow_remnants_to_void(self):
        """A cleansing spell to delete temporary graphical and CSV artifacts."""
        remnants = [
            os.path.join(self.bone_lair, "constellation_map.svg"),
            os.path.join(self.bone_lair, "ouroboros_mirror.svg"),
            os.path.join(self.bone_lair, "alchemical_treasures.csv"),
            os.path.join(self.bone_lair, "ultimate_parchment.html")
        ]
        banished_count = 0
        for artifact in remnants:
            if os.path.exists(artifact):
                os.remove(artifact)
                banished_count += 1
        return banished_count
