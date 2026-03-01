import sys
import argparse
from Logic import GrandAlchemistOfDoubleStrands
from DB import GuardianDragonOfMatrices

def embark_on_terminal_crusade():
    """CLI mode explicitly focused on the dark lords of the terminal."""
    
    help_epilogue = """
-----------------------------------------------------------------------------------
📖 PALINDROME (GRIMOIRE OF HELP)
-----------------------------------------------------------------------------------
Welcome to the Palindrome. This tool fuses dark fantasy with 
professional bioinformatics to analyze DNA sequences at a profound level.

NEW MAGICS CONJURED:
  - --batch: Invokes the Hellfire Forge to process a full directory of FASTA tomes.
             Example: python Palindrome.py --batch ./my_fasta_folder
  - --crispr: Awakens the Homunculus Cas9 to identify NGG PAM targets.
  - --mask: Casts the Dust Shield, transmuting tandem repeats into 'N's.
  - --trim: Wields the Blacksmith Scissors to strip terminal 'N's.
  - --report: Scribes the Ultimate HTML Parchment (Pergaminho Definitivo).
  - --clean: Banishes all shadow remnants (SVG, CSV, HTML) to the void.
  - Massive Scroll Reading: Handles multi-FASTA files perfectly with the -f flag.
  - PDF Demon Exorcism: Purifies messy PDF copied texts automatically.
-----------------------------------------------------------------------------------
"""

    war_herald = argparse.ArgumentParser(
        description="⚔️ Palindrome CLI Grandmaster - The insane explorer of genomic spells.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=help_epilogue
    )
    
    war_herald.add_argument("-s", "--sequence", type=str, default=None, 
                                  help="The raw DNA potion to be scanned manually.")
    war_herald.add_argument("-f", "--file", type=str, default=None, 
                                  help="Path to an ancient genomic tome (handles multi-FASTA).")
    war_herald.add_argument("-n", "--ncbi", type=str, default=None, 
                                  help="Accession ID to abduct from the NCBI dimension.")
    war_herald.add_argument("--mock", type=int, default=None, 
                                  help="Length of random mock DNA homunculus to conjure for testing.")
    war_herald.add_argument("--batch", type=str, default=None, 
                                  help="Directory path to ignite the Hellfire Batch Forge.")
    war_herald.add_argument("-min", "--minimum", type=int, default=4, 
                                  help="Minimum power (length) of the palindrome.")
    war_herald.add_argument("-max", "--maximum", type=int, default=12, 
                                  help="Maximum power (length) of the palindrome.")
    war_herald.add_argument("-mis", "--mismatches", type=int, default=0, 
                                  help="Tolerated curses (mismatches within the palindrome).")
    war_herald.add_argument("--regex", type=str, default=None, 
                                  help="Search for a regex beast pattern within the sequence.")
    war_herald.add_argument("--crispr", action="store_true", 
                                  help="Invoke Homunculus Cas9 to map CRISPR cleavage sites.")
    war_herald.add_argument("--mask", action="store_true", 
                                  help="Raise the Dust Shield to mask repetitive regions.")
    war_herald.add_argument("--trim", action="store_true", 
                                  help="Wield Blacksmith Scissors to trim poor quality edges.")
    war_herald.add_argument("--report", action="store_true", 
                                  help="Forge the Pergaminho Definitivo (Interactive HTML report).")
    war_herald.add_argument("--clean", action="store_true", 
                                  help="Clean up SVG, CSV, and HTML shadows left in the database cavern.")
    
    scroll_args = war_herald.parse_args()

    print("\n🧙‍♂️ The Hermit Mage lights the bonfire and begins chanting...\n")
    
    supreme_mage = GrandAlchemistOfDoubleStrands()
    db_dragon = GuardianDragonOfMatrices()
    
    if scroll_args.clean:
        cleansed = db_dragon.banish_shadow_remnants_to_void()
        print(f"🧹 Void Cleansing Complete. {cleansed} shadows banished from the physical realm.")
        # If the user only wanted to clean, exit gracefully.
        if not any([scroll_args.file, scroll_args.sequence, scroll_args.ncbi, scroll_args.mock, scroll_args.batch]):
            sys.exit(0)

    # Triggering the Batch Forge bypasses the single-file pipeline
    if scroll_args.batch:
        try:
            print(f"🏭 Igniting the Hellfire Batch Forge in realm: {scroll_args.batch}")
            scrolls_forged, summary_tome = supreme_mage.ignite_the_hellfire_batch_forge(scroll_args.batch, db_dragon)
            print(f"✅ Glorious Success! Forged {scrolls_forged} ancient scrolls.")
            print(f"📄 The Grand Summary Tome has been mapped at: {summary_tome}\n")
            sys.exit(0)
        except Exception as batch_disaster:
            print(f"☠️ The Batch Forge collapsed: {str(batch_disaster)}")
            sys.exit(1)

    try:
        raw_sequence = ""
        
        # Determine the source of the magical runes
        if scroll_args.mock:
            print(f"🎲 Conjuring a Homunculus of size {scroll_args.mock}...")
            raw_sequence = supreme_mage.conjure_homunculus_mock_dna(scroll_args.mock)
        elif scroll_args.ncbi:
            print(f"🛸 Opening dimensional portal to NCBI for {scroll_args.ncbi}...")
            raw_sequence = supreme_mage.abduct_alien_dna_from_ncbi_dimension(scroll_args.ncbi)
        elif scroll_args.file:
            print(f"📜 Unrolling the ancient tome (supports multi-FASTA): {scroll_args.file}...")
            raw_sequence = supreme_mage.decipher_ancient_genomic_tome(scroll_args.file)
        elif scroll_args.sequence:
            raw_sequence = scroll_args.sequence
        else:
            print("☠️ Catastrophe! You must provide a source: -s, -f, --ncbi, --batch, or --mock!")
            sys.exit(1)

        # Step 1: Purify the text (Removes PDF spaces, numbers, formatting)
        purified_seq = supreme_mage.exorcise_pdf_demons_and_purify(raw_sequence)
        
        # Apply pre-processing spells if invoked
        if scroll_args.mask:
            print("🛡️ Raising the Dust Shield... Masking low complexity regions.")
            purified_seq = supreme_mage.raise_dust_shield_masking(purified_seq)
            
        if scroll_args.trim:
            print("✂️ Wielding the Blacksmith Scissors... Trimming loose edges.")
            purified_seq = supreme_mage.wield_blacksmith_scissors_trimming(purified_seq)
        
        # Regex search if requested
        if scroll_args.regex:
            hunted = supreme_mage.track_aberrations_via_regex(purified_seq, scroll_args.regex)
            print(f"\n🔍 Regex Hunt Results for '{scroll_args.regex}':")
            for h in hunted:
                print(f"   - Found {h['found_beast']} at coord {h['start']}-{h['end']}")
            print("-" * 50)
            
        # Homunculus Cas9 search if requested
        if scroll_args.crispr:
            cas9_targets = supreme_mage.summon_homunculus_cas9_cleavage(purified_seq)
            print(f"\n🧪 Homunculus Cas9 Awakened! Found {len(cas9_targets)} targets:")
            for target in cas9_targets[:15]: # Show up to 15 to avoid flooding the terminal
                print(f"   - Target {target['target']} at {target['start']} (Strand: {target['strand']})")
            if len(cas9_targets) > 15:
                print(f"   ... and {len(cas9_targets) - 15} more hidden in the dark.")
            print("-" * 50)
            
        # Core analysis
        total_gc = supreme_mage.summon_dragon_breath_gc(purified_seq)
        tally = supreme_mage.tally_the_magical_runes(purified_seq)
        treasures = supreme_mage.conjure_search_for_palindromic_treasures(
            purified_seq, scroll_args.minimum, scroll_args.maximum, scroll_args.mismatches
        )
        
        # Save output artifacts
        csv_file = db_dragon.devour_csv_tribute_advanced("cli_alchemical_treasures.csv", treasures)
        svg_file = supreme_mage.paint_runic_constellations_on_parchment(len(purified_seq), treasures, db_dragon.bone_lair)
        
        # Espelho do Ouroboros map integration
        mystical_features = [{"start": t["start_idx"], "end": t["start_idx"] + t["size"]} for t in treasures]
        ouroboros_file = supreme_mage.gaze_into_ouroboros_mirror(len(purified_seq), mystical_features, db_dragon.bone_lair)
        
        # Report integration
        if scroll_args.report:
            spoils = {
                'global_gc': total_gc,
                'tally': tally,
                'palindromes': treasures,
            }
            html_file = db_dragon.inscribe_ultimate_parchment_html("cli_ultimate_parchment.html", spoils)
            print(f"📄 Pergaminho Definitivo (HTML) forged at: {html_file}")
        
        # Output Results to Terminal
        print(f"🐉 Dragon's Breath (Total GC%): {total_gc:.2f}%")
        print(f"🧮 Nucleotide Tally: {tally}")
        
        print(f"\n✨ Glory! {len(treasures)} Palindromic Relics unearthed:")
        print(f"📂 CSV Saved to: {csv_file}")
        print(f"🎨 Linear SVG Map drawn at: {svg_file}")
        if ouroboros_file:
            print(f"🌀 Ouroboros Circular Map drawn at: {ouroboros_file}")
        print("\n" + "=" * 170)
        print(f"{'COORD.':<10} | {'SEQ':<14} | {'MIS':<3} | {'GC%':<5} | {'TM(°C)':<6} | {'ΔG (kcal)':<9} | {'ENTROPY':<7} | {'METHYL.':<12} | {'PATHOGEN SIG.':<15}")
        print("-" * 170)
        for t in treasures:
            print(f"{t['position']:<10} | {t['sequence']:<14} | {t['mismatches']:<3} | {t['gc']:<5} | {t['tm_celsius']:<6} | {t['delta_g']:<9} | {t['entropy']:<7} | {t['methylation']:<12} | {t['pathogen']:<15}")
            
    except Exception as magical_disaster:
        print(f"\n☠️ A catastrophic failure corroded the system: {str(magical_disaster)}")
        sys.exit(1)

if __name__ == "__main__":
    embark_on_terminal_crusade()
