import os

SCROLL_OF_COLORS = {
    'A': '#3498db', 
    'T': '#e74c3c', 
    'C': '#2ecc71', 
    'G': '#f39c12', 
    'U': '#9b59b6', 
    'N': '#95a5a6'
}

NINJA_SESSION_FILE = os.path.join(os.path.expanduser("~"), ".bio_ninja_session.json")

NARUTO_HOKAGE_STYLE = """
QMainWindow { 
    background-color: #121212; 
}
QWidget { 
    color: #e0e0e0; 
    font-family: 'Segoe UI', Arial, sans-serif; 
    font-size: 10pt; 
}

/* Top Menu Bar Fix */
QMenuBar {
    background-color: #1a1a1a;
    color: #e0e0e0;
    border-bottom: 1px solid #333333;
}
QMenuBar::item {
    background-color: transparent;
    padding: 6px 12px;
}
QMenuBar::item:selected {
    background-color: #ffd900;
    color: #121212;
    font-weight: bold;
}
QMenu {
    background-color: #1a1a1a;
    color: #e0e0e0;
    border: 1px solid #333333;
}
QMenu::item:selected {
    background-color: #ffd900;
    color: #121212;
}

/* Sidebar Styling */
QListWidget#sidebar_ninja {
    background-color: #1a1a1a;
    border: none;
    border-right: 2px solid #0a0a0a;
    outline: none;
}
QListWidget#sidebar_ninja::item {
    padding: 15px 15px;
    color: #888888;
    border-left: 4px solid transparent;
    border-bottom: 1px solid #222222;
    font-size: 10pt;
}
QListWidget#sidebar_ninja::item:selected {
    background-color: #2a2a2a;
    color: #ffd900;
    border-left: 4px solid #ffd900;
    font-weight: bold;
}
QListWidget#sidebar_ninja::item:hover:!selected { 
    background-color: #222222; 
    color: #ffffff;
}

/* GroupBoxes and Panels */
QGroupBox { 
    border: 1px solid #333333; 
    border-radius: 6px; 
    margin-top: 18px; 
    font-weight: bold; 
    color: #ffd900; 
    background-color: #1a1a1a;
}
QGroupBox::title { 
    subcontrol-origin: margin; 
    left: 15px; 
    padding: 0 10px; 
    background-color: #1a1a1a;
    border-radius: 4px;
}

/* Buttons */
QPushButton { 
    background-color: #ffd900; 
    color: #121212; 
    border-radius: 4px; 
    padding: 8px 14px; 
    font-weight: bold; 
    border: none;
}
QPushButton:hover { background-color: #ffea4d; }
QPushButton:pressed { background-color: #ccad00; }
QPushButton:disabled { background-color: #333333; color: #666666; }

/* Inputs and Texts */
QComboBox, QSpinBox, QLineEdit, QTextEdit { 
    background-color: #121212; 
    color: #ffffff; 
    border: 1px solid #444444; 
    border-radius: 4px; 
    padding: 6px; 
}
QComboBox:drop-down { border: none; }
QComboBox QAbstractItemView {
    background-color: #121212;
    color: #ffffff;
    selection-background-color: #ffd900;
    selection-color: #121212;
}

/* Tables (Neji VCF Viewer) */
QTableWidget {
    background-color: #121212;
    color: #e0e0e0;
    border: 1px solid #333333;
    gridline-color: #444444;
}
QHeaderView::section {
    background-color: #1a1a1a;
    color: #ffd900;
    padding: 4px;
    border: 1px solid #333333;
    font-weight: bold;
}
QTableWidget::item:selected {
    background-color: #ffd900;
    color: #121212;
}

/* Inner Tabs (Plots vs Stats) */
QTabWidget::pane { 
    border: 1px solid #333333; 
    border-radius: 6px; 
    background-color: #1a1a1a;
}
QTabBar::tab { 
    background: #121212; 
    color: #888888; 
    padding: 10px 25px; 
    border-top-left-radius: 6px; 
    border-top-right-radius: 6px; 
    margin-right: 2px;
    border: 1px solid #333333;
    border-bottom: none;
}
QTabBar::tab:selected { 
    background: #1a1a1a; 
    color: #ffd900; 
    font-weight: bold; 
}
QTabBar::tab:hover:!selected {
    background: #222222;
    color: #ffffff;
}

/* Checkboxes & Scrollbars */
QCheckBox { color: #e0e0e0; font-weight: bold; }
QCheckBox::indicator { 
    width: 16px; height: 16px; 
    background-color: #121212; 
    border: 2px solid #444444; 
    border-radius: 3px;
}
QCheckBox::indicator:checked { 
    background-color: #ffd900; 
    border: 2px solid #ffd900;
}
QScrollBar:vertical, QScrollBar:horizontal {
    border: none; background: #121212; width: 14px; height: 14px; border-radius: 7px;
}
QScrollBar::handle:vertical, QScrollBar::handle:horizontal {
    background: #444444; border-radius: 7px;
}
QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover { background: #ffd900; }

/* Docks - For Floating Anbu Terminal */
QDockWidget {
    color: #ffd900;
    font-weight: bold;
}
QDockWidget::title {
    background: #1a1a1a;
    text-align: center;
    border: 1px solid #333333;
    padding: 4px;
}
"""

DOCUMENTATION_HTML = """
<div style="color: #e0e0e0; font-family: 'Segoe UI', Arial, sans-serif; line-height: 1.6; max-width: 900px; margin: auto;">
<h1 style="color: #ffd900; border-bottom: 2px solid #333; padding-bottom: 10px; text-align: center;">📜 Instruction Scroll: LenghtCount (Six Paths Ninja Edition)</h1>
<p style="text-align: justify;">Welcome to <b>LenghtCount</b>. This tool was forged for high-performance genomic and proteomic analysis. Each sidebar tab represents a class of jutsus. Read this scroll carefully to master all features.</p>

<h2 style="color: #ffd900; margin-top: 30px;">👁️ Assembly Statistics (Sharingan)</h2>
<p>This tab is the core for statistical visualization and primary file manipulation.</p>
<ul>
    <li><b>Genomic File Management:</b> Allows loading files (FASTA, FASTQ, GBK, VCF, GFF, BAM). You can load entire directories. The <i>Chibaku Tensei</i> button concatenates all selected files into a single FASTA.</li>
    <li><b>Pain Bansho Tenin (NCBI Downloader):</b> Enter an <i>Accession Number</i> (e.g., NC_045512) for the system to pull the FASTA directly from NCBI servers via Entrez.</li>
    <li><b>Henge Format Shifter:</b> Utility file converter. Transforms BAM to SAM (requires <i>samtools</i> in the system PATH) and FASTQ/GBK to FASTA instantly.</li>
    <li><b>Virtual Trimming Engine:</b> Trims sequences based on size (Min/Max bp) and plots the results.</li>
    <li><b>Graph Renderers:</b> Choose between Histograms, Boxplots, Violins, KDE, and ECDF. You can render HTML reports or export the raw matrix to CSV.</li>
</ul>

<h2 style="color: #ffd900; margin-top: 30px;">🧬 Base Composition (Mokuton)</h2>
<p>Dedicated to thermodynamic analysis, transcription, and exact nucleotide composition.</p>
<ul>
    <li><b>In-Silico Translation & Reverse Complement:</b> Quick action buttons to convert your target sequence into Protein, RNA, or Reverse Complement.</li>
    <li><b>Sharingan Regex Eye:</b> Search for specific patterns or motifs using Regular Expressions (e.g., <code>AT[GC]T</code>). Matches will be highlighted in bright red.</li>
    <li><b>K-mer Swarm Matrix:</b> Counts and plots the absolute frequency of K-mers in your sequence (e.g., Hexamers). Excellent for building genomic signatures.</li>
    <li><b>Inuzuka Dog Taxonomy (Kraken-lite):</b> Attempts to infer the sample's organism based on GC balance against a built-in mini-database.</li>
    <li><b>Shannon Entropy & Gaara Dust Complexity:</b> Measures sequence complexity. DUST masks (transforms to 'N') low complexity (repetitive) areas.</li>
    <li><b>Zabuza Restriction Map & Gel:</b> Scans for popular restriction enzyme cleavages (EcoRI, BamHI, HindIII) and draws a virtual agarose gel (log scale) with the fragments.</li>
    <li><b>Orochimaru Cursed PAIs:</b> Identifies "Pathogenicity Islands" (alien DNA/plasmids) by tracking abrupt GC drops in the middle of the contig.</li>
    <li><b>Telomere Mapper & Rho Odds:</b> Identifies TTAGGG repeats at chromosome ends and calculates the statistical representation of dinucleotide pairs (intrinsic Ti/Tv).</li>
</ul>

<h2 style="color: #ffd900; margin-top: 30px;">⚡ Protein & Primers (Chidori)</h2>
<p>Biophysical tools focused on the final product of translation (Amino Acids) and molecular biology (PCR and CRISPR).</p>
<ul>
    <li><b>Sage Mode ORF Scanner:</b> Finds valid <i>Open Reading Frames</i> (ORFs) (Started by ATG and ended with Stop) in all 6 reading frames.</li>
    <li><b>Gai Eight Gates Properties:</b> Calculates pI (Isoelectric Point), Molecular Weight, Aliphatic Index, GRAVY Score (global hydrophobicity), molar extinction coefficient, and 2D secondary propensity (Alpha-helix and Beta-sheet).</li>
    <li><b>Hydropathy Plot (Kyte-Doolittle):</b> Plots hydrophobicity along the protein chain, identifying transmembrane zones.</li>
    <li><b>Kakashi Copy CRISPR PAM Scan:</b> Scans the provided DNA sequence for PAM sites (NGG) for the SpCas9 nuclease, generates the 20bp sgRNA, its GC%, and cross-references the genome to predict simple Off-Target risk.</li>
    <li><b>Oligo Thermodynamics:</b> Provides the Melting Temperature (Tm) and Gibbs free energy (ΔG) of a primer, evaluating hairpin risk.</li>
</ul>

<h2 style="color: #ffd900; margin-top: 30px;">🌀 FastQ Metrics (Rasengan)</h2>
<p>Area dedicated to primary Quality Control (QC) of HTS sequencing (Illumina, Nanopore, etc).</p>
<ul>
    <li><b>Gaara Sand Trimmer:</b> Trims the ends of FASTQ Reads that have a Phred Score below a specified threshold, saving a new clean FASTQ file.</li>
    <li><b>Phred Assessment Algorithm:</b> Reads the FASTQ in blocks (generator) to save memory and plots the average quality per cycle (base). Ideal for viewing the 3' drop-off.</li>
    <li><b>Kisame Saturation Chart:</b> Plots a curve showing coverage evolution (X) against the number of reads, allowing you to see if the sequencing "saturated" the genome.</li>
    <li><b>Zetsu PCR Duplication Estimator:</b> Analyzes the first 50bp of each read in the FASTQ. High incidence of 100% identical reads at the beginning indicates a PCR duplication artifact.</li>
</ul>

<h2 style="color: #ffd900; margin-top: 30px;">🌌 Advanced Sequence (Rinnegan)</h2>
<p>Comparative Biology Jutsus for alignments and advanced extraction.</p>
<ul>
    <li><b>Global Align & Kabuto Predictor:</b> Performs a Smith-Waterman global alignment. The predictor deduces the Ti/Tv ratio (Transitions vs Transversions) to check if the mutation is biological or an error, and infers if the SNP mutation generated an Amino Acid change (Synonymous/Non-Synonymous).</li>
    <li><b>Deidara Explosive Local BLAST+:</b> Allows calling a local NCBI BLAST+ installation. Insert the binary path (<code>blastn</code> on Linux or <code>C:\\...\\blastn.exe</code> on Windows) via the browse button to run against a Subject.</li>
    <li><b>Hashirama Mokuton UPGMA Tree:</b> Requires a FASTA file containing a Multiple Sequence Alignment (MSA). Draws the phylogenetic dendrogram tree by calculating matrix distances.</li>
    <li><b>Yamato Wood Extractor & Density:</b> Load a GFF and your reference FASTA. Enter a gene name and the tool surgically cuts it from the sequence. The Density button will plot gene density by 10kb windows along the chromosome.</li>
</ul>

<h2 style="color: #ffd900; margin-top: 30px;">🪬 Variant Viewer (Neji Eye)</h2>
<p>Dedicated to VCF (Variant Call Format) files.</p>
<ul>
    <li>Load a VCF or VCF.GZ file. The tool will build a searchable DataFrame allowing you to search by Chromosome, Position, ID, or Mutation type (REF/ALT).</li>
</ul>

<h2 style="color: #ffd900; margin-top: 30px;">⚔️ Anbu Floating HQ (Embedded Shell)</h2>
<p>In the bottom corner of the screen, you have a dockable tab called <b>Anbu Floating HQ</b>. You can "pull it" out of the program and use it as a real auxiliary terminal. It's useful for firing commands in Docker containers (e.g., <code>docker run ...</code>) or running third-party python scripts while you tweak the graphs in parallel.</p>

<hr style="border-color: #333; margin-top: 30px;">
<p style="text-align: center; color: #888;">Developed by <b>VictorSC</b>. The Will of Fire applied to Bioinformatics.</p>
</div>
"""
