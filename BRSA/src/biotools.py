import subprocess
import os
import sys
import requests
import networkx as nx

# Submits the raw fastq files to the FastQC interrogation chamber before anything else
def summon_fastqc_inquisitor(fastq_paths: list, output_dir: str, sample_id: str, threads: int):
    fastqc_out_dir = os.path.join(output_dir, "fastqc_interrogation_room")
    os.makedirs(fastqc_out_dir, exist_ok=True)
    
    fastqc_cmd = ["fastqc", "-o", fastqc_out_dir, "-t", str(threads)]
    fastqc_cmd.extend(fastq_paths)
    
    try:
        subprocess.run(fastqc_cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] FastQC inquisitor stumbled for {sample_id}:\n{e.stderr}")

# Summons the all-seeing eye of MultiQC to gaze upon all generated reports
def conjure_multiqc_all_seeing_eye(bunker_dir: str):
    try:
        subprocess.run(["multiqc", bunker_dir, "-o", bunker_dir, "-n", "MultiQC_All_Seeing_Eye_Report.html"], 
                       check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] The all-seeing eye of MultiQC was blinded:\n{e.stderr}")

# Chops off the nasty adapter tails and low quality junk using the fastp guillotine
def execute_fastp_guillotine(fastq_paths: list, output_dir: str, sample_id: str, threads: int) -> list:
    clean_reads = []
    html_report = os.path.join(output_dir, f"{sample_id}_fastp_report.html")
    json_report = os.path.join(output_dir, f"{sample_id}_fastp_report.json")
    
    fastp_cmd = ["fastp", "-w", str(threads), "-h", html_report, "-j", json_report]
    
    if len(fastq_paths) == 2:
        out1 = os.path.join(output_dir, f"{sample_id}_clean_R1.fastq.gz")
        out2 = os.path.join(output_dir, f"{sample_id}_clean_R2.fastq.gz")
        fastp_cmd.extend(["-i", fastq_paths[0], "-I", fastq_paths[1], "-o", out1, "-O", out2])
        clean_reads.extend([out1, out2])
    else:
        out1 = os.path.join(output_dir, f"{sample_id}_clean.fastq.gz")
        fastp_cmd.extend(["-i", fastq_paths[0], "-o", out1])
        clean_reads.append(out1)
        
    try:
        subprocess.run(fastp_cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Fastp execution failed for sample {sample_id}:\n{e.stderr}")
        
    return clean_reads

# Throws the reads into the HISAT2 vortex to find their genomic coordinates
def cast_hisat2_alignment_spell(fastq_paths: list, ref_index: str, sam_output: str, threads: int):
    hisat2_cmd = ["hisat2", "-x", ref_index, "-p", str(threads), "-S", sam_output]
    if len(fastq_paths) == 2:
        hisat2_cmd.extend(["-1", fastq_paths[0], "-2", fastq_paths[1]])
    else:
        hisat2_cmd.extend(["-U", fastq_paths[0]])
        
    try:
        subprocess.run(hisat2_cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"HISAT2 alignment error:\n{e.stderr}")

# Forces samtools to do the tedious job of translating SAM to BAM and sorting the mess
def invoke_samtools_sorting_ritual(sam_input: str, bam_output: str, sorted_bam_output: str):
    try:
        with open(bam_output, "wb") as b_out:
            subprocess.run(["samtools", "view", "-bS", sam_input], stdout=b_out, check=True, stderr=subprocess.PIPE)
        subprocess.run(["samtools", "sort", "-o", sorted_bam_output, bam_output], check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Samtools sorting process failed:\n{e.stderr}")

# Harvests the reads falling into the exons using featureCounts like a ruthless farmer
def reap_featurecounts_harvest(bam_input: str, annotation_file: str, counts_output: str, threads: int, is_gff3: bool = False):
    fc_cmd = [
        "featureCounts",
        "-a", annotation_file,
        "-o", counts_output,
        "-T", str(threads)
    ]
    if is_gff3:
        fc_cmd.extend(["-t", "exon", "-g", "ID"])
        
    fc_cmd.append(bam_input)
    
    try:
        subprocess.run(fc_cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"FeatureCounts quantification failed:\n{e.stderr}")

# Unleashes the ultra-fast pseudo-alignment ghost (Kallisto)
def unleash_kallisto_quantification(fastq_paths: list, index_path: str, output_dir: str, sample_id: str, threads: int):
    sample_out_dir = os.path.join(output_dir, f"{sample_id}_kallisto")
    os.makedirs(sample_out_dir, exist_ok=True)
    
    kal_cmd = ["kallisto", "quant", "-i", index_path, "-o", sample_out_dir, "-t", str(threads)]
    
    if len(fastq_paths) == 2:
        kal_cmd.extend([fastq_paths[0], fastq_paths[1]])
    else:
        kal_cmd.extend(["--single", "-l", "200", "-s", "20", fastq_paths[0]])
        
    try:
        subprocess.run(kal_cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Kallisto quantification failed for {sample_id}:\n{e.stderr}")
        
    return os.path.join(sample_out_dir, "abundance.tsv")

# Reaches into the void of STRING DB to pull out physical protein bonds
def summon_string_protein_web(significant_genes_list: list) -> list:
    print("[RITUAL] Summoning STRING DB spirits for protein interactomes...")
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(significant_genes_list[:100]), 
        "species": 9606, 
        "caller_identity": "BRSA_Necromancer"
    }
    try:
        response = requests.post(string_api_url, data=params, timeout=15)
        if response.status_code == 200:
            return response.json()
    except Exception as e:
        print(f"[WARNING] The interactome spirits ignored your call: {e}")
    return []

# Binds the summoned STRING entities into a tangible NetworkX graph
def invoke_interactome_spirits(string_json: list, confidence_score: float = 0.4) -> nx.Graph:
    G = nx.Graph()
    if not string_json: return G
    for edge in string_json:
        if edge.get('score', 0) >= confidence_score:
            G.add_edge(edge['preferredName_A'], edge['preferredName_B'], weight=edge['score'])
    return G

# Carves the exact computational environment into a stone tablet for reproducibility
def scribe_reproducibility_chronicles(output_path: str, parameters: dict):
    try:
        with open(output_path, "w") as f:
            f.write("="*50 + "\n")
            f.write("BRSA BIOINFORMATICS REPRODUCIBILITY SCROLL\n")
            f.write("="*50 + "\n\n")
            
            f.write("--- APPLIED EXPERIMENTAL PARAMETERS ---\n")
            for key, value in parameters.items():
                if key != 'harvested_treasures':
                    f.write(f"{key.upper()}: {value}\n")
            
            f.write("\n--- PYTHON ALCHEMICAL ENVIRONMENT (pip freeze) ---\n")
            pip_freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True)
            f.write(pip_freeze.stdout)
            
            f.write("\n" + "="*50 + "\n")
    except Exception as e:
        print(f"Failed to forge the reproducibility scroll: {e}")
