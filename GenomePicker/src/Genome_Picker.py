#!/usr/bin/env python3
import argparse
import multiprocessing
import pandas as pd
import re
import random
import string
import gzip
import os
import urllib.request
import urllib.error
import shutil
import time
from datetime import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import plotly.express as px

def extract_bacterial_strain_from_table(csv_table_row):
    if isinstance(csv_table_row.get("Organism Infraspecific Names Strain", None), str):
        formatted_strain = re.sub(r"((?![\.A-z0-9_-]).)", "_", csv_table_row["Organism Infraspecific Names Strain"])
    elif isinstance(csv_table_row.get("Organism Infraspecific Names Isolate", None), str):
        formatted_strain = re.sub(r"((?![\.A-z0-9_-]).)", "_", csv_table_row["Organism Infraspecific Names Isolate"])
    elif isinstance(csv_table_row.get("Organism Infraspecific Names Breed", None), str):
        formatted_strain = re.sub(r"((?![\.A-z0-9_-]).)", "_", csv_table_row["Organism Infraspecific Names Breed"])
    elif isinstance(csv_table_row.get("Organism Infraspecific Names Cultivar", None), str):
        formatted_strain = re.sub(r"((?![\.A-z0-9_-]).)", "_", csv_table_row["Organism Infraspecific Names Cultivar"])
    elif isinstance(csv_table_row.get("Organism Infraspecific Names Ecotype", None), str):
        formatted_strain = re.sub(r"((?![\.A-z0-9_-]).)", "_", csv_table_row["Organism Infraspecific Names Ecotype"])
    else:
        formatted_strain = "Unknown"
    return formatted_strain

def download_ncbi_genome_files(download_parameters):
    https_link, input_file_name, output_file_name, bacteria_genus, bacteria_species, bacteria_strain, locustag_identifier, genome_accession, destination_directory, execute_prokka = download_parameters
    
    if os.path.exists(output_file_name):
        print(f"Arquivo {output_file_name} ja existe. Pulando download.")
        return

    download_attempts = 1
    prokka_commands = []
    
    while download_attempts < 6:
        print(f"Iniciando tentativa de download: {download_attempts} - {genome_accession} {bacteria_genus} {bacteria_species} {bacteria_strain}")
        try:
            http_request = urllib.request.Request(https_link, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
            with urllib.request.urlopen(http_request, timeout=45) as server_response:
                temporary_file = f"temp_{input_file_name}_{''.join(random.choices(string.ascii_letters, k=5))}"
                with open(temporary_file, 'wb') as temp_write_file:
                    shutil.copyfileobj(server_response, temp_write_file)
            
            with gzip.open(temporary_file, 'rb') as compressed_read_file:
                with open(output_file_name, 'wb') as final_write_file:
                    shutil.copyfileobj(compressed_read_file, final_write_file)
            
            print(f"Arquivo {output_file_name} foi baixado com sucesso...")
            os.remove(temporary_file)
            
            prokka_execution_command = f"prokka --addgenes --force --species {bacteria_species} --genus {bacteria_genus} --strain {bacteria_strain} {output_file_name} --prefix {locustag_identifier} --outdir {locustag_identifier} --locustag {locustag_identifier}\n"
            prokka_commands.append(prokka_execution_command)
            break
        except urllib.error.URLError as url_error:
            download_attempts += 1
            time.sleep(3)
            if download_attempts == 6:
                print(f"Nao foi possivel baixar o arquivo {output_file_name} do servidor NCBI. Erro: {url_error}")
        except Exception as general_error:
            download_attempts += 1
            time.sleep(3)
            if download_attempts == 6:
                print(f"Erro inesperado ao baixar {output_file_name}: {general_error}")

    if execute_prokka:
        prokka_script_path = os.path.join(destination_directory, "PROKKA.sh")
        with open(prokka_script_path, "a") as prokka_script_file:
            for individual_command in prokka_commands:
                prokka_script_file.write(individual_command)

def process_ncbi_genomes_table(csv_table_file, file_extension_type, available_cpus, user_output_directory):
    if file_extension_type == "gbff":
        final_extension = "_genomic.gbff.gz"
    elif file_extension_type == "gbk":
        final_extension = "_genomic.gbff.gz"
    elif file_extension_type == "fna":
        final_extension = "_genomic.fna.gz"
    elif file_extension_type == "fasta":
        final_extension = "_genomic.fna.gz"
    elif file_extension_type == "faa":
        final_extension = "_protein.faa.gz"
    elif file_extension_type == "gff":
        final_extension = "_genomic.gff.gz"
        
    file_separator = "," if str(csv_table_file).lower().endswith(".csv") else "\t"
    
    try:
        pandas_table = pd.read_csv(csv_table_file, sep=file_separator, index_col="Assembly Accession", low_memory=False, on_bad_lines='skip')
    except Exception as error_reading:
        print(f"Erro ao ler o arquivo {csv_table_file}. Verifique a integridade do arquivo. Erro: {error_reading}")
        return

    duplicated_strains_list = []
    multiprocessing_arguments_list = []
    
    if user_output_directory:
        working_directory = user_output_directory
    else:
        working_directory = "Genomas_" + datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
        
    for row_index in pandas_table.index:
        try:
            organism_genus = str(pandas_table["Organism Name"][row_index]).split(" ")[0]
            organism_species = str(pandas_table["Organism Name"][row_index]).split(" ")[1]
        except IndexError:
            organism_genus = "Unknown"
            organism_species = "sp"
            
        organism_strain = extract_bacterial_strain_from_table(pandas_table.loc[row_index])
        if organism_strain == "Unknown":
            continue
            
        assembly_name = str(pandas_table["Assembly Name"][row_index]).strip()
        
        if organism_strain not in duplicated_strains_list:
            duplicated_strains_list.append(organism_strain)
        else:
            random_code = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
            organism_strain = f"{organism_strain}_dup_{random_code}"
            duplicated_strains_list.append(organism_strain)
            
        output_file_path = os.path.join(working_directory, f"{organism_genus[0]}{organism_species}_{organism_strain}.{file_extension_type}")
        locustag_identifier = f"{organism_genus[0]}{organism_species[0]}_{organism_strain}"
        
        accession_trios = re.findall("...", str(row_index).replace("_", ""))
        base_https_link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        
        for characters_trio in accession_trios:
            base_https_link = f"{base_https_link}{characters_trio}/"
            
        complete_https_link = f"{base_https_link}{row_index}_{assembly_name}/{row_index}_{assembly_name}{final_extension}"
        server_input_file_name = f"{row_index}_{assembly_name}{final_extension}"
        
        execute_prokka_flag = True if (file_extension_type == "fna" or file_extension_type == "fasta") else False
        function_arguments = (complete_https_link, server_input_file_name, output_file_path, organism_genus, organism_species, organism_strain, locustag_identifier, row_index, working_directory, execute_prokka_flag)
        multiprocessing_arguments_list.append(function_arguments)
        
    existing_prokka_path = os.path.join(working_directory, "PROKKA.sh")
    if os.path.exists(existing_prokka_path) and (file_extension_type == "fna" or file_extension_type == "fasta"):
        os.remove(existing_prokka_path)
        
    real_cpu_count = multiprocessing.cpu_count()
    if available_cpus >= real_cpu_count:
        available_cpus = max(1, real_cpu_count - 1)
        
    print(f"Utilizando {available_cpus} threads para download. Total de genomas: {len(multiprocessing_arguments_list)}\n")
    
    if len(multiprocessing_arguments_list) > 0:
        chunk_size = max(1, len(multiprocessing_arguments_list) // (available_cpus * 4))
        process_pool = multiprocessing.Pool(processes=available_cpus)
        process_pool.map(download_ncbi_genome_files, multiprocessing_arguments_list, chunksize=chunk_size)
        process_pool.close()
        process_pool.join()

def download_single_metadata_file(metadata_parameters):
    https_biosample_link, biosample_accession, assembly_accession, bacteria_genus, bacteria_species, bacteria_strain, genome_accession, output_file_name = metadata_parameters
    print(f"Avaliando metadados de {biosample_accession} - {bacteria_genus} {bacteria_species} {bacteria_strain}")
    
    metadata_attempts = 1
    complete_html_text = ""
    
    while metadata_attempts < 6:
        try:
            http_meta_request = urllib.request.Request(https_biosample_link, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(http_meta_request, timeout=45) as html_response:
                html_lines = html_response.readlines()
                for single_line in html_lines:
                    complete_html_text += single_line.decode('utf-8').strip().replace("\n", "")
            break
        except urllib.error.URLError as meta_connection_error:
            metadata_attempts += 1
            time.sleep(3)
            if metadata_attempts == 6:
                print(f"Falha ao acessar metadados de {biosample_accession}: {meta_connection_error}")
                return
        except Exception as meta_general_error:
            metadata_attempts += 1
            time.sleep(3)
            if metadata_attempts == 6:
                print(f"Erro inesperado nos metadados de {biosample_accession}: {meta_general_error}")
                return

    found_host = "NA"
    found_disease = "NA"
    geographic_location = "NA"
    isolation_source = "NA"

    host_index = complete_html_text.find("<th>host</th><td>")
    if host_index != -1:
        found_host = complete_html_text[host_index+17 : complete_html_text.find("<", host_index+17)].strip().capitalize()

    disease_index = complete_html_text.find("host disease</th><td>")
    if disease_index != -1:
        found_disease = complete_html_text[disease_index+21 : complete_html_text.find("<", disease_index+21)].strip().capitalize()

    geography_index = complete_html_text.find("geo_loc_name=")
    if geography_index != -1:
        geographic_location = complete_html_text[geography_index+13 : complete_html_text.find("&", geography_index+13)].strip()

    source_index = complete_html_text.find(">isolation source</th><td>")
    if source_index != -1:
        isolation_source = complete_html_text[source_index+26 : complete_html_text.find("<", source_index+26)].strip().capitalize()

    metadata_record_line = f"{genome_accession}\t{bacteria_genus} {bacteria_species}\t{bacteria_strain}\t{biosample_accession}\t{found_host}\t{found_disease}\t{isolation_source}\t{geographic_location}\n"
    
    with open(output_file_name, "a", encoding="utf-8") as metadata_write_file:
        metadata_write_file.write(metadata_record_line)

def extract_all_metadata_from_table(csv_table_file, available_cpus, user_output_directory):
    file_separator = "," if str(csv_table_file).lower().endswith(".csv") else "\t"
    
    try:
        pandas_meta_table = pd.read_csv(csv_table_file, sep=file_separator, index_col="Assembly Accession", low_memory=False, on_bad_lines='skip')
    except Exception as error_reading_meta:
        print(f"Erro ao ler o arquivo {csv_table_file}. Erro: {error_reading_meta}")
        return None

    metadata_parameters_list = []
    
    if user_output_directory:
        base_directory = user_output_directory
        if not os.path.exists(base_directory):
            os.makedirs(base_directory)
    else:
        base_directory = "."
        
    final_metadata_file_name = os.path.join(base_directory, f"Metadados_{datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}.tsv")
    
    with open(final_metadata_file_name, "w", encoding="utf-8") as initialization_file:
        initialization_file.write("Genome Accession\tSpecies\tStrain\tBioSample\tHost\tDisease\tIsolation Source\tGeographic Localization\n")
        
    for accession_index in pandas_meta_table.index:
        try:
            org_genus = str(pandas_meta_table["Organism Name"][accession_index]).split(" ")[0]
            org_species = str(pandas_meta_table["Organism Name"][accession_index]).split(" ")[1]
        except IndexError:
            org_genus = "Unknown"
            org_species = "sp"
            
        org_strain = extract_bacterial_strain_from_table(pandas_meta_table.loc[accession_index])
        
        if org_strain == "Unknown":
            continue
            
        meta_assembly_name = str(pandas_meta_table["Assembly Name"][accession_index]).strip()
        meta_biosample_accession = str(pandas_meta_table["Assembly BioSample Accession"][accession_index]).strip()
        ncbi_biosample_link = f"https://www.ncbi.nlm.nih.gov/biosample/{meta_biosample_accession}"
        
        metadata_parameters_list.append((ncbi_biosample_link, meta_biosample_accession, meta_assembly_name, org_genus, org_species, org_strain, accession_index, final_metadata_file_name))
        
    real_cpu_count_meta = multiprocessing.cpu_count()
    if available_cpus >= real_cpu_count_meta:
        available_cpus = max(1, real_cpu_count_meta - 1)
        
    print(f"Utilizando {available_cpus} threads para baixar metadados. Total: {len(metadata_parameters_list)}\n")
    
    if len(metadata_parameters_list) > 0:
        chunk_size_meta = max(1, len(metadata_parameters_list) // (available_cpus * 4))
        meta_process_pool = multiprocessing.Pool(processes=available_cpus)
        meta_process_pool.map(download_single_metadata_file, metadata_parameters_list, chunksize=chunk_size_meta)
        meta_process_pool.close()
        meta_process_pool.join()
    
    return final_metadata_file_name

def generate_geographic_maps_basemap_plotly(tsv_metadata_file):
    if not tsv_metadata_file or not os.path.exists(tsv_metadata_file):
        print("Arquivo de metadados invalido ou inexistente para gerar mapa.")
        return

    coordinates_url = "https://raw.githubusercontent.com/VictorCaricatte/DataBase-for-Bioinformatics/main/Additional/Latlon/latlon.csv"
    coordinates_file_name = "latlon.csv"
    
    try:
        coord_request = urllib.request.Request(coordinates_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(coord_request, timeout=45) as coord_response:
            with open(coordinates_file_name, 'wb') as local_coord_file:
                shutil.copyfileobj(coord_response, local_coord_file)
    except Exception as coord_error:
        print(f"Erro ao baixar tabela de coordenadas: {coord_error}")
        return

    metadata_table = pd.read_csv(tsv_metadata_file, sep="\t", low_memory=False, on_bad_lines='skip')
    
    if "Geographic Localization" not in metadata_table.columns:
        print("Coluna Geographic Localization ausente nos metadados.")
        return
        
    extracted_countries_list = metadata_table["Geographic Localization"].tolist()
    
    try:
        countries_keys_table = pd.read_csv(coordinates_file_name, sep=";", index_col="homecontinent")
    except Exception as coord_read_error:
        print(f"Nao foi possivel abrir o arquivo latlon.csv: {coord_read_error}")
        return
        
    clean_unique_countries = []
    all_countries_list = []
    
    for country_name in extracted_countries_list:
        if isinstance(country_name, str) and country_name != "NA":
            if ":" not in country_name:
                formatted_country = country_name
            else:
                formatted_country = country_name.split(":")[0]
                
            all_countries_list.append(formatted_country)
            if formatted_country not in clean_unique_countries:
                clean_unique_countries.append(formatted_country)
                
    countries_latitudes = []
    countries_longitudes = []
    occurrences_count = []
    successfully_identified_countries = []
    
    for target_country in clean_unique_countries:
        if target_country in countries_keys_table.index.values.tolist():
            countries_latitudes.append(countries_keys_table["homelat"][target_country])
            countries_longitudes.append(countries_keys_table["homelon"][target_country])
            occurrences_count.append(all_countries_list.count(target_country))
            successfully_identified_countries.append(target_country)
        else:
            print(f"Pais nao identificado na base de coordenadas: {target_country}")
            
    if not successfully_identified_countries:
        print("Nenhum pais valido encontrado para plotar no mapa.")
        if os.path.exists(coordinates_file_name):
            os.remove(coordinates_file_name)
        return

    map_data_dictionary = {
        "Country": successfully_identified_countries,
        "Latitude": countries_latitudes,
        "Longitude": countries_longitudes,
        "Quantity": occurrences_count
    }
    map_dataframe = pd.DataFrame.from_dict(map_data_dictionary)
    
    current_directory = os.path.dirname(tsv_metadata_file)
    if not current_directory:
        current_directory = "."
        
    count_file_name = os.path.join(current_directory, f"Metadados_contagem_paises_{datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}.tsv")
    map_dataframe.to_csv(count_file_name, sep="\t", index=False)
    
    try:
        plt.figure(figsize=(20, 15))
        basemap_map = Basemap(llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=80)
        basemap_map.drawmapboundary(fill_color='#A6CAE0', linewidth=0)
        basemap_map.fillcontinents(color='green', alpha=0.3)
        basemap_map.drawcoastlines(linewidth=0.1, color="white")
        map_dataframe['Numeric_Labels'] = pd.factorize(map_dataframe['Country'])[0]
        
        basemap_map.scatter(
            x=map_dataframe['Longitude'].values, 
            y=map_dataframe['Latitude'].values, 
            s=map_dataframe['Quantity'].values * 10, 
            alpha=0.4, 
            c=map_dataframe['Numeric_Labels'].values, 
            cmap="plasma"
        )
        png_image_name = os.path.join(current_directory, f"Mapa_Basemap_{datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}.png")
        plt.savefig(png_image_name, dpi=600, bbox_inches="tight")
        plt.close()
    except Exception as basemap_error:
        print(f"Erro ao gerar mapa com Basemap: {basemap_error}")
    
    try:
        plotly_figure = px.scatter_geo(
            map_dataframe,
            lat="Latitude",
            lon="Longitude",
            size="Quantity",
            hover_name="Country",
            projection="natural earth",
            title="Geographic Distribution of Isolates",
            color="Quantity",
            color_continuous_scale=px.colors.sequential.Plasma
        )
        html_image_name = os.path.join(current_directory, f"Mapa_Interativo_Plotly_{datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}.html")
        plotly_figure.write_html(html_image_name)
    except Exception as plotly_error:
        print(f"Erro ao gerar mapa com Plotly: {plotly_error}")
        
    if os.path.exists(coordinates_file_name):
        os.remove(coordinates_file_name)

def download_protein_structures_ncbi_uniprot(ids_list_file, bank_type, file_format, user_output_directory):
    if user_output_directory:
        structures_directory = user_output_directory
    else:
        structures_directory = f"Proteinas_{bank_type}_{datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}"
        
    if not os.path.exists(structures_directory):
        os.makedirs(structures_directory)
        
    try:
        with open(ids_list_file, "r") as ids_file:
            protein_ids_list = [line.strip() for line in ids_file.readlines() if line.strip()]
    except Exception as ids_read_error:
        print(f"Erro ao ler arquivo de IDs: {ids_read_error}")
        return

    for protein_identifier in protein_ids_list:
        save_path = os.path.join(structures_directory, f"{protein_identifier}.{file_format}")
        
        if os.path.exists(save_path):
            print(f"Arquivo {save_path} ja existe. Pulando.")
            continue
            
        if bank_type == "uniprot":
            if file_format == "fasta":
                protein_download_url = f"https://rest.uniprot.org/uniprotkb/{protein_identifier}.fasta"
            elif file_format == "pdb":
                protein_download_url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_identifier}-F1-model_v4.pdb"
            elif file_format == "cif":
                protein_download_url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_identifier}-F1-model_v4.cif"
        elif bank_type == "ncbi":
            protein_download_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_identifier}&rettype={file_format}&retmode=text"
            
        protein_attempts = 1
        while protein_attempts < 6:
            try:
                http_prot_request = urllib.request.Request(protein_download_url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(http_prot_request, timeout=45) as prot_response:
                    with open(save_path, "wb") as prot_output_file:
                        shutil.copyfileobj(prot_response, prot_output_file)
                print(f"Download concluido: {save_path}")
                break
            except urllib.error.URLError as prot_error:
                protein_attempts += 1
                time.sleep(3)
                if protein_attempts == 6:
                    print(f"Falha ao baixar {protein_identifier} de {bank_type}: {prot_error}")
            except Exception as prot_general_error:
                protein_attempts += 1
                time.sleep(3)
                if protein_attempts == 6:
                    print(f"Erro inesperado ao baixar {protein_identifier}: {prot_general_error}")

if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser(
        description="Script to download genomes, proteomes, metadata from NCBI and protein structures from UniProt/NCBI.",
        epilog="Example: python genome_picker.py --input table.tsv --genome --gff --cpu 4"
    )
    
    argument_parser.add_argument("--input", type=str, help="TSV or CSV file obtained from NCBI Genome.")
    argument_parser.add_argument("--genome", action="store_true", help="Download genome file (fna).")
    argument_parser.add_argument("--proteome", action="store_true", help="Download proteome file (faa).")
    argument_parser.add_argument("--genbank", action="store_true", help="Download genbank file (gbff).")
    argument_parser.add_argument("--gff", action="store_true", help="Download annotation file (gff3).")
    argument_parser.add_argument("--fasta", action="store_true", help="Download genome file and save as .fasta.")
    argument_parser.add_argument("--gbk", action="store_true", help="Download genbank file and save as .gbk.")
    argument_parser.add_argument("--metadata", action="store_true", help="Download metadata and generate geographic maps (Basemap and Plotly).")
    argument_parser.add_argument("--cpu", type=int, default=1, help="Number of threads for concurrent download.")
    argument_parser.add_argument("--outdir", type=str, help="Standardized output directory to save all files.")
    
    argument_parser.add_argument("--uniprot_fasta", type=str, help="Text file with UniProt IDs list to download FASTA.")
    argument_parser.add_argument("--uniprot_pdb", type=str, help="Text file with UniProt IDs list to download PDB via AlphaFold.")
    argument_parser.add_argument("--uniprot_cif", type=str, help="Text file with UniProt IDs list to download CIF via AlphaFold.")
    argument_parser.add_argument("--ncbi_protein", type=str, help="Text file with NCBI Protein IDs list to download FASTA.")

    provided_arguments = argument_parser.parse_args()

    if not any(vars(provided_arguments).values()):
        argument_parser.print_help()
        exit()

    if provided_arguments.input:
        input_file_extension = os.path.splitext(provided_arguments.input)[1].lower()
        if input_file_extension not in [".tsv", ".csv"]:
            print("AVISO: A extensao do arquivo input DEVE ser .tsv ou .csv!")
            exit()
            
        if provided_arguments.genome:
            process_ncbi_genomes_table(provided_arguments.input, "fna", provided_arguments.cpu, provided_arguments.outdir)
        if provided_arguments.proteome:
            process_ncbi_genomes_table(provided_arguments.input, "faa", provided_arguments.cpu, provided_arguments.outdir)
        if provided_arguments.genbank:
            process_ncbi_genomes_table(provided_arguments.input, "gbff", provided_arguments.cpu, provided_arguments.outdir)
        if provided_arguments.gff:
            process_ncbi_genomes_table(provided_arguments.input, "gff", provided_arguments.cpu, provided_arguments.outdir)
        if provided_arguments.fasta:
            process_ncbi_genomes_table(provided_arguments.input, "fasta", provided_arguments.cpu, provided_arguments.outdir)
        if provided_arguments.gbk:
            process_ncbi_genomes_table(provided_arguments.input, "gbk", provided_arguments.cpu, provided_arguments.outdir)
        if provided_arguments.metadata:
            metadata_file_path = extract_all_metadata_from_table(provided_arguments.input, provided_arguments.cpu, provided_arguments.outdir)
            if metadata_file_path:
                generate_geographic_maps_basemap_plotly(metadata_file_path)

    if provided_arguments.uniprot_fasta:
        download_protein_structures_ncbi_uniprot(provided_arguments.uniprot_fasta, "uniprot", "fasta", provided_arguments.outdir)
    if provided_arguments.uniprot_pdb:
        download_protein_structures_ncbi_uniprot(provided_arguments.uniprot_pdb, "uniprot", "pdb", provided_arguments.outdir)
    if provided_arguments.uniprot_cif:
        download_protein_structures_ncbi_uniprot(provided_arguments.uniprot_cif, "uniprot", "cif", provided_arguments.outdir)
    if provided_arguments.ncbi_protein:
        download_protein_structures_ncbi_uniprot(provided_arguments.ncbi_protein, "ncbi", "fasta", provided_arguments.outdir)
