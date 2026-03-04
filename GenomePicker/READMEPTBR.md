#  GenomePicker

**GenomePicker** é uma ferramenta de linha de comando escrita em Python para download e processamento automatizado de dados genômicos do NCBI e estruturas proteicas do UniProt/AlphaFold. Suporta downloads paralelizados, extração de metadados e visualização geográfica de isolados.

---

##  Índice

- [Visão Geral](#visão-geral)
- [Requisitos](#requisitos)
- [Instalação](#instalação)
- [Arquivos de Entrada](#arquivos-de-entrada)
- [Uso](#uso)
- [Argumentos](#argumentos)
- [Funcionalidades](#funcionalidades)
- [Arquivos de Saída](#arquivos-de-saída)
- [Exemplos](#exemplos)
- [Observações e Limitações](#observações-e-limitações)

---

## Visão Geral

O GenomePicker automatiza três fluxos de trabalho principais:

1. **Download de Genomas/Proteomas do NCBI** — Baixa arquivos genômicos (`.fna`, `.faa`, `.gbff`, `.gff`) diretamente do servidor FTP do NCBI, utilizando uma tabela de metadados (TSV/CSV) exportada do [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/).
2. **Extração de Metadados e Mapeamento Geográfico** — Acessa as páginas BioSample do NCBI para extrair metadados epidemiológicos (hospedeiro, doença, fonte de isolamento, localização geográfica) e gera mapas de distribuição geográfica.
3. **Download de Estruturas Proteicas** — Baixa sequências proteicas (FASTA) e estruturas 3D (PDB/CIF) do UniProt/AlphaFold ou do banco de dados NCBI Protein.

---

## Requisitos

| Dependência | Finalidade |
|---|---|
| `Python >= 3.8` | Ambiente de execução |
| `pandas` | Leitura e processamento de tabelas |
| `matplotlib` | Geração de mapas estáticos |
| `mpl_toolkits.basemap` | Renderização do mapa-múndi |
| `plotly` | Geração de mapas interativos |

Instale todas as dependências com:

```bash
pip install pandas matplotlib plotly basemap
```

> **Atenção:** O pacote `basemap` pode exigir dependências adicionais a nível de sistema. Consulte o [guia de instalação do Basemap](https://matplotlib.org/basemap/stable/users/installation.html) para instruções específicas de cada plataforma.

---

## Instalação

Nenhuma instalação é necessária. Basta clonar ou baixar o script e executá-lo diretamente:

```bash
git clone <url-do-repositório>
cd genomepicker
python genome_picker.py --help
```

---

## Arquivos de Entrada

### Tabela do NCBI Genome (TSV/CSV)

A entrada principal é uma tabela de metadados baixada do [portal NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/). O arquivo deve conter as seguintes colunas:

| Coluna | Descrição |
|---|---|
| `Assembly Accession` | Número de acesso único do genoma (usado como índice) |
| `Assembly Name` | Nome do assembly |
| `Organism Name` | Nome completo do organismo (gênero + espécie) |
| `Organism Infraspecific Names Strain` | Identificador de cepa (primário) |
| `Organism Infraspecific Names Isolate` | Identificador de isolado (alternativo) |
| `Organism Infraspecific Names Breed` | Identificador de raça (alternativo) |
| `Organism Infraspecific Names Cultivar` | Identificador de cultivar (alternativo) |
| `Organism Infraspecific Names Ecotype` | Identificador de ecótipo (alternativo) |
| `Assembly BioSample Accession` | ID do BioSample (necessário para `--metadata`) |

> Formatos suportados: `.tsv` (separado por tabulação) e `.csv` (separado por vírgula).

### Listas de IDs de Proteínas (TXT)

Para downloads de proteínas/estruturas, forneça um arquivo de texto simples com um identificador por linha:

```
P12345
Q9UNQ0
A0A000
```

---

## Uso

```bash
python genome_picker.py [OPÇÕES]
```

Se nenhum argumento for fornecido, o menu de ajuda é exibido automaticamente.

---

## Argumentos

### Download de Genomas (requer `--input`)

| Argumento | Descrição |
|---|---|
| `--input ARQUIVO` | Caminho para a tabela TSV ou CSV do NCBI Genome |
| `--genome` | Baixa sequências de genoma em formato FASTA nucleotídico (`.fna`) |
| `--fasta` | Igual a `--genome`, mas salva os arquivos com extensão `.fasta` |
| `--proteome` | Baixa sequências de proteínas (`.faa`) |
| `--genbank` | Baixa arquivos no formato GenBank (`.gbff`) |
| `--gbk` | Igual a `--genbank`, mas salva os arquivos com extensão `.gbk` |
| `--gff` | Baixa arquivos de anotação genômica (`.gff`) |
| `--metadata` | Extrai metadados do BioSample e gera mapas geográficos |
| `--cpu N` | Número de threads paralelas para download (padrão: `1`) |
| `--outdir DIR` | Diretório de saída personalizado (gerado automaticamente se não informado) |

### Download de Proteínas/Estruturas

| Argumento | Descrição |
|---|---|
| `--uniprot_fasta ARQUIVO` | Arquivo de texto com IDs UniProt para baixar sequências FASTA |
| `--uniprot_pdb ARQUIVO` | Arquivo de texto com IDs UniProt para baixar estruturas PDB via AlphaFold |
| `--uniprot_cif ARQUIVO` | Arquivo de texto com IDs UniProt para baixar estruturas CIF via AlphaFold |
| `--ncbi_protein ARQUIVO` | Arquivo de texto com IDs NCBI Protein para baixar sequências FASTA |

---

## Funcionalidades

###  Downloads Paralelos
Os downloads são paralelizados com `multiprocessing.Pool`, reduzindo significativamente o tempo de execução para grandes conjuntos de dados. O número de workers é limitado automaticamente a `total_CPUs - 1` para evitar sobrecarga do sistema.

###  Tentativas Automáticas de Reconexão
Cada download é tentado até **5 vezes**, com intervalo de 3 segundos entre as tentativas, garantindo robustez contra falhas transitórias de rede.

###  Tratamento de Cepas Duplicadas
Se duas entradas compartilharem o mesmo identificador de cepa, um sufixo alfanumérico aleatório de 5 caracteres (`_dup_XXXXX`) é adicionado automaticamente ao nome do arquivo para evitar conflitos.

###  Geração Automática de Script PROKKA
Ao baixar arquivos de genoma (`--genome` ou `--fasta`), um script shell `PROKKA.sh` é gerado automaticamente no diretório de saída. Ele contém comandos prontos do [Prokka](https://github.com/tseemann/prokka) para cada genoma baixado, já preenchidos com gênero, espécie, cepa e locus tag.

###  Mapeamento Geográfico
Com o uso de `--metadata`, o GenomePicker:
1. Acessa cada página BioSample no NCBI para extrair hospedeiro, doença, fonte de isolamento e localização geográfica.
2. Baixa uma tabela de referência de coordenadas (`latlon.csv`) para converter nomes de países em latitude/longitude.
3. Gera dois mapas:
   - **PNG estático** — Mapa-múndi em alta resolução (600 DPI) renderizado com Basemap (`matplotlib`).
   - **HTML interativo** — Mapa com zoom e tooltips renderizado com Plotly Express.

### Fontes para Estruturas Proteicas
- **UniProt FASTA**: obtido de `rest.uniprot.org`
- **AlphaFold PDB/CIF**: obtido de `alphafold.ebi.ac.uk` (modelo v4)
- **NCBI Protein FASTA**: obtido via E-utilities do NCBI (`efetch`)

---

## Arquivos de Saída

### Downloads de Genomas

Os arquivos são salvos no diretório de saída com nomes padronizados:

```
{Inicial_Gênero}{espécie}_{cepa}.{extensão}
```

Exemplo: `Ecoli_O157-H7.fna`

### Metadados

| Arquivo | Descrição |
|---|---|
| `Metadados_DD-MM-AAAA_HH-MM-SS.tsv` | Tabela de metadados com hospedeiro, doença, fonte de isolamento e localização geográfica por genoma |
| `Metadados_contagem_paises_*.tsv` | Contagem de isolados por país, utilizada para geração dos mapas |
| `Mapa_Basemap_*.png` | Mapa-múndi estático em alta resolução (600 DPI) |
| `Mapa_Interativo_Plotly_*.html` | Mapa de distribuição geográfica interativo |

### Prokka

| Arquivo | Descrição |
|---|---|
| `PROKKA.sh` | Script shell com comandos Prokka pré-configurados para cada genoma baixado |

### Proteínas e Estruturas

Os arquivos são salvos como `{ID}.{formato}` (ex.: `P12345.fasta`, `Q9UNQ0.pdb`).

---

## Exemplos

### Baixar sequências de genoma usando 8 threads
```bash
python genome_picker.py --input tabela_assembly.tsv --genome --cpu 8 --outdir MeusGenomas/
```

### Baixar genomas e arquivos de anotação GFF simultaneamente
```bash
python genome_picker.py --input tabela_assembly.tsv --genome --gff --cpu 4
```

### Extrair metadados e gerar mapas geográficos
```bash
python genome_picker.py --input tabela_assembly.tsv --metadata --cpu 4 --outdir MeusMetadados/
```

### Baixar proteomas e arquivos GenBank
```bash
python genome_picker.py --input tabela_assembly.tsv --proteome --genbank --cpu 6
```

### Baixar FASTA de proteínas do UniProt
```bash
python genome_picker.py --uniprot_fasta ids_uniprot.txt --outdir UniProtFASTA/
```

### Baixar estruturas PDB via AlphaFold
```bash
python genome_picker.py --uniprot_pdb ids_uniprot.txt --outdir Estruturas/
```

### Baixar FASTA de proteínas do NCBI Protein
```bash
python genome_picker.py --ncbi_protein ids_ncbi_proteina.txt --outdir ProteínasNCBI/
```

### Execução completa combinada
```bash
python genome_picker.py \
  --input tabela_assembly.tsv \
  --genome --gff --metadata \
  --cpu 8 \
  --outdir AnaliseCompleta/
```

---

## Observações e Limitações

- **Entradas sem identificador de cepa/isolado/raça/cultivar/ecótipo são ignoradas automaticamente.** Apenas registros com pelo menos um nome infraespecífico são processados.
- O flag `--metadata` acessa páginas HTML do NCBI BioSample e pode ser lento para grandes conjuntos de dados. O uso de `--cpu` com múltiplas threads é fortemente recomendado.
- A funcionalidade de mapa geográfico depende da disponibilidade do arquivo de referência de coordenadas (`latlon.csv`) hospedado no GitHub. Países não presentes nesse arquivo não aparecerão no mapa.
- Downloads são ignorados caso o arquivo de saída já exista, permitindo a retomada segura de execuções interrompidas.
- O número de CPUs é automaticamente limitado a `total_CPUs - 1`, independentemente do valor passado em `--cpu`, para preservar a estabilidade do sistema.

---
