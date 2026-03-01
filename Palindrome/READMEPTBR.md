# ⚔️ PALINDROME — O Grimório da Feitiçaria Genômica

> *"Na dupla hélice, cada fita espelha sua gêmea — assim como todo herói é definido por sua sombra."*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![Licença](https://img.shields.io/badge/Licença-MIT-green)](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)
[![Status](https://img.shields.io/badge/Status-Ativo-brightgreen)]()
[![Bioinformática](https://img.shields.io/badge/Domínio-Bioinformática-purple)]()

---

## 📜 A Lenda — O Que É Este Artefato?

**Palindrome** é uma ferramenta bioinformática de linha de comando desenvolvida em Python para a análise profunda de sequências palindrômicas de DNA. Sequências palindrômicas — regiões onde o DNA de fita dupla é lido de forma idêntica em ambas as direções 5'→3' — possuem importância biológica fundamental: atuam como sítios de reconhecimento para endonucleases de restrição, sítios de ligação para fatores de transcrição e reguladores estruturais da expressão gênica.

Esta ferramenta oferece um pipeline unificado para entrada e purificação de sequências, detecção de palíndromos, caracterização termodinâmica, mapeamento de sítios CRISPR, identificação de ilhas CpG, mascaramento de repetições, predição de ORFs e geração de relatórios completos — tudo a partir de um único comando.

---

## 🏛️ Os Cinco Tomos Sagrados — Arquitetura do Projeto

| Módulo | Codinome | Função |
|---|---|---|
| `Palindrome.py` | *O Pergaminho do Arauto* | Ponto de entrada da CLI e despachante de argumentos |
| `Logic.py` | *O Grande Alquimista* | Algoritmos centrais de bioinformática e análises |
| `DB.py` | *O Dragão Guardião* | Gerenciamento de saídas: CSV, SVG, relatórios HTML |
| `Config.py` | *O Pergaminho das Verdades Eternas* | Gerenciamento persistente de configurações JSON |
| `Interface.py` | *A Câmara do Oráculo* | Interface gráfica do usuário (GUI) |
| `Dependences.py` | *O Convocador do Arsenal* | Gerenciamento de ferramentas externas (BLAST, DIAMOND) |

---

## ⚗️ Arsenal Mágico — Funcionalidades

### Análises Principais
- **Detecção de Palíndromos** com faixa de comprimento e tolerância a mismatches configuráveis
- **Cálculo de Conteúdo GC** para sequências completas e palíndromos individuais
- **Contagem de Frequência de Nucleotídeos** (A, T, C, G, N)
- **Entropia de Shannon** para avaliação da complexidade de sequências
- **Perfil Termodinâmico**: Temperatura de Fusão (Tm) e Energia Livre de Gibbs (ΔG) via modelo de Vizinhos Mais Próximos

### Molecular e Estrutural
- **Reconhecimento de Enzimas de Restrição**: Compara palíndromos com 15 sítios de restrição canônicos (EcoRI, BamHI, HindIII, NotI e outros)
- **Tradução de Códons**: Tradução in silico de proteínas nos três quadros de leitura
- **Detecção de Motivos de Metilação**: Identificação de sítios Dam (GATC), Dcm (CCWGG) e CpG
- **Detecção de ORFs**: Predição de Fases de Leitura Aberta nos seis quadros (ATG → códon de parada) com tamanho mínimo configurável
- **Mapeamento de Ilhas CpG**: Análise por janela deslizante com critérios GC% ≥ 50% e razão O/E ≥ 0,6
- **Detecção de Repetições em Tandem**: Identifica ecos de motivos sequenciais (ex.: ATATAT...)
- **Identificação de Shine-Dalgarno**: Predição do sítio de ligação ribossomal a montante de portais ATG
- **Predição de Terminadores de Transcrição**: Detecção de hairpins intrínsecos ricos em GC + cauda poli-T

### Utilitários Avançados
- **Mapeamento de Alvos CRISPR/Cas9**: Identificação de sítios PAM NGG em ambas as fitas
- **Busca por Padrões Regex**: Caça de motivos customizados dentro das sequências
- **Detecção de Blocos de Sintenia**: Identifica blocos k-mer conservados entre duas sequências
- **Cálculo da Razão Ti/Tv**: Razão transição/transversão para análise evolutiva
- **Árvore Filogenética UPGMA**: Geração de árvore simples em formato Newick a partir de múltiplas sequências
- **Design de Primers PCR**: Design heurístico de primers flanqueando cada palíndromo (GC%, Tm, verificação de hairpin)
- **Varredura de Assinaturas Patogênicas**: Compara sequências com banco de motivos de virulência
- **Cálculo de Peso Molecular**: Peso exato do DNA de fita dupla em Daltons
- **Probabilidade de Ocorrência Aleatória**: Expectativa estatística de um palíndromo ocorrer ao acaso

### Pré-Processamento
- **Exorcismo de Demônios de PDF**: Remove artefatos de formatação de sequências copiadas de PDFs
- **Mascaramento de Baixa Complexidade** (`--mask`): Transmuta regiões de repetição em tandem para 'N'
- **Trimagem Terminal** (`--trim`): Remove bases 'N' soltas nas extremidades e adaptadores

### Formatos de Saída
- **Relatório CSV**: Saída tabular completa com todos os palíndromos e atributos
- **Mapa SVG Linear**: Visualização vetorial das posições dos palíndromos ao longo da sequência
- **Mapa SVG Circular (Ouroboros)**: Mapa genômico em coordenadas polares destacando os loci palindrômicos
- **Relatório HTML Interativo**: Relatório completo de análise genômica com tema escuro

### Fontes de Entrada
- String de sequência bruta (`-s`)
- Arquivos FASTA / multi-FASTA / FASTQ / GenBank (`-f`)
- Busca remota via NCBI E-utilities por ID de acesso (`--ncbi`)
- Geração de sequência mock aleatória para testes (`--mock`)
- Processamento em lote de diretório com arquivos FASTA (`--batch`)

---

## 🗡️ Instalação — Forjando a Arma

### Pré-Requisitos

```bash
Python >= 3.8
```

### Clonar o Repositório

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Palindrome
```

### Instalar Dependências

```bash
pip install matplotlib numpy biopython
```

> **Nota:** `matplotlib` e `numpy` são necessários para geração dos mapas SVG e circular. `biopython` é exigido apenas para leitura de arquivos GenBank (`.gbk`). Todas as demais funcionalidades rodam com a biblioteca padrão do Python.

---

## 📖 Uso — Lançando os Feitiços

### Varredura Básica de Palíndromos

```bash
python Palindrome.py -s "ATGCGCAT"
```

### A Partir de um Arquivo FASTA

```bash
python Palindrome.py -f genoma.fasta
```

### Arquivo Multi-FASTA

```bash
python Palindrome.py -f multi_registros.fasta
```

### A Partir do NCBI por Acesso

```bash
python Palindrome.py --ncbi NC_000913.3
```

### Gerar Sequência Mock Aleatória

```bash
python Palindrome.py --mock 5000
```

### Execução Completa com Todas as Funcionalidades

```bash
python Palindrome.py -f genoma.fasta -min 4 -max 12 -mis 1 --crispr --mask --trim --report
```

### Processamento em Lote de Diretório

```bash
python Palindrome.py --batch ./minha_pasta_fasta/
```

### Limpar Todos os Artefatos de Saída

```bash
python Palindrome.py --clean
```

---

## 🧙 Referência de Argumentos da CLI

| Argumento | Tipo | Padrão | Descrição |
|---|---|---|---|
| `-s`, `--sequence` | `str` | `None` | Entrada de sequência DNA bruta |
| `-f`, `--file` | `str` | `None` | Caminho para arquivo genômico (FASTA, FASTQ, GenBank) |
| `-n`, `--ncbi` | `str` | `None` | ID de acesso NCBI para busca remota |
| `--mock` | `int` | `None` | Comprimento da sequência DNA aleatória para testes |
| `--batch` | `str` | `None` | Caminho do diretório para processamento em lote de FASTAs |
| `-min`, `--minimum` | `int` | `4` | Comprimento mínimo do palíndromo |
| `-max`, `--maximum` | `int` | `12` | Comprimento máximo do palíndromo |
| `-mis`, `--mismatches` | `int` | `0` | Máximo de mismatches tolerados no palíndromo |
| `--regex` | `str` | `None` | Padrão regex para buscar dentro da sequência |
| `--crispr` | flag | `False` | Ativar detecção de sítios PAM NGG do CRISPR/Cas9 |
| `--mask` | flag | `False` | Mascarar regiões de baixa complexidade com 'N' |
| `--trim` | flag | `False` | Remover bases 'N' terminais das extremidades da sequência |
| `--report` | flag | `False` | Gerar relatório HTML interativo completo |
| `--clean` | flag | `False` | Remover todos os artefatos SVG, CSV e HTML de saída |

---

## 📊 Saída — As Colunas Sagradas

Cada palíndromo detectado é reportado com os seguintes atributos:

| Coluna | Descrição |
|---|---|
| `COORD.` | Coordenadas genômicas (início–fim) |
| `SEQ` | Sequência palindrômica |
| `MIS` | Número de mismatches |
| `GC%` | Percentual de conteúdo GC |
| `TM(°C)` | Temperatura de fusão |
| `ΔG (kcal)` | Energia livre de Gibbs do dúplex |
| `ENTROPY` | Entropia de Shannon (complexidade da sequência) |
| `METHYL.` | Motivos de metilação detectados |
| `PATHOGEN SIG.` | Correspondências com motivos de virulência |

---

## 🔬 Métodos Científicos e Referências

Os algoritmos implementados nesta ferramenta são fundamentados em métodos estabelecidos de bioinformática e biologia molecular:

### Detecção de Palíndromos
Comparação por complemento reverso adaptada dos princípios fundamentais de análise de sequências.
> Watson, J.D., Crick, F.H.C. (1953). Molecular structure of nucleic acids. *Nature*, 171, 737–738. https://doi.org/10.1038/171737a0

### Cálculo do Conteúdo GC
> Marmur, J., Doty, P. (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. *Journal of Molecular Biology*, 5(1), 109–118. https://doi.org/10.1016/S0022-2836(62)80066-7

### Temperatura de Fusão — Regra de Wallace (sequências curtas)
> Wallace, R.B., Shaffer, J., Murphy, R.F., et al. (1979). Hybridization of synthetic oligodeoxyribonucleotides to φX174 DNA: the effect of single base pair mismatch. *Nucleic Acids Research*, 6(11), 3543–3557. https://doi.org/10.1093/nar/6.11.3543

### Temperatura de Fusão — Fórmula para Sequências Longas
> Baldino, F. Jr., Chesselet, M.F., Lewis, M.E. (1989). High-resolution in situ hybridization histochemistry. *Methods in Enzymology*, 168, 761–777. https://doi.org/10.1016/0076-6879(89)68055-3

### Termodinâmica — Modelo ΔG de Vizinhos Mais Próximos
> SantaLucia, J. Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the National Academy of Sciences*, 95(4), 1460–1465. https://doi.org/10.1073/pnas.95.4.1460

### Entropia de Shannon para Complexidade de Sequências
> Shannon, C.E. (1948). A mathematical theory of communication. *The Bell System Technical Journal*, 27(3), 379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x

### Identificação de Ilhas CpG
Algoritmo de janela deslizante com GC% ≥ 50% e razão CpG observada/esperada ≥ 0,6.
> Gardiner-Garden, M., Frommer, M. (1987). CpG islands in vertebrate genomes. *Journal of Molecular Biology*, 196(2), 261–282. https://doi.org/10.1016/0022-2836(87)90689-9

### Sítios de Reconhecimento de Enzimas de Restrição
> Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2015). REBASE — a database for DNA restriction and modification: enzymes, genes and genomes. *Nucleic Acids Research*, 43(D1), D298–D299. https://doi.org/10.1093/nar/gku1046

### Tabela de Tradução de Códons
> Nomenclature Committee of the International Union of Biochemistry (NC-IUB). (1985). Nomenclature for incompletely specified bases in nucleic acid sequences. *European Journal of Biochemistry*, 150(1), 1–5. https://doi.org/10.1111/j.1432-1033.1985.tb08977.x

### Motivos de Metilação do DNA (Dam, Dcm, CpG)
> Casadesús, J., Low, D.A. (2013). Programmed heterogeneity: epigenetic mechanisms in bacteria. *Journal of Biological Chemistry*, 288(20), 13929–13935. https://doi.org/10.1074/jbc.R113.472274

### Detecção de Sítios PAM CRISPR/Cas9 (NGG)
> Jinek, M., Chylinski, K., Fonfara, I., et al. (2012). A programmable dual-RNA–guided DNA endonuclease in adaptive bacterial immunity. *Science*, 337(6096), 816–821. https://doi.org/10.1126/science.1225829

### Predição de Fases de Leitura Aberta (ORF)
> Zhu, H., Hu, G.Q., Yang, Y.F., Wang, J., She, Z.S. (2007). MED: a new non-supervised gene prediction algorithm for bacterial and archaeal genomes. *BMC Bioinformatics*, 8, 97. https://doi.org/10.1186/1471-2105-8-97

### Reconhecimento da Sequência Shine-Dalgarno
> Shine, J., Dalgarno, L. (1974). The 3'-terminal sequence of Escherichia coli 16S ribosomal RNA: complementarity to nonsense triplets and ribosome binding sites. *Proceedings of the National Academy of Sciences*, 71(4), 1342–1346. https://doi.org/10.1073/pnas.71.4.1342

### Predição de Terminadores Intrínsecos de Transcrição
> d'Aubenton Carafa, Y., Brody, E., Thermes, C. (1990). Prediction of rho-independent Escherichia coli transcription terminators. *Journal of Molecular Biology*, 216(4), 835–858. https://doi.org/10.1016/S0022-2836(99)80005-9

### Razão Transição/Transversão (Ti/Tv)
> Freese, E. (1962). On the evolution of the base composition of DNA. *Journal of Theoretical Biology*, 3(1), 82–101. https://doi.org/10.1016/0022-5193(62)90071-8

### Construção de Árvore Filogenética UPGMA
> Sokal, R.R., Michener, C.D. (1958). A statistical method for evaluating systematic relationships. *University of Kansas Science Bulletin*, 38, 1409–1438.

### Detecção de Repetições em Tandem
> Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research*, 27(2), 573–580. https://doi.org/10.1093/nar/27.2.573

### Cálculo de Peso Molecular
> Thermo Fisher Scientific. (2020). Calculating the molecular weight of double-stranded DNA. *Technical Reference*. Disponível em: https://www.thermofisher.com

### Alinhamento Local Smith-Waterman
> Smith, T.F., Waterman, M.S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195–197. https://doi.org/10.1016/0022-2836(81)90087-5

### Busca de Sequências via NCBI E-Utilities
> Sayers, E.W., Beck, J., Bolton, E.E., et al. (2022). Database resources of the National Center for Biotechnology Information. *Nucleic Acids Research*, 50(D1), D20–D26. https://doi.org/10.1093/nar/gkab1112

---

## 📂 Estrutura de Arquivos de Saída

Após a execução, os seguintes arquivos são gerados no diretório `local_db_dungeon/`:

```
local_db_dungeon/
├── cli_alchemical_treasures.csv    # Tabela completa de palíndromos (CSV)
├── constellation_map.svg           # Mapa linear de DNA (SVG)
├── ouroboros_mirror.svg            # Mapa circular do genoma (SVG)
└── cli_ultimate_parchment.html     # Relatório HTML interativo (com --report)

batch_results_dungeon/
└── hellfire_batch_summary.csv      # Sumário do processamento em lote (com --batch)
```

---

## 🔧 Integração com Ferramentas Externas

O Palindrome suporta integração opcional com:

- **BLAST** (`makeblastdb`): Para construção de bancos de dados locais de nucleotídeos
- **DIAMOND** (`diamond makedb`): Para construção de bibliotecas de comparação de proteínas/DNA

Configure os caminhos via GUI (`Interface.py`) ou editando diretamente o arquivo `secret_configurations.json`.

---

## 🧪 Teste Rápido

```bash
# Gera uma sequência mock de 2000 pb e executa análise completa com mapeamento CRISPR e relatório HTML
python Palindrome.py --mock 2000 -min 4 -max 12 --crispr --report
```

---

## 📜 Licença

Este projeto é distribuído sob a **Licença MIT**.
Consulte o texto completo da licença em: [LICENSE](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)

---

## 🧙‍♂️ Autor

**Victor Caricatte**
GitHub: [@VictorCaricatte](https://github.com/VictorCaricatte)
Repositório: [BasicBioinfo/Palindrome](https://github.com/VictorCaricatte/BasicBioinfo/tree/main/Palindrome)

---

> *"A sequência fala. O palíndromo sussurra de volta. O alquimista escuta."*
