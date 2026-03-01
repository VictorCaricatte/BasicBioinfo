<div align="center">

# 🍃 LenghtCount — Six Paths Ninja Suite 🍃

> *"A Vontade do Fogo aplicada à Bioinformática."*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python)](https://www.python.org/)
[![BioPython](https://img.shields.io/badge/BioPython-1.80%2B-green?style=for-the-badge)](https://biopython.org/)
[![Platform](https://img.shields.io/badge/Plataforma-Linux%20%7C%20Windows%20%7C%20macOS-lightgrey?style=for-the-badge)]()
[![Repository](https://img.shields.io/badge/GitHub-LengthCounter-orange?style=for-the-badge&logo=github)](https://github.com/VictorCaricatte/BasicBioinfo/tree/main/LengthCounter)

**Uma suíte completa de análise genômica e proteômica com interface gráfica (PyQt6) e uma poderosa CLI — forjada na Vila Oculta da Bioinformática.**

</div>

---

## 📜 Sumário

1. [Visão Geral](#-visão-geral)
2. [Arquitetura](#-arquitetura)
3. [Requisitos & Instalação](#-requisitos--instalação)
4. [Início Rápido](#-início-rápido)
5. [Módulo I — Estatísticas de Assembly (Sharingan)](#-módulo-i--estatísticas-de-assembly-sharingan)
6. [Módulo II — Composição de Bases (Mokuton)](#-módulo-ii--composição-de-bases-mokuton)
7. [Módulo III — Proteínas & Primers (Chidori)](#-módulo-iii--proteínas--primers-chidori)
8. [Módulo IV — Métricas FastQ (Rasengan)](#-módulo-iv--métricas-fastq-rasengan)
9. [Módulo V — Análise Avançada de Sequências (Rinnegan)](#-módulo-v--análise-avançada-de-sequências-rinnegan)
10. [Módulo VI — Visualizador de Variantes (Olho do Neji)](#-módulo-vi--visualizador-de-variantes-olho-do-neji)
11. [Módulo VII — QG Flutuante Anbu (Terminal Integrado)](#-módulo-vii--qg-flutuante-anbu-terminal-integrado)
12. [Referência da CLI](#-referência-da-cli)
13. [Metodologias Científicas & Referências](#-metodologias-científicas--referências)
14. [Formatos de Arquivo Suportados](#-formatos-de-arquivo-suportados)
15. [Contribuindo](#-contribuindo)

---

## 🌀 Visão Geral

**LenghtCount** é uma plataforma bioinformática de código aberto e multimodal, projetada para caracterização genômica e proteômica abrangente. Ela integra análise estatística de assembly, perfilamento de composição nucleotídica, caracterização físico-química de proteínas, controle de qualidade de sequenciamento, filogenética, chamada de variantes e predição de alvos CRISPR — tudo em uma única ferramenta coesa.

O projeto é estruturado em quatro arquivos principais:

| Arquivo | Função |
|---|---|
| `Lenght.py` | Ponto de entrada — parsing de argumentos CLI e orquestração da execução |
| `logic.py` | Motor computacional central — todos os algoritmos biológicos |
| `Interface.py` | Camada de interface gráfica em PyQt6 |
| `config.py` | Constantes de estilização e documentação HTML |

A ferramenta é completamente **autossuficiente**: não requer serviços web externos além da API Entrez do NCBI (para o downloader de accession) e de uma instalação local opcional do BLAST+.

---

## 🏗 Arquitetura

```
LenghtCount/
│
├── Lenght.py             # Ponto de entrada CLI (argparse)
├── logic.py              # Motor central de bioinformática
├── Interface.py          # Interface gráfica PyQt6
├── config.py             # Temas, cores, strings de documentação
│
└── Classes Analisadoras (logic.py):
    ├── UchihaSharinganAnalyzer     → Estatísticas de Assembly (FASTA/GBK)
    ├── SenjuMokutonAnalyzer        → Composição de Bases & Motivos
    ├── UzumakiChakraFastqAnalyzer  → Métricas de QC para FASTQ
    ├── RinneganAdvancedBioTools    → Alinhamento, ORF, Filogenética, CRISPR
    ├── KatsuyuSlugFormatParser     → Parsing de VCF / GFF / SAM
    ├── NejiVariantEyeParser        → Extração de Variantes de VCF
    ├── YamatoWoodFeatureExtractor  → Extração de genes guiada por GFF
    ├── PainBanshoTeninNCBIFetcher  → Recuperação de sequências via NCBI Entrez
    ├── GaaraSandTrimmer            → Trimagem de qualidade de FASTQ
    ├── HengeFormatShifter          → Conversão de formatos (BAM/SAM/FASTQ/GBK)
    └── AnbuBlackOpsExternalTools   → Wrapper de subprocess para BLAST+
```

---

## ⚙️ Requisitos & Instalação

### Dependências do Sistema

- Python **3.8** ou superior
- `samtools` (opcional — necessário apenas para conversão BAM → SAM)
- NCBI BLAST+ (opcional — necessário apenas para alinhamento BLAST local)

### Dependências Python

```bash
pip install biopython numpy pandas matplotlib seaborn scipy PyQt6
```

Lista completa de dependências:

| Pacote | Finalidade |
|---|---|
| `biopython >= 1.80` | I/O de sequências, alinhamento, análise de restrição, filogenética, Entrez |
| `numpy` | Operações vetorizadas (GC skew, detecção de CpG, dot plots) |
| `pandas` | Manipulação tabular de VCF, distribuições de comprimento, cobertura |
| `matplotlib` | Renderização de figuras (todas as visualizações) |
| `seaborn` | Plots estatísticos (boxplot, violin, KDE, ECDF, histplot) |
| `scipy` | Funções estatísticas auxiliares |
| `PyQt6` | Framework de interface gráfica |

### Instalação

```bash
# Clonar o repositório
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/LengthCounter

# Instalar dependências Python
pip install -r requirements.txt

# Iniciar a interface gráfica
python Lenght.py --gui

# Ou usar diretamente pela linha de comando
python Lenght.py genome.fasta --histogram --stats relatorio.txt
```

---

## 🚀 Início Rápido

```bash
# Estatísticas básicas de assembly
python Lenght.py assembly.fasta --stats saida.txt

# Gerar histograma e boxplot
python Lenght.py assembly.fasta --histogram --boxplot --out-prefix minha_analise

# Baixar sequência do NCBI e analisar
python Lenght.py --ncbi NC_045512

# Perfilamento de qualidade FASTQ
python Lenght.py leituras.fastq.gz --phred --rarefaction

# Análise físico-química de proteínas
python Lenght.py --protein-file proteinas.faa

# Iniciar a interface gráfica
python Lenght.py --gui
```

---

## 🔴 Módulo I — Estatísticas de Assembly (Sharingan)

> *"O Sharingan enxerga tudo. Este módulo também."*

Este módulo fornece uma estrutura estatística abrangente para avaliação da qualidade de montagens genômicas, por meio de múltiplas visualizações de distribuição e métricas-resumo.

### Métricas de Qualidade de Assembly

O motor estatístico central (`shikamaru_shadow_stats`) computa as seguintes métricas sobre um conjunto de comprimentos de contigs/scaffolds:

| Métrica | Descrição |
|---|---|
| **N50** | Comprimento tal que contigs ≥ N50 cobrem ao menos 50% do comprimento total do assembly |
| **N90** | Comprimento tal que contigs ≥ N90 cobrem ao menos 90% do comprimento total do assembly |
| **L50 / L90** | Número de contigs necessários para atingir N50/N90, respectivamente |
| **NG50** | N50 referenciado ao genoma esperado — calculado contra um tamanho *G* fornecido |
| **LG50** | L50 referenciado ao genoma esperado |
| **Mín / Máx / Média / Mediana / DP** | Estatísticas descritivas da distribuição |
| **Q1 / Q3 (IQR)** | Primeiro e terceiro quartis |

**Largura dos bins** no histograma é determinada pela regra de Freedman-Diaconis:

```
largura_bin = 2 × IQR × n^(-1/3)
```

Esta regra é preferida à regra de Sturges para distribuições genômicas grandes e assimétricas (Freedman & Diaconis, 1981).

### Visualizações

| Gráfico | Renderizador | Caso de Uso |
|---|---|---|
| **Histograma** | Binagem de Freedman-Diaconis | Distribuição geral de comprimentos por dataset |
| **Boxplot** | Baseado em IQR, outliers suprimidos | Comparação entre datasets |
| **Violin Plot** | KDE espelhado | Forma completa da distribuição |
| **KDE (Estimativa de Densidade por Kernel)** | Kernel gaussiano | Função de densidade suavizada por dataset |
| **ECDF (Função de Distribuição Empírica)** | Função degrau, `seaborn.ecdfplot` | Análise de proporção cumulativa |

### Operações com Arquivos

- **Carregamento de múltiplos arquivos e diretórios**: Todos os formatos genômicos comuns (FASTA, FASTQ, GBK, VCF, GFF, BAM) são auto-detectados pelo sniffer de formato (`neji_byakugan_format_sniffer`).
- **Concatenação (Chibaku Tensei)**: Mescla múltiplos arquivos FASTA selecionados em um único arquivo de saída.
- **Trimagem Virtual**: Filtra sequências por uma janela Min/Max bp configurável antes da plotagem.
- **Exportação CSV**: Exporta os dados brutos de comprimento por contig para análise posterior em R ou Excel.
- **Relatório HTML**: Incorpora todas as figuras geradas como PNGs codificados em base64 em um arquivo HTML autônomo.

### Downloader NCBI Entrez

O `PainBanshoTeninNCBIFetcher` encapsula o `Entrez.efetch` do BioPython para recuperar sequências em formato FASTA diretamente do Nucleotídeo do NCBI através do número de acesso (ex.: `NC_045512` para o SARS-CoV-2).

```bash
python Lenght.py --ncbi NC_045512 --out-prefix sarscov2
```

---

## 🟢 Módulo II — Composição de Bases (Mokuton)

> *"De um único nucleotídeo, uma floresta cresce."*

Este módulo foca nas propriedades composicionais e termodinâmicas intrínsecas das sequências nucleotídicas.

### Contagem de Nucleotídeos & Conteúdo GC

Frequências brutas das bases (A, T, C, G, U, N) são contadas e reportadas com composição percentual. O conteúdo GC é utilizado como principal característica discriminatória ao longo deste módulo.

### Transcrição, Tradução & Complemento Reverso

Implementado via `Bio.Seq`:

- **Transcrição**: DNA → RNA (substituição T→U)
- **Tradução**: DNA/RNA → Proteína, usando o código genético padrão do NCBI (`translate(to_stop=False)`)
- **Complemento Reverso**: Complementaridade Watson-Crick com reversão de fita

### Busca de Motivos por Expressão Regular

O **Olho Sharingan Regex** permite matching completo com padrões `re` do Python sobre a sequência carregada. Códigos de ambiguidade IUPAC e padrões complexos (ex.: `AT[GC]T`, `CG{3,}`) são suportados. Os matches são destacados posicionalmente na interface.

### Análise de Frequência de K-mers

A contagem de K-mers utiliza um gerador de janela deslizante com `collections.Counter` do Python, reportando os 15 K-mers mais frequentes. Esta abordagem é equivalente ao método de assinatura genômica baseada em frequência descrito por Karlin et al. (1997), usado para perfilamento taxonômico e detecção de transferência horizontal de genes.

```bash
python Lenght.py genome.fasta --kmer 6
```

### Entropia de Shannon & Mascaramento de Complexidade

A **Entropia de Shannon** é calculada sobre uma janela deslizante de tamanho configurável:

```
H(w) = - Σ p_i × log2(p_i)
```

onde *p_i* é a frequência da base *i* na janela *w*. Alta entropia indica sequência diversa e rica em informação; baixa entropia sinaliza regiões repetitivas ou de baixa complexidade (Shannon, 1948).

**Mascaramento DUST** (`gaara_dust_complexity_shield`): Janelas com H < 1,0 bit são substituídas por `N`, análogo ao algoritmo DUST implementado no BLAST+ (Morgulis et al., 2006).

```bash
python Lenght.py genome.fasta --entropy 50
```

### Detecção de Ilhas CpG (Vetorizada)

Ilhas CpG são identificadas usando os critérios de Gardiner-Garden & Frommer (1987):

- Razão Observado/Esperado de CpG ≥ 0,60
- Conteúdo GC ≥ 50% (implícito pelo enriquecimento na janela deslizante)

O cálculo é **completamente vetorizado com NumPy** via `np.convolve` sobre máscaras binárias de bases, fornecendo desempenho O(N) em sequências de escala cromossômica.

```bash
python Lenght.py genome.fasta --cpg
```

### GC Skew

O GC skew é calculado com uma janela deslizante de 1.000 bp (padrão):

```
GC Skew = (G - C) / (G + C)
```

Mudanças de sinal no GC skew são indicadores clássicos da **origem de replicação (oriC)** e do **término de replicação (ter)** em cromossomos bacterianos (Grigoriev, 1998).

```bash
python Lenght.py genome.fasta --gc-skew
```

### Mapeamento de Enzimas de Restrição & Gel Virtual

O `zabuza_executioner_restriction_map` utiliza `Bio.Restriction` para identificar sítios de clivagem das enzimas EcoRI, BamHI, HindIII, XhoI e NotI. Os tamanhos dos fragmentos são plotados em um **strip plot em escala logarítmica** para simular uma imagem de gel de agarose.

### Detecção de Ilhas de Patogenicidade (PAI)

O `orochimaru_cursed_islands_pai` detecta janelas genômicas com conteúdo GC significativamente abaixo da média global do genoma (queda > 5%). Esta heurística imita o método de desvio de GC para identificar DNA adquirido horizontalmente, como ilhas de patogenicidade (Hacker & Kaper, 2000).

### Mapeamento de Telômeros

Repetições do hexâmero TTAGGG/CCCTAA são buscadas dentro dos primeiros e últimos 5.000 bp de cada contig usando regex que requer ≥ 3 cópias em tandem consecutivas — consistente com a estrutura da unidade de repetição telomérica canônica em vertebrados (Meyne et al., 1989).

### Razão Rho de Dinucleotídeos

A estatística Rho (ρ) mede a sobre- ou sub-representação de dinucleotídeos:

```
ρ(XY) = f(XY) / (f(X) × f(Y))
```

Valores > 1,2 indicam sobre-representação; valores < 0,8 indicam supressão. A supressão de CpG (ρ(CpG) << 1) é uma característica marcante de genomas de vertebrados (Karlin & Ladunga, 1994).

### Inferência Taxonômica Kraken-lite

Um classificador heurístico de organismos baseado em conteúdo GC, comparado a um pequeno banco de referência. Destinado à demonstração educacional do princípio de perfilamento por GC; não substitui o KRAKEN2 (Wood & Salzberg, 2014) em fluxos de trabalho de produção.

---

## 🟡 Módulo III — Proteínas & Primers (Chidori)

> *"Mil pássaros cantando verdades físico-químicas."*

### Localizador de ORFs — Todos os 6 Quadros de Leitura

O `jiraiya_sage_mode_orf_finder` varre ambas as fitas nos três quadros de leitura (6 no total) em busca de ORFs iniciados por ATG e terminados por códon de parada, utilizando `Bio.Seq.translate(to_stop=True)`. Apenas ORFs ≥ 100 aminoácidos (padrão) são reportados, ordenados por comprimento decrescente.

```bash
python Lenght.py genome.fasta --orf
```

### Propriedades Físico-Químicas de Proteínas (Oito Portões do Gai)

Implementado via `Bio.SeqUtils.ProtParam.ProteinAnalysis`:

| Propriedade | Algoritmo / Referência |
|---|---|
| **Peso Molecular** | Soma das massas dos resíduos + água (Gasteiger et al., 2005) |
| **Ponto Isoelétrico (pI)** | Busca iterativa de pH pelas equações de Henderson-Hasselbalch |
| **Índice de Instabilidade** | Método da matriz DIWV (Guruprasad et al., 1990). Escore < 40 = estável |
| **Escore GRAVY** | Média Grand da Hidropatiticidade (Kyte & Doolittle, 1982) |
| **Índice Alifático** | `100 × (Ala + 2,9Val + 3,9(Ile + Leu)) / N` (Ikai, 1980) |
| **Coeficiente de Extinção Molar** | Pace et al. (1995); formas reduzida e oxidada (dissulfeto Cys-Cys) |
| **Fração de Estrutura Secundária 2D** | Propensidades estatísticas de Chou-Fasman: α-hélice, β-folha, volta (Chou & Fasman, 1978) |

```bash
python Lenght.py --protein-file proteinas.faa
```

### Gráfico de Hidropatiticidade de Kyte-Doolittle

Uma média de hidrofobicidade em janela deslizante (padrão: 9 resíduos) é computada usando a escala de Kyte-Doolittle (Kyte & Doolittle, 1982). Segmentos acima de **+1,6** são anotados como putativos domínios transmembrana, consistente com o critério TMbase (Hofmann & Stoffel, 1993).

```bash
python Lenght.py --hydro-file proteinas.faa
```

### Termodinâmica de Oligonucleotídeos — Wizard de Primers

O `minato_rasengan_primer_wizard` calcula para qualquer sequência de primer:

| Parâmetro | Método |
|---|---|
| **Tm (Temperatura de Desnaturação)** | Modelo termodinâmico de Vizinhos Mais Próximos (`Bio.SeqUtils.MeltingTemp.Tm_NN`) — SantaLucia, 1998 |
| **ΔG (Energia Livre de Gibbs)** | `ΔG = ΔH - T × ΔS` a 37°C (310,15 K) |
| **GC%** | Contagem direta; faixa ótima: 40–60% |
| **Risco de Auto-Dímero** | Detecção por sobreposição de prefixo do complemento reverso |
| **Risco de Hairpin** | Detecção de match do complemento reverso intra-sequência (haste ≥ 4 pb) |

```bash
python Lenght.py --primer-file primers.fasta
```

### Scanner de PAM para CRISPR SpCas9

O `kakashi_crispr_copy_sgrna` varre a sequência em busca de todos os sítios PAM `NGG` usando uma regex com lookahead. Para cada sítio válido, extrai:

- O protoespaçador de 20 nt (sgRNA)
- GC% do guia (faixa ótima: 40–70%)
- Estimativa de risco off-target baseada na frequência exata do 20-mer na sequência fornecida

O algoritmo segue as regras canônicas de targeting do CRISPR-Cas9 de Doench et al. (2014) e Hsu et al. (2013).

```bash
python Lenght.py --crispr genome.fasta
```

---

## 🔵 Módulo IV — Métricas FastQ (Rasengan)

> *"Qualidade é tudo. Até o Naruto aprendeu isso."*

### Parsing de FASTQ (Baseado em Gerador/Blocos)

O `naruto_rasengan_read_fastq_blocks` usa o `Bio.SeqIO.QualityIO.FastqGeneralIterator` — um leitor baseado em gerador que processa arquivos FASTQ linha a linha sem carregar o arquivo inteiro na memória. Suporta `.fastq` simples e `.fastq.gz` comprimidos com gzip.

### Perfil de Qualidade Phred por Ciclo

Os escores de qualidade Phred (Q) são decodificados a partir da codificação ASCII:

```
Q = ord(char) - 33  [Phred+33, Illumina 1.8+]
```

A qualidade média por posição de ciclo é computada e plotada, revelando a característica **queda de qualidade na extremidade 3'** comum no sequenciamento Illumina de leituras curtas (Ewing & Green, 1998).

```bash
python Lenght.py leituras.fastq.gz --phred
```

### Estimativa de Cobertura do Genoma

Equação de cobertura de Lander-Waterman (Lander & Waterman, 1988):

```
Cobertura (C) = (N_leituras × C_leitura × multiplicador) / G
```

Onde:
- `N_leituras` = contagem total de leituras
- `C_leitura` = comprimento médio das leituras (configurável, padrão 150 pb)
- `multiplicador` = 2 para paired-end, 1 para single-end
- `G` = tamanho esperado do genoma (configurável, padrão 1 Mb)

```bash
python Lenght.py leituras.fastq --coverage --gsize 3000000000 --rlen 150
```

### Curva de Rarefação / Saturação

Simula amostragem progressiva de leituras e plota a cobertura cumulativa de k-mers únicos (representativos) contra o número de leituras. Um platô na curva indica saturação do sequenciamento — ou seja, leituras adicionais trarão retornos decrescentes de informação nova (Lander & Waterman, 1988; Good, 1953).

```bash
python Lenght.py leituras.fastq --rarefaction
```

### Estimador de Duplicatas de PCR

Analisa os primeiros 50 pb de cada leitura. Uma alta proporção (> 30%) de fragmentos de 50-mers 100% idênticos é um forte indicador de **viés de amplificação por PCR** — um artefato comum na preparação de bibliotecas (Kozarewa et al., 2009). Esta heurística é análoga ao flag de deduplicação do FastQC.

### Trimador de FASTQ (Gaara Sand Trimmer)

Trimagem baseada em qualidade da extremidade 3': leituras são escritas em um novo FASTQ somente se seu escore Phred médio exceder um limiar mínimo configurável (padrão Q20, equivalente a 99% de acurácia na chamada de base).

```bash
python Lenght.py --trim-fastq leituras.fastq --min-phred 20
```

---

## 🟣 Módulo V — Análise Avançada de Sequências (Rinnegan)

> *"Os olhos que enxergam todos os caminhos."*

### Alinhamento Global — Needleman-Wunsch

O `sasuke_sharingan_snp_deducer_global` implementa alinhamento global par-a-par usando o `Align.PairwiseAligner` do BioPython em **modo global** com o seguinte esquema de pontuação:

| Parâmetro | Valor |
|---|---|
| Pontuação de match | +1 |
| Penalidade de mismatch | -2 |
| Penalidade de abertura de gap | -5 |
| Penalidade de extensão de gap | -1 |

Equivalente ao algoritmo de Needleman-Wunsch (Needleman & Wunsch, 1970) com penalidades afins de gap.

**Análise de SNPs pós-alinhamento:**

- **Razão Ti/Tv**: Transições (A↔G, C↔T) vs Transversões (A/G↔C/T). Uma razão de ~2,0 é esperada para mutação biológica; valores menores sugerem artefatos de erro de sequenciamento (Collins & Jukes, 1994).
- **Classificação Sinônima/Não-Sinônima**: Cada SNP é classificado comparando o aminoácido codificado pelo códon de referência vs o códon mutante.

```bash
python Lenght.py --snp-files referencia.fasta query.fasta
```

### Matriz de Dot Plot

O `madara_rinnegan_dotplot_matrix_numpy` gera uma **matriz de similaridade de sequências** usando broadcasting vetorizado do NumPy para comparação O(N × M), seguido de convolução diagonal com `scipy.signal.convolve2d` para filtrar matches por threshold de janela deslizante. Limitado a sequências ≤ 4.000 pb para evitar estouro de memória. Esta abordagem espelha o método clássico de dotplot de Gibbs & McIntyre (1970).

```bash
python Lenght.py --dotplot-files seq1.fasta seq2.fasta
```

### Detecção de Repetições em Tandem

O `choji_expansion_jutsu_tandem_repeats` utiliza padrões regex para detectar microssatélites e repetições em tandem de comprimento de unidade 2–9 pb com ≥ 3 cópias, equivalente à abordagem do TRF (Tandem Repeats Finder; Benson, 1999) em seu modo básico de detecção de motivos.

```bash
python Lenght.py genome.fasta --tandem
```

### Árvore Filogenética UPGMA

O `hashirama_phylo_tree_builder` constrói um dendrograma filogenético a partir de um **alinhamento múltiplo de sequências (FASTA)** utilizando:

1. **Matriz de distância por identidade** (`Bio.Phylo.TreeConstruction.DistanceCalculator('identity')`)
2. **Clusterização UPGMA (Método de Grupo Par Não-Ponderado com Média Aritmética)** (`DistanceTreeConstructor.upgma()`)

O UPGMA pressupõe um relógio molecular e taxas de substituição uniformes (Sokal & Michener, 1958). É adequado para sequências proximamente relacionadas ou demonstrações educacionais.

```bash
python Lenght.py --phylo alinhamento.fasta
```

### Integração com BLAST+ Local

O `deidara_explosive_blast_art` encapsula o binário `blastn` do NCBI BLAST+ via subprocess, escrevendo arquivos FASTA temporários, executando o alinhamento e parseando o XML resultante via `Bio.Blast.NCBIXML`. Requer instalação local do BLAST+ (Altschul et al., 1990).

```bash
python Lenght.py --blast-path /usr/bin/blastn --blast-query query.fasta --blast-subject subject.fasta
```

### Extração de Genes a partir de GFF

O `YamatoWoodFeatureExtractor.mokuton_extract_feature` parseia um arquivo de anotação GFF3 para recuperar as coordenadas cromossômicas de uma feature nomeada, extraindo cirurgicamente a sequência correspondente do FASTA de referência via `SeqIO` do BioPython. A orientação da fita é respeitada (complemento reverso aplicado para features na fita `-`).

```bash
python Lenght.py --gff anotacao.gff3 --fasta-ref genome.fasta --gene "BRCA1"
```

### Gráfico de Densidade Gênica

Utiliza `seaborn.histplot` com bins de 10.000 pb e overlay de KDE para visualizar a distribuição das posições de início dos genes ao longo do cromossomo — uma métrica comum de qualidade de anotação genômica.

---

## 🔷 Módulo VI — Visualizador de Variantes (Olho do Neji)

> *"O Byakugan enxerga cada SNP."*

O `NejiVariantEyeParser.eight_trigrams_parse_vcf` parseia arquivos **VCF (Variant Call Format)** padrão (v4.x), incluindo `.vcf.gz` comprimidos com gzip. Extrai as colunas obrigatórias:

```
CHROM | POS | ID | REF | ALT | QUAL | INFO
```

A interface gráfica constrói um `QTableWidget` completamente interativo e pesquisável a partir dos registros resultantes, com filtragem por cromossomo, posição, ID da variante e alelos REF/ALT. O parser ignora corretamente as linhas de metadados `##` e a linha de cabeçalho `#CHROM`.

---

## ⚔️ Módulo VII — QG Flutuante Anbu (Terminal Integrado)

O terminal integrado (`AnbuBlackOpsExternalTools.itachi_mangekyou_run_command`) fornece um widget de shell acoplável dentro da interface gráfica. Executa comandos shell arbitrários via `subprocess.run` do Python com `shell=True`, capturando tanto os fluxos stdout quanto stderr.

**Usos recomendados:**
- Iniciar containers Docker com ferramentas bioinformáticas (ex.: `docker run broadinstitute/gatk ...`)
- Executar scripts de terceiros enquanto a GUI principal permanece ativa
- Disparar samtools, BWA, STAR ou outras ferramentas externas em paralelo

> ⚠️ **Aviso:** O shell executa com as mesmas permissões do processo Python em execução. Use com responsabilidade.

---

## 📟 Referência da CLI

```
uso: Lenght.py [-h] [--gui] [--out-prefix OUT_PREFIX] [--label LABEL]
               [--histogram] [--boxplot] [--stats STATS] [--csv CSV]
               [--kmer K] [--entropy WINDOW] [--cpg] [--gc-skew]
               [--orf] [--tandem]
               [--coverage] [--gsize GSIZE] [--rlen RLEN]
               [--phred] [--rarefaction]
               [--primer-file PRIMER_FILE] [--protein-file PROTEIN_FILE]
               [--hydro-file HYDRO_FILE]
               [--snp-files REF QUERY] [--dotplot-files SEQ1 SEQ2]
               [--ncbi NCBI] [--phylo PHYLO]
               [--restriction RESTRICTION] [--crispr CRISPR]
               [--vcf-table VCF_TABLE]
               [--blast-path BLAST_PATH] [--blast-query BLAST_QUERY]
               [--blast-subject BLAST_SUBJECT]
               [--gff GFF] [--fasta-ref FASTA_REF] [--gene GENE]
               [--trim-fastq TRIM_FASTQ] [--min-phred MIN_PHRED]
               [--convert-in CONVERT_IN] [--convert-out CONVERT_OUT]
               [--fmt-in FMT_IN] [--fmt-out FMT_OUT]
               [arquivos ...]
```

### Argumentos Principais

| Flag | Descrição |
|---|---|
| `--gui` | Iniciar a interface gráfica PyQt6 |
| `--histogram` | Gerar histograma de comprimento (Freedman-Diaconis) |
| `--boxplot` | Gerar boxplot de comprimento de sequências |
| `--stats ARQUIVO` | Escrever relatório de estatísticas de assembly em ARQUIVO |
| `--csv ARQUIVO` | Exportar matriz de comprimentos brutos para CSV |
| `--kmer K` | Calcular e plotar frequência de K-mer (ex.: `--kmer 6`) |
| `--entropy W` | Plotar entropia de Shannon com janela de W pb |
| `--cpg` | Detectar e plotar ilhas CpG |
| `--gc-skew` | Plotar GC skew deslizante (predição de oriC) |
| `--orf` | Localizar ORFs em todos os 6 quadros de leitura |
| `--tandem` | Detectar microssatélites e repetições em tandem |
| `--coverage` | Estimar cobertura do genoma (usar com `--gsize`, `--rlen`) |
| `--phred` | Plotar perfil de qualidade Phred por ciclo |
| `--rarefaction` | Plotar curva de saturação/rarefação do sequenciamento |
| `--primer-file` | Analisar termodinâmica de primers a partir de FASTA |
| `--protein-file` | Propriedades físico-químicas de proteínas a partir de FASTA |
| `--hydro-file` | Gráfico de hidropatiticidade de Kyte-Doolittle |
| `--snp-files REF QUERY` | Alinhamento global + Ti/Tv + análise de sinonímia |
| `--dotplot-files S1 S2` | Matriz de dot plot de sintenia |
| `--ncbi ACESSO` | Baixar sequência do NCBI por número de acesso |
| `--phylo MSA.fasta` | Construir árvore filogenética UPGMA a partir de MSA |
| `--restriction SEQ` | Digestão de restrição + gel de agarose virtual |
| `--crispr SEQ` | Varredura de PAM SpCas9 + lista de sgRNAs |
| `--vcf-table VCF` | Parsear e exibir registros de variantes VCF |
| `--blast-query / --blast-subject` | Executar alinhamento BLAST+ local |
| `--gff + --fasta-ref + --gene` | Extrair sequência de gene a partir de anotação GFF |
| `--trim-fastq + --min-phred` | Trimar arquivo FASTQ por qualidade |
| `--convert-in/out --fmt-in/out` | Converter entre BAM/SAM/FASTQ/GBK/FASTA |

---

## 📚 Metodologias Científicas & Referências

| Método | Implementação | Referência |
|---|---|---|
| Métricas N50 / NG50 de assembly | `shikamaru_shadow_stats` | Schatz et al. (2010); Miller et al. (2010) |
| Binagem de Freedman-Diaconis | `sasuke_amaterasu_histogram` | Freedman & Diaconis (1981) |
| Entropia de Shannon | `orochimaru_curse_mark_entropy_plot` | Shannon (1948) |
| Mascaramento DUST de baixa complexidade | `gaara_dust_complexity_shield` | Morgulis et al. (2006) |
| Detecção de Ilhas CpG (razão O/E) | `hidan_jashin_cpg_islands_vectorized` | Gardiner-Garden & Frommer (1987) |
| GC Skew / Predição de oriC | `kakashi_kamui_gc_skew_rolling` | Grigoriev (1998) |
| Razão Rho de Dinucleotídeos | `neji_trigram_rho_odds` | Karlin & Ladunga (1994) |
| Assinatura genômica por K-mer | `shino_kikaichu_kmer_swarm` | Karlin et al. (1997) |
| Detecção de Ilhas de Patogenicidade | `orochimaru_cursed_islands_pai` | Hacker & Kaper (2000) |
| Mapeamento de Telômeros | `kimimaro_bone_telomere_mapper` | Meyne et al. (1989) |
| Análise de enzimas de restrição | `zabuza_executioner_restriction_map` | BioPython: Cock et al. (2009) |
| Predição de ORF (6 quadros) | `jiraiya_sage_mode_orf_finder` | Código genético padrão; BioPython |
| pI, PM e instabilidade de proteínas | `raikage_lightning_armor_protein_properties` | Gasteiger et al. (2005); Guruprasad et al. (1990) |
| GRAVY / Índice Alifático | `raikage_lightning_armor_protein_properties` | Kyte & Doolittle (1982); Ikai (1980) |
| Coeficiente de extinção molar | `raikage_lightning_armor_protein_properties` | Pace et al. (1995) |
| Estrutura secundária 2D de Chou-Fasman | `raikage_lightning_armor_protein_properties` | Chou & Fasman (1978) |
| Hidropatiticidade de Kyte-Doolittle | `kisame_water_prison_hydrophobicity_plot` | Kyte & Doolittle (1982) |
| Tm por Vizinhos Mais Próximos (ΔG, ΔH, ΔS) | `minato_rasengan_primer_wizard` | SantaLucia (1998) |
| Varredura de PAM CRISPR | `kakashi_crispr_copy_sgrna` | Hsu et al. (2013); Doench et al. (2014) |
| Pontuação de qualidade Phred | `naruto_rasengan_read_fastq_blocks` | Ewing & Green (1998) |
| Cobertura de Lander-Waterman | `naruto_oodama_rasengan_calculate_genomes` | Lander & Waterman (1988) |
| Estimativa de duplicatas de PCR | Análise de fingerprint 50-mer em FASTQ | Kozarewa et al. (2009) |
| Alinhamento global Needleman-Wunsch | `sasuke_sharingan_snp_deducer_global` | Needleman & Wunsch (1970) |
| Razão Ti/Tv | Parser de SNPs pós-alinhamento | Collins & Jukes (1994) |
| Dot plot / Matriz de sintenia | `madara_rinnegan_dotplot_matrix_numpy` | Gibbs & McIntyre (1970) |
| Detecção de repetições em tandem | `choji_expansion_jutsu_tandem_repeats` | Benson (1999) |
| Árvore filogenética UPGMA | `hashirama_phylo_tree_builder` | Sokal & Michener (1958) |
| BLAST+ local | `deidara_explosive_blast_art` | Altschul et al. (1990) |

---

## 📂 Formatos de Arquivo Suportados

| Formato | Extensão(ões) | Operações |
|---|---|---|
| FASTA | `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz` | Leitura, escrita, mesclagem, extração, fatiamento |
| FASTQ | `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` | Leitura, QC, trimagem, conversão |
| GenBank | `.gb`, `.gbk`, `.gbff` | Leitura, parsing de features, conversão para FASTA |
| VCF | `.vcf`, `.vcf.gz` | Parsing, contagem, visualização em tabela |
| GFF3 | `.gff`, `.gff3` | Parsing, extração de features, gráfico de densidade |
| BED | `.bed` | Contagem de features |
| SAM | `.sam` | Contagem de alinhamentos |
| BAM | `.bam` | Detecção; conversão para SAM via `samtools` |

---

## 🤝 Contribuindo

Contribuições, reportes de bugs e solicitações de features são bem-vindos. Por favor, abra uma issue ou envie um pull request.

1. Faça um fork do repositório
2. Crie sua branch de feature (`git checkout -b feature/novo-jutsu`)
3. Faça o commit das suas mudanças (`git commit -m 'Add: novo jutsu bioinformático'`)
4. Envie para a branch (`git push origin feature/novo-jutsu`)
5. Abra um Pull Request

---

<div align="center">

**Desenvolvido por VictorSC**  
*"A Vontade do Fogo aplicada à Bioinformática."*

🍃 *Acredite nisso!* 🍃

</div>
