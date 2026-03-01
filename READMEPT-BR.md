# Ferramentas básicas para bioinformatas

Olá amigos,

Esse repositório são algumas ferramentas que eu fiz enquanto trabalho com bioinformática, elas foram feitas devido demandas de pesquisas que eu tive e ainda tenho. Porém acabei criando temas divertidos para elas, e melhorei algumas coisas para elas ficarem mais completas, porém ainda são ferramentas de análises básicas para bioinformática. Sinta-se a vontade para usar, dar fork, e também fazer a sua própria a partir dessas. Sintam-se livres para reportar erros, dar sujestões e também pedir melhorias ou implementações (novas análises e funcionalidades), pois posso adicionar ao longo do tempo.


Atenciosamente,
VictorSC

---

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/Licença-MIT-green)](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)
[![Status](https://img.shields.io/badge/Status-Ativo-brightgreen)]()
[![Institution](https://img.shields.io/badge/Instituição-UFMG-blue)]()

**Repositório:** [https://github.com/VictorCaricatte/BasicBioinfo](https://github.com/VictorCaricatte/BasicBioinfo)

---

## Sumário

- [Visão Geral das Ferramentas](#visão-geral-das-ferramentas)
- [BRSA — Análise Básica de RNA-Seq](#brsa--análise-básica-de-rna-seq)
- [DELIMITER — Filtro de Arquivos Delimitados](#delimiter--filtro-de-arquivos-delimitados)
- [LenghtCount — Suíte de Análise Genômica e Proteômica](#lenghtcount--suíte-de-análise-genômica-e-proteômica)
- [PALINDROME — Analisador de Palíndromos no DNA](#palindrome--analisador-de-palíndromos-no-dna)
- [Instalação](#instalação)
- [Licença](#licença)
- [Autor](#autor)

---

## Visão Geral das Ferramentas

| Ferramenta | Descrição |
|------------|-----------|
| [BRSA](#brsa--análise-básica-de-rna-seq) | Pipeline completo de análise de RNA-Seq |
| [DELIMITER](#delimiter--filtro-de-arquivos-delimitados) | Aplicativo desktop para filtragem e transformação de arquivos tabulares e genômicos |
| [LenghtCount](#lenghtcount--suíte-de-análise-genômica-e-proteômica) | Suíte multi-módulo para caracterização genômica e proteômica |
| [PALINDROME](#palindrome--analisador-de-palíndromos-no-dna) | Detecção e análise de sequências palindrômicas no DNA |

---

## BRSA — Análise Básica de RNA-Seq

**Localização:** `BRSA/`

O BRSA é um pipeline completo de análise de sequenciamento de RNA desenvolvido na Universidade Federal de Minas Gerais (UFMG). Ele guia os pesquisadores por todo o fluxo de trabalho transcriptômico — desde matrizes de contagens brutas até figuras prontas para publicação, relatórios HTML interativos e logs de reprodutibilidade.

O pipeline integra métodos estatísticos clássicos com aprendizado de máquina e biologia de redes, expondo todos os módulos analíticos por meio de uma CLI unificada ou de uma interface gráfica baseada em PyQt5.

### Principais Funcionalidades

- **Análise de Expressão Diferencial** usando um modelo binomial negativo inspirado no DESeq2 com correção FDR de Benjamini-Hochberg
- **Normalização** via CPM, TPM ou estabilização de variância do DESeq2
- **Redução de Dimensionalidade** com PCA, t-SNE e UMAP
- **Aprendizado de Máquina** com classificadores Random Forest e SVM para descoberta de biomarcadores
- **Análise de Co-expressão** baseada em clustering hierárquico no estilo WGCNA
- **Interatoma Proteína-Proteína** consultado no banco de dados STRING
- **Análise de Sobrevivência** com curvas de Kaplan-Meier estratificadas por expressão gênica
- **Análise de Enriquecimento Gênico (GSEA e ORA)** via GSEApy, com suporte a GO, KEGG, Reactome e MSigDB Hallmarks
- **Gráficos Interativos** (artefatos HTML em Plotly) para visualizações de vulcão, MA e PCA
- **Log de Reprodutibilidade** com registro completo de parâmetros e snapshot do `pip freeze`
- Suporte à **correção de efeito de lote** via covariáveis de metadados

### Início Rápido

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/BRSA
pip install -r requirements.txt

# Uso via CLI
python BRSA.py --counts data/counts_matrix.csv --meta data/metadata.csv \
  --cond1 Tumor --cond2 Normal --norm DESeq2 --dim_red UMAP --ml --ppi

# Uso via interface gráfica
python interface.py
```

### Dependências Externas

Requer as seguintes ferramentas disponíveis no `PATH` para as etapas de pré-processamento: `fastqc`, `multiqc`, `fastp`, `hisat2`, `samtools`, `featureCounts` e `kallisto`.

---

## DELIMITER — Filtro de Arquivos Delimitados

**Localização:** `Delimiter/`

O DELIMITER é um aplicativo desktop para inspeção, filtragem, transformação e exportação de arquivos de dados delimitados. Construído com PyQt5 e Pandas, oferece uma interface de processamento de dados multi-threaded projetada principalmente para dados tabulares de bioinformática, mas aplicável a qualquer arquivo estruturado.

### Principais Funcionalidades

- **Carregamento de Arquivos** CSV, TSV e Excel (`.xlsx`) com detecção automática de delimitador; a leitura ocorre em threads dedicadas para evitar travamentos na interface
- **Filtragem Multi-condição** com pods de filtros empilhados, suportando operadores como igualdade, contém, começa/termina com, regex e comparações numéricas
- **Gerenciamento de Presets de Filtros** — salve e recarregue configurações de filtros em formato JSON
- **Operações de Colunas** — crie colunas derivadas via concatenação, adição, subtração, multiplicação ou divisão
- **Tratamento de NaN** — remova, preencha com média/mediana ou com um valor fixo personalizado
- **Edição de Células** com histórico completo de desfazer/refazer
- **Destaque Condicional de Células** — defina regras visuais para sinalizar valores na tabela
- **Exportação** para CSV, TSV, Excel e FASTA, todos gravados em threads com relatório de progresso
- **Painel de Estatísticas Descritivas** para qualquer coluna selecionada
- **Internacionalização** com suporte nativo a inglês, português brasileiro e espanhol
- **Temas Claro e Escuro** alternáveis em tempo de execução

### Início Rápido

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Delimiter
pip install PyQt5 pandas openpyxl psutil
python delimeter.py
```

---

## LenghtCount — Suíte de Análise Genômica e Proteômica

**Localização:** `LengthCounter/`

O LenghtCount é uma plataforma bioinformática multi-modal para caracterização genômica e proteômica abrangente. Ele consolida uma ampla gama de tarefas de análise de sequências em uma única ferramenta acessível via CLI ou interface gráfica PyQt6.

### Principais Funcionalidades

**Estatísticas de Montagem e Sequências**
- N50, N90, L50, L90, NG50 e estatísticas descritivas completas para montagens FASTA
- Histograma, boxplot, gráfico de violino, ECDF e ridgeline dos comprimentos de contigs

**Análise Nucleotídica e de Composição**
- Conteúdo GC, desvio GC (predição de oriC), detecção de ilhas CpG, perfil de frequência de k-mers
- Entropia de Shannon, mascaramento de baixa complexidade (DUST), razão Rho de dinucleotídeos
- Detecção de repetições em tandem e microssatélites, assinaturas de ilhas de patogenicidade, mapeamento de telômeros

**Controle de Qualidade de FASTQ**
- Perfis de qualidade Phred por ciclo, composição de bases, distribuição de GC
- Estimativa de profundidade de sequenciamento (Lander-Waterman), curvas de rarefação/saturação, estimativa de duplicatas de PCR
- Trimagem de FASTQ baseada em qualidade

**Análise de Proteínas**
- Propriedades físico-químicas: peso molecular, ponto isoelétrico, índice de instabilidade, escore GRAVY, índice alifático, coeficiente de extinção molar
- Gráficos de hidrofobicidade de Kyte-Doolittle
- Predição de estrutura secundária (Chou-Fasman)

**Genômica Comparativa e Estrutural**
- Alinhamento global par a par (Needleman-Wunsch), razão Ti/Tv, análise de substituições sinônimas/não-sinônimas
- Dot plot / matriz de sintenia
- Mapeamento de enzimas de restrição e simulação de gel de agarose virtual
- Varredura de sítios PAM NGG (CRISPR/SpCas9)
- Predição de ORF (6 fases de leitura)
- Construção de árvore filogenética UPGMA a partir de alinhamentos múltiplos de sequências
- Integração com BLAST+ local

**Análise de Variantes e Anotação**
- Leitura de arquivos VCF e visualizador interativo de variantes em tabela
- Extração de sequências gênicas guiada por anotação GFF3

**Conversão de Formatos**
- Interconversão entre FASTA, FASTQ, GenBank, SAM e BAM (via samtools)

**Integração com NCBI**
- Download de sequências por número de acesso via NCBI Entrez

### Início Rápido

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/LengthCounter
pip install biopython numpy pandas matplotlib seaborn scipy PyQt6

# Estatísticas de montagem
python Lenght.py assembly.fasta --stats output.txt --histogram

# Abrir interface gráfica
python Lenght.py --gui
```

---

## PALINDROME — Analisador de Palíndromos no DNA

**Localização:** `Palindrome/`

O PALINDROME é uma ferramenta de linha de comando em Python para detecção e caracterização de sequências palindrômicas no DNA. Regiões palindrômicas — onde o DNA de fita dupla é lido de forma idêntica nas direções 5'→3' em ambas as fitas — são biologicamente relevantes como sítios de reconhecimento de enzimas de restrição, sítios de ligação de fatores de transcrição e reguladores estruturais da expressão gênica.

### Principais Funcionalidades

- **Detecção de Palíndromos** com intervalo de comprimento e tolerância a mismatches configuráveis
- **Perfil Termodinâmico** — temperatura de fusão (Tm) e energia livre de Gibbs (ΔG) pelo modelo de vizinhos mais próximos
- **Conteúdo GC** e **Entropia de Shannon** por palíndromo
- **Correspondência com Enzimas de Restrição** em 15 sítios canônicos (EcoRI, BamHI, HindIII, NotI, entre outros)
- **Detecção de Motivos de Metilação** — Dam (GATC), Dcm (CCWGG) e sítios CpG
- **Mapeamento de Sítios PAM do CRISPR/Cas9** (NGG em ambas as fitas)
- **Identificação de Ilhas CpG** via análise por janela deslizante de GC%/razão O/E
- **Predição de ORF** em 6 fases de leitura
- **Detecção de Repetições em Tandem** e **Mascaramento de Baixa Complexidade**
- **Design de Primers de PCR** com verificação de GC%, Tm e hairpin flanqueando cada palíndromo
- **Tradução de Códons** nas três fases de leitura
- **Predição de Shine-Dalgarno e de Terminadores de Transcrição**
- **Razão Ti/Tv** e **Árvore Filogenética UPGMA** para análise comparativa
- **Múltiplas Fontes de Entrada** — sequência bruta, arquivos FASTA/FASTQ/GenBank, número de acesso NCBI, geração de sequência simulada ou processamento em lote de diretórios
- **Formatos de Saída** — tabela CSV, mapa linear SVG, mapa circular SVG e relatório HTML interativo

### Início Rápido

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Palindrome
pip install matplotlib numpy biopython

# Varredura básica
python Palindrome.py -s "ATGCGCAT"

# Análise completa a partir de um arquivo FASTA
python Palindrome.py -f genome.fasta -min 4 -max 12 -mis 1 --crispr --report

# Gerar sequência simulada e testar
python Palindrome.py --mock 2000 --crispr --report
```

---

## Instalação

Cada ferramenta possui suas próprias dependências. Consulte a seção correspondente acima para instruções específicas de instalação. Todos os projetos compartilham o seguinte requisito base:

```
Python >= 3.8
```

Clone o repositório completo para acessar todas as ferramentas:

```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
```

---

## Licença

Todas as ferramentas deste repositório são distribuídas sob a **Licença MIT**.  
Veja o texto completo da licença em: [LICENSE](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)

---

## Autor

**Victor S. Caricatte De Araújo**  
Universidade Federal de Minas Gerais (UFMG)  
GitHub: [@VictorCaricatte](https://github.com/VictorCaricatte)  
Repositório: [BasicBioinfo](https://github.com/VictorCaricatte/BasicBioinfo)
