<div align="center">

```
██████╗ ██████╗ ███████╗ █████╗ 
██╔══██╗██╔══██╗██╔════╝██╔══██╗
██████╔╝██████╔╝███████╗███████║
██╔══██╗██╔══██╗╚════██║██╔══██║
██████╔╝██║  ██║███████║██║  ██║
╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝
```

# 🕷️ BRSA — Análise Básica de RNA-Seq 🕷️
### *"Das Catacumbas das Reads Brutas, a Vida Há de Ressurgir."*

[![Python](https://img.shields.io/badge/Python-3.9%2B-8a2be2?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![Licença](https://img.shields.io/badge/Licença-MIT-ff6600?style=for-the-badge)](https://github.com/VictorCaricatte/BasicBioinfo/blob/main/LICENSE)
[![Status](https://img.shields.io/badge/Status-Vivo%20(por%20pouco)-darkred?style=for-the-badge)]()
[![Plataforma](https://img.shields.io/badge/Plataforma-Windows%20%7C%20Linux%20%7C%20macOS-8a2be2?style=for-the-badge)]()

*Um pipeline completo e profano de análise de RNA-Seq, invocado das profundezas da arcana bioinformática — com expressão diferencial, aprendizado de máquina, interatoma proteico e análise de sobrevivência.*

---

</div>

## 🪦 Sumário

- [A Maldição (Visão Geral)](#-a-maldição-visão-geral)
- [O Aquelarre (Funcionalidades)](#-o-aquelarre-funcionalidades)
- [Ingredientes da Invocação (Instalação)](#-ingredientes-da-invocação-instalação)
- [Ressuscitando os Mortos (Uso)](#-ressuscitando-os-mortos-uso)
  - [O Ritual Negro via Terminal (CLI)](#o-ritual-negro-via-terminal-cli)
  - [A Sessão Espírita Gráfica (GUI)](#a-sessão-espírita-gráfica-gui)
- [A Arquitetura Sagrada](#-a-arquitetura-sagrada)
- [O Grimório dos Métodos](#-o-grimório-dos-métodos)
  - [Controle de Qualidade e Pré-processamento](#i-controle-de-qualidade-e-pré-processamento)
  - [Alinhamento e Quantificação](#ii-alinhamento-e-quantificação)
  - [Normalização](#iii-normalização)
  - [Análise de Expressão Diferencial](#iv-análise-de-expressão-diferencial)
  - [Redução de Dimensionalidade](#v-redução-de-dimensionalidade)
  - [Golens de Aprendizado de Máquina](#vi-golens-de-aprendizado-de-máquina)
  - [Análise de Redes de Co-expressão](#vii-análise-de-redes-de-co-expressão)
  - [Interatoma Proteína-Proteína](#viii-interatoma-proteína-proteína)
  - [Análise de Sobrevivência](#ix-análise-de-sobrevivência)
  - [Análise de Enriquecimento de Conjuntos Gênicos](#x-análise-de-enriquecimento-de-conjuntos-gênicos)
- [Os Parâmetros Eldritchianos (Referência CLI)](#-os-parâmetros-eldritchianos-referência-cli)
- [Os Artefatos do Caldeirão](#-os-artefatos-do-caldeirão)
- [Pergaminhos de Reprodutibilidade](#-pergaminhos-de-reprodutibilidade)
- [Referências e Tomos Ancestrais](#-referências-e-tomos-ancestrais)
- [Licença e Pacto das Trevas](#-licença-e-pacto-das-trevas)
- [O Necromante (Autor)](#-o-necromante-autor)

---

## 🕸️ A Maldição (Visão Geral)

O **BRSA** é um pipeline abrangente e completo de análise de sequenciamento de RNA, desenvolvido na **Universidade Federal de Minas Gerais (UFMG)**. Foi conjurado para guiar pesquisadores por todo o ritual transcriptômico — desde os arquivos FASTQ brutos até figuras prontas para publicação, artefatos HTML interativos e registros de reprodutibilidade.

O pipeline integra métodos estatísticos clássicos com aprendizado de máquina moderno e biologia de redes, permitindo que bioinformatas mortais empunhem um arsenal analítico por meio de uma única interface unificada — seja por uma invocação no terminal ou por uma sessão espírita gráfica completa.

> *"Aquilo que não mata o pipeline apenas o torna mais estatisticamente significativo."*

---

## 🦇 O Aquelarre (Funcionalidades)

| Módulo | Poder Concedido |
|--------|----------------|
| 🧟 **Controle de Qualidade** | Interrogatório FastQC + olho onividente MultiQC |
| ⚗️ **Trimagem de Adaptadores** | Guilhotina fastp para purgar adaptadores e baixa qualidade |
| 🗺️ **Alinhamento** | Vórtice de alinhamento spliced-aware do HISAT2 |
| 💀 **Pseudo-alinhamento** | Fantasma ultra-veloz de quantificação Kallisto |
| 🔮 **Normalização** | CPM, TPM e estabilização de variância DESeq2 |
| 🌋 **Expressão Diferencial** | Teste binomial negativo inspirado no DESeq2 |
| 🕳️ **Redução de Dimensionalidade** | Projeções PCA, t-SNE e UMAP |
| 🤖 **Aprendizado de Máquina** | Random Forest + SVM para descoberta de biomarcadores |
| 🕸️ **Interatoma PPI** | Invocação de rede proteica via STRING DB |
| 🧬 **Análise de Isoformas** | Detecção de eventos de splicing alternativo |
| 💉 **Análise de Sobrevivência** | Curvas de doom mortal de Kaplan-Meier |
| 📜 **GSEA** | Enriquecimento gênico com Enrichr/GSEApy |
| 🗺️ **Gráficos Interativos** | Artefatos HTML forjados com Plotly |
| 📖 **Reprodutibilidade** | Registro completo de parâmetros e ambiente computacional |

---

## 🧪 Ingredientes da Invocação (Instalação)

### Pré-requisitos — Os Demônios Externos

Os seguintes espíritos de linha de comando devem estar instalados e acessíveis no `PATH` do seu sistema antes de invocar o BRSA:

```bash
# Demônios de Controle de Qualidade
fastqc --version
multiqc --version

# Espíritos de Trimagem e Alinhamento
fastp --version
hisat2 --version
samtools --version
featureCounts -v      # Pacote Subread

# Fantasma de Pseudo-alinhamento
kallisto version
```

### Ambiente Alquímico Python

```bash
# Clone o repositório proibido
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BRSA

# Conjure um santuário virtual (recomendado)
python -m venv brsa_aquelarre
source brsa_aquelarre/bin/activate          # Linux/macOS
# brsa_aquelarre\Scripts\activate.bat       # Windows (se você ousar)

# Instale as bibliotecas arcanas necessárias
pip install -r requirements.txt
```

### `requirements.txt` — Os Ingredientes Proibidos

```
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.6.0
seaborn>=0.12.0
scipy>=1.9.0
scikit-learn>=1.1.0
pydeseq2>=0.3.0
gseapy>=1.0.0
plotly>=5.10.0
mplcursors>=0.5.0
networkx>=2.8.0
requests>=2.28.0
lifelines>=0.27.0
umap-learn>=0.5.3
PyQt5>=5.15.0        # Apenas para a sessão espírita gráfica (GUI)
```

---

## 💀 Ressuscitando os Mortos (Uso)

### O Ritual Negro via Terminal (CLI)

Para os que preferem trabalhar nas sombras do terminal, o BRSA pode ser invocado diretamente sem carregar nenhum espírito gráfico:

```bash
python BRSA.py \
  --counts dados/matriz_contagens.csv \
  --meta dados/metadados.csv \
  --cond1 Tumor \
  --cond2 Normal \
  --norm DESeq2 \
  --dim_red UMAP \
  --pval 0.05 \
  --min_counts 10 \
  --batch coluna_batch \
  --ml \
  --ppi
```

Ao término, dois artefatos sombrios se materializarão:
- `CLI_differential_results.csv` — O tomo dos genes significativos
- `CLI_reproducibility_log.txt` — O pergaminho eterno de reprodutibilidade

### A Sessão Espírita Gráfica (GUI)

Para os que preferem seus rituais com botões e controles deslizantes:

```bash
python interface.py
```

O portal gráfico se materializará, permitindo carregar suas matrizes, configurar parâmetros e invocar todos os módulos analíticos com um clique.

---

## 🏚️ A Arquitetura Sagrada

```
BRSA/
│
├── 🧠 BRSA.py              — Ponto de entrada principal: CLI vs GUI
├── 🎛️ interface.py         — Sessão espírita gráfica PyQt5 (GUI)
├── ⚙️ job.py               — FrankensteinBioinformaticsOverlordThread
│                              (Thread de trabalho: executa todas as análises)
├── 🔬 biotools.py          — Wrappers de ferramentas externas e APIs
│                              (FastQC, HISAT2, Kallisto, STRING DB...)
├── 📊 plot.py              — Todas as funções de visualização
│                              (Vulcão, MA, UMAP, Clustermap, GSEA...)
└── 📜 args.py              — Parser CLI e dataclass SacredSessionScroll
```

**Fluxo de dados pelas catacumbas:**

```
Arquivos FASTQ ──► FastQC/fastp ──► HISAT2/Kallisto ──► Matriz de Contagens
                                                               │
Metadados CSV ─────────────────────────────────────────────────┤
                                                               ▼
                                                        Normalização
                                                       (CPM/TPM/DESeq2)
                                                               │
                                            ┌──────────────────┤
                                            ▼                  ▼
                                     Redução Dim.         Teste DE DESeq2
                                  (PCA/t-SNE/UMAP)              │
                                            │                  ▼
                                            │          Vulcão / MA Plot
                                            │          GSEA / ORA
                                            │          Golens de ML
                                            │          Rede PPI
                                            │          Análise de Sobrevivência
                                            └──────► Artefatos HTML Interativos
```

---

## 📚 O Grimório dos Métodos

> *Todo ritual sombrio tem sua justificativa. Aqui jazem os tomos acadêmicos que deram origem a esses métodos.*

---

### I. Controle de Qualidade e Pré-processamento

O **FastQC** é empregado como a primeira inquisição sobre as reads de sequenciamento bruto, avaliando escores de qualidade por base, distribuição de conteúdo GC, contaminação por adaptadores e níveis de duplicação. Os resultados de múltiplas amostras são agregados pelo **MultiQC** em um relatório de vigilância unificado.

A trimagem de adaptadores e a filtragem de qualidade são executadas via **fastp**, que realiza simultaneamente a trimagem por janela deslizante (Phred ≥ 20), detecção automática de adaptadores e geração de relatórios de QC em HTML e JSON.

> Andrews, S. (2010). *FastQC: A Quality Control Tool for High Throughput Sequence Data*. Bioinformatics Group, Babraham Institute. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
>
> Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
>
> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

---

### II. Alinhamento e Quantificação

As reads são alinhadas ao genoma de referência usando o **HISAT2**, um algoritmo de alinhamento baseado em grafos que lida eficientemente com alinhamentos spliced através de junções éxon-éxon, utilizando um Índice FM Hierárquico por Grafos (HGFM). As reads alinhadas são convertidas e ordenadas por coordenadas via **SAMtools**.

As contagens de reads em nível de gene são colhidas pelo **featureCounts** (pacote Subread), que atribui reads alinhadas a features genômicas com base em anotação GTF/GFF3, utilizando uma estratégia de sobreposição por união.

Como alternativa, o **Kallisto** oferece quantificação ultra-rápida de abundância transcricional por pseudo-alinhamento — mapeando reads a um índice de k-mers do transcriptoma de referência sem alinhamento completo, produzindo estimativas de TPM com quantificação de incerteza via bootstrap.

> Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. *Nature Methods*, 12(4), 357–360. https://doi.org/10.1038/nmeth.3317
>
> Li, H., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
>
> Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656
>
> Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

---

### III. Normalização

O BRSA suporta três exorcismos de normalização para remover os demônios de profundidade de sequenciamento e vieses de composição:

- **CPM (Contagens Por Milhão):** Escala as contagens brutas pelo tamanho total da biblioteca, corrigindo para profundidade de sequenciamento. Adequado para análise exploratória e visualização.
- **TPM (Transcritos Por Milhão):** Normaliza primeiro pelo comprimento do gene (RPK), depois pela soma total de RPK. Preserva comparabilidade entre amostras para abundância transcricional.
- **Transformação Estabilizadora de Variância DESeq2 (VST):** Aplica estabilização de variância baseada em modelo binomial negativo, contraindo a dependência média-variância. Recomendado para testes estatísticos subsequentes e agrupamentos.

> Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
>
> Wagner, G. P., Kin, K., & Lynch, V. J. (2012). Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory in Biosciences*, 131(4), 281–285. https://doi.org/10.1007/s12064-012-0162-3

---

### IV. Análise de Expressão Diferencial

A expressão diferencial entre duas condições é testada usando o framework de **modelo linear generalizado binomial negativo**, conforme implementado no paradigma DESeq2. O pipeline:

1. Filtra genes com baixa expressão (padrão: < 10 reads em todas as amostras)
2. Estima fatores de tamanho via normalização pela **mediana das razões**
3. Estima a dispersão gene a gene usando *shrinkage* empírico bayesiano
4. Testa a hipótese nula (Log₂FC = 0) por meio do **teste de Wald**
5. Ajusta os p-valores pelo procedimento de **Benjamini-Hochberg** para controlar a Taxa de Falsa Descoberta (FDR)

Os resultados são visualizados como **Volcano Plots** (Log₂FC vs −log₁₀ FDR) e **Gráficos de Constelação MA** (expressão média vs Log₂FC), com anotações interativas ao passar o cursor, via `mplcursors`.

> Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
>
> Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

---

### V. Redução de Dimensionalidade

Três armas espectrais estão disponíveis para projetar o espaço de expressão de alta dimensionalidade em um plano 2D visível:

- **PCA (Análise de Componentes Principais):** Decomposição linear via decomposição em valores singulares (SVD), preservando a estrutura global de variância. A variância explicada por componente é reportada.
- **t-SNE (t-Distributed Stochastic Neighbor Embedding):** Redução de dimensionalidade não-linear que preserva a estrutura de vizinhança local minimizando a divergência KL entre distribuições de probabilidade de alta e baixa dimensão.
- **UMAP (Uniform Manifold Approximation and Projection):** Aprendizado de variedades que preserva topologia, baseado em geometria Riemanniana e teoria de conjuntos fuzzy simpliciais, equilibrando estrutura local e global com velocidade computacional superior.

> Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. *The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science*, 2(11), 559–572.
>
> van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE. *Journal of Machine Learning Research*, 9, 2579–2605. https://jmlr.org/papers/v9/vandermaaten08a.html
>
> McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv preprint*. https://arxiv.org/abs/1802.03426

---

### VI. Golens de Aprendizado de Máquina

Quando invocado com `--ml`, o BRSA treina dois golens preditivos sobre a matriz de expressão normalizada:

- **Random Forest:** Um ensemble de árvores de decisão descorrelacionadas, treinadas via bootstrap aggregation (bagging). As importâncias de features são calculadas como a diminuição média de impureza (importância Gini) em todas as árvores, ranqueando genes pelo seu poder preditivo para classificação de condições.
- **Support Vector Machine (SVM):** Um classificador de margem máxima usando kernel de função de base radial (RBF), projetando os dados em um espaço de features de maior dimensão para encontrar o hiperplano separador ótimo.

Os escores de importância de features do Random Forest são visualizados como gráficos de barras horizontais, permitindo a descoberta de biomarcadores.

> Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
>
> Cortes, C., & Vapnik, V. (1995). Support-vector networks. *Machine Learning*, 20(3), 273–297. https://doi.org/10.1007/BF00994018

---

### VII. Análise de Redes de Co-expressão

O agrupamento hierárquico dos genes mais variáveis é realizado e visualizado como um **clustermap** (mapa de calor com dendrogramas duplos), usando ligação de Ward sobre distâncias Euclidianas. O agrupamento em nível de amostras revela efeitos de batch e agrupamentos biológicos, enquanto o agrupamento em nível de genes expõe módulos transcricionais co-regulados.

Esta abordagem segue os princípios da **Análise de Redes de Co-expressão Gênica Ponderada (WGCNA)**, em que genes com perfis de expressão correlacionados entre amostras são agrupados em módulos funcionais.

> Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559
>
> Ward, J. H. (1963). Hierarchical Grouping to Optimize an Objective Function. *Journal of the American Statistical Association*, 58(301), 236–244. https://doi.org/10.2307/2282967

---

### VIII. Interatoma Proteína-Proteína

Quando `--ppi` é ativado, o BRSA consulta a API REST do banco de dados **STRING** para recuperar dados de interação proteína-proteína para a lista de genes diferencialmente expressos significativos. As interações são filtradas por um limiar de escore de confiança (padrão: 0,4) e montadas em um grafo não-direcionado **NetworkX**, renderizado como uma rede com layout de força dirigida (*spring layout*) — a Bola de Cabelo do Caos.

> Szklarczyk, D., et al. (2023). The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any of 14000+ organisms. *Nucleic Acids Research*, 51(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000
>
> Hagberg, A. A., Schult, D. A., & Swart, P. J. (2008). Exploring Network Structure, Dynamics, and Function using NetworkX. In *Proceedings of the 7th Python in Science Conference (SciPy2008)*, pp. 11–15.

---

### IX. Análise de Sobrevivência

O BRSA implementa a **estimativa da curva de sobrevivência de Kaplan-Meier** para correlacionar os níveis de expressão gênica individual com dados de desfecho clínico dos pacientes. As amostras são estratificadas em grupos de expressão "Alta" e "Baixa" com base na mediana de expressão. As funções de sobrevivência são estimadas de forma não-paramétrica para cada grupo e plotadas como funções em degrau com intervalos de confiança.

Esta análise requer colunas de metadados clínicos para tempo de sobrevivência e indicador de ocorrência do evento (indicador de censura).

> Kaplan, E. L., & Meier, P. (1958). Nonparametric Estimation from Incomplete Observations. *Journal of the American Statistical Association*, 53(282), 457–481. https://doi.org/10.2307/2281868
>
> Davidson-Pilon, C. (2019). lifelines: survival analysis in Python. *Journal of Open Source Software*, 4(40), 1317. https://doi.org/10.21105/joss.01317

---

### X. Análise de Enriquecimento de Conjuntos Gênicos

O BRSA emprega o **GSEApy** para realizar tanto:

- **GSEA (Análise de Enriquecimento de Conjuntos Gênicos):** Ranqueia todos os genes detectados pela sua mudança de fold (Log₂FC) com sinal e testa o enriquecimento coordenado de conjuntos gênicos curados nas extremidades da lista ranqueada, calculando um Escore de Enriquecimento Progressivo (RES) e escore de enriquecimento normalizado (NES).
- **Análise de Sobre-Representação (ORA):** Testa a sobreposição significativa entre a lista de genes significativos e conjuntos gênicos anotados usando o teste exato de Fisher com correção FDR.

Múltiplos bancos de dados de conjuntos gênicos são suportados, incluindo GO Processo Biológico, KEGG, Reactome e MSigDB Hallmarks.

> Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102
>
> Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*, 39(1), btac757. https://doi.org/10.1093/bioinformatics/btac757

---

## ⚰️ Os Parâmetros Eldritchianos (Referência CLI)

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--counts` | `str` | *obrigatório* | Caminho para a matriz de contagens CSV (genes × amostras) |
| `--meta` | `str` | *obrigatório* | Caminho para o CSV de metadados (amostras × covariáveis) |
| `--cond1` | `str` | *obrigatório* | Condição de referência (Grupo A / Alpha) |
| `--cond2` | `str` | *obrigatório* | Condição de teste (Grupo B / Omega) |
| `--min_counts` | `int` | `10` | Limiar mínimo de reads para filtragem de genes |
| `--pval` | `float` | `0.05` | Limiar de FDR para significância |
| `--norm` | `str` | `DESeq2` | Normalização: `CPM`, `TPM` ou `DESeq2` |
| `--batch` | `str` | `None` | Coluna de metadados para correção de efeito de batch |
| `--dim_red` | `str` | `PCA` | Redução de dimensionalidade: `PCA`, `t-SNE`, `UMAP` |
| `--ml` | flag | `False` | Invoca golens de ML (Random Forest + SVM) |
| `--ppi` | flag | `False` | Invoca a rede proteica do STRING DB |
| `--isoforms` | flag | `False` | Caça eventos de splicing alternativo |

### Formatos de Entrada Esperados

**Matriz de Contagens (`--counts`):**
```
gene_id,Amostra_1,Amostra_2,Amostra_3,...
ENSG00000000003,1500,2300,800,...
ENSG00000000005,45,12,67,...
```

**Metadados (`--meta`):**
```
sample,condition,batch
Amostra_1,Tumor,batch_A
Amostra_2,Normal,batch_A
Amostra_3,Tumor,batch_B
```

---

## 🏺 Os Artefatos do Caldeirão

Após a conclusão bem-sucedida, os seguintes artefatos sombrios serão forjados:

```
outputs/
├── 📄 CLI_differential_results.csv         — Tabela completa de genes DE (LFC, p-adj, significância)
├── 📄 CLI_reproducibility_log.txt          — Registro de parâmetros + snapshot do pip freeze
│
└── 📁 interactive_artifacts/
    ├── 🌐 interactive_pca.html             — PCA/t-SNE/UMAP interativo com Plotly
    ├── 🌐 interactive_volcano.html         — Volcano interativo com anotações gênicas
    └── 🌐 interactive_ma_constellation.html — Gráfico MA interativo
```

**Colunas da tabela de resultados DE:**

| Coluna | Descrição |
|--------|-----------|
| `gene` | Identificador do gene |
| `log2fc` | Log₂ Fold Change (cond2 vs cond1) |
| `p-adj` | P-valor ajustado por Benjamini-Hochberg |
| `mean_1` | Expressão normalizada média em cond1 |
| `mean_2` | Expressão normalizada média em cond2 |
| `significant` | Booleano: True se FDR < limiar |

---

## 📜 Pergaminhos de Reprodutibilidade

Todo ritual do BRSA automaticamente esculpe um registro de reprodutibilidade contendo:

- Todos os parâmetros e limiares aplicados (estado completo do `SacredSessionScroll`)
- A saída completa do `pip freeze`, capturando versões exatas das bibliotecas
- Contexto de execução e informações do ambiente computacional

Isso garante que as condições analíticas exatas possam ser reconstruídas por qualquer necromante futuro — ou revisor de periódico — que tente replicar os resultados.

---

## 📖 Referências e Tomos Ancestrais

1. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Ewels, P., et al. (2016). MultiQC. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
3. Chen, S., et al. (2018). fastp. *Bioinformatics*, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560
4. Kim, D., et al. (2015). HISAT. *Nature Methods*, 12(4), 357–360. https://doi.org/10.1038/nmeth.3317
5. Li, H., et al. (2009). SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
6. Liao, Y., et al. (2014). featureCounts. *Bioinformatics*, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656
7. Bray, N. L., et al. (2016). Kallisto. *Nature Biotechnology*, 34(5), 525–527. https://doi.org/10.1038/nbt.3519
8. Love, M. I., et al. (2014). DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
9. Benjamini, Y., & Hochberg, Y. (1995). FDR. *JRSS-B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
10. van der Maaten, L., & Hinton, G. (2008). t-SNE. *JMLR*, 9, 2579–2605. https://jmlr.org/papers/v9/vandermaaten08a.html
11. McInnes, L., et al. (2018). UMAP. *arXiv*. https://arxiv.org/abs/1802.03426
12. Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
13. Cortes, C., & Vapnik, V. (1995). SVM. *Machine Learning*, 20(3), 273–297. https://doi.org/10.1007/BF00994018
14. Langfelder, P., & Horvath, S. (2008). WGCNA. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559
15. Szklarczyk, D., et al. (2023). STRING 2023. *Nucleic Acids Research*, 51(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000
16. Kaplan, E. L., & Meier, P. (1958). Kaplan-Meier. *JASA*, 53(282), 457–481. https://doi.org/10.2307/2281868
17. Davidson-Pilon, C. (2019). lifelines. *JOSS*, 4(40), 1317. https://doi.org/10.21105/joss.01317
18. Subramanian, A., et al. (2005). GSEA. *PNAS*, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102
19. Fang, Z., et al. (2023). GSEApy. *Bioinformatics*, 39(1), btac757. https://doi.org/10.1093/bioinformatics/btac757

---

## ⚖️ Licença e Pacto das Trevas

Este software é distribuído sob a **Licença MIT**. Consulte o arquivo `LICENSE` para detalhes.

Ao utilizar o BRSA, você implicitamente concorda em citar as referências metodológicas relevantes em quaisquer publicações resultantes — que seus manuscritos não sejam amaldiçoados pelo Revisor #3.

---

## 🧙 O Necromante (Autor)

**Victor S. Caricatte De Araújo (Frankestein)**  
*Universidade Federal de Minas Gerais (UFMG)*  
Versão: `0.9.0` 

---

<div align="center">

*"No cemitério dos erros de sequenciamento e efeitos de batch, o BRSA acende a lanterna da expressão diferencial."*

🕯️ *Que seus p-valores sejam baixos e seus tamanhos amostrais sejam altos.* 🕯️

</div>
