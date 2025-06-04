# BRSA - Basic RNA-Seq Analysis

**Autor**: Victor Silveira Caricatte  
**Instituição**: Universidade Federal de Minas Gerais (UFMG)

## Visão Geral

O **BRSA (Basic RNA-Seq Analysis)** é uma aplicação gráfica desenvolvida com PyQt6 que permite a análise de dados de RNA-seq de forma acessível, reprodutível e automatizada. Ele realiza desde a importação de dados, normalização, análise de expressão diferencial, até visualizações gráficas e análises de enriquecimento funcional.

---

## Funcionalidades

### 📂 Importação de Arquivos
- Arquivo de contagem de genes (formato CSV): genes nas linhas, amostras nas colunas.
- Arquivo de metadados (formato CSV): deve conter colunas `sample` e `condition`.

### ⚙️ Parâmetros de Análise
- Seleção de duas condições experimentais para comparação.
- Filtro de genes por contagem mínima.
- Corte de valor de p ajustado por FDR.

### 🔬 Etapas da Análise
1. **Pré-processamento**:
   - Filtragem de genes de baixa expressão.
   - Verificação de integridade dos dados.

2. **Normalização**:
   - Usando CPM (Counts Per Million).

3. **Análise de Expressão Diferencial**:
   - Teste t de Student.
   - Correção de múltiplas comparações (FDR - Benjamini-Hochberg).
   - Cálculo de log2 Fold Change.

4. **Visualizações**:
   - PCA (Análise de Componentes Principais).
   - Heatmap dos 50 genes mais variáveis.
   - Volcano plot interativo com `mplcursors`.

5. **Análise de Enriquecimento Funcional**:
   - GO (BP, MF, CC)
   - KEGG
   - Reactome
   - Disease Ontology
   - Utiliza `gprofiler`.

6. **Exportação de Resultados**:
   - 📄 PDF com gráficos e tabelas (via `reportlab`).
   - 📊 Planilha Excel com abas para resultados DE, contagens normalizadas e enriquecimento (via `openpyxl`).

7. **Gerenciamento de Sessão**:
   - Salvar/Carregar sessões de análise com arquivos `.json`.

---

## Interface Gráfica

O BRSA conta com uma GUI dividida em abas:

- **Input**: seleção de arquivos e parâmetros.
- **Results**:
  - `Differential Expression`: resumo dos genes significativos.
  - `Visualization`: seleção e exibição dos gráficos.
  - `Enrichment Analysis`: execução e exibição da análise funcional.
- **Help**: instruções e informações sobre o uso da aplicação.

---

## Requisitos

### 🐍 Python
- Versão: Python 3.7 ou superior

### 📦 Bibliotecas

Instale todas as dependências com:

```bash
pip install -r requirements.txt
