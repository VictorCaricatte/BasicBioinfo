# LengthCounter.py - Analisador de Sequências FASTA/FASTQ

**Autor:** Victor S. Caricatte De Araújo
**Email:** victorleniwys@gmail.com ou victorsc@ufmg.br
**Instituição:** Universidade Federal de Minas Gerais (UFMG)
**Versão:** 1.1.9
**Data da Versão:** 28 de Abril

## Visão Geral

`LengthAnalay.py` é uma ferramenta gráfica (GUI) e de linha de comando (CLI) desenvolvida em Python para realizar análises em arquivos de sequência biológica. Suas principais funcionalidades incluem:

1.  **Análise de Comprimento de Sequências (Arquivos FASTA):**
    * Carrega múltiplos arquivos FASTA.
    * Calcula estatísticas descritivas do comprimento das sequências (mínimo, máximo, média, mediana, desvio padrão, quartis).
    * Gera visualizações como histogramas (com opção de KDE) e boxplots para comparar a distribuição dos comprimentos.
2.  **Análise de Comprimento de Sequências (Arquivos FASTA):**
    * Permite a inserção manual de uma sequência de DNA/RNA.
    * Calcula a contagem e frequência de cada base (A, T, C, G, U, N).
    * Gera um gráfico de barras da composição de bases.
    * Permite a personalização das cores de cada base no gráfico.
3.  **Estimativa de Cobertura Genômica (Arquivos FASTQ):**
    * Carrega múltiplos arquivos FASTQ.
    * Calcula a cobertura genômica estimada com base no número de leituras, tamanho do genoma, comprimento da leitura e tipo de sequenciamento (single-end ou paired-end).

A aplicação utiliza bibliotecas como PyQt5 para a interface gráfica, Matplotlib e Seaborn para a geração de gráficos, Biopython para o parsing de arquivos FASTA/FASTQ, e NumPy para cálculos estatísticos.

## Funcionalidades Principais

* **Interface Gráfica Intuitiva:** Fácil de usar, com abas para diferentes tipos de análise.
* **Análise de Comprimento (FASTA):**
    * Adicionar e remover múltiplos arquivos FASTA.
    * Atribuir rótulos personalizados para cada conjunto de dados.
    * Escolher entre histograma e boxplot para visualização.
    * Ajustar o número de bins e a exibição de KDE para histogramas.
    * Visualizar estatísticas detalhadas.
    * Salvar gráficos (PNG, JPG, PDF, SVG) e estatísticas (TXT).
* **Análise de Comprimento de Sequências (FASTA):**
    * Entrada de sequência de DNA/RNA diretamente na interface.
    * Visualização gráfica da composição de bases.
    * Personalização das cores das bases no gráfico.
    * Exibição de estatísticas de composição (contagens e percentuais).
    * Salvar gráfico de composição.
* **Estimativa de Cobertura (FASTQ):**
    * Adicionar e remover múltiplos arquivos FASTQ.
    * Configurar o tamanho do genoma, comprimento da leitura e tipo de sequenciamento (paired-end/single-end).
    * Visualizar os resultados da estimativa de cobertura.
    * Salvar resultados (TXT).
* **Aba de Ajuda:** Documentação integrada explicando o uso de cada funcionalidade e a fórmula de cálculo da cobertura.
* **Interface de Linha de Comando (CLI):** Para análises rápidas de comprimento de sequências FASTA sem a necessidade de abrir a GUI.

## Pré-requisitos

* Python 3.x
* Bibliotecas Python:
    * PyQt5
    * Matplotlib
    * Seaborn
    * NumPy
    * Biopython

## Instalação

1.  Clone o repositório ou baixe o arquivo `LengthCounter.py`.
2.  Instale as dependências necessárias. Você pode usar o pip:

    ```bash
    pip install PyQt5 matplotlib seaborn numpy biopython
    ```

## Uso

### Interface Gráfica (GUI)

Para iniciar a aplicação no modo GUI, execute o script sem argumentos:

```bash
python lenghtanaly.py
