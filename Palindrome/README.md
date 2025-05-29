# DNA Palindrome Finder

**Versão:** 1.1
**Autor:** Victor S Caricatte De Araújo
**Email:** victorleniwys@gmail.com ou victorsc@ufmg.br
**Instituição:** Universidade Federal de Minas Gerais (UFMG)
**Data:** Abr, 14

## Sumário

1.  [Descrição](#descrição)
2.  [O que são Palíndromos de DNA?](#o-que-são-palíndromos-de-dna)
3.  [Funcionalidades](#funcionalidades)
4.  [Como Usar](#como-usar)
5.  [Interface do Usuário](#interface-do-usuário)
    * [Aba Principal (Análise de DNA)](#aba-principal-análise-de-dna)
    * [Aba de Ajuda](#aba-de-ajuda)
    * [Menu](#menu)
6.  [Requisitos](#requisitos)
7.  [Instalação de Dependências](#instalação-de-dependências)
8.  [Como Executar](#como-executar)
9.  [Detalhes Técnicos](#detalhes-técnicos)
10. [Sobre o Autor](#sobre-o-autor)

## Descrição

O **DNA Palindrome Finder** é uma ferramenta de bioinformática desenvolvida para identificar sequências palindrômicas em moléculas de DNA. Palíndromos de DNA são sequências que leem o mesmo em uma fita e em sua fita complementar reversa. Essas sequências são biologicamente significativas, pois frequentemente constituem sítios de reconhecimento e clivagem para enzimas de restrição.

Esta aplicação fornece uma interface gráfica amigável para facilitar a análise de sequências de DNA, permitindo que usuários insiram sequências manualmente, carreguem-nas de arquivos FASTA e configurem parâmetros para a busca de palíndromos.

## O que são Palíndromos de DNA?

Em biologia molecular, um palíndromo de DNA é uma sequência de nucleotídeos que é idêntica à sua sequência complementar reversa.
Por exemplo, a sequência 5'-GAATTC-3' é um palíndromo porque sua fita complementar é 3'-CTTAAG-5'. Se lermos esta fita complementar reversa (da direita para a esquerda, no sentido 5' para 3'), obtemos 5'-GAATTC-3', que é idêntica à sequência original.
Este exemplo específico (GAATTC) é o sítio de reconhecimento para a enzima de restrição EcoRI.

## Funcionalidades

* **Entrada de Sequência Flexível:**
    * Entrada manual de sequências de DNA na área de texto.
    * Carregamento de sequências a partir de arquivos no formato FASTA (`.fasta`, `.fa`, `.faa`, `.fna`).
* **Parâmetros de Busca Configuráveis:**
    * Definição do comprimento mínimo e máximo dos palíndromos a serem identificados.
* **Validação de Sequência:**
    * Verifica se a sequência de entrada contém apenas bases válidas (A, T, C, G), ignorando espaços e quebras de linha.
* **Análise Eficiente:**
    * Algoritmo otimizado para busca de palíndromos.
    * Uso de threading para evitar que a interface gráfica congele durante análises longas.
* **Visualização de Resultados:**
    * Apresentação dos palíndromos encontrados em uma tabela com colunas: Posição (início-fim), Sequência Palindrômica e Comprimento.
* **Gerenciamento de Resultados:**
    * Copiar resultados para a área de transferência.
    * Salvar resultados em arquivos de texto (`.txt`) ou CSV (`.csv`).
    * Limpar campos de entrada e resultados para nova análise.
* **Interface Amigável:**
    * Interface gráfica intuitiva construída com Tkinter.
    * Aba de "Ajuda" com documentação e informações sobre o programa.
    * Barra de status para feedback ao usuário.

## Como Usar

1.  **Inserir Sequência de DNA:**
    * Digite ou cole a sequência de DNA diretamente na área de texto "DNA Input".
    * Alternativamente, vá em `Menu > File > Open FASTA file...` para carregar uma sequência de um arquivo FASTA. Apenas caracteres A, T, C, G serão considerados (maiúsculas ou minúsculas).
2.  **Definir Parâmetros:**
    * Ajuste os campos "Min length" (comprimento mínimo) e "Max length" (comprimento máximo) para os palíndromos que deseja encontrar. O padrão é 4 e 12, respectivamente.
3.  **Analisar:**
    * Clique no botão "Analyze DNA".
    * A barra de status indicará o progresso da análise.
4.  **Visualizar Resultados:**
    * Os palíndromos encontrados serão listados na tabela "Results", mostrando sua posição (1-baseado), a sequência e o comprimento.
5.  **Gerenciar Resultados (Opcional):**
    * **Copy Results:** Copia os resultados da tabela para a área de transferência.
    * **Save Results:** Salva os resultados em um arquivo de texto ou CSV.
    * **Clear All:** Limpa a área de entrada de DNA, a tabela de resultados e redefine os comprimentos para os valores padrão.

## Interface do Usuário

A interface principal é dividida em abas e um menu.

### Aba Principal (Análise de DNA)

* **DNA Input:** Uma área de texto para inserir ou visualizar a sequência de DNA.
* **Settings Frame:**
    * **Min length:** Campo para definir o comprimento mínimo do palíndromo.
    * **Max length:** Campo para definir o comprimento máximo do palíndromo.
    * **Analyze DNA:** Botão para iniciar a busca por palíndromos.
* **Results Frame:**
    * **Tabela de Resultados:** Exibe os palíndromos encontrados com colunas "Position", "Sequence", e "Length".
    * **Copy Results:** Botão para copiar os dados da tabela.
    * **Save Results:** Botão para salvar os dados da tabela em um arquivo.
    * **Clear All:** Botão para limpar todos os campos e resultados.
* **Status Bar:** Localizada na parte inferior da janela, exibe o status atual da aplicação.

### Aba de Ajuda

Contém informações sobre o programa, como usá-lo, a definição técnica de palíndromos de DNA e detalhes sobre o desenvolvedor.

### Menu

* **File:**
    * **Open FASTA file...:** Abre um diálogo para selecionar e carregar um arquivo FASTA.
    * **Save results...:** Abre um diálogo para salvar os resultados atuais.
    * **Exit:** Fecha a aplicação.
* **Help:**
    * **Documentation:** Seleciona a aba de Ajuda.
    * **About:** Exibe uma janela com informações sobre a versão e o propósito do software.

## Requisitos

* Python 3.x
* Tkinter (geralmente incluído na instalação padrão do Python)
* Biopython (a biblioteca `Bio.Seq` é importada, embora a lógica de verificação de palíndromo seja customizada no script)

## Instalação de Dependências

Se você não tiver o Biopython instalado, pode instalá-lo usando pip:

```bash
pip install biopython
