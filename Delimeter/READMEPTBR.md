# Filtro de Arquivos Delimitados para Arquivos Genéticos

**Versão:** 1.1.2
**Data:** 11 de Abril
**Autor:** Victor S Caricatte De Araújo
**Email:** victorleniwys@gmail.com ou victorsc@ufmg.br
**Instituição:** Universidade Federal de Minas Gerais (UFMG)
**GitHub:** [VictorCaricatte](https://github.com/VictorCaricatte)

## Descrição

Uma aplicação gráfica que lê diversos formatos de arquivo (CSV, TSV, Excel) e permite filtrar linhas com base em valores específicos em uma determinada coluna. Possui uma interface com abas, controles de filtragem e uma seção de ajuda abrangente. Ideal para manipulação e filtragem de dados, especialmente útil em contextos de análise de dados genéticos ou tabulares em geral.

## Funcionalidades Principais

* **Suporte a Múltiplos Formatos:** Lê arquivos CSV, TSV e Excel (.xlsx, .xls).
* **Interface com Abas:** Organiza as funcionalidades em "Filtro" e "Ajuda".
* **Múltiplas Condições de Filtragem:**
    * Igual a (Equals)
    * Contém (Contains)
    * Começa com (Starts With)
    * Termina com (Ends With)
    * Expressão Regular (Regex)
* **Opções de Sensibilidade:** Permite filtragem sensível ou insensível a maiúsculas/minúsculas.
* **Inverter Correspondência:** Exclui as linhas que correspondem aos critérios de filtro.
* **Manipulação de Cabeçalho:** Opção para tratar a primeira linha como cabeçalho.
* **Pré-visualização dos Resultados:** Exibe uma prévia dos dados filtrados (até 1000 linhas).
* **Exportação Flexível:** Exporta os resultados filtrados para CSV, TSV ou Excel.
* **Alternância de Tema:** Suporte a tema claro (light) e escuro (dark) para melhor visualização.

## Pré-requisitos

Para executar este software, você precisará ter o Python 3 instalado, juntamente com as seguintes bibliotecas:

* pandas
* PyQt5
* openpyxl (para suporte a arquivos Excel)

## Instalação

1.  Clone ou baixe este repositório.
2.  Navegue até o diretório do projeto.
3.  Instale as dependências usando pip:

    ```bash
    pip install pandas PyQt5 openpyxl
    ```

## Como Usar

1.  Execute o script `delimitergenfile.py`:

    ```bash
    python delimitergenfile.py
    ```

2.  **Na aba "Filtro":**
    * **Seleção de Arquivo:**
        * Clique em "Procurar..." para selecionar seu arquivo de dados.
        * O tipo de arquivo pode ser detectado automaticamente ("Auto Detect") ou especificado (CSV, TSV, Excel).
        * Clique em "Carregar Arquivo" para carregar os dados.
    * **Configurações do Filtro:**
        * **Coluna:** Selecione a coluna na qual o filtro será aplicado (habilitado após carregar o arquivo).
        * **Operação:** Escolha o tipo de operação de filtragem (Igual a, Contém, etc.).
        * **Valor:** Digite o valor a ser procurado/comparado.
        * **Opções:**
            * `Sensível a Maiúsculas`: Marque para diferenciar maiúsculas de minúsculas.
            * `Inverter Correspondência`: Marque para obter as linhas que *não* correspondem ao filtro.
            * `Possui Cabeçalho`: Marque se o seu arquivo tem uma linha de cabeçalho (marcado por padrão).
    * **Ações:**
        * Clique em "Aplicar Filtro" para visualizar os resultados na tabela "Pré-visualização dos Resultados".
        * Clique em "Resetar Filtro" para limpar o filtro aplicado e mostrar todos os dados originais.
    * **Exportar Resultados:**
        * Após filtrar, clique em "Exportar para CSV", "Exportar para TSV" ou "Exportar para Excel" para salvar os dados filtrados.

3.  **Na aba "Ajuda":**
    * Encontre um guia rápido, descrição das opções de filtro, formatos suportados, dicas e informações de contato.

4.  **Menu "View" (Visualizar):**
    * Permite alternar entre o tema escuro (padrão) e o tema claro.

## Dicas

* Para arquivos grandes, a aplicação de filtros complexos pode levar algum tempo.
* As operações de Regex utilizam a sintaxe completa de expressões regulares do Python.
* Você pode pré-visualizar os resultados antes de exportar.

## Contribuições

Contribuições são bem-vindas! Sinta-se à vontade para abrir uma *issue* ou enviar um *pull request*.

## Licença

Este projeto possui uma licença MIT até o momento.
