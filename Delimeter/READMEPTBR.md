# 🧬 DELIMITER — Filtro de Arquivos Delimitados Quântico
### *Um Acelerador de Partículas para Dados Genômicos. Um Reator. Uma Fera.*

> **"Aqui não processamos arquivos. Nós os detonamos."**
> — *Supremo Overlord Genômico, Sessão Iniciada*

---

```
╔══════════════════════════════════════════════════════════════════════╗
║   ██████  ███████ ██      ██ ███    ███ ██ ████████ ███████ ██████  ║
║   ██   ██ ██      ██      ██ ████  ████ ██    ██    ██      ██   ██ ║
║   ██   ██ █████   ██      ██ ██ ████ ██ ██    ██    █████   ██████  ║
║   ██   ██ ██      ██      ██ ██  ██  ██ ██    ██    ██      ██   ██ ║
║   ██████  ███████ ███████ ██ ██      ██ ██    ██    ███████ ██   ██ ║
╚══════════════════════════════════════════════════════════════════════╝
          REATOR DE FILTRAGEM QUÂNTICA — CONTENÇÃO ATIVA
```

[![Versão](https://img.shields.io/badge/versão-0.3.1-brightgreen?style=for-the-badge)](#)
[![Python](https://img.shields.io/badge/python-3.8+-blue?style=for-the-badge&logo=python)](#)
[![PyQt5](https://img.shields.io/badge/UI-PyQt5-orange?style=for-the-badge)](#)
[![Licença](https://img.shields.io/badge/licença-MIT-purple?style=for-the-badge)](#)
[![Reator](https://img.shields.io/badge/reator-ESTÁVEL-success?style=for-the-badge)](#)
[![Biossegurança](https://img.shields.io/badge/NSF%2FANSI%2049-CONFORME-red?style=for-the-badge)](#)

---

## ☢️ AVISO: ÁREA DO REATOR — SOMENTE PESSOAL AUTORIZADO

**Delimiter** é uma aplicação desktop para **inspeção, filtragem, transformação e síntese** de arquivos de dados delimitados (CSV, TSV, Excel, FASTA). Construído sobre PyQt5 e Pandas, oferece um pipeline de processamento de dados de alta performance e multi-thread, embrulhado em uma interface que trata cada clique de botão como um evento nuclear.

Desenvolvido na **Universidade Federal de Minas Gerais (UFMG)** por **Victor S. Caricatte De Araújo** (`victorsc@ufmg.br`), o Delimiter foi concebido com pipelines de bioinformática em mente — mas seu reator é poderoso o suficiente para processar *qualquer* dado tabular estruturado.

---

## 🗺️ Sumário

1. [Arquitetura do Núcleo](#-arquitetura-do-núcleo)
2. [Pré-requisitos](#%EF%B8%8F-pré-requisitos)
3. [Instalação](#-instalação--sequência-de-ignição-do-reator)
4. [Estrutura do Projeto](#-estrutura-do-projeto)
5. [Funcionalidades — Subsistemas do Reator](#-funcionalidades--subsistemas-do-reator)
6. [Protocolos de Síntese para Exportação](#-protocolos-de-síntese-para-exportação)
7. [Internacionalização — O Ribossomo Pedra de Rosetta](#-internacionalização--o-ribossomo-pedra-de-rosetta)
8. [Temas — Mudança do Espectro Quântico](#-temas--mudança-do-espectro-quântico)
9. [Logs — Runas de Biocontenção](#-logs--runas-de-biocontenção)
10. [Anomalias Conhecidas e Solução de Problemas](#-anomalias-conhecidas-e-solução-de-problemas)
11. [Diretrizes de Biossegurança](#-diretrizes-de-biossegurança)
12. [Autor e Atribuições](#-autor-e-atribuições)

---

## 🔬 Arquitetura do Núcleo

A aplicação é organizada em torno de quatro módulos — cada um cumprindo um papel específico no pipeline do reator:

| Módulo | Codinome | Função |
|---|---|---|
| `delimeter.py` | **Acelerador de Partículas** | Ponto de entrada principal. Orquestra todas as threads, sinais e eventos de UI. Instancia o `Supreme_Genomic_Overlord`. |
| `Interface.py` | **Holodeck Genômico** | Interface gráfica completa em PyQt5. Todos os widgets, diálogos e definições de layout. |
| `logic.py` | **Motor de Dados Esterilizado** | Toda a lógica de processamento de dados. Filtragem, merge, transmutação de colunas, irradiação de NaN, deduplicação e exportação. |
| `xenoglossy_codex.py` | **Ribossomo Pedra de Rosetta** | Motor de internacionalização. Armazena e recupera strings da UI em EN, PT-BR e ES. |

A classe `Supreme_Genomic_Overlord` em `delimeter.py` atua como a **barra de controle central** do reator — conecta todos os sinais PyQt5 da UI aos seus métodos correspondentes no motor de dados, gerencia workers `QThread` e impede a aplicação de ir ao estado crítico.

---

## 🛠️ Pré-requisitos

Antes de ignitar o reator, certifique-se de que as seguintes partículas estão disponíveis no seu ambiente:

- **Python** ≥ 3.8
- **PyQt5** — Framework de UI
- **pandas** — Motor de dados tabulares
- **openpyxl** — Síntese de arquivos Excel
- **psutil** *(opcional, mas recomendado)* — Detecção de sobrecarga de RAM antes da injeção

---

## 🚀 Instalação — Sequência de Ignição do Reator

**Passo 1 — Clone o repositório genético:**
```bash
git clone https://github.com/VictorCaricatte/BasicBioinfo.git
cd BasicBioinfo/Delimeter/src
```

**Passo 2 — Instale as dependências:**
```bash
pip install PyQt5 pandas openpyxl psutil
```

**Passo 3 — Verifique se o diretório `themes/`** existe com os arquivos `dark_theme.qss` e `light_theme.qss`. O reator inicia no modo **Matéria Escura** por padrão.

**Passo 4 — Inicie o Overlord Supremo:**
```bash
python delimeter.py
```

Se uma falha catastrófica ocorrer na inicialização, um arquivo `CRASH_REPORT_BIOLAB.txt` será sintetizado no diretório de trabalho. Consulte-o para obter o traceback completo antes de registrar um relatório de brecha de contenção.

---

## 📁 Estrutura do Projeto

```
Delimeter/src/
│
├── delimeter.py              # ⚛️  Núcleo principal do reator — ponto de entrada
├── Interface.py              # 🖥️  Holodeck Genômico (PyQt5)
├── logic.py                  # 🧠  Motor de Dados Esterilizado & threads trabalhadoras
├── xenoglossy_codex.py       # 🌐  Internacionalização — Ribossomo Pedra de Rosetta
│
├── themes/
│   ├── dark_theme.qss        # 🌙  Stylesheet Matéria Escura
│   └── light_theme.qss       # ☀️  Stylesheet Fóton Solar
│
├── panvita_quarantine.log    # 📋  Log de runtime — Runas de Biocontenção
└── CRASH_REPORT_BIOLAB.txt   # 💥  Gerado somente em falha catastrófica
```

---

## ⚡ Funcionalidades — Subsistemas do Reator

### 🧬 Injeção de Biomassa (Carregamento de Arquivos)
- Suporta os formatos **CSV**, **TSV** e **Excel** (`.xlsx`).
- **Detecção automática** de delimitadores via `csv.Sniffer` — o `DarkMatter_Sniffer` sonda os primeiros 2048 bytes do arquivo e determina o separador correto automaticamente.
- Arquivos podem ser carregados arrastando para o campo de caminho ou usando o botão **🧬 Injetar**.
- A detecção de linha de cabeçalho é configurável por injeção.
- Antes do carregamento, um **Alerta de Sobrecarga de Massa** é acionado se o `psutil` detectar que o tamanho do arquivo × 3 excede a RAM disponível — dando ao operador a chance de abortar antes que o reator colapse.
- O carregamento de arquivos roda em uma **`QThread` dedicada** (`Hyperdrive_Loader_Thread`) — a UI nunca trava durante a injeção, independentemente do tamanho do material biológico.

---

### 🎛️ Reator de Filtros Matryoshka
A funcionalidade mais marcante do Delimiter. Os filtros são estruturados como **pods de condição empilhados e aninhados** — como uma boneca russa Matryoshka — cada um estreitando ainda mais o conjunto de dados.

Cada pod de filtro expõe a seguinte configuração:

| Parâmetro | Descrição |
|---|---|
| **Coluna** | A fita genômica (coluna) a ser visada |
| **Operação** | O operador de comparação (ver tabela abaixo) |
| **Valor** | O valor-alvo biológico |
| **Sensível a Maiúsculas** | Alterna a correspondência de string com distinção de maiúsculas/minúsculas |
| **Ignorar NaN** | Exclui automaticamente linhas nulas/NaN do passe de filtragem |
| **Somente NaN** | Isola e retorna apenas as linhas onde a coluna é nula |

**Operadores Quânticos Suportados:**

| Operador | Descrição |
|---|---|
| `Equals` | Igualdade estrita de string |
| `Contains` | Verificação de presença de substring |
| `Starts With` | Correspondência de prefixo |
| `Ends With` | Correspondência de sufixo |
| `Regex` | Correspondência por expressão regular completa (validada antes da detonação) |
| `Greater Than (>)` | Comparação numérica |
| `Less Than (<)` | Comparação numérica |
| `Numeric Equals (=)` | Igualdade numérica exata |

Múltiplos pods podem ser adicionados dinamicamente via **➕ Adicionar Regra**. Todos os filtros ativos são aplicados sequencialmente quando **⚡ Detonar Colisor** é pressionado. A colisão roda em uma `QThread` separada (`Quantum_Collider_Thread`) e reporta estatísticas de sobrevivência no painel Raio-X Biológico ao término.

---

### 💾 Cristalização e Ressurreição de Presets
Configurações de filtros podem ser salvas e recarregadas entre sessões:

- **💾 Cristalizar Memória** — Serializa todos os pods Matryoshka ativos em um arquivo `.json` via `Elephant_Brain_Storage`.
- **📂 Ressuscitar Filtros** — Desserializa e repopula a UI com um conjunto de filtros salvo anteriormente. Todos os pods existentes são vaporizados antes da ressurreição para evitar contaminação.

---

### 💀 Purgar Clones (Deduplicação)
A operação **💀 Purgar Clones** remove todas as linhas duplicadas do conjunto de dados ativo em todas as colunas. Executada de forma síncrona contra o `Sterilized_Data_Engine` e dispara uma atualização completa da UI.

---

### 🧬 Splicing Alienígena (Merge de Datasets)
Injete material genético estrangeiro — carregue um arquivo secundário e faça o merge com o dataset ativo usando uma coluna-chave. Esta é uma operação padrão de `inner join` do Pandas, encapsulada em um diálogo dedicado (`Xenomorph_Splicer_Dialog`). A abominação resultante substitui o dataset de trabalho atual. Um `gc.collect()` é executado após o merge para recuperar memória.

---

### ☢️ Radiação N-Terminal (Imputação de NaN)
Preencha valores ausentes (`NaN`) em uma coluna selecionada usando um de três métodos de irradiação:

| Método | Descrição |
|---|---|
| `Mean` | Substitui NaN pela média aritmética da coluna |
| `Median` | Substitui NaN pelo valor mediano da coluna |
| `Fixed Value` | Substitui NaN por uma constante definida pelo usuário |

---

### ⚗️ Transcrição de Coluna (Transmutação de Coluna)
Crie uma **nova coluna derivada** a partir de duas existentes usando operações aritméticas ou de string:

| Feitiço | Descrição |
|---|---|
| `Concat` | Concatenação de string de duas colunas |
| `Add (+)` | Adição numérica |
| `Subtract (-)` | Subtração numérica |
| `Multiply (*)` | Multiplicação numérica |
| `Divide (/)` | Divisão numérica |

A coluna resultante é anexada ao dataset com um nome definido pelo usuário.

---

### 📊 Espectrometria (Estatísticas de Coluna)
Abre um diálogo exibindo estatísticas descritivas (`pandas .describe()`) para a coluna atualmente selecionada no painel de varredura profunda. Inclui contagem, valores únicos, média, desvio padrão, min/max e distribuição por quartis.

---

### 🚨 Iluminar Anomalias (Destaque Condicional de Células)
Defina regras visuais para destacar células na tabela de dados com base em seus valores. Cada regra visa uma coluna (ou `ALL_COLUMNS`), especifica um operador e um valor de limiar, e atribui uma cor de destaque (Vermelho, Verde ou Amarelo). Essas regras são aplicadas ao vivo no `Infinite_Hologram_Model` — o `QAbstractTableModel` customizado que alimenta a view de tabela.

---

### ⏪ Desfazer / ⏩ Refazer (Dilatação Temporal)
O `Sterilized_Data_Engine` mantém uma pilha de histórico de estados. Cada operação mutante (detonação de filtro, edição CRISPR de célula, purga de clones, splice alienígena, radiação, transmutação) salva um snapshot do dataset filtrado atual. As operações de desfazer (`rewind_time_dilation`) e refazer (`temporal_redo`) percorrem esse histórico, permitindo que operadores revertam ou reapliquem transformações sem recarregar do disco.

---

### ✏️ Edição CRISPR de Células
Células individuais na tabela de dados são diretamente editáveis. Um duplo clique em uma célula permite modificação in-place. As alterações são roteadas através de `signal_crispr_edit` e confirmadas no `Sterilized_Data_Engine`, que atualiza o dataset filtrado e empurra um novo estado de desfazer.

---

### 🔬 Painel Raio-X Biológico
Um dashboard de estatísticas ao vivo, atualizado após cada colisão de filtro e mutação de dados:

- **Leituras Brutas** — Número total de linhas no arquivo injetado originalmente.
- **Sobreviventes** — Número de linhas restantes após a configuração de filtro atual.
- **Taxa** — Taxa de sobrevivência como porcentagem, calculada como `(sobreviventes / leituras_brutas) × 100`.

---

## 📤 Protocolos de Síntese para Exportação

Quando pronto para sintetizar a saída, quatro formatos estão disponíveis. Todas as exportações rodam em uma `QThread` dedicada (`Chunk_Devourer_Exporter` ou `Fasta_Devourer_Exporter`) com relatório de progresso ao vivo e um **Aviso de Paradoxo Temporal** caso o arquivo de destino já exista.

| Formato | Thread | Observações |
|---|---|---|
| **CSV** | `Chunk_Devourer_Exporter` | Escrito em blocos de 50.000 linhas para eficiência de memória |
| **TSV** | `Chunk_Devourer_Exporter` | Mesmo pipeline em blocos, separado por tabulação |
| **Excel** | `Chunk_Devourer_Exporter` | Passagem única via `pandas.to_excel` (`openpyxl`) |
| **FASTA** | `Fasta_Devourer_Exporter` | Requer mapeamento de uma coluna de cabeçalho e uma coluna de sequência. Exporta no formato padrão `>cabeçalho\nsequência`, reportando progresso a cada 1.000 linhas. |

Todos os arquivos exportados têm como nome padrão `{nome_original}_mutated.{ext}`.

---

## 🌐 Internacionalização — O Ribossomo Pedra de Rosetta

Todas as strings da UI são gerenciadas pelo `Rosetta_Stone_Ribosome` em `xenoglossy_codex.py`. O codex mapeia chaves de string para traduções em três dialetos suportados:

| Código | Idioma |
|---|---|
| `EN` | Inglês |
| `PT-BR` | Português Brasileiro |
| `ES` | Espanhol |

**Uso em Runtime:**
```python
# Troca o dialeto ativo
Rosetta_Stone_Ribosome.mutate_dialect("PT-BR")

# Recupera uma string localizada
label = Rosetta_Stone_Ribosome.extract_peptide("btn_inject")
# → "🧬 Injetar"
```

A troca de dialeto está exposta na UI. Ao selecionar um idioma no combo box, o sinal `signal_dialect_mutation` é disparado, que chama `mutate_dialect()` e depois aciona `retranscribe_ui()` no Holodeck — fazendo com que todos os rótulos, botões e caixas de grupo visíveis sejam re-renderizados com o novo dialeto instantaneamente, sem reiniciar a aplicação.

Para **adicionar um novo idioma**, estenda cada entrada em `_genetic_codex` com o novo código de idioma e sua string de tradução correspondente.

---

## 🎨 Temas — Mudança do Espectro Quântico

O Delimiter suporta dois modos visuais, alternados em runtime pelo botão de alternância **☀️ / 🌙**:

| Tema | Arquivo | Descrição |
|---|---|---|
| **Matéria Escura** | `themes/dark_theme.qss` | Padrão. Fundo escuro, cores de destaque de alto contraste. |
| **Fóton Solar** | `themes/light_theme.qss` | Modo claro para operadores que temem o vazio. |

Os temas são carregados por `Chromatic_Mutator.force_mutation()`, que lê o arquivo `.qss` e o aplica à instância `QApplication`. Se um arquivo de tema não for encontrado em `themes/`, o sistema retorna à busca no diretório de trabalho e registra um aviso de **Anomalia Estética** no log.

Para criar um tema personalizado, crie um novo arquivo `.qss` seguindo a sintaxe de stylesheet Qt e carregue-o modificando a chamada `Chromatic_Mutator.force_mutation()`.

---

## 📋 Logs — Runas de Biocontenção

Todas as operações significativas são gravadas em `panvita_quarantine.log` via `Runes_Of_Biocontainment.engrave_biocontainment_runes()`. Isso inclui:

- Início de sessão
- Colisões de filtro (quantos pods Matryoshka foram aplicados)
- Conclusões de síntese de FASTA / CSV / TSV / Excel
- Eventos de splice alienígena (quais colunas-chave foram usadas)
- Eventos de transmutação de coluna (qual nova fita foi criada)
- Rastreamentos de erro em brecha de contenção

Formato do log:
```
AAAA-MM-DD HH:MM:SS,ms - BIOLOGICAL ALERT - [mensagem]
```

Em uma **falha catastrófica de inicialização**, um `CRASH_REPORT_BIOLAB.txt` separado é gravado no diretório de trabalho com o traceback completo do Python antes de o processo encerrar com código `1`.

---

## 🐛 Anomalias Conhecidas e Solução de Problemas

| Sintoma | Causa Provável | Resolução |
|---|---|---|
| Aplicação não inicia | Dependência ausente | Execute `pip install PyQt5 pandas openpyxl` |
| Tema não aplicado | Diretório `themes/` ausente | Crie a pasta `themes/` e adicione os arquivos `.qss` |
| `Alerta de Sobrecarga de Massa` em arquivos pequenos | `psutil` não instalado | `pip install psutil` — ou prossiga e aceite o risco |
| Filtro retorna 0 linhas inesperadamente | Incompatibilidade de tipo de dado ou sensibilidade a maiúsculas | Verifique o tipo de operador. Use `Numeric Equals` para números, não `Equals`. |
| Exportação FASTA falha | Colunas selecionadas contêm NaN | Aplique a Radiação N-Terminal antes de exportar |
| Relatório de crash gerado | Falha crítica de importação ou exceção não tratada | Leia o `CRASH_REPORT_BIOLAB.txt` — verifique o traceback |
| UI trava durante o carregamento | Problema de threading — não deveria ocorrer | Registre um relatório de brecha de contenção com o log |

---

## ⚠️ Diretrizes de Biossegurança

> *O Brasil não possui norma ABNT para Cabines de Segurança Biológica. Na ausência de regulamentação nacional, a norma aplicável é a **NSF/ANSI 49**.*

Esta aplicação processa dados reais com operações reais. Os seguintes protocolos são obrigatórios:

- **Sempre mantenha uma cópia de backup** dos arquivos de origem antes da injeção e do processamento.
- **Verifique as taxas de sobrevivência** no painel Raio-X após cada colisão de filtro antes de exportar.
- **Não force o encerramento** durante uma thread de síntese ativa — arquivos de saída parciais serão produzidos.
- **Audite o `panvita_quarantine.log`** após processar grandes datasets para verificar que todas as operações foram concluídas conforme esperado.
- Quando um **Aviso de Paradoxo Temporal** for levantado para sobrescritas de arquivo, confirme deliberadamente — esta ação é irreversível dentro da aplicação.

---

## 👤 Autor e Atribuições

| Campo | Detalhe |
|---|---|
| **Autor** | Victor S. Caricatte De Araújo |
| **E-mail** | victorsc@ufmg.br |
| **Instituição** | Universidade Federal de Minas Gerais (UFMG) |
| **Versão** | 0.3.1 |
| **Codinome** | Acelerador de Partículas |
| **Repositório** | [github.com/VictorCaricatte/BasicBioinfo](https://github.com/VictorCaricatte/BasicBioinfo) |

---

<p align="center">
  <b>🧬 DELIMITER — Filtragem concluída. Sobreviventes reportados. Reator esfriando. 🧬</b><br>
  <i>Construído na UFMG. Alimentado por Python, Pandas e a vontade de detonar dados.</i>
</p>
