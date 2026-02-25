# CoDHy — Computational Drug Hypothesis Generation System

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/LICENSE-CC%20BY--NC--SA%204.0-blue.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![NetworkX 3.4.2](https://img.shields.io/badge/networkx-3.4.2-0086B3.svg)](https://networkx.org/)
[![Node2Vec 0.4.6](https://img.shields.io/badge/node2vec-0.4.6-8E44AD.svg)](https://github.com/eliorc/node2vec)
[![Neo4j 5.26.0](https://img.shields.io/badge/neo4j-5.26.0-008CC1.svg)](https://neo4j.com/)
[![Reactome2Py 3.0.1](https://img.shields.io/badge/reactome2py-3.0.1-E35F21.svg)](https://github.com/reactome/reactome2py)
[![CIViCPy 1.3.1](https://img.shields.io/badge/civicpy-1.3.1-224488.svg)](https://civicpy.readthedocs.io/)
[![ChEMBL Client 0.10.9](https://img.shields.io/badge/chembl_client-0.10.9-2E8B57.svg)](https://github.com/chembl/chembl_webresource_client)
[![RDKit 2025.09.1](https://img.shields.io/badge/rdkit-2025.09.1-713B2A.svg)](https://www.rdkit.org/)
[![Biopython 1.84](https://img.shields.io/badge/biopython-1.84-4B8BBE.svg)](https://biopython.org/)
[![SpaCy 3.8.0](https://img.shields.io/badge/spacy-3.8.0-09A3D5.svg)](https://spacy.io/)
[![SciSpaCy 0.5.4](https://img.shields.io/badge/scispacy-0.5.4-4A90E2.svg)](https://allenai.github.io/scispacy/)
[![TheFuzz 0.22.1](https://img.shields.io/badge/thefuzz-0.22.1-A9A9A9.svg)](https://github.com/seatgeek/thefuzz)
[![Sentence Transformers 3.4.1](https://img.shields.io/badge/sentence_transformers-3.4.1-FFD21E.svg)](https://www.sbert.net/)

**CoDHy** is an interactive, human-in-the-loop decision support system for **biomarker-guided drug discovery and hypothesis generation** in translational oncology. It integrates multi-source biomedical knowledge graphs, literature mining, and large language models (LLMs) to generate, validate, and rank combination therapy hypotheses for cancer treatment.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [System Architecture](#system-architecture)
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
- [Data Sources](#data-sources)
- [Project Structure](#project-structure)
- [License](#license)

---

## Overview

CoDHy automates the early-stage drug hypothesis generation pipeline by:

1. Fetching and integrating genomic, clinical, and pharmacological data from multiple public databases.
2. Building a **Neo4j-backed biomedical knowledge graph (KG)** connecting genes, drugs, diseases, pathways, and clinical trials.
3. Mining the scientific literature (PubMed) to extract novel relationship triples and enrich the KG.
4. Using **graph embeddings** (Node2Vec) to learn latent representations of biomedical entities.
5. Generating combination therapy hypotheses via a hybrid LLM + KG approach.
6. Validating, ranking, and reporting hypotheses with supporting evidence and synergy data.

---

## Key Features

- **Multi-source Knowledge Graph** — integrates CIViC, TCGA/GDC, DepMap, Reactome, ChEMBL, STRING, SynLethDB, SIDER, and clinical trials data.
- **Literature Mining** — automated PubMed abstract retrieval and biomedical NLP triple extraction using SciSpaCy and sentence transformers.
- **Graph Embedding** — Node2Vec-based latent representation learning for semantic similarity scoring.
- **LLM-driven Hypothesis Generation** — leverages locally deployed Llama 3.1 (via Ollama) for context-aware drug combination proposal.
- **Automated Validation** — PubMed evidence retrieval, novelty checks, plausibility assessment, and toxicity risk profiling.
- **Evidence-based Ranking** — composite scoring incorporating graph evidence, AI safety scores, and synergy database lookups.
- **Human-in-the-Loop Design** — transparent reasoning with full provenance (database URLs, PubMed links, graph paths).

---

## System Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    Data Ingestion Layer                 │
│  CIViC · GDC/TCGA · DepMap · Reactome · ChEMBL          │
│  STRING · SynLethDB · SIDER · ClinicalTrials.gov        │
└────────────────────────┬────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────┐
│              Knowledge Graph (Neo4j)                    │
│  Nodes: Gene · Drug · Disease · Pathway · Trial         │
│  Edges: TARGETS · DRIVES · ASSOCIATED_WITH · ...        │
└──────┬─────────────────┬───────────────────────┬────────┘
       │                 │                       │
┌──────▼──────┐  ┌───────▼──────┐   ┌────────────▼──────┐
│ Literature  │  │    Graph     │   │   LLM Agent       │
│  Mining     │  │  Embedding   │   │ (Llama 3.1/Ollama)│
│ (SciSpaCy)  │  │  (Node2Vec)  │   │                   │
└──────┬──────┘  └───────┬──────┘   └───────────┬───────┘
       │                 │                       │
┌──────▼─────────────────▼───────────────────────▼───────┐
│           Hypothesis Generator (DrugScreeningAgent)    │
│         Hybrid: KG-Only · LLM-Only · Full              │
└────────────────────────┬───────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────┐
│         Validation & Ranking Pipeline                   │
│  ValidationAgent · RankingAgent · SynergyChecker        │
└────────────────────────┬────────────────────────────────┘
                         │
┌────────────────────────▼────────────────────────────────┐
│                  Final Report (JSON)                    │
│  Ranked drug combinations with scores & evidence links  │
└─────────────────────────────────────────────────────────┘
```

---

## Installation

Find below the instructions to run CoDHy locally. Follow the steps below to set up the environment:

### 1. Clone the repository or copy the project files

Make sure to keep all `.py` module files in the same file as the `main.py` file. Be sure to change the path for `folder_path` (main path where all `.py` files are saved), `SAVE_PATH` (where data caches from biological databases and embeddings are saved) and `REPORT_DIR` (where final generated reports are saved):
```
folder_path = '/content/drive/MyDrive/Colab Notebooks/AIcoS/FINAL'
SAVE_PATH = "/content/drive/MyDrive/Colab Notebooks/AIcoS/FINAL/data_cache/"
REPORT_DIR = "/content/drive/MyDrive/Colab Notebooks/AIcoS/FINAL/reports/"
```

### 2. Install dependencies

Make sure to install all the library dependencies. Run:

```python
!pip install reactome2py civicpy chembl-webresource-client rdkit neo4j biopython spacy scispacy -q
!pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz -q
!pip install thefuzz sentence-transformers node2vec gensim networkx pandas requests -q
!pip install ollama -q
```

### 3. Set up Ollama and the LLM

```bash
!apt-get install -y zstd
!curl -fsSL https://ollama.com/install.sh | sh
!ollama serve &
!ollama pull llama3.1
```

### 4. Configure Neo4j credentials

Store the following secrets in your environment (via `userdata`):
- `NEO4J_URI`
- `NEO4J_USER`
- `NEO4J_PASSWORD`

---

## Requirements

See [`requirements.txt`](requirements.txt) for the full list of Python dependencies. Key packages include:

| Package | Purpose |
|---|---|
| `neo4j` | Knowledge graph database driver |
| `biopython` | PubMed/Entrez API access |
| `civicpy` | CIViC clinical evidence database |
| `chembl-webresource-client` | ChEMBL drug data API |
| `reactome2py` | Reactome pathway analysis |
| `rdkit` | Cheminformatics |
| `spacy` + `scispacy` | Biomedical NLP pipeline |
| `sentence-transformers` | Semantic similarity for entity linking |
| `node2vec` + `gensim` | Graph embedding |
| `networkx` | Graph construction and analysis |
| `pandas` | Data manipulation |
| `thefuzz` | Fuzzy string matching |
| `requests` | HTTP API calls |
| `ollama` | Local LLM interface (Llama 3.1) |

---

## Usage

1. Go to your project folder and add it to the Python path.
2. Set your gene(s) of interest and cancer type(s) in `main.py`.
3. Configure your Neo4j connection credentials.
4. Run `main.py`:

```python
INPUT_GENE = "EGFR"
INPUT_DISEASE = "Breast Invasive Carcinoma"
CHOSEN_METHOD = "Full"   # Options: "Full", "KG-Only", "LLM-only"
LIT_TO_FETCH = 50

results = run_full_experiment(kg, embedder, INPUT_GENE, INPUT_DISEASE,
                              method=CHOSEN_METHOD, num_lit=LIT_TO_FETCH)
```

Results are printed to the console and saved as a timestamped JSON report in the configured `REPORT_DIR`.

---

## Data Sources

| Database | Description |
|---|---|
| [CIViC](https://civicdb.org/) | Clinical interpretation of variants in cancer |
| [TCGA/GDC](https://portal.gdc.cancer.gov/) | Genomic data commons, cancer genomics |
| [DepMap](https://depmap.org/) | Cancer dependency map, gene essentiality |
| [Reactome](https://reactome.org/) | Pathway analysis |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | Drug and bioactivity data |
| [STRING](https://string-db.org/) | Protein-protein interaction network |
| [SynLethDB](http://synlethdb.sist.shanghaitech.edu.cn/) | Synthetic lethality database |
| [SIDER](http://sideeffects.embl.de/) | Drug side effects |
| [ClinicalTrials.gov](https://clinicaltrials.gov/) | Clinical trial data |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/) | Biomedical literature |

---

## Project Structure

```
CoDHy/
├── main.py               # Main orchestration script (entry point)
├── build_kg.py           # Knowledge graph construction and Neo4j integration
├── fetch_data_util.py    # Data fetching utilities for all external databases
├── nlp_processor.py      # Literature mining and NLP triple extraction
├── graph_embedding.py    # Node2Vec graph embedding agent
├── hgenerator.py         # Drug combination hypothesis generator (LLM agent)
├── hvalidator.py         # Hypothesis validation agent
├── hranker.py            # Hypothesis ranking and scoring
├── synergycheck.py       # Drug synergy verification (DrugComb database)
├── requirements.txt      # Python package dependencies
├── data_cache/           # Cached CSV data from external databases
├── reports/              # Generated JSON reports
└── README.md             # This file
```

---

## License

This project is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License**.

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/LICENSE-CC%20BY--NC--SA%204.0-blue.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

You are free to share and adapt this work for **non-commercial purposes**, provided you give appropriate credit and distribute your contributions under the same license.
