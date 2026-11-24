"""
login_node_preprocess.py  —  Run on the LOGIN NODE (internet access, no GPU)
=============================================================================
Steps performed here:
  1. Download the chosen biomedical LLM from HuggingFace to a local directory.
  2. Pre-fetch PubMed articles for all gene/disease combinations.

Everything saved to disk is then used by the compute node without any internet.

Usage
-----
  python login_node_preprocess.py                               # uses secrets.env
  python login_node_preprocess.py --list-models                 # print model registry and exit
  python login_node_preprocess.py --model biom_mistral_7b       # override model key
  python login_node_preprocess.py --pubmed-only                 # skip model download
  python login_node_preprocess.py --model-only                  # skip PubMed prefetch
"""

import os
import sys
import json
import argparse
from dotenv import load_dotenv

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="CoDHy login-node preprocessing")
parser.add_argument("--list-models",  action="store_true",
                    help="Print the model registry and exit.")
parser.add_argument("--model",        type=str, default=None,
                    help="Override LLM_MODEL_KEY from secrets.env.")
parser.add_argument("--pubmed-only",  action="store_true",
                    help="Skip model download, only prefetch PubMed.")
parser.add_argument("--model-only",   action="store_true",
                    help="Skip PubMed prefetch, only download model.")
parser.add_argument("--env",          type=str, default="secrets.env",
                    help="Path to secrets.env (default: secrets.env).")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Import registry (must come after argparse so --list-models works without deps)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from llm_backend import MODEL_REGISTRY, list_models

if args.list_models:
    list_models()
    sys.exit(0)

# ---------------------------------------------------------------------------
# Load config
# ---------------------------------------------------------------------------
load_dotenv(args.env)

MODEL_KEY  = args.model or os.getenv("LLM_MODEL_KEY", "jsl_medllama_3_8b")
HF_TOKEN   = os.getenv("HF_TOKEN", "")
LOCAL_DIR  = os.getenv("LLM_LOCAL_DIR", "/user/suvinava.basak/u19797/CoDHy_v2/models")
SAVE_PATH  = os.getenv("PUBMED_CACHE_DIR",
                        "/user/suvinava.basak/u19797/CoDHy_v2/data_cache/pubmed/")
EMAIL      = os.getenv("EMAIL", "")

os.makedirs(LOCAL_DIR, exist_ok=True)
os.makedirs(SAVE_PATH, exist_ok=True)

# ===========================================================================
# STEP 1 — Download LLM
# ===========================================================================
if not args.pubmed_only:
    if MODEL_KEY not in MODEL_REGISTRY:
        print(f"ERROR: Unknown model key '{MODEL_KEY}'.")
        list_models()
        sys.exit(1)

    repo_id, gated, vram, notes = MODEL_REGISTRY[MODEL_KEY]
    model_dir = os.path.join(LOCAL_DIR, MODEL_KEY)

    print("=" * 70)
    print(f"  STEP 1: Downloading model  ->  {MODEL_KEY}")
    print(f"  Repo : {repo_id}")
    print(f"  Dest : {model_dir}")
    print(f"  Notes: {notes}")
    if gated:
        print(f"  WARNING: This model is gated. Make sure HF_TOKEN has access.")
    print("=" * 70)

    if os.path.isdir(model_dir) and len(os.listdir(model_dir)) > 3:
        print(f"  Model directory already populated. Skipping download.")
        print(f"  Delete {model_dir} manually to force re-download.")
    else:
        from huggingface_hub import snapshot_download
        try:
            snapshot_download(
                repo_id=repo_id,
                local_dir=model_dir,
                token=HF_TOKEN if HF_TOKEN else None,
                # Skip framework weights we won't use (saves disk space)
                ignore_patterns=["*.msgpack", "flax_model*", "tf_model*", "rust_model*"],
            )
            print(f"  Model downloaded to {model_dir}")
        except Exception as e:
            print(f"  Download failed: {e}")
            sys.exit(1)

# ===========================================================================
# STEP 2 — Pre-fetch PubMed articles
# ===========================================================================
if not args.model_only:
    print()
    print("=" * 70)
    print("  STEP 2: Pre-fetching PubMed abstracts")
    print("=" * 70)

    # These must match main.py exactly
    GENES_OF_INTEREST = [
        "MYC", "EGFR", "CCNE1", "MDM2", "PIK3CA", "PTEN", "TP53", "KRAS",
        "CHK1", "PARP1", "WEE1", "ATR", "AURKA", "CDK4", "CDK6", "ATM"
    ]
    DISEASES = [
        "Glioblastoma Multiforme", "Brain Lower Grade Glioma",
        "Breast Invasive Carcinoma", "Lung Adenocarcinoma",
        "Lung Squamous Cell Carcinoma", "Colon Adenocarcinoma",
        "Rectum Adenocarcinoma", "Skin Cutaneous Melanoma",
        "Head and Neck Squamous Cell Carcinoma",
        "Kidney Renal Clear Cell Carcinoma",
        "Liver Hepatocellular Carcinoma",
        "Ovarian Serous Cystadenocarcinoma",
        "Pancreatic Adenocarcinoma", "Prostate Adenocarcinoma",
        "Bladder Urothelial Carcinoma", "Stomach Adenocarcinoma",
    ]

    try:
        from Bio import Entrez
        from nlp_processor import LiteratureProcessor
    except ImportError as e:
        print(f"  Could not import LiteratureProcessor: {e}")
        print("  Make sure nlp_processor.py is in the same directory.")
        sys.exit(1)

    Entrez.email = EMAIL
    nlp = LiteratureProcessor(existing_nodes=[], existing_relations=[], email=EMAIL)

    total   = len(GENES_OF_INTEREST) * len(DISEASES)
    cached  = 0
    fetched = 0
    failed  = 0

    for gene in GENES_OF_INTEREST:
        for disease in DISEASES:
            cache_key  = f"pubmed_{gene}_{disease.replace(' ', '')}.json"
            cache_file = os.path.join(SAVE_PATH, cache_key)

            if os.path.exists(cache_file):
                cached += 1
                continue

            try:
                papers = nlp.fetch_papers(gene, disease, max_results=50)
                with open(cache_file, 'w') as f:
                    json.dump(papers, f)
                fetched += 1
                print(f"  [{fetched + cached}/{total}] Cached {len(papers):>3} papers "
                      f"for  {gene:<8} + {disease}")
            except Exception as exc:
                failed += 1
                print(f"  FAILED: {gene} + {disease}  ->  {exc}")

    print()
    print(f"  PubMed prefetch complete.")
    print(f"   Already cached : {cached}")
    print(f"   Newly fetched  : {fetched}")
    print(f"   Failed         : {failed}")
    print(f"   Total combos   : {total}")

print()
print("Login-node preprocessing finished.")
print("You can now submit the SLURM job on the compute node.")