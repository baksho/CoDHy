import io
import os
import json
import requests
import zipfile
import urllib.parse
import pandas as pd
from rdkit import Chem
from civicpy import civic
from reactome2py import analysis
from chembl_webresource_client.new_client import new_client

"""### **Reactome Pathways Database**"""

def fetch_reactome_pathways(genes):
    """
    Fetches Reactome pathways using the Analysis Service.
    FIXED: Decodes the token and uses the correct CSV download endpoint.
    """
    print(f"--- Reactome: Analysis Service for {len(genes)} genes ---")

    # 1. Prepare Payload
    gene_str = ",".join(genes)

    try:
        # 2. Submit Analysis
        result = analysis.identifiers(ids=gene_str)

        # 3. Get and DECODE the Token
        # The token often comes as '...%3D'. We need '...='
        raw_token = result.get("summary", {}).get("token")
        if not raw_token:
            print("⚠️ Failed to get Analysis Token.")
            return pd.DataFrame()

        token = urllib.parse.unquote(raw_token)
        print(f"   Analysis Token: {token} (Decoded)")

        # 4. Download the Mapping CSV
        # Endpoint: /download/{token}/entities/found/{resource}/{filename}.csv
        # Resource = TOTAL (all resources), Filename = result.csv (arbitrary)
        url = f"https://reactome.org/AnalysisService/download/{token}/pathways/TOTAL/result.csv"

        r = requests.get(url)
        if r.status_code != 200:
            print(f"⚠️ Download failed (Status {r.status_code}). URL: {url}")
            return pd.DataFrame()

        # 5. Parse CSV
        # The Reactome CSV format usually has columns:
        # [Pathway Identifier, Pathway Name, ..., Found Identifiers, ...]
        df_raw = pd.read_csv(io.StringIO(r.text))
        # print(df_raw.head())

        data = []
        pmid_cache = {}
        # 5. Iterate and Parse
        for _, row in df_raw.iterrows():
            p_id = row.get("Pathway identifier")
            p_name = row.get("Pathway name")

            # Extract the Gene String (e.g., "CHK1;ATM;TP53")
            found_ids = row.get("Submitted entities found")

            if pd.isna(found_ids):
                continue

            # --- Fetch PubMed IDs from Content Service ---
            if p_id not in pmid_cache:
                pmids = []
                try:
                    # Query the Content Service for pathway metadata
                    content_url = f"https://reactome.org/ContentService/data/query/{p_id}"
                    content_resp = requests.get(content_url, timeout=10)
                    if content_resp.status_code == 200:
                        content_data = content_resp.json()
                        # Extract IDs from the literatureReference list
                        refs = content_data.get("literatureReference", [])
                        pmids = [str(ref.get("pubMedIdentifier")) for ref in refs if ref.get("pubMedIdentifier")]
                except:
                    pass # Fallback to empty list if API fails
                pmid_cache[p_id] = "; ".join(pmids)


            # Split found genes and build rows
            found_list = str(found_ids).split(";")

            for found_gene in found_list:
                found_gene = found_gene.strip() # Remove any extra whitespace

                # Double-check it's in our target list (optional but safe)
                if found_gene in genes:
                    data.append({
                        "gene": found_gene,
                        "pathway_id": p_id,
                        "pathway_name": p_name,
                        "pmids": pmid_cache[p_id],
                        "source_url": f"https://reactome.org/content/detail/{p_id}"
                    })

        if not data:
            print("   No pathways found for these genes.")
            return pd.DataFrame()

        df = pd.DataFrame(data).drop_duplicates()
        print(f"✅ Found {len(df)} gene-pathway associations from reactome pathways database.")
        return df

    except Exception as e:
        print(f"❌ Reactome Error: {e}")
        return pd.DataFrame()

"""### **CIViC - Clinical Interpretation of Variants in Cancer**"""

def fetch_civic_database(genes):
    """
    Fetches rich clinical data from CIViC to populate the refined schema:
    Gene -> Variant -> Disease -> Drug -> Evidence (PMID)
    """
    print(f"--- CIViC (Rich Data): Fetching for {len(genes)} genes ---")
    data = []

    # CIViCpy caches data, so the first run might take a few seconds to initialize
    for gene_symbol in genes:
        try:
            # 1. Get Gene Object
            gene = civic.get_gene_by_name(gene_symbol)
            if not gene:
                print(f"❌ Gene '{gene_symbol}' not found in CIViC.")
                continue

            # 2. Iterate through Variants (e.g., V600E)
            for variant in gene.variants:
                variant_name = variant.name

                # 3. Iterate through Evidence Items (Papers) via Molecular Profiles
                # CIViC V2 structure: Variant -> Molecular Profile -> Evidence
                for mp in variant.molecular_profiles:

                    # --- A. PROCESS EVIDENCE ITEMS (Research Papers) ---
                    for evidence in mp.evidence_items:

                        # Filter for Predictive Evidence (Drug Response)
                        if evidence.evidence_type == 'PREDICTIVE' and evidence.status.upper() != 'REJECTED' and evidence.therapies:

                            # A. Extract Drugs (Handle Combinations)
                            # We keep them as a list to UNWIND them later in Neo4j
                            drug_list = [t.name for t in evidence.therapies]

                            # B. Extract Disease
                            disease_name = evidence.disease.name if evidence.disease else "Cancer"

                            # C. Extract PubMed ID & Source URL
                            pmid = evidence.source.citation_id if evidence.source else "N/A"
                            source_url = evidence.site_link # The specific CIViC URL

                            # D. Significance (Sensitivity vs Resistance)
                            significance = evidence.significance.lower() if evidence.significance else "unknown"

                            # E. Evidence Level (A, B, C, D, E)
                            level = evidence.evidence_level if evidence.evidence_level else "N/A"

                            data.append({
                                "gene": gene_symbol,
                                "variant": variant_name,
                                "drugs": drug_list,        # List: ['Cisplatin', 'Olaparib']
                                "source_type": "Evidence",
                                "disease": disease_name,
                                "pmid": str(pmid),
                                "significance": significance,
                                "level": level,
                                "confidence": 0.8,         # Base confidence for evidence items
                                "source_url": source_url   # Stored on Edge
                            })

                    # --- B. PROCESS ASSERTIONS (Clinical Guidelines / FDA) ---
                    # These are critical for "Ground Truth" (Confidence = 1.0)
                    for assertion in mp.assertions:
                        if assertion.status.upper() == 'ACCEPTED' and assertion.therapies:

                            # Extract Fields
                            drug_list = [t.name for t in assertion.therapies]
                            disease_name = assertion.disease.name if assertion.disease else "Cancer"

                            # Assertions might summarize multiple papers, so PMID is often N/A or a list.
                            # We use the Assertion ID as the primary reference if no single PMID exists.
                            # However, sometimes they link to an NCCN guideline ID.
                            pmid = "N/A" # Default for assertions unless specific
                            source_url = assertion.site_link

                            significance = assertion.significance.lower() if assertion.significance else "unknown"

                            # Map AMP Level (Tier I/II) to our Schema
                            amp_level = assertion.amp_level if hasattr(assertion, 'amp_level') else "N/A"

                            # Set High Confidence for FDA/Tier I
                            conf = 1.0 if "TIER_I" in str(amp_level) else 0.9

                            data.append({
                                "gene": gene_symbol,
                                "variant": variant_name,
                                "drugs": drug_list,         # List of strings
                                "source_type": "Assertion",
                                "disease": disease_name,
                                "pmid": str(pmid),          # Often N/A for guidelines
                                "significance": significance,
                                "level": amp_level,         # e.g., "TIER_I_LEVEL_A"
                                "confidence": conf,         # HIGHER CONFIDENCE
                                "source_url": source_url
                            })

        except Exception as e:
            # print(f"Error for {gene_symbol}: {e}")
            continue

    if not data: return pd.DataFrame()
    df = pd.DataFrame(data)
    print(f"✅ Found {len(df)} data from CIViC database.")
    return df

"""### **The Cancer Genomic Atlas (TCGA) data from GDC portal**"""

def fetch_gdc_data(genes, projects_meta, size=1000):
    """Fetches patient mutation data from NCI GDC."""
    print(f"--- GDC: Fetching mutations for {len(genes)} genes in {len(list(projects_meta.keys()))} projects ---")
    url = "https://api.gdc.cancer.gov/ssm_occurrences"
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "case.project.project_id", "value": list(projects_meta.keys())}},
            {"op": "in", "content": {"field": "ssm.consequence.transcript.gene.symbol", "value": genes}}
        ]
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "ssm.consequence.transcript.gene.symbol,ssm.consequence.transcript.consequence_type,case.project.project_id,case.submitter_id",
        "size": size,
        "format": "JSON"
    }

    try:
        r = requests.get(url, params=params, timeout=60)
        data = []
        if r.status_code == 200:
            hits = r.json().get("data", {}).get("hits", [])
            for hit in hits:
                # Extract nested fields safely
                case = hit.get("case", {})
                project_id = case.get("project", {}).get("project_id")
                ssm = hit.get("ssm", {}).get("consequence", [{}])[0] # Take first transcript
                gene = ssm.get("transcript", {}).get("gene", {}).get("symbol")

                if gene and case.get("submitter_id"):
                    # Get metadata from our input dictionary
                    meta = projects_meta.get(project_id, {})

                    data.append({
                        "sample_id": case.get("submitter_id"),
                        "project_id": project_id,
                        "disease_name": meta.get("full_name", "Unknown"),
                        "broad_disease": meta.get("broad_disease", "Unknown"),
                        "gene_symbol": gene,
                        "mutation_type": ssm.get("transcript", {}).get("consequence_type", "Unknown"),
                        "source_url": f"https://portal.gdc.cancer.gov/cases/{case.get('submitter_id')}"
                    })
        df = pd.DataFrame(data)
        print(f"✅ Successfully enriched {len(df)} mutations from TCGA-GDC data portal for {len(genes)} genes in {len(list(projects_meta.keys()))} projects.")
        return df
    except Exception as e:
        print(f"GDC Error: {e}")
        return pd.DataFrame()

"""### **ClinicalTrials.gov**"""

def fetch_clinical_trials(genes, disease_keyword="Cancer"):
    """
    Fetches active clinical trials via ClinicalTrials.gov API v2.
    """
    print(f"--- ClinicalTrials.gov: Searching for {len(genes)} genes ---")
    base_url = "https://clinicaltrials.gov/api/v2/studies"
    data = []

    # We batch requests or do a broad search to save time, here we loop for precision
    for gene in genes: # Limit to first 5 for speed in demo, remove slice for full run
        params = {
            "query.term": f"{gene} AND {disease_keyword}",
            "filter.overallStatus": "RECRUITING",
            "pageSize": 5,
            "format": "json"
        }

        try:
            r = requests.get(base_url, params=params, timeout=60)
            if r.status_code == 200:
                studies = r.json().get("studies", [])
                for study in studies:
                    protocol = study.get("protocolSection", {})
                    id_module = protocol.get("identificationModule", {})
                    status_module = protocol.get("statusModule", {})
                    design_module = protocol.get("designModule", {})
                    conditions = protocol.get("conditionsModule", {}).get("conditions", [])

                    # Status: e.g., RECRUITING, COMPLETED, TERMINATED
                    overall_status = status_module.get("overallStatus", "UNKNOWN")
                    # Phase: e.g., ["PHASE1", "PHASE2"]
                    phases_list = design_module.get("enrollmentInfo", {}).get("phase", []) # Fallback location
                    # Usually found here in v2:
                    phases_list = design_module.get("phases", [])
                    phase_str = ", ".join(phases_list) if phases_list else "N/A"

                    # Extract Interventions (Drugs)
                    arms = protocol.get("armsInterventionsModule", {}).get("interventions", [])
                    drugs = [i["name"] for i in arms if i.get("type") == "DRUG"]

                    if drugs:
                        nct_id = id_module.get("nctId")
                        data.append({
                            "nct_id": nct_id,
                            "title": id_module.get("briefTitle"),
                            "gene": gene,
                            "drugs": ", ".join(drugs),
                            "condition": conditions[0] if conditions else disease_keyword,
                            "status": overall_status,
                            "phase": phase_str,
                            "source_url": f"https://clinicaltrials.gov/study/{nct_id}"
                        })
        except Exception:
            continue

    if not data:
        print("⚠️ No clinical trials found for these genes.")
        return pd.DataFrame()
    else:
        df = pd.DataFrame(data)
        print(f"✅ Found {len(df)} clinical trials for {len(genes)} genes from ClinicalTrials.gov.")
        return df

"""### **Gene Metadata: mygene.info**"""

def fetch_gene_metadata(genes):
    """
    Fetches Gene Descriptions and Metadata from MyGene.info.
    This replaces the broken DepMap/CCLE code.
    """
    print(f"--- MyGene.info: Fetching metadata for {len(genes)} genes ---")

    # MyGene.info allows batch queries, which is much faster than one-by-one
    url = "https://mygene.info/v3/query"
    params = {
        "q": ",".join(genes),
        "scopes": "symbol",
        "fields": "name,summary,entrezgene,ensembl.gene",
        "species": "human"
    }

    try:
        r = requests.post(url, data=params, timeout=15)
        if r.status_code != 200:
            print(f"❌ API Error: {r.status_code}")
            return pd.DataFrame()

        results = r.json()
        data = []

        for item in results:
            gene_symbol = item.get("query")
            # If the gene wasn't found, skip it
            if item.get("notfound"):
                continue

            data.append({
                "gene": gene_symbol,
                "common_name": item.get("name"),
                "description": item.get("summary", "No summary available."),
                "entrez_id": item.get("entrezgene"),
                "ensembl_id": item.get("ensembl", {}).get("gene") if isinstance(item.get("ensembl"), dict) else None,
                "source_url": f"https://www.ncbi.nlm.nih.gov/gene/{item.get('entrezgene')}"
            })

        df = pd.DataFrame(data)
        print(f"✅ Successfully enriched {len(df)} genes from mygene.info database.")
        return df

    except Exception as e:
        print(f"⚠️ Error: {e}")
        return pd.DataFrame()

"""### **ChEMBL**"""

def fetch_chembl_data(genes):
    """
    Fetches drug candidates and properties from ChEMBL.
    Includes: Mechanism of Action, Toxicity flags (Max Phase).
    """
    print(f"--- ChEMBL: Fetching compounds and Trial IDs for {len(genes)} targets ---")
    data = []

    # 1. Get Target ChEMBL IDs for our Genes
    Target = new_client.target
    Mechanism = new_client.mechanism
    Molecule = new_client.molecule
    Indication = new_client.drug_indication

    valid_targets = {}

    for gene in genes:
        res = Target.filter(target_synonym__icontains=gene).filter(target_type="SINGLE PROTEIN").only(["target_chembl_id"])
        if res:
            valid_targets[gene] = res[0]["target_chembl_id"]

    for gene, target_id in valid_targets.items():
        mechs = Mechanism.filter(target_chembl_id=target_id).only(["molecule_chembl_id", "action_type"])

        for m in mechs[:20]:
            mol_id = m["molecule_chembl_id"]
            mol = Molecule.get(mol_id)
            synonyms = list(set([s['molecule_synonym'] for s in mol.get('molecule_synonyms', [])]))
            if not mol: continue

            data.append({
                "gene": gene,
                "drug": mol.get("pref_name") or mol_id,
                "synonyms": synonyms,
                "chembl_id": mol_id,
                "action": m.get("action_type"),
                "max_phase": mol.get("max_phase"),
                "is_withdrawn": mol.get("withdrawn_flag"),
                "source_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol_id}/"
            })

    if not data:
        print("⚠️ No ChEMBL data found for these genes.")
    else:
        df = pd.DataFrame(data)
        print(f"✅ Found {len(df)} ChEMBL compounds for {len(genes)} genes from ChEMBL database.")
        return df

"""### **STRING: functional protein association networks**"""

def fetch_string_interactions(genes, species=9606, score_cutoff=700):
    """
    Fetches high-confidence Protein-Protein Interactions from STRING.
    - species: 9606 is the NCBI taxon ID for Humans.
    - score_cutoff: 700 is the threshold for 'High Confidence'.
    """
    print(f"--- STRING: Fetching interactions for {len(genes)} genes ---")

    # 1. Map Gene Symbols to STRING Identifiers
    method = "get_string_ids"
    params = {
        "identifiers": "\r".join(genes), # API expects carriage returns
        "species": species,
        "caller_identity": "my_knowledge_graph_project"
    }

    base_url = "https://string-db.org/api/json/"
    response = requests.post(base_url + method, data=params)

    if response.status_code != 200:
        print("❌ Error mapping IDs to STRING.")
        return pd.DataFrame()

    # Create a mapping of STRING ID -> Original Symbol
    mapping_df = pd.DataFrame(response.json())
    string_ids = mapping_df['stringId'].tolist()

    # 2. Fetch the Network Interactions
    method = "network"
    params = {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "required_score": score_cutoff, # High confidence only
        "caller_identity": "my_knowledge_graph_project"
    }

    response = requests.post(base_url + method, data=params)

    if response.status_code == 200:
        interactions = response.json()
        df = pd.DataFrame(interactions)

        # Clean up column names to match your KG style
        if not df.empty:
            df = df[['preferredName_A', 'preferredName_B', 'score', 'nscore', 'escore', 'dscore', 'tscore']]
            # score: Combined confidence
            # escore: Experimental evidence
            # dscore: Database evidence
            # tscore: Text-mining evidence
            df_filtered = df[df['score'] >= 0.700]

        print(f"✅ Found {len(df_filtered)} interactions for {len(genes)} genes from STRING database.")
        return df_filtered
    else:
        print("⚠️ No interactions found or API error.")
        return pd.DataFrame()

"""### **SynlethDB**"""

def fetch_synlethdb_local(genes_of_interest):
    """
    Processes SynlethDB CSVs via direct URLs.
    Filters for specific genes to keep memory low.
    """
    sl_url = "https://synlethdb.sist.shanghaitech.edu.cn/v2/static/download/SL/Human_SL.csv"
    nonsl_url = "https://synlethdb.sist.shanghaitech.edu.cn/v2/static/download/nonSL/Human_nonSL.csv"

    results = []

    # Process SL (Positive) pairs
    print("--- SynlethDB: Processing Positive SL Pairs ---")
    try:
        # We read the whole thing because it's usually small enough (~10-20MB)
        # Columns: gene_a, gene_b, score, source, pubmed_id
        df_sl = pd.read_csv(sl_url)
        # print(df_sl.head())

        # Filter: Either gene_a OR gene_b must be in our 16 genes
        mask = (df_sl['n1.name'].isin(genes_of_interest)) | (df_sl['n2.name'].isin(genes_of_interest))
        df_sl_filtered = df_sl[mask].copy()
        df_sl_filtered['is_lethal'] = True
        results.append(df_sl_filtered)
        print(f"   Found {len(df_sl_filtered)} lethal interactions.")
    except Exception as e:
        print(f"❌ Error reading SL file: {e}")

    # Process non-SL (Negative) pairs
    print("--- SynlethDB: Processing non-SL Pairs ---")
    try:
        df_nonsl = pd.read_csv(nonsl_url)
        # print(df_nonsl.head())
        mask_non = (df_nonsl['n1.name'].isin(genes_of_interest)) | (df_nonsl['n2.name'].isin(genes_of_interest))
        df_nonsl_filtered = df_nonsl[mask_non].copy()
        df_nonsl_filtered['is_lethal'] = False
        results.append(df_nonsl_filtered)
        print(f"   Found {len(df_nonsl_filtered)} non-lethal interactions.")
    except Exception as e:
        print(f"❌ Error reading non-SL file: {e}")

    if not results:
        return pd.DataFrame()

    df = pd.concat(results, ignore_index=True)
    print(f"✅ Successfully fetched {len(df)} lethalities for {len(genes_of_interest)} genes from SynlethDB.")
    return df

"""### **DepMap: The Cancer Dependency Map Project**"""

def get_depmap_direct_links(release_name="DepMap Public 25Q3"):
    # This API endpoint returns a table of ALL available files and their direct links
    index_url = "https://depmap.org/portal/api/download/files"
    r = requests.get(index_url)

    if r.status_code == 200:
        df_files = pd.read_csv(index_url)
        print(f"Available columns in DepMap files index: {df_files.columns.tolist()}") # Added for debugging
        # Filter for the release and the specific files you want
        my_release = df_files[df_files['release'] == release_name] # Changed to use 'release' column

        if my_release.empty:
            print(f"❌ Could not find release '{release_name}' or appropriate column for release name.")
            return None

        # We search for our three targets
        targets = ["Model.csv", "CRISPRGeneDependency.csv", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"]
        links = {}

        for target in targets:
            match = my_release[my_release['filename'] == target]
            if not match.empty:
                links[target] = match.iloc[0]['url'] # Changed 'downloadUrl' to 'url' based on columns
                print(f"✅ Found direct link for {target}")
            else:
                print(f"❌ Could not find {target} in release {release_name}")
        return links
    else:
        print("Could not access DepMap File Index.")
        return None

def fetch_depmap_essentials(genes_of_interest):
    """
    Processes Model, CRISPR, and Omics data using direct links.
    Filters for specific genes to keep memory usage low.
    """
    direct_links = get_depmap_direct_links()

    results = {}

    # 1. Process Model.csv (Relatively small)
    print("--- Processing Model Metadata ---")

    model_headers = pd.read_csv(direct_links['Model.csv'], nrows=0).columns.tolist()
    # print(f"Model.csv headers: {model_headers}")

    df_models = pd.read_csv(direct_links['Model.csv'], usecols=['ModelID', 'CellLineName', 'OncotreeLineage'])
    results['models'] = df_models
    # print(df_models.head())

    # 2. Process CRISPRGeneDependency.csv (Large)
    print("--- Processing CRISPR Dependency (Filtering Genes) ---")
    crispr_url = direct_links['CRISPRGeneDependency.csv']

    # Temporarily read headers to inspect column names
    crispr_headers = pd.read_csv(direct_links['CRISPRGeneDependency.csv'], nrows=0).columns.tolist()
    # print(f"CRISPRGeneDependency.csv headers: {crispr_headers}")
    model_col = crispr_headers[0]

    target_cols = [model_col]
    for col in crispr_headers[1:]:
        symbol = col.split(' ')[0]
        if symbol in genes_of_interest:
            target_cols.append(col)

    # 2. Load and Transform
    df_crispr = pd.read_csv(crispr_url, usecols=target_cols)
    df_crispr = df_crispr.rename(columns={model_col: 'ModelID'})
    new_column_names = {
        col: col.split(' ')[0] for col in df_crispr.columns if col != 'ModelID'
    }
    df_crispr = df_crispr.rename(columns=new_column_names)
    # print(df_crispr.head())

    # 3. Process OmicsExpression (Extremely Large)
    print("--- Processing Omics Expression (Filtering Genes) ---")
    omics_url = direct_links['OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv']
    omics_headers = pd.read_csv(omics_url, nrows=0).columns.tolist()
    # print(f"OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv headers: {omics_headers}")

    # 1. Define the metadata columns we need to keep
    metadata_cols = ['ModelID', 'IsDefaultEntryForModel']

    # 2. Identify target gene columns
    gene_cols = [c for c in omics_headers if c.split(' ')[0] in genes_of_interest]

    # 3. Load only what we need
    df_omics = pd.read_csv(omics_url, usecols=metadata_cols + gene_cols)

    # 4. FILTER: Keep only the default sequencing entries
    df_omics = df_omics[df_omics['IsDefaultEntryForModel'] == "Yes"].copy()

    # 5. Clean Column Names (as we did with CRISPR)
    df_omics.columns = [c.split(' ')[0] if '(' in c else c for c in df_omics.columns]

    # 6. Drop the flag column now that we've used it
    df_omics = df_omics.drop(columns=['IsDefaultEntryForModel'])

    # print(df_omics.head())

    print("--- Merging DepMap Datasets ---")

    # 1. Melt CRISPR to Long Format
    df_crispr_long = df_crispr.melt(id_vars='ModelID', var_name='gene', value_name='dependency')

    # 2. Melt Omics to Long Format
    df_omics_long = df_omics.melt(id_vars='ModelID', var_name='gene', value_name='expression')

    # 3. Merge the two long dataframes on ModelID and Gene
    # This aligns the dependency score and expression value for every gene-model pair
    df_merged = pd.merge(df_crispr_long, df_omics_long, on=['ModelID', 'gene'], how='inner')

    # 4. (Optional) Final merge with Model metadata to get Cell Line names and Lineage
    df_final = pd.merge(df_merged, df_models, on='ModelID', how='left')

    print(f"✅ Successfully fetched {len(df_final)} CRISPR and Omics dependencies from DepMap: The Cancer Dependency Map Project.")
    return df_final

"""### **SIDER Side Effect Resource**"""

def fetch_sider_data(drug_names):
    """
    Fetches Side Effects from SIDER (via MyChem.info API).
    Now includes a direct 'source_url' to the API result for transparency.
    """
    # Filter out empty or short names
    clean_drugs = [d for d in drug_names if len(d) > 2]
    if not clean_drugs: return pd.DataFrame()

    print(f"--- SIDER (via MyChem): Fetching side effects for {len(clean_drugs)} drugs ---")

    url = "http://mychem.info/v1/query"
    data = []

    # MyChem accepts batch queries
    payload = {
        "q": ",".join(clean_drugs),
        "scopes": "name,alias",
        "fields": "sider",
        "species": "human"
    }

    try:
        r = requests.post(url, data=payload, timeout=15)
        if r.status_code != 200:
            print(f"❌ API Error: {r.status_code}")
            return pd.DataFrame()

        results = r.json()

        for item in results:
            drug_query = item.get("query", "Unknown")

            if "sider" not in item:
                continue

            # Create a clickable URL that users can check
            # We encode the drug name to make it a valid URL
            safe_drug_name = urllib.parse.quote(drug_query)
            source_link = f"https://mychem.info/v1/query?q={safe_drug_name}"

            sider_data = item["sider"]
            if isinstance(sider_data, dict):
                sider_data = [sider_data]

            for entry in sider_data:
                effect_name = entry.get("side_effect", {}).get("name")
                if not effect_name: continue

                freq = entry.get("frequency", "unknown")

                data.append({
                    "drug": drug_query,
                    "side_effect": effect_name.lower(),
                    "frequency": freq,
                    "source_url": source_link  # <--- NEW: Direct Verification Link
                })

    except Exception as e:
        print(f"⚠️ SIDER Fetch Error: {e}")
        return pd.DataFrame()

    if not data:
        print("⚠️ No SIDER data found for these drugs.")
        return pd.DataFrame()

    df = pd.DataFrame(data).drop_duplicates()
    print(f"✅ Found {len(df)} side effect associations from SIDER database.")
    return df

"""### **DrugBank and Drug Central**"""

def fetch_drug_data_integrated(genes_of_interest):
    """
    1. Downloads DrugBank Open Vocabulary & Structures via direct URLs.
    2. Fetches clinical targets from DrugCentral.
    3. Merges all three into a single 'Consensus' DataFrame.
    """

    def download_and_extract(url, zip_name, target_folder):
        print(f"--- Downloading {zip_name} ---")
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(zip_name, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            with zipfile.ZipFile(zip_name, 'r') as zip_ref:
                zip_ref.extractall(target_folder)
            print(f"✅ Extracted to {target_folder}")
            return True
        else:
            print(f"❌ Failed to download {zip_name}. Status: {r.status_code}")
            return False

    # URLs provided
    vocab_url = "https://go.drugbank.com/releases/5-1-14/downloads/all-drugbank-vocabulary"
    struct_url = "https://go.drugbank.com/releases/5-1-14/downloads/all-open-structures"

    # Step 1: Download and Load DrugBank Data
    success_v = download_and_extract(vocab_url, "vocab.zip", "vocab_data")
    success_s = download_and_extract(struct_url, "struct.zip", "struct_data")

    if not (success_v and success_s):
        print("Required DrugBank files could not be downloaded.")
        return pd.DataFrame()

    # Step 2: Load DataFrames
    # Note: File names inside zips are usually 'drugbank vocabulary.csv' and 'open structures.sdf'
    df_vocab = pd.read_csv("vocab_data/drugbank vocabulary.csv")

    # Correctly read SDF file using RDKit
    sdf_file_path = "struct_data/open structures.sdf"
    supplier = Chem.SDMolSupplier(sdf_file_path)
    sdf_data = []
    for mol in supplier:
        if mol is not None:
            # Extract DRUGBANK_ID, SMILES, and InChI from the molecule properties
            db_id = mol.GetProp("DRUGBANK_ID") if mol.HasProp("DRUGBANK_ID") else None
            smiles = Chem.MolToSmiles(mol) if mol.HasProp("DRUGBANK_ID") else None # Ensure SMILES if ID exists
            inchi = Chem.MolToInchi(mol) if mol.HasProp("DRUGBANK_ID") else None # Ensure InChI if ID exists
            if db_id:
                sdf_data.append({'drugbank_id': db_id, 'smiles': smiles, 'inchi': inchi})
    df_struct = pd.DataFrame(sdf_data)

    # Step 3: Fetch DrugCentral Targets
    print("--- Fetching DrugCentral Targets ---")
    # Updated URL for DrugCentral drug-target data
    dc_url = "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"

    try:
        # Read gzipped TSV directly
        df_dc = pd.read_csv(dc_url, sep='\t', compression='gzip')
        # print(f"DrugCentral columns: {df_dc.columns.tolist()}") # Debugging line to show columns
        # print(df_dc.head(10))
        df_dc = df_dc[df_dc['GENE'].isin(genes_of_interest)]
    except Exception as e:
        print(f"❌ DrugCentral data loading failed: {e}")
        print(f"Please verify the URL: {dc_url} and its contents.")
        return pd.DataFrame() # Return empty DataFrame if API fails

    # Step 4: Merge Everything (The "Consensus" Logic)
    # Standardize names for matching
    df_dc['match_name'] = df_dc['DRUG_NAME'].str.lower()
    df_vocab['match_name'] = df_vocab['Common name'].str.lower()

    # Merge Targets with IDs (Vocabulary)
    merged = pd.merge(df_dc, df_vocab[['match_name', 'DrugBank ID', 'CAS', 'UNII', 'Synonyms']],
                      on='match_name', how='left')

    # Merge with Structures using DrugBank ID
    final_df = pd.merge(merged, df_struct[['drugbank_id', 'smiles', 'inchi']],
                        left_on='DrugBank ID', right_on='drugbank_id', how='left')

    # Cleanup
    final_df.rename(columns={'GENE': 'gene', 'DRUG_NAME': 'drug_name', 'TARGET_NAME': 'target_name', 'Synonyms': 'synonyms',
                             'TARGET_CLASS': 'target_class', 'DrugBank ID': 'db_id', 'ACT_TYPE': 'action'}, inplace=True)

    print(f"✅ Found {len(final_df)} structural data for {len(genes_of_interest)} genes from DrugBank and Drug Central.")
    return final_df[['gene', 'drug_name', 'db_id', 'synonyms', 'action', 'CAS', 'UNII', 'smiles', 'inchi']]