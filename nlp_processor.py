import spacy
import scispacy
from Bio import Entrez
from thefuzz import process
from scispacy.linking import EntityLinker
from sentence_transformers import SentenceTransformer, util

class LiteratureProcessor:
    def __init__(self, existing_nodes, existing_relations, email="you@example.com"):
        print("Initializing NLP Pipeline...")
        Entrez.email = email

        # 1. Load NER Model
        try:
            self.nlp = spacy.load("en_ner_bionlp13cg_md")
        except OSError:
            print("BioNLP missing. Fallback to generic.")
            if not spacy.util.is_package("en_core_web_sm"):
                spacy.cli.download("en_core_web_sm")
            self.nlp = spacy.load("en_core_web_sm")

        # 2. Load Semantic Model (Sentence Transformer)
        # 'all-MiniLM-L6-v2' is lightweight and perfect for short phrases
        print("Loading Sentence Transformer...")
        self.encoder = SentenceTransformer('all-MiniLM-L6-v2')

        # 3. Store Graph Ontology
        self.graph_nodes = existing_nodes
        self.graph_relations = existing_relations

        # Pre-compute relationship embeddings for fast semantic matching
        self.rel_embeddings = self.encoder.encode(self.graph_relations, convert_to_tensor=True)

        # Mapping spaCy biological entity types to our Graph Labels
        self.type_map = {
            "GENE_OR_GENE_PRODUCT": "Gene",
            "SIMPLE_CHEMICAL": "Drug",
            "CANCER": "Disease"
        }

    def fetch_papers(self, target_gene, diseases, max_results=5):
        """Fetches papers enforcing Gene + Disease presence."""
        disease_query = " OR ".join([f'"{d}"' for d in diseases])
        query = f'("{target_gene}") AND ({disease_query})'

        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            ids = record["IdList"]

            if not ids: return []

            handle = Entrez.efetch(db="pubmed", id=ids, rettype="fetch", retmode="xml")
            papers = Entrez.read(handle)
            results = []

            for paper in papers.get('PubmedArticle', []):
                try:
                    pmid = str(paper['MedlineCitation']['PMID'])
                    article = paper['MedlineCitation']['Article']
                    title = article.get('ArticleTitle', '')
                    abs_raw = article.get('Abstract', {}).get('AbstractText', [])
                    abstract = " ".join(abs_raw) if isinstance(abs_raw, list) else str(abs_raw)
                    results.append({"pmid": pmid, "text": f"{title}. {abstract}", "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
                except: continue
            return results
        except Exception as e:
            print(f"⚠️ PubMed API Error: {e}")
            return []

    def map_entity(self, entity_name, entity_type):
        """Fuzzy matches extracted entities to existing Graph nodes."""
        graph_label = self.type_map.get(entity_type)
        if not graph_label:
            return None, None # Ignore entity types we don't care about (e.g., 'Organism')

        existing_names = self.graph_nodes.get(graph_label, [])
        if not existing_names:
            return entity_name, graph_label # Graph is empty for this label, accept new name

        # Fuzzy Match
        best_match, score = process.extractOne(entity_name, existing_names)

        # If confidence is > 85%, use the existing graph node name
        if score > 85:
            return best_match, graph_label
        else:
            # It's a brand new entity! Return it so we can add it to the graph.
            return entity_name, graph_label

    def map_relation(self, extracted_verb):
        """Uses Sentence Transformer to find the closest semantic graph edge."""
        verb_embedding = self.encoder.encode(extracted_verb, convert_to_tensor=True)

        # Compare verb against all existing graph relations
        cosine_scores = util.cos_sim(verb_embedding, self.rel_embeddings)[0]
        best_match_idx = cosine_scores.argmax().item()
        best_score = cosine_scores[best_match_idx].item()

        # If semantic similarity is decent, map it. Otherwise, use a safe fallback.
        if best_score > 0.45:
            return self.graph_relations[best_match_idx]
        else:
            return "LINKED_WITH"

    def extract_triples(self, papers):
        """Runs NER, finds co-occurring entities, and maps to ontology."""
        triples = []

        for paper in papers:
            doc = self.nlp(paper["text"])

            # Sentence-level Relation Extraction
            for sent in doc.sents:

                # 1. Get entities natively from the sentence
                mapped_sent_ents = []
                for ent in sent.ents:
                    clean_name, graph_label = self.map_entity(ent.text, ent.label_)
                    if clean_name:
                        mapped_sent_ents.append({
                            "name": clean_name,
                            "label": graph_label
                        })

                # 2. If we have at least 2 entities in the same sentence, look for a verb
                if len(mapped_sent_ents) >= 2:
                    verbs = [token.lemma_ for token in sent if token.pos_ == "VERB"]
                    if not verbs:
                        continue

                    # Map the main verb to our graph ontology
                    mapped_rel = self.map_relation(verbs[0])

                    # 3. Create triples (Link the first entity to the others)
                    subject = mapped_sent_ents[0]
                    for obj in mapped_sent_ents[1:]:
                        # Prevent self-loops (e.g., MYC -> MYC)
                        if subject["name"] == obj["name"]:
                            continue

                        triples.append({
                            "subject_name": subject["name"],
                            "subject_label": subject["label"],
                            "relation": mapped_rel,
                            "object_name": obj["name"],
                            "object_label": obj["label"],
                            "pmid": paper["pmid"],
                            "evidence": paper["text"][:200] + "...",
                            "source_url": paper["source_url"]
                        })

        return triples