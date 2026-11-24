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
            print("Loaded BioNLP model for NER.")
        except OSError:
            print("BioNLP missing. Fallback to generic.")
            if not spacy.util.is_package("en_core_web_sm"):
                spacy.cli.download("en_core_web_sm")
            self.nlp = spacy.load("en_core_web_sm")

        # Add UMLS Entity Linker
        # self.nlp.add_pipe(
        #     "scispacy_linker",
        #     config={
        #         "linker_name": "mesh",
        #         "resolve_abbreviations": True
        #     }
        # )

        # self.linker = self.nlp.get_pipe("scispacy_linker")

        # 2. Load Semantic Model (Sentence Transformer)
        print("Loading Sentence Transformer...")
        # self.encoder = SentenceTransformer('all-MiniLM-L6-v2')
        self.encoder = SentenceTransformer("BAAI/bge-base-en-v1.5")
        # self.encoder = SentenceTransformer("BAAI/bge-large-en-v1.5")

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
        """
        Fetches papers enforcing Gene + Disease presence.

        `diseases` may be a single disease string (e.g. "Breast Invasive
        Carcinoma") or an iterable of disease strings. A bare string is
        wrapped in a list — iterating over it directly would otherwise walk
        it character-by-character and silently build a garbage query.
        """
        if isinstance(diseases, str):
            diseases = [diseases]
        # disease_query = " OR ".join([f'"{d}"' for d in diseases])
        # query = f'("{target_gene}") AND ({disease_query})'
        disease_query = " OR ".join(diseases)
        query = f'({target_gene}) AND ({disease_query})'

        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            ids = record["IdList"]
            # print(f"Query: {query}")
            # print(f"PubMed interpreted query as: {record.get('QueryTranslation', 'N/A')}")
            # print(f"Found {len(ids)} PubMed IDs: {ids}")

            if not ids: return []

            handle = Entrez.efetch(db="pubmed", id=ids, rettype="fetch", retmode="xml")
            papers = Entrez.read(handle)
            print(f"Articles fetched: {len(papers.get('PubmedArticle', []))}")
            results = []

            for paper in papers.get('PubmedArticle', []):
                try:
                    pmid = str(paper['MedlineCitation']['PMID'])
                    article = paper['MedlineCitation']['Article']
                    title = article.get('ArticleTitle', '')
                    abs_raw = article.get('Abstract', {}).get('AbstractText', [])
                    abstract = " ".join(abs_raw) if isinstance(abs_raw, list) else str(abs_raw)
                    results.append({"pmid": pmid, "text": f"{title}. {abstract}", "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
                except Exception as e:
                    print(f"Error parsing paper: {e}")
                    continue
            
            print(f"Results extracted: {len(results)}")
            return results
            
        except Exception as e:
            print(f"⚠️ PubMed API Error: {e}")
            return []

    # def map_entity(self, entity_name, entity_type):
    #     """Fuzzy matches extracted entities to existing Graph nodes."""
    #     graph_label = self.type_map.get(entity_type)
    #     if not graph_label:
    #         return None, None # Ignore entity types we don't care about (e.g., 'Organism')

    #     existing_names = self.graph_nodes.get(graph_label, [])
    #     if not existing_names:
    #         return entity_name, graph_label # Graph is empty for this label, accept new name

    #     # Fuzzy Match
    #     best_match, score = process.extractOne(entity_name, existing_names)

    #     # If confidence is > 85%, use the existing graph node name
    #     if score > 85:
    #         return best_match, graph_label
    #     else:
    #         # It's a brand new entity! Return it so we can add it to the graph.
    #         return entity_name, graph_label

    def map_entity(self, ent):
        """
        Fuzzy matches a spaCy entity to an existing graph node using
        thefuzz, instead of the UMLS linker (removed — it required
        loading ~2GB of KB/index data into memory, which OOMs on the
        standard free-tier Colab runtime).
        """
        graph_label = self.type_map.get(ent.label_)
        if not graph_label:
            return None, None       # Ignore entity types we don't care about (e.g., 'Organism')

        # # If linker couldn't map it
        # if not ent._.kb_ents:
        #     return ent.text, graph_label

        # # Highest-confidence UMLS candidate
        # umls_id, confidence = ent._.kb_ents[0]

        existing_names = self.graph_nodes.get(graph_label, [])
        if not existing_names:
            return ent.text, graph_label # Graph is empty for this label, accept new name

        # Fuzzy match against existing graph nodes
        best_match, score = process.extractOne(ent.text, existing_names)
        # canonical_name = linker.kb.cui_to_entity[umls_id].canonical_name

        if score > 85:
            return best_match, graph_label
        else:
            # It's a brand new entity! Return it so we can add it to the graph.
            return ent.text, graph_label

        # for node in existing_names:
        #     if node.lower() == canonical_name.lower():
        #         return node, graph_label

        # # New entity
        # return canonical_name, graph_label

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

    # ------------------------------------------------------------------ #
    #  Dependency-path relation extraction (agenda item 4)               #
    # ------------------------------------------------------------------ #
    # Replaces the previous heuristic (first-verb-in-sentence linked to
    # every entity pair, regardless of syntactic role) with a proper
    # shortest-dependency-path approach: for a pair of entities, walk up
    # each one's ancestor chain in the dependency tree and look for a
    # lowest common ancestor that is a VERB. This is a standard relation-
    # extraction technique and only requires spaCy/scispaCy's regular
    # dependency parse (token.head / token.dep_ / token.pos_) — no new
    # dependency needed. Verified against a hand-built dependency tree
    # (see Phase 2 changelog) since no biomedical parser model could be
    # downloaded in the dev sandbox; please spot-check a few extracted
    # triples once this runs against real abstracts in Colab.
    #
    # This directly fixes a real correctness gap in the old heuristic:
    # it had no way to tell "Gefitinib inhibits EGFR" apart from
    # "Gefitinib does NOT inhibit EGFR" — both produced an identical
    # INHIBITS-mapped triple. Negated relations are now detected and
    # dropped instead of silently extracted as if affirmative.

    @staticmethod
    def _ancestors_with_self(token):
        chain = [token]
        seen = {token.i}
        while token.head.i != token.i and token.head.i not in seen:
            token = token.head
            chain.append(token)
            seen.add(token.i)
        return chain

    def _find_connecting_verb(self, ent_a_root, ent_b_root):
        """Lowest common ancestor of two entity heads, if it's a VERB."""
        anc_a = self._ancestors_with_self(ent_a_root)
        anc_b_by_i = {t.i: t for t in self._ancestors_with_self(ent_b_root)}
        for tok in anc_a:
            if tok.i in anc_b_by_i and tok.pos_ == "VERB":
                return tok
        return None

    @staticmethod
    def _is_negated(verb_token):
        return any(child.dep_ == "neg" for child in verb_token.children)

    def extract_triples(self, papers):
        """Runs NER, then dependency-path relation extraction, then maps to ontology."""
        triples = []

        for paper in papers:
            doc = self.nlp(paper["text"])

            for sent in doc.sents:
                # 1. Get entities natively from the sentence, keeping the
                #    spaCy Span object (not just its text) so we can reach
                #    its syntactic head (.root) for dependency-tree walks.
                mapped_sent_ents = []
                for ent in sent.ents:
                    # clean_name, graph_label = self.map_entity(ent.text, ent.label_)
                    clean_name, graph_label = self.map_entity(ent)
                    if clean_name:
                        mapped_sent_ents.append({
                            "name": clean_name,
                            "label": graph_label,
                            "root": ent.root,
                        })

                if len(mapped_sent_ents) < 2:
                    continue

                sentence_had_dependency_triple = False

                # 2. Try dependency-path extraction for every UNORDERED pair.
                #    (Bug fix: iterating ordered pairs (i,j) and (j,i) both
                #    resolve to the same lowest-common-ancestor verb and the
                #    same correct subject/object direction via
                #    _is_syntactic_subject — so trying both orderings just
                #    emitted every real relation twice. One pass over i<j is
                #    correct and sufficient; direction is decided inside the
                #    loop body, not by iteration order.)
                for i in range(len(mapped_sent_ents)):
                    for j in range(i + 1, len(mapped_sent_ents)):
                        ent_a, ent_b = mapped_sent_ents[i], mapped_sent_ents[j]
                        if ent_a["name"] == ent_b["name"]:
                            continue

                        verb_tok = self._find_connecting_verb(ent_a["root"], ent_b["root"])
                        if verb_tok is None:
                            continue
                        if self._is_negated(verb_tok):
                            # Correctness fix: a negated relation ("does NOT
                            # inhibit") must not be extracted as if it were
                            # the corresponding affirmative INHIBITS edge.
                            continue

                        # Direction: prefer the entity that is the verb's
                        # syntactic subject as the triple's subject.
                        subj_is_a = self._is_syntactic_subject(ent_a["root"], verb_tok)
                        subject, obj = (ent_a, ent_b) if subj_is_a else (ent_b, ent_a)

                        mapped_rel = self.map_relation(verb_tok.lemma_)
                        triples.append({
                            "subject_name": subject["name"],
                            "subject_label": subject["label"],
                            "relation": mapped_rel,
                            "object_name": obj["name"],
                            "object_label": obj["label"],
                            "pmid": paper["pmid"],
                            "evidence": paper["text"][:200] + "...",
                            "source_url": paper["source_url"],
                            "extraction_method": "dependency-parsed",
                        })
                        sentence_had_dependency_triple = True

                # 3. Fallback: if no pair in this sentence had a clean
                #    dependency-verb link (e.g. parser limitations, elided
                #    verbs, appositives), fall back to the old co-occurrence
                #    heuristic — but tag it as lower-confidence so downstream
                #    consumers (and the paper's methodology section) can
                #    distinguish the two.
                if not sentence_had_dependency_triple:
                    verbs = [token for token in sent if token.pos_ == "VERB"]
                    if not verbs:
                        continue
                    if self._is_negated(verbs[0]):
                        continue
                    mapped_rel = self.map_relation(verbs[0].lemma_)
                    subject = mapped_sent_ents[0]
                    for obj in mapped_sent_ents[1:]:
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
                            "source_url": paper["source_url"],
                            "extraction_method": "cooccurrence-fallback",
                        })

        return triples

    @staticmethod
    def _is_syntactic_subject(ent_root, verb_tok):
        """
        True if ent_root is the *semantic* subject (the actor) of verb_tok,
        False if it's the semantic object (the thing acted upon).

        This distinguishes active from passive voice:
          - "Gefitinib inhibits EGFR"        -> Gefitinib: nsubj      -> True
          - "EGFR is inhibited by Gefitinib" -> EGFR: nsubjpass       -> False (!)
                                                 Gefitinib: pobj of the
                                                 "agent" (by-phrase)   -> True
        Bug fixed: nsubjpass looks like a "subject" dependency label, but
        it marks the passive SUBJECT position, whose semantic role is the
        OBJECT of the action. Treating it as the semantic subject (the
        original naive version of this function did) silently reverses the
        relation direction for any passive-voice sentence — and passive
        voice ("X is inhibited by Y", "Y is associated with X") is very
        common in biomedical abstracts, so this would have been a
        systematic, not rare, source of backwards triples.
        """
        tok = ent_root
        seen = {tok.i}
        while True:
            if tok.head.i == verb_tok.i:
                if tok.dep_ == "nsubj":
                    return True                  # active subject
                if tok.dep_ == "nsubjpass":
                    return False                 # passive subject -> semantic object
                if tok.dep_ in ("dobj", "attr", "oprd"):
                    return False                  # active object
            # Passive agent phrase: "by Gefitinib" — Gefitinib is pobj of
            # "by", and "by" itself is attached to the verb with dep_="agent".
            if (tok.dep_ == "pobj" and tok.head.dep_ == "agent"
                    and tok.head.head.i == verb_tok.i):
                return True                       # the actor in a by-phrase
            if tok.head.i == tok.i or tok.head.i in seen:
                break
            seen.add(tok.head.i)
            tok = tok.head
        # Ambiguous / couldn't determine — default to treating ent_a as subject
        # (preserves old behaviour's bias rather than guessing differently).
        return True