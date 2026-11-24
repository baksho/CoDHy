# import networkx as nx
# from node2vec import Node2Vec
# from gensim.models import KeyedVectors

import os
import torch
import numpy as np
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory
from pykeen.models import RotatE, TransE

# Registry of supported KGE model classes, used both for training (via the
# pykeen pipeline's `model=` string) and for reconstructing the right class
# when reloading a checkpoint from disk.
KGE_MODEL_REGISTRY = {
    "RotatE": RotatE,
    "TransE": TransE,
}

class GraphEmbeddingAgent:
    """
    Layer 2 of the three-layer reasoning stack.
    Learns relational KG embeddings (RotatE by default, TransE optionally)
    and predicts missing (drug → gene → pathway → cancer) links.
    """

    def __init__(self, kg, model_name="RotatE"):
        if model_name not in KGE_MODEL_REGISTRY:
            raise ValueError(
                f"Unknown model_name='{model_name}'. "
                f"Valid options: {list(KGE_MODEL_REGISTRY.keys())}"
            )
        self.kg = kg
        self.model_name = model_name
        self.model      = None   # trained pykeen model (RotatE or TransE)
        self.tf         = None   # TriplesFactory (holds entity/relation index)
        self.embeddings = None   # kept for backward-compat with find_similar_nodes

    # ------------------------------------------------------------------ #
    #  Training                                                          #
    # ------------------------------------------------------------------ #

    def generate_embeddings(self, file_path="graph_embeddings_rotate"):
        print(f"\nEmbedding Agent ({self.model_name}) is learning relational structure...")

        # 1. Pull all triples from Neo4j (subject, relation_type, object)
        query = """
        MATCH (s)-[r]->(o)
        WHERE s.name IS NOT NULL AND o.name IS NOT NULL
          AND type(r) <> 'MENTIONED_IN'
        RETURN s.name AS h, type(r) AS rel, o.name AS t
        """
        with self.kg.driver.session() as session:
            rows = session.run(query).data()

        if not rows:
            print("Graph has no valid edges. Cannot train RotatE.")
            return

        # Build numpy array of string triples  [N, 3]
        triples_np = np.array([[r["h"], r["rel"], r["t"]] for r in rows], dtype=str)
        print(f"   -> {len(triples_np)} triples loaded for training.")

        self.tf = TriplesFactory.from_labeled_triples(triples_np)

        # 2. Train via pykeen pipeline (model name is a plain string key,
        #    e.g. "RotatE" or "TransE" — see KGE_MODEL_REGISTRY above)
        result = pipeline(
            training=self.tf,
            testing=self.tf,
            model=self.model_name,
            model_kwargs={"embedding_dim": 64},
            training_kwargs={"num_epochs": 100, "batch_size": 512},
            optimizer="Adam",
            optimizer_kwargs={"lr": 0.001},
            random_seed=42,
            device="cpu",
        )
        self.model = result.model
        self.model.eval()
        print(f"✅ {self.model_name} embeddings trained successfully.")

        # Save for reuse — model_name is recorded so load_embeddings() can
        # reconstruct the correct pykeen model class later.
        save_path = file_path + ".pt"
        torch.save({
            "model_name": self.model_name,
            "model_state": self.model.state_dict(),
            "entity_to_id": self.tf.entity_to_id,
            "relation_to_id": self.tf.relation_to_id,
        }, save_path)
        print(f"   -> Saved to {save_path}")


    def load_embeddings(self, file_path="graph_embeddings_rotate"):
        """Reload a previously trained model without retraining."""
        pt_path = file_path + ".pt"
        if not os.path.exists(pt_path):
            print(f"No saved model found at {pt_path}. Train first.")
            return
        checkpoint = torch.load(pt_path, map_location="cpu")
        entity_to_id   = checkpoint["entity_to_id"]
        relation_to_id = checkpoint["relation_to_id"]

        # Older checkpoints (saved before TransE support was added) didn't
        # record model_name — assume RotatE for those, matching the original
        # hardcoded behaviour.
        saved_model_name = checkpoint.get("model_name", "RotatE")
        self.model_name = saved_model_name

        model_cls = KGE_MODEL_REGISTRY.get(saved_model_name)
        if model_cls is None:
            raise ValueError(
                f"Checkpoint was trained with unknown model '{saved_model_name}'. "
                f"Valid options: {list(KGE_MODEL_REGISTRY.keys())}"
            )

        # # Rebuild a minimal TriplesFactory stub so pykeen can reconstruct the model
        # n_entities   = len(entity_to_id)
        # n_relations  = len(relation_to_id)
        # self.model = RotatE(
        #     triples_factory=None,
        #     embedding_dim=64,
        #     num_entities=n_entities,
        #     num_relations=n_relations,
        # )

        mapped_triples = torch.empty((0, 3), dtype=torch.long)
        tf_stub = TriplesFactory(
            mapped_triples=mapped_triples,
            entity_to_id=entity_to_id,
            relation_to_id=relation_to_id,
        )
        self.model = model_cls(triples_factory=tf_stub, embedding_dim=64)
        self.model.load_state_dict(checkpoint["model_state"])
        self.model.eval()

        # Recreate lookup maps
        self.tf = type("_TF", (), {
            "entity_to_id":   entity_to_id,
            "relation_to_id": relation_to_id,
            "entity_id_to_label":   {v: k for k, v in entity_to_id.items()},
            "relation_id_to_label": {v: k for k, v in relation_to_id.items()},
        })()
        print(f"✅ {saved_model_name} model loaded from {pt_path}.")


    # ------------------------------------------------------------------ #
    #  Link prediction (Layer 2)                                         #
    # ------------------------------------------------------------------ #

    def predict_missing_links(self, head_entity, relation_type, top_k=5):
        """
        Given a head entity and relation type, returns the top_k most likely
        tail entities that RotatE predicts as plausible completions.

        Example:
            predict_missing_links("Gefitinib", "INHIBITS", top_k=5)
            → [{"entity": "EGFR", "score": 0.92}, ...]
        """
        if self.model is None or self.tf is None:
            return []

        eid = self.tf.entity_to_id.get(head_entity)
        rid = self.tf.relation_to_id.get(relation_type)
        if eid is None or rid is None:
            return []

        with torch.no_grad():
            # Score all possible tail entities.
            hr_batch = torch.tensor([[eid, rid]])
            scores = self.model.score_t(hr_batch).squeeze()
            scores = torch.sigmoid(scores).numpy()


        id_to_entity = {v: k for k, v in self.tf.entity_to_id.items()}
        top_indices  = scores.argsort()[::-1][:top_k + 1]

        results = []
        for idx in top_indices:
            entity = id_to_entity.get(int(idx), "")
            if entity and entity != head_entity:
                results.append({"entity": entity, "score": float(scores[idx])})
            if len(results) >= top_k:
                break
        return results

    def build_predicted_paths(self, focus_gene, cancer_type, top_k=5):
        """
        Uses RotatE to discover hidden drug→gene→pathway→cancer chains
        and returns them in the same dict format as Layer 1 paths,
        tagged with layer='PREDICTED' and an embedding_score.

        The prediction chain tried is:
          Drug -[INHIBITS]→ focus_gene  (drug targets the gene)
          focus_gene -[IN_PATHWAY]→ Pathway
          Pathway -[ASSOCIATED_WITH]→ cancer_type
        """
        if self.model is None:
            return []

        predicted_paths = []

        # Step A: which drugs are predicted to inhibit focus_gene?
        candidate_drugs = self.predict_missing_links(focus_gene, "INHIBITS", top_k=top_k)

        for drug_hit in candidate_drugs:
            drug_name   = drug_hit["entity"]
            drug_score  = drug_hit["score"]

            # Step B: which pathways does focus_gene likely activate?
            pathway_hits = self.predict_missing_links(focus_gene, "IN_PATHWAY", top_k=3)

            for pw_hit in pathway_hits:
                pw_name  = pw_hit["entity"]
                pw_score = pw_hit["score"]

                combined_score = (drug_score + pw_score) / 2.0

                predicted_paths.append({
                    "subject":    drug_name,
                    "relation":   "INHIBITS",
                    "object":     focus_gene,
                    "layer":      "PREDICTED",
                    "confidence": drug_score,
                    "source":     self.model_name,
                    "source_url": "",
                    # Extra fields that capture the multi-hop chain
                    "intermediate": focus_gene,
                    "pathway":      pw_name,
                    "cancer":       cancer_type,
                    "embedding_score": combined_score,
                })

        return predicted_paths

    
    # ------------------------------------------------------------------ #
    #  Backward compatibility                                            #
    # ------------------------------------------------------------------ #

    def find_similar_nodes(self, target_node, top_k=5):
        """
        Kept for backward compatibility with hgenerator.py implicit context block.
        Returns cosine-nearest entities in embedding space.
        """
        if self.model is None or self.tf is None:
            return []
        eid = self.tf.entity_to_id.get(target_node)
        if eid is None:
            return []

        with torch.no_grad():
            all_ids     = torch.arange(len(self.tf.entity_to_id))
            entity_embs = self.model.entity_representations[0](all_ids)
            if entity_embs.is_complex():
                entity_embs = torch.view_as_real(entity_embs).flatten(start_dim=1)
            query_emb = entity_embs[eid].unsqueeze(0)
            sims = torch.nn.functional.cosine_similarity(query_emb, entity_embs)
            sims[eid] = -1.0  # exclude self
            top = sims.topk(top_k).indices.tolist()

        id_to_entity = {v: k for k, v in self.tf.entity_to_id.items()}
        return [{"name": id_to_entity[i], "score": float(sims[i])} for i in top]