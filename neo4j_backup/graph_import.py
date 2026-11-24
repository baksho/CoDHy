"""
This script reads a JSON file containing nodes and relationships exported from a Neo4j database and imports them into another Neo4j instance.
Make sure to update the connection details (URI, username, password) in the script before running it.
The JSON file should be in the format produced by the corresponding export script (see graph_export_from_neo4j_local.py).
"""
from neo4j import GraphDatabase
import json

driver = GraphDatabase.driver("bolt://127.0.0.1:7687", auth=None)

with open("/user/suvinava.basak/u19797/CoDHy_v2/neo4j_backup/graph_export.json") as f:
    data = json.load(f)

with driver.session() as session:
    # Import nodes
    print("Importing nodes...")
    for node in data["nodes"]:
        labels = ":".join(node["labels"])
        session.run(f"""
            CREATE (n:{labels})
            SET n = $props
            SET n._import_id = $id
        """, props=node["props"], id=node["id"])

    # Import relationships
    print("Importing relationships...")
    for rel in data["relationships"]:
        session.run(f"""
            MATCH (a {{_import_id: $from_id}})
            MATCH (b {{_import_id: $to_id}})
            CREATE (a)-[r:{rel['type']}]->(b)
            SET r = $props
        """, from_id=rel["from_id"], to_id=rel["to_id"], props=rel["props"])

    # Clean up temp import IDs
    session.run("MATCH (n) REMOVE n._import_id")

print("Import complete.")
driver.close()