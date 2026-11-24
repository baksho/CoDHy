"""
This script connects to a Neo4j database, retrieves all nodes and relationships, and exports them to a JSON file.
The exported JSON file can then be imported into another Neo4j instance using the corresponding import script (see graph_import.py).
Make sure to update the connection details (URI, username, password) in the script before running it.

NOTE: This file is intended to be run in a local environment where you have access to the Neo4j Aura database.
"""

from neo4j import GraphDatabase
from neo4j.time import DateTime, Date, Time, Duration
import json

def neo4j_serializer(obj):
    if isinstance(obj, (DateTime, Date, Time, Duration)):
        return str(obj)
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

driver = GraphDatabase.driver(
    "neo4j+s://61abe14d.databases.neo4j.io",
    auth=("neo4j", "68Jl4KfoqSiBpvhZ9aalzwvgHqm298ud535muH1Tk8E")
)

with driver.session() as session:
    nodes = session.run("""
        MATCH (n) 
        RETURN id(n) as id, labels(n) as labels, properties(n) as props
    """).data()
    
    rels = session.run("""
        MATCH (a)-[r]->(b) 
        RETURN id(a) as from_id, id(b) as to_id, type(r) as type, properties(r) as props
    """).data()

with open("graph_export.json", "w") as f:
    json.dump({"nodes": nodes, "relationships": rels}, f, indent=2, default=neo4j_serializer)

print(f"Exported {len(nodes)} nodes and {len(rels)} relationships")
driver.close()