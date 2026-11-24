#!/bin/bash
#SBATCH -p grete:interactive
#SBATCH -t 2:00:00
#SBATCH -G 2g.10gb
#SBATCH --job-name=codhy-gpu
#SBATCH -o CoDHy_output_%j.log                  # Output saved to this file (%j = job ID)

echo "Job started on node: ${SLURM_NODELIST}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "--------------------------------------"

# Load required modules
module load gcc
module load cuda

# Activate Python environment (packages already installed on login node)
source $HOME/python_env_baksho/bin/activate

# Java + Neo4j
export JAVA_HOME=/user/suvinava.basak/u19797/jdk-17.0.19+10
export PATH=$JAVA_HOME/bin:$PATH
export NEO4J_HOME=/user/suvinava.basak/u19797/neo4j
export PATH=$NEO4J_HOME/bin:$PATH

# Start Neo4j (graph already built; this just starts the local server)
$NEO4J_HOME/bin/neo4j start
sleep 20

# Run CoDHy — no internet access needed here because:
#   • Neo4j graph is already built and stored on disk
#   • LLM weights are pre-downloaded by login_node_preprocess.py
#   • PubMed abstracts are pre-cached by login_node_preprocess.py
srun python -u main.py

$NEO4J_HOME/bin/neo4j stop

echo "--------------------------------------"
echo "Job finished."