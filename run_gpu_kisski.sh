#!/bin/bash
#SBATCH -p kisski               # KISSKI partition (exclusive, not shared)
#SBATCH -t 2:00:00
#SBATCH --qos=2h
#SBATCH -G A100:4               # Request all 4 GPUs (required for kisski)
#SBATCH --job-name=gpu-test
#SBATCH -o CoDHy_output_%j.log   # Output saved to this file (%j = job ID)

echo "Job started on node: ${SLURM_NODELIST}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "--------------------------------------"

# Load required modules
# module load miniforge3
module load gcc
module load cuda

# Activate your conda environment
source $HOME/python_env_baksho/bin/activate       # ← Replace with your environment name

# Install dependencies (if not already installed in the environment)
pip install -r requirements.txt
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz

# Run the test script
srun python -u main.py                  # ← srun launches it on the compute node 

echo "--------------------------------------"
echo "Job finished."
