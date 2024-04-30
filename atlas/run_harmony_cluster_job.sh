#!/bin/bash
#SBATCH -J harmony
#SBATCH --partition=amd
#SBATCH -t 36:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB

ml any/python/3.9.9
source /home/scvi_venv/bin/activate
cd /home/immune_ageing

python harmony_cluster.py -i /home/immune_ageing_out/harmony_integrated.h5ad -o /home/immune_ageing_out