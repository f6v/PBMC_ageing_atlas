#!/bin/bash
#SBATCH -J scanpy
#SBATCH --partition=amd
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB

ml any/python/3.9.9
source /home/scvi_venv/bin/activate
cd /home/immune_ageing/

python no_integration.py -i /home/immune_ageing_in/prepared.h5ad -o //home/immune_ageing_out