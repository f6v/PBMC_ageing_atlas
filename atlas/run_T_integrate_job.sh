#!/bin/bash
#SBATCH -J scvi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tesla:1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB

ml any/python/3.9.9
source /home/scvi_venv/bin/activate
cd /home/ihor/immune_ageing/

python T_integrate.py -i /home/immune_ageing_in/prepared.h5ad -m /home/immune_ageing_in/annotated_metadata.tsv -o /home/immune_ageing_out